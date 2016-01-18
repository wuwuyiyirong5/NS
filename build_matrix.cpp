#include "ISOP2P1.h"
#define DIM 2

void ISOP2P1::buildMatrix()
{
	/// 计算一下各空间自由度和总自由度.
	int n_dof_v = fem_space_v.n_dof();
	int n_dof_p = fem_space_p.n_dof();
	int n_total_dof =  2 * n_dof_v + n_dof_p;
	if (n_total_dof != sp_stokes.n_rows())
	{
		std::cerr << "ERROR: the demision of matrix is not correct!" << std::endl;
		exit(-1);
	}
	else
		std::cout << "dof no. of v: " << n_dof_v << ", "
		<< "dof no. of p: " << n_dof_p << ", "
		<< "total dof no.: " << n_total_dof << std::endl;

	/// 构建系数矩阵和右端项.
	mat_v_stiff.reinit(sp_vxvx);
	mat_v_mass.reinit(sp_vxvx);
	mat_pvx_divT.reinit(sp_pvx);
	mat_pvy_divT.reinit(sp_pvy);

	FEMSpace<double, DIM>::ElementIterator the_element_v = fem_space_v.beginElement();
	FEMSpace<double, DIM>::ElementIterator end_element_v = fem_space_v.endElement();

	/// 遍历速度单元, 拼装相关系数矩阵和右端项.
	for (the_element_v = fem_space_v.beginElement();
			the_element_v != end_element_v; ++the_element_v)
	{
		/// 当前单元信息.
		double volume = the_element_v->templateElement().volume();
		/// 积分精度, u 和 p 都是 1 次, 梯度和散度 u 都是常数. 因此矩阵拼
		/// 装时积分精度不用超过 1 次. (验证一下!)
		const QuadratureInfo<DIM>& quad_info = the_element_v->findQuadratureInfo(3);
		std::vector<double> jacobian = the_element_v->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element_v->local_to_global(quad_info.quadraturePoint());
		/// 速度单元信息.
		std::vector<std::vector<std::vector<double> > > basis_gradient_v = the_element_v->basis_function_gradient(q_point);
		std::vector<std::vector<double> >  basis_value_v = the_element_v->basis_function_value(q_point);
		const std::vector<int>& element_dof_v = the_element_v->dof();
		int n_element_dof_v = the_element_v->n_dof();
		// std::cout << the_element_v->index() << std::endl;
		/// 压力单元信息.
		Element<double, DIM> &p_element = fem_space_p.element(index_ele_v2p[the_element_v->index()]);

		const std::vector<int>& element_dof_p = p_element.dof();
		std::vector<std::vector<double> >  basis_value_p = p_element.basis_function_value(q_point);
		int n_element_dof_p = p_element.n_dof();
		/// 实际拼装.
		for (int l = 0; l < n_quadrature_point; ++l)
		{
			double Jxw = quad_info.weight(l) * jacobian[l] * volume;
			for (int i = 0; i < n_element_dof_v; ++i)
			{
				for (int j = 0; j < n_element_dof_v; ++j)
				{
					double cont = Jxw * innerProduct(basis_gradient_v[i][l], basis_gradient_v[j][l]);
					/// V Stiff
					mat_v_stiff.add(element_dof_v[i], element_dof_v[j], cont);
					cont = Jxw * basis_value_v[i][l] * basis_value_v[j][l];
					/// V Mass
					mat_v_mass.add(element_dof_v[i], element_dof_v[j], cont);
				}
				for (int j = 0; j < n_element_dof_p; ++j)
				{
					/// DivT x
					double cont = -Jxw * (basis_gradient_v[i][l][0] * basis_value_p[j][l]);
					mat_pvx_divT.add(element_dof_v[i], element_dof_p[j], cont);

					/// DivT y
					cont = -Jxw * (basis_gradient_v[i][l][1] * basis_value_p[j][l]);
					mat_pvy_divT.add(element_dof_v[i], element_dof_p[j], cont);
				}
			}
		}
	}

	/// 构建系数矩阵和右端项.
	mat_p_mass.reinit(sp_mass_p);
	mat_p_stiff.reinit(sp_mass_p);
	mat_vxp_div.reinit(sp_vxp);
	mat_vyp_div.reinit(sp_vyp);
	FEMSpace<double, DIM>::ElementIterator the_element_p = fem_space_p.beginElement();
	FEMSpace<double, DIM>::ElementIterator end_element_p = fem_space_p.endElement();

	/// 遍历压力单元. 拼装矩阵和右端项.
	for (the_element_p = fem_space_p.beginElement();
	     the_element_p != end_element_p; ++the_element_p)
	{
		/// 当前单元信息.
		double volume_p = the_element_p->templateElement().volume();
		const QuadratureInfo<DIM>& quad_info_p = the_element_p->findQuadratureInfo(3);
		std::vector<double> jacobian_p = the_element_p->local_to_global_jacobian(quad_info_p.quadraturePoint());
		int n_quadrature_point = quad_info_p.n_quadraturePoint();
		std::vector<Point<DIM> > q_point_p = the_element_p->local_to_global(quad_info_p.quadraturePoint());
		/// 压力单元信息.
		std::vector<std::vector<double> >  basis_value_p = the_element_p->basis_function_value(q_point_p);
		std::vector<std::vector<std::vector<double> > >  basis_gradient_p = the_element_p->basis_function_gradient(q_point_p);
		const std::vector<int>& element_dof_p = the_element_p->dof();
		int n_element_dof_p = the_element_p->n_dof();
		for (int i = 0; i < n_element_dof_p; ++i)
		{
			for (int l = 0; l < n_quadrature_point; ++l)
			{
				for (int j = 0; j < n_element_dof_p; ++j)
				{
					double Jxw = quad_info_p.weight(l) * jacobian_p[l] * volume_p;
					double cont = Jxw * basis_value_p[i][l] * basis_value_p[j][l];
					mat_p_mass.add(element_dof_p[i], element_dof_p[j], cont);
					/// mat_p_stiff.
					cont = Jxw * innerProduct(basis_gradient_p[i][l], basis_gradient_p[j][l]);
					mat_p_stiff.add(element_dof_p[i], element_dof_p[j], cont);
				}
			}
			int idx_p = the_element_p->index();
			int n_chi = index_ele_p2v[idx_p].size();
			for (int k = 0; k < n_chi; k++)
			{
				/// 速度单元信息.
				Element<double, DIM> &v_element = fem_space_v.element(index_ele_p2v[idx_p][k]);

				/// 几何信息.
				double volume = v_element.templateElement().volume();
				const QuadratureInfo<DIM>& quad_info = v_element.findQuadratureInfo(3);
				std::vector<double> jacobian = v_element.local_to_global_jacobian(quad_info.quadraturePoint());
				int n_quadrature_point = quad_info.n_quadraturePoint();
				std::vector<Point<DIM> > q_point = v_element.local_to_global(quad_info.quadraturePoint());

				const std::vector<int>& element_dof_v = v_element.dof();
				std::vector<std::vector<std::vector<double> > > basis_gradient_v = v_element.basis_function_gradient(q_point);
				int n_element_dof_v = v_element.n_dof();

				/// 压力单元信息.
				std::vector<std::vector<std::vector<double> > > basis_gradient_p = the_element_p->basis_function_gradient(q_point);
				std::vector<std::vector<double> >  basis_value_p = the_element_p->basis_function_value(q_point);
				/// 具体拼装.
				for (int l = 0; l < n_quadrature_point; ++l)
				{
					double Jxw = quad_info.weight(l) * jacobian[l] * volume;

					for (int j = 0; j < n_element_dof_v; ++j)
					{
						/// Div x
						double cont = -Jxw * basis_value_p[i][l] * basis_gradient_v[j][l][0];
						mat_vxp_div.add(element_dof_p[i], element_dof_v[j], cont);
						/// Div y
						cont = -Jxw * basis_value_p[i][l] * basis_gradient_v[j][l][1];
						mat_vyp_div.add(element_dof_p[i], element_dof_v[j], cont);
					}
				}
			}
		}
	}

	// 输出测试矩阵, 只限小规模矩阵.
//	outputMat("Axx",mat_v_stiff);
//	outputMat("Mxx",mat_v_mass);
//	outputMat("Apx",mat_pvx_divT);
//	outputMat("Apy",mat_pvy_divT);
//	outputMat("Axp",mat_vxp_div);
//	outputMat("Ayp",mat_vyp_div);

	std::cout << "Basic matrixes builded." << std::endl;
};

void ISOP2P1::updateNonlinearMatrix()
{
	int n_dof_v = fem_space_v.n_dof();
	int n_dof_p = fem_space_p.n_dof();
	int n_total_dof = 2 * n_dof_v + n_dof_p;

	/// 更新对流矩阵块.
	mat_v_convection.reinit(sp_vxvx);
	/// 更新 Jacobi 矩阵诸块.
	mat_v_Jacobi_xx.reinit(sp_vxvx);
	mat_v_Jacobi_xy.reinit(sp_vyvx);
	mat_v_Jacobi_yx.reinit(sp_vxvy);
	mat_v_Jacobi_yy.reinit(sp_vxvx);

	FEMSpace<double,2>::ElementIterator the_element_v = fem_space_v.beginElement();
	FEMSpace<double,2>::ElementIterator end_element_v = fem_space_v.endElement();
	FEMSpace<double,2>::ElementIterator the_element_p = fem_space_p.beginElement();
	FEMSpace<double,2>::ElementIterator end_element_p = fem_space_p.endElement();
	/// 遍历速度单元, 拼装相关系数矩阵和右端项.
	for (the_element_v = fem_space_v.beginElement();
	     the_element_v != end_element_v; ++the_element_v)
	{
		/// 当前单元信息.
		double volume = the_element_v->templateElement().volume();
		/// 积分精度, u 和 p 都是 1 次, 梯度和散度 u 都是常数. 因此矩阵拼
		/// 装时积分精度不用超过 1 次. (验证一下!)
		const QuadratureInfo<DIM>& quad_info = the_element_v->findQuadratureInfo(4);
		std::vector<double> jacobian
		= the_element_v->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element_v->local_to_global(quad_info.quadraturePoint());
		/// 速度单元信息.
		std::vector<std::vector<std::vector<double> > > basis_gradient_v
		= the_element_v->basis_function_gradient(q_point);
		std::vector<std::vector<double> >  basis_value_v
		= the_element_v->basis_function_value(q_point);
		const std::vector<int>& element_dof_v = the_element_v->dof();
		int n_element_dof_v = the_element_v->n_dof();
		std::vector<double> vx_value = v_h[0].value(q_point, *the_element_v);
		std::vector<double> vy_value = v_h[1].value(q_point, *the_element_v);
		std::vector<std::vector<double> > vx_gradient = v_h[0].gradient(q_point, *the_element_v);
		std::vector<std::vector<double> > vy_gradient = v_h[1].gradient(q_point, *the_element_v);
		/// 压力单元信息.
		Element<double, DIM> &p_element = fem_space_p.element(index_ele_v2p[the_element_v->index()]);
		const std::vector<int>& element_dof_p = p_element.dof();
		std::vector<std::vector<std::vector<double> > > basis_gradient_p
		= p_element.basis_function_gradient(q_point);
		std::vector<std::vector<double> >  basis_value_p = p_element.basis_function_value(q_point);
		std::vector<double> p_value = p_h.value(q_point, p_element);
		int n_element_dof_p = p_element.n_dof();
		/// 实际拼装.
		for (int l = 0; l < n_quadrature_point; ++l)
		{
			double Jxw = quad_info.weight(l) * jacobian[l] * volume;
			for (int i = 0; i < n_element_dof_v; ++i)
			{
				for (int j = 0; j < n_element_dof_v; ++j)
				{
					double cont = Jxw * (vx_value[l] * basis_gradient_v[j][l][0] +
							vy_value[l] * basis_gradient_v[j][l][1])
					* basis_value_v[i][l];
					mat_v_convection.add(element_dof_v[i], element_dof_v[j], cont);
					cont = Jxw * vx_gradient[l][0] * basis_value_v[j][l] * basis_value_v[i][l];
					mat_v_Jacobi_xx.add(element_dof_v[i], element_dof_v[j], cont);
					/// (0, 1).
					cont = Jxw * vx_gradient[l][1] * basis_value_v[j][l] * basis_value_v[i][l];
					mat_v_Jacobi_xy.add(element_dof_v[i], element_dof_v[j], cont);
					/// (1, 0).
					cont = Jxw * vy_gradient[l][0] * basis_value_v[j][l] * basis_value_v[i][l];
					mat_v_Jacobi_yx.add(element_dof_v[i], element_dof_v[j], cont);
					/// (1, 1).
					cont = Jxw * vy_gradient[l][1] * basis_value_v[j][l] * basis_value_v[i][l];
					mat_v_Jacobi_yy.add(element_dof_v[i], element_dof_v[j], cont);
				}
			}
		}
	}
};
#undef DIM
