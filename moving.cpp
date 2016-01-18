#include "ISOP2P1.h"
#include "preconditioner.h"
#include "functions.h"
#define DIM 2
void ISOP2P1::syncMesh()
{

	RegularMesh<DIM> &mesh_v = irregular_mesh_v->regularMesh();
	RegularMesh<DIM> &mesh_p = irregular_mesh_p->regularMesh();
        /// 更新P网格的三个顶点.
	for (int i = 0; i < mesh_p.n_geometry(0); ++i)
	{
		(*mesh_p.h_geometry<0>(i))[0] = point(i)[0];
		(*mesh_p.h_geometry<0>(i))[1] = point(i)[1];
	}
	/// 更新p网格的中点. 不能替换单元对应.
	for (int j = 0; j < mesh_p.n_geometry(1); ++j)
	{
		GeometryBM &bnd = mesh_p.geometry(1, j);
		(*mesh_p.h_geometry<1>(bnd.index())->child[1]->vertex[0])[0]
			= 0.5 * ((*mesh_p.h_geometry<0>(bnd.vertex(0)))[0] +
				 (*mesh_p.h_geometry<0>(bnd.vertex(1)))[0]);
		(*mesh_p.h_geometry<1>(bnd.index())->child[1]->vertex[0])[1]
			= 0.5 * ((*mesh_p.h_geometry<0>(bnd.vertex(0)))[1] +
				 (*mesh_p.h_geometry<0>(bnd.vertex(1)))[1]);
	}

        /// 非正则网格做半正则化, 不做的话，irregular_mesh下的regularMesh上点的坐标不会变。
	irregular_mesh_p->semiregularize();
	/// 正则化, 但不用重新编号.
	irregular_mesh_p->regularize(false);
	irregular_mesh_v->semiregularize();
	irregular_mesh_v->regularize(false);
	/// 顺便把有限元空间的插值点更新一下。
	fem_space_v.updateDofInterpPoint();
	fem_space_p.updateDofInterpPoint();
};
void ISOP2P1::getMonitor()
{
	/// 每一次计算monitor都需要将网格同步一下.
	syncMesh();
	/// 网格点变化，插值点发生变化，意味着基函数也发生了变化，因此所有的矩阵需要重新拼装。
	/// 虽然网格变化了，但是网格逻辑关系没有发生变化，因此不需要重新生成矩阵结构。
//	buildMatrixStruct();
	buildMatrix();
	/// 将速度有限元空间的数值解插值到压力有限元空间中，
	/// 为了方便计算压力网格单元上的monitor的值.
	FEMFunction<double, DIM> _u_h(fem_space_p);
	FEMFunction<double, DIM> _v_h(fem_space_p);
	Operator::L2Interpolate(v_h[0], _u_h);
	Operator::L2Interpolate(v_h[1], _v_h);

       /// 注释的这部分代码是选用文章中的后验误差(4.3) - (4.4)做monitor.       

//	 RegularMesh<DIM> &mesh_p = irregular_mesh_p->regularMesh();
//
//	 int n_face = mesh_p.n_geometry(1);
//	 std::vector<double> side_length(n_face);
//         std::vector<bool> flag(n_face, false);
//         std::vector<double> jump_ux(n_face);
//         std::vector<double> jump_uy(n_face);
//
//         FEMSpace<double,DIM>::ElementIterator the_element_p = fem_space_p.beginElement();
//         FEMSpace<double,DIM>::ElementIterator end_element_p = fem_space_p.endElement();
//         /// 遍历所有速度单元,计算每条边上的 jump / ||v||_L2.
//         for (; the_element_p != end_element_p; ++the_element_p)
//         {
//         	/// 几何信息.
//         	double volume = the_element_p->templateElement().volume();
//         	const QuadratureInfo<2>& quad_info = the_element_p->findQuadratureInfo(4);
//         	std::vector<double> jacobian = the_element_p->local_to_global_jacobian(quad_info.quadraturePoint());
//         	int n_quadrature_point = quad_info.n_quadraturePoint();
//         	std::vector<Point<2> > q_point = the_element_p->local_to_global(quad_info.quadraturePoint());
//
//         	const GeometryBM& geometry = the_element_p->geometry();
//
//         	for (int i = 0; i < geometry.n_boundary(); ++i)
//         	{
//         		int j = geometry.boundary(i);
//         		const GeometryBM& side = mesh_p.geometry(1, j);
//         		const Point<DIM>& p0 = mesh_p.point(mesh_p.geometry(0, side.boundary(0)).vertex(0));
//         		const Point<DIM>& p1 = mesh_p.point(mesh_p.geometry(0, side.boundary(1)).vertex(0));
//         		std::vector<double> vx_gradient = _u_h.gradient(midpoint(p0, p1), *the_element_p);
//         		std::vector<double> vy_gradient = _v_h.gradient(midpoint(p0, p1), *the_element_p);
//         		double vx_value = _u_h.value(midpoint(p0, p1), *the_element_p);
//         		double vy_value = _v_h.value(midpoint(p0, p1), *the_element_p);
//
//         		double v_L2norm = sqrt(vx_value * vx_value + vy_value * vy_value + eps);
//         		side_length[j] = distance(p0, p1);
//         		if (flag[j])
//         		{
//         			jump_ux[j] -= (vx_gradient[0] * (p0[1] - p1[1]) + vx_gradient[1] * (p1[0] - p0[0]));
//         			jump_uy[j] -= (vy_gradient[0] * (p0[1] - p1[1]) + vy_gradient[1] * (p1[0] - p0[0]));
//         			flag[j] = false;
//         		}
//         		else
//         		{
//         			jump_ux[j] = (vx_gradient[0] * (p0[1] - p1[1]) + vx_gradient[1] * (p1[0] - p0[0]));
//         			jump_uy[j] = (vy_gradient[0] * (p0[1] - p1[1]) + vy_gradient[1] * (p1[0] - p0[0]));
//         			flag[j] = true;
//         		}
//         	}
//         }
//         the_element_p = fem_space_p.beginElement();
//         for(int i = 0; the_element_p != end_element_p; ++the_element_p, ++i)
//         {
//         	const GeometryBM &geometry = the_element_p->geometry();
//         	double cont = 0.0, total_ele_side_length = 0.0;
//         	for (int l = 0; l < geometry.n_boundary(); ++l)
//         	{
//         		int j = geometry.boundary(l);
//         		/// 如果编号是j的边是边界.
//         		if (flag[j])
//         			continue;
//         		/// 将DIM个梯度沿法向方向上的跳跃平方值累加.
//         		cont += (jump_ux[j] * jump_ux[j] + jump_uy[j] * jump_uy[j]) * side_length[j];
//         		total_ele_side_length += side_length[j];
//         	}
//         	monitor(i) = sqrt(cont / total_ele_side_length);
//         }
	
	FEMSpace<double, DIM>::ElementIterator the_element = fem_space_p.beginElement();
	FEMSpace<double, DIM>::ElementIterator end_element = fem_space_p.endElement();
	for(; the_element != end_element; ++the_element)
	{
		double volume = the_element->templateElement().volume();
		const QuadratureInfo<DIM> &quad_info = the_element->findQuadratureInfo(3);
		std::vector<double> jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
		std::vector<std::vector<double> >  basis_value = the_element->basis_function_value(q_point);
		std::vector<double> u_h_value = _u_h.value(q_point, *the_element);
		std::vector<double> v_h_value = _v_h.value(q_point, *the_element);
		std::vector<std::vector<double> > u_h_gradient = _u_h.gradient(q_point, *the_element);
		std::vector<std::vector<double> > v_h_gradient = _v_h.gradient(q_point, *the_element);
		float d = 0.0, area = 0.0, cont = 0.0;
		for(int l = 0; l < n_quadrature_point; ++l)
		{
			double Jxw = quad_info.weight(l) * jacobian[l] * volume;
			area += Jxw;
			d += Jxw * sqrt(innerProduct(u_h_gradient[l], u_h_gradient[l]) + innerProduct(v_h_gradient[l], v_h_gradient[l]));
		}
		monitor(the_element->index()) = d / area;
	}
	std::cout << "max monitor=" << *std::max_element(monitor().begin(), monitor().end())
		  << "\t min monitor=" << *std::min_element(monitor().begin(), monitor().end())
		  << std::endl;
	/// 对monitor的光滑。
	smoothMonitor(4);
	double max_monitor = *std::max_element(monitor().begin(), monitor().end());
	for (int i = 0; i < n_geometry(DIM); ++i)
		         monitor(i) = 1.0 / sqrt(1.0 + alpha * pow(monitor(i), beta));
}

void ISOP2P1::updateSolution()
{
	/// 网格移动步长.
	const double &msl = moveStepLength();

	/// 虚拟时间步长，取值为1，主要是由于网格移动的速度是O(h)的量，步长可以取的大些。
	double _tau = 1.0;
	double tau_step = _tau / 1;

	int n_dof_v = fem_space_v.n_dof();
	int n_dof_p = fem_space_p.n_dof();
	int n_total_dof_v = DIM * n_dof_v;
	int n_total_dof = DIM * n_dof_v + n_dof_p;
	/// 备份数值解。
	FEMFunction<double, DIM> _u_h(v_h[0]);
	FEMFunction<double, DIM> _v_h(v_h[1]);
	/// 开始更新数值解，论文上是用三阶显示Runge-Kutta方法，我们这里为了调试的便利，
	/// 仅仅使用一步显示Euler方法.
	for (int m = 1; m > 0; --m)
	{
		/// 首先声名系数矩阵，这个系数矩阵跟Stokes方程离散完的系数矩阵结构基本相似，
		/// 只是在两个对角块上为两个质量阵，所以稀疏矩阵结构就采用sp_stokes。
		matrix.reinit(sp_stokes);
		/// 我们需要用到的分块矩阵是mat_v_mass, Bx^T,By^T, Bx, By,我们将上述矩阵
		/// 加到对应的系数矩阵matrix上去.
		/// (0, 0)块
		for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
			matrix.global_entry(index_vxvx[i]) = mat_v_mass.global_entry(i);
		/// (1, 1)块
		for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
			matrix.global_entry(index_vyvy[i]) = mat_v_mass.global_entry(i);
		/// (0, 2)块
		for (int i = 0; i < sp_pvx.n_nonzero_elements(); ++i)
			matrix.global_entry(index_pvx[i]) = tau_step / m * mat_pvx_divT.global_entry(i);
		/// (1, 2)块
		for (int i = 0; i < sp_pvy.n_nonzero_elements(); ++i)
			matrix.global_entry(index_pvy[i]) = tau_step / m * mat_pvy_divT.global_entry(i);
		/// (2, 0)块.
		for (int i = 0; i < sp_vxp.n_nonzero_elements(); ++i)
			matrix.global_entry(index_vxp[i]) = tau_step / m * mat_vxp_div.global_entry(i);
		/// (2, 1)块
		for (int i = 0; i < sp_vyp.n_nonzero_elements(); ++i)
			matrix.global_entry(index_vyp[i]) = tau_step / m * mat_vyp_div.global_entry(i);
		/// 右端项在这里声名.
		rhs.reinit(n_total_dof);
		/// 开始拼装右端项.
		FEMSpace<double, DIM>::ElementIterator the_element = fem_space_v.beginElement();
		FEMSpace<double, DIM>::ElementIterator end_element = fem_space_v.endElement();
		for(; the_element != end_element; ++the_element)
		{
			/// 参考单元单元面积.
			double volume = the_element->templateElement().volume();
			/// 根据积分精度寻找积分信息，包括积分点。
			const QuadratureInfo<DIM> &quad_info = the_element->findQuadratureInfo(3);
			/// 参考单元到实际单元的jacobian矩阵变换。
			std::vector<double>  jacobian = the_element->local_to_global_jacobian(quad_info.quadraturePoint());
			/// 当前单元上的积分点的个数.
			int n_quadrature_point = quad_info.n_quadraturePoint();
			/// 实际单元上的积分点.
			std::vector<Point<DIM> > q_point = the_element->local_to_global(quad_info.quadraturePoint());
			/// 速度单元信息.
			std::vector<std::vector<double> > basis_value = the_element->basis_function_value(q_point);
			std::vector<std::vector<std::vector<double> > > basis_gradient = the_element->basis_function_gradient(q_point);
			std::vector<double> u_h_value = _u_h.value(q_point, *the_element);
			std::vector<double> v_h_value = _v_h.value(q_point, *the_element);
			/// 速度的梯度.
			std::vector<std::vector<double> > u_h_gradient = v_h[0].gradient(q_point, *the_element);
			std::vector<std::vector<double> > v_h_gradient = v_h[1].gradient(q_point, *the_element);
			/// 当前单元上的自由度编号.
			const std::vector<int> &ele_dof_v = the_element->dof();
			/// 当前单元自由度的个数.
			int n_element_dof_v = the_element->n_dof();

			/// 对应的压力单元.
			Element<double, DIM> &p_element = fem_space_p.element(index_ele_v2p[the_element->index()]);
			const std::vector<int>& element_dof_p = p_element.dof();
			std::vector<std::vector<std::vector<double> > > basis_gradient_p = p_element.basis_function_gradient(q_point);
			std::vector<std::vector<double> >  basis_value_p = p_element.basis_function_value(q_point);
			int n_element_dof_p = p_element.n_dof();

			/// 速度积分点的移动方向, 因为我们移动的是P网格,把速度网格上的积分点看作是大的P网格上的一个点.
			/// 先要根据但前速度单元编号找到对应的压力单元编号,再获取压力单元上的移动方向。
			/// 已经仔细测试过这个move_vector.
			std::vector<std::vector<double> > move_vector = moveDirection(q_point, p_element.index());
			for (int l = 0; l < n_quadrature_point; ++l)
			{
				double Jxw = quad_info.weight(l) * jacobian[l] * volume;
				for (int i = 0; i < n_element_dof_v; ++i)
				{
					// for(int j = 0; j < n_element_dof_v; ++j)
					// {
					// 	double cont = -Jxw * tau_step / m * msl * basis_value[i][l] * innerProduct(basis_gradient[j][l], move_vector[l]);
					// 	matrix.add(ele_dof_v[i], ele_dof_v[j], cont);
					// 	matrix.add(ele_dof_v[i] + n_dof_v, ele_dof_v[j] + n_dof_v, cont);
					// }
					/// x方向.
					double rhs_cont = basis_value[i][l] * (u_h_value[l] + tau_step / m * msl * innerProduct(u_h_gradient[l], move_vector[l]));
					rhs_cont *= Jxw;
					rhs(ele_dof_v[i]) += rhs_cont;
					/// y方向.
					rhs_cont = basis_value[i][l] * (v_h_value[l] + tau_step / m * msl * innerProduct(v_h_gradient[l], move_vector[l]));
					rhs_cont *= Jxw;
					rhs(ele_dof_v[i] + n_dof_v) += rhs_cont;
				}
			}
		}
		/// 整体数值解.
		Vector<double> x(n_total_dof);
		/// 边界条件处理, 跟Stokes方程求解的边界处理相同.
		if (fabs(t) < 1e-12)
			boundaryValueStokes(x, t);
		else
			boundaryValueStokes(x, t + dt);
		std::cout << "boundary condition for updateSolution OK!" << std::endl;

		/// 拼装完矩阵和右端项，做完边界条件处理，开始矩阵求解。
		dealii::SolverControl solver_control (n_dof_v, l_Euler_tol * rhs.l2_norm(), 0);

		SolverGMRES<Vector<double> >::AdditionalData para(2000, false, true);
		/// 不用para算不动，残差不下降.
		SolverGMRES<Vector<double> > gmres (solver_control, para);

		/// 不完全LU分解.
		dealii::SparseILU <double> preconditioner;
		preconditioner.initialize(matrix);

		/// 不用预处理，求解速度非常慢，不能接受.
		gmres.solve(matrix, x, rhs, preconditioner);

		// /// 矩阵求解.
		// SparseMatrix<double> mat_BTx(sp_pvx);
		// SparseMatrix<double> mat_BTy(sp_pvy);
		// SparseMatrix<double> mat_Bx(sp_vxp);
		// SparseMatrix<double> mat_By(sp_vyp);
		// SparseMatrix<double> mat_Ax(sp_vxvx);
		// SparseMatrix<double> mat_Ay(sp_vyvy);

		// for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
		// 	mat_Ax.global_entry(i) = matrix.global_entry(index_vxvx[i]);
		// for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
		// 	mat_Ay.global_entry(i) = matrix.global_entry(index_vyvy[i]);
		// for (int i = 0; i < sp_pvx.n_nonzero_elements(); ++i)
		// 	mat_BTx.global_entry(i) = matrix.global_entry(index_pvx[i]);
		// for (int i = 0; i < sp_pvy.n_nonzero_elements(); ++i)
		// 	mat_BTy.global_entry(i) = matrix.global_entry(index_pvy[i]);
		// for (int i = 0; i < sp_vxp.n_nonzero_elements(); ++i)
		// 	mat_Bx.global_entry(i) = matrix.global_entry(index_vxp[i]);
		// for (int i = 0; i < sp_vyp.n_nonzero_elements(); ++i)
		// 	mat_By.global_entry(i) = matrix.global_entry(index_vyp[i]);
	
		// /// alp对AMGSolver的初始化影响比较大, 如果取得很小，初始化很快.
		// double alp = dt * viscosity;
		// AMGSolver solverQ(mat_Ax, 1.0e-12, 3, 100, 0.382, alp);
                // //        AMGSolver solverQ(mat_Ax);
		// InverseMatrix AInv(mat_Ax, solverQ);
		// /// 这里没有对速度质量阵进行边界条件处理.
		// InverseMatrix QInv(mat_v_mass, solverQ);
		// SchurComplement schur_complement(mat_BTx, mat_BTy, mat_Bx, mat_By, mat_v_mass, QInv, QInv);
		// /// 压力块的AMG 预处理用压力刚度矩阵效果比质量阵好.
		// AMGSolver solverP(mat_p_stiff);
		// ApproxSchurComplement asc(mat_p_stiff, solverP);
		
		// /// updateSolution 预处理.
		// updateSolutionPreconditioner update_solution_preconditioner(mat_BTx, mat_BTy, AInv, AInv, schur_complement, asc);
		// dealii::SolverControl solver_control (n_dof_p, l_Euler_tol * rhs.l2_norm(), 1);
		// SolverGMRES<Vector<double> >::AdditionalData para(100, false, true);
		// SolverGMRES<Vector<double> > gmres (solver_control, para);

		// clock_t t_cost = clock();
		// gmres.solve(matrix, x, rhs, update_solution_preconditioner);
		// t_cost = clock() - t_cost;
		// std::cout << "time cost: " << (((float)t_cost) / CLOCKS_PER_SEC) << std::endl;

		/// 将数值解分成v_h[0], v_h[1], 和 p_h.
		for (int i = 0; i < n_dof_v; ++i)
		{
			v_h[0](i) = x(i);
			v_h[1](i) = x(i + n_dof_v);
		}
		for (int i = 0; i < n_dof_p; ++i)
			p_h(i) = x(i + 2 * n_dof_v);
		Vector<double> res(n_total_dof);
		matrix.vmult(res, x);
		res *= -1;
		res += rhs;
		/// 输出残量.
		std::ofstream output;
		output.open("residual", std::ofstream::out);
		output.setf(std::ios::fixed);
		output.precision(20);
	
		output << "residual = [" << std::endl;
		for (int i = 0; i < n_total_dof; ++i)
		{
			output << res(i) << std::endl;
		}
		output << "];" << std::endl;

		std::cout << "res_l2norm =" << res.l2_norm() << std::endl;
	}
};

void ISOP2P1::outputSolution()
{
};
void ISOP2P1::movingMesh()
{
	/// 网格移动部分。
	moveMesh();
	/// 移动完同步网格。
	syncMesh();
	/// 顺便把有限元空间的插值点更新一下。
	fem_space_v.updateDofInterpPoint();
	fem_space_p.updateDofInterpPoint();
};

