#include "ISOP2P1.h"
#include "preconditioner.h"
#include "functions.h"
#define DIM 2

void ISOP2P1::solveStokes()
{
	buildStokesSys();
	std::cout << "Stokes system builded." << std::endl;
	int n_dof_v = fem_space_v.n_dof();
	int n_dof_p = fem_space_p.n_dof();
	int n_total_dof = 2 * n_dof_v + n_dof_p;

	/// 构建系数矩阵和右端项.
	/// 这个存放整体的数值解. 没有分割成 u_h[0], u_h[1] 和 p_h.
	Vector<double> x(n_total_dof);
	/// 将数值解合并一个向量便于边界处理.
	for (int i = 0; i < n_dof_v; ++i)
	{
		x(i) = v_h[0](i);
		x(n_dof_v + i) = v_h[1](i);
	}

	for (int i = 0; i < n_dof_p; ++i)
		x(2 * n_dof_v + i) = p_h(i);

	rhs.reinit(n_total_dof);

	/// 边界条件一起处理了. 修改了matrix, rhs和x.
	boundaryValueStokes(x, t);
	std::cout << "Stokes boundary applied." << std::endl;

// 	/// 矩阵求解.
// 	dealii::SolverControl solver_control (400000, l_tol * rhs.l2_norm(), 1);
// 	/// 不完全LU分解.
// 	dealii::SparseILU<double> preconditioner;
// 	preconditioner.initialize(matrix);

// 	/// 求解Stokes方程, MinRes要比GMRES求解器要快很多,
// 	/// 当矩阵规模稍微大点的时候,GMRES会出现不收敛的情况.
// 	SolverMinRes<Vector<double> > minres (solver_control);
// 	SolverGMRES<Vector<double> >::AdditionalData para(1000, false, true);
// 	/// 不用para算不动, 但是均匀网格下可以不用para,算其它的例子也不用para,是不是跟计算区域有关系？
// 	SolverGMRES<Vector<double> > gmres(solver_control, para);

//         // /// 移动网格和时间发展中，这个预处理失效.
// 	//  StokesPreconditioner preconditioner;
// 	 // /// 预处理矩阵.
// 	 // SparseMatrix<double> matrix_vxvx(sp_vxvx);
// 	 // SparseMatrix<double> matrix_vyvy(sp_vyvy);
// 	 // /// 这里从 Stokes 取是因为加了边界条件.
// 	 // for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
// 	 // 	matrix_vxvx.global_entry(i) = matrix.global_entry(index_vxvx[i]);
// 	 // for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
// 	 // 	matrix_vyvy.global_entry(i) = matrix.global_entry(index_vyvy[i]);

// 	 // preconditioner.initialize(mat_v_stiff, mat_v_stiff, mat_p_mass);


// 	clock_t t_cost = clock();
//         minres.solve (matrix, x, rhs, PreconditionIdentity());
// //	gmres.solve(matrix, x, rhs, preconditioner);
// 	t_cost = clock() - t_cost;

// 	std::cout << "time cost: " << (((float)t_cost) / CLOCKS_PER_SEC) << std::endl;

// 	/// 将整体数值解分割成速度和压力.
// 	for (int i = 0; i < n_dof_v; ++i)
// 	{
// 		v_h[0](i) = x(i);
// 		v_h[1](i) = x(i + n_dof_v);
// 	}
// 	for (int i = 0; i < n_dof_p; ++i)
// 		p_h(i) =  x(i + 2 * n_dof_v);
//        /// 计算误差, t为时间参数.
// 	computeError(t);
// 	/// debug, 计算惨量的L2 norm。
// 	Vector<double> res(n_total_dof);
// 	matrix.vmult(res, x);
// 	res *= -1;
// 	res += rhs;
// 	std::cout << "res_l2norm =" << res.l2_norm() << std::endl;

	/// 矩阵求解.
	SparseMatrix<double> mat_BTx(sp_pvx);
	SparseMatrix<double> mat_BTy(sp_pvy);
	SparseMatrix<double> mat_Bx(sp_vxp);
	SparseMatrix<double> mat_By(sp_vyp);
	SparseMatrix<double> mat_Ax(sp_vxvx);
	SparseMatrix<double> mat_Ay(sp_vyvy);

	for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
		mat_Ax.global_entry(i) = matrix.global_entry(index_vxvx[i]);
	for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
		mat_Ay.global_entry(i) = matrix.global_entry(index_vyvy[i]);
	for (int i = 0; i < sp_pvx.n_nonzero_elements(); ++i)
		mat_BTx.global_entry(i) = matrix.global_entry(index_pvx[i]);
	for (int i = 0; i < sp_pvy.n_nonzero_elements(); ++i)
		mat_BTy.global_entry(i) = matrix.global_entry(index_pvy[i]);
	for (int i = 0; i < sp_vxp.n_nonzero_elements(); ++i)
		mat_Bx.global_entry(i) = matrix.global_entry(index_vxp[i]);
	for (int i = 0; i < sp_vyp.n_nonzero_elements(); ++i)
		mat_By.global_entry(i) = matrix.global_entry(index_vyp[i]);

	/// alp对AMGSolver的初始化影响比较大, 如果取得很小，初始化很快.
        double alp = dt * viscosity;
	AMGSolver solverQ(mat_Ax, 1.0e-12, 3, 100, 0.382, alp);
//	AMGSolver solverQ(mat_Ax);
	InverseMatrix AInv(mat_Ax, solverQ);
	/// 这里没有对速度质量阵进行边界条件处理.
	InverseMatrix QInv(mat_v_mass, solverQ);
	SchurComplement schur_complement(mat_BTx, mat_BTy, mat_Bx, mat_By, mat_v_mass, QInv, QInv);
	AMGSolver solverP(mat_p_mass);
	ApproxSchurComplement asc(mat_p_mass, solverQ);
	
	LSCPreconditioner lsc_preconditioner(mat_BTx, mat_BTy, mat_Bx, mat_By, mat_Ax, mat_Ax, mat_v_mass, schur_complement, asc, QInv,
					     AInv, AInv);
	/// 矩阵求解.
	dealii::SolverControl solver_control (400000, l_Euler_tol * rhs.l2_norm(), 0);
	SolverMinRes<Vector<double> > minres(solver_control);
	
	minres.solve(matrix, x, rhs, lsc_preconditioner);
	/// 将整体数值解分割成速度和压力.
	for (int i = 0; i < n_dof_v; ++i)
	{
		v_h[0](i) = x(i);
		v_h[1](i) = x(i + n_dof_v);
	}
	for (int i = 0; i < n_dof_p; ++i)
		p_h(i) =  x(i + 2 * n_dof_v);
	
	/// 计算误差, t为时间参数.
	computeError(t);


// 	/// 矩阵求解.
// 	SparseMatrix<double> mat_BTx(sp_pvx);
// 	SparseMatrix<double> mat_BTy(sp_pvy);
// 	SparseMatrix<double> mat_Bx(sp_vxp);
// 	SparseMatrix<double> mat_By(sp_vyp);
// 	SparseMatrix<double> mat_Ax(sp_vxvx);
// 	SparseMatrix<double> mat_Ay(sp_vyvy);

// 	for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
// 		mat_Ax.global_entry(i) = matrix.global_entry(index_vxvx[i]);
// 	for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
// 		mat_Ay.global_entry(i) = matrix.global_entry(index_vyvy[i]);
// 	for (int i = 0; i < sp_pvx.n_nonzero_elements(); ++i)
// 		mat_BTx.global_entry(i) = matrix.global_entry(index_pvx[i]);
// 	for (int i = 0; i < sp_pvy.n_nonzero_elements(); ++i)
// 		mat_BTy.global_entry(i) = matrix.global_entry(index_pvy[i]);
// 	for (int i = 0; i < sp_vxp.n_nonzero_elements(); ++i)
// 		mat_Bx.global_entry(i) = matrix.global_entry(index_vxp[i]);
// 	for (int i = 0; i < sp_vyp.n_nonzero_elements(); ++i)
// 		mat_By.global_entry(i) = matrix.global_entry(index_vyp[i]);
	
// 	Vector<double> tmp1(n_dof_v);
// 	Vector<double> tmp2(n_dof_v);
// 	Vector<double> rhs_vx(n_dof_v);
// 	Vector<double> rhs_vy(n_dof_v);
// 	Vector<double> rhs_p(n_dof_p);

// 	for (int i = 0; i < n_dof_v; ++i)
// 	{
// 		rhs_vx(i) = rhs(i);
// 		v_h[0](i) = x(i);
// 		rhs_vy(i) = rhs(n_dof_v + i);
// 		v_h[1](i) = x(n_dof_v + i);
// 	}
// 	for (int i = 0; i < n_dof_p; ++i)
// 	{
// 		rhs_p(i) = rhs(2 * n_dof_v + i);
// 		p_h(i) = x(2 * n_dof_v + i);
// 	}

// 	Vector<double> schur_rhs (n_dof_p);
// 	AMGSolver solverQ(mat_Ax);
// 	InverseMatrix M(mat_Ax, solverQ);
// 	M.vmult (tmp1, rhs_vx);
// 	M.vmult (tmp2, rhs_vy);
// 	mat_Bx.vmult(schur_rhs, tmp1);
// 	mat_By.vmult_add(schur_rhs, tmp2);
// 	schur_rhs -= rhs_p;

// 	SchurComplement schur_complement(mat_BTx, mat_BTy, mat_Bx, mat_By, mat_v_mass, M, M, dt);

// 	SolverControl solver_control_cg (n_dof_p * 2,
// 					 1e-12*schur_rhs.l2_norm());
// 	SolverCG<>    cg (solver_control_cg);

// 	AMGSolver AQ(mat_p_mass);
// 	ApproxSchurComplement asc(mat_p_mass, AQ);
// 	cg.solve (schur_complement, p_h, schur_rhs, asc);

// //	cg.solve (schur_complement, p_h, schur_rhs, PreconditionIdentity());
// 	std::cout << solver_control_cg.last_step()
//             		  << " CG Schur complement iterations to obtain convergence."
//             		  << std::endl;

// 	mat_BTx.vmult(tmp1, *dynamic_cast<const Vector<double>* >(&p_h));
// 	mat_BTy.vmult(tmp2, *dynamic_cast<const Vector<double>* >(&p_h));
// 	tmp1 *= -1;
// 	tmp2 *= -1;
// 	tmp1 += rhs_vx;
// 	tmp2 += rhs_vy;

// 	M.vmult(v_h[0], tmp1);
// 	M.vmult(v_h[1], tmp2);
// 	std::cout << "Stokes system solved." << std::endl;
// 	/// 计算误差, t为时间参数.
// 	computeError(t);
// 	/// debug, 计算惨量的L2 norm。
// 	Vector<double> res(n_total_dof);
// 	matrix.vmult(res, x);
// 	res *= -1;
// 	res += rhs;
// 	std::cout << "res_l2norm =" << res.l2_norm() << std::endl;
};

void ISOP2P1::solveNS(int method)
{
	int n_dof_v = fem_space_v.n_dof();
	int n_dof_p = fem_space_p.n_dof();
	int n_total_dof = n_dof_v * 2 + n_dof_p;

	/// 开始迭代.
	double error_N = 1.0;
	int iteration_times = 0;
	while (error_N > n_tol)
	{
		/// Newton 迭代或 Picard 迭代.
		/// 先更新和速度场有关的矩阵块.
		updateNonlinearMatrix();
		/// 构建迭代矩阵.
		if (method == 1)
			buildNewtonSys4NS();
		else if (method == 2)
			buildPicardSys4NS();
		else if (method == 3)
			if (iteration_times < 2)
				buildPicardSys4NS();
			else
				buildNewtonSys4NS();
		else
		{
			std::cout << "Newton: 1, Picard: 2, Hybrid: 3." << std::endl;
			exit(1);
		}

		/// 建立右端项.
		rhs.reinit(n_total_dof);

		FEMSpace<double, DIM>::ElementIterator the_element_v = fem_space_v.beginElement();
		FEMSpace<double, DIM>::ElementIterator end_element_v = fem_space_v.endElement();
		FEMSpace<double, DIM>::ElementIterator the_element_p = fem_space_p.beginElement();
		FEMSpace<double, DIM>::ElementIterator end_element_p = fem_space_p.endElement();
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
			std::vector<Point<DIM> > q_point
				= the_element_v->local_to_global(quad_info.quadraturePoint());
			/// 速度单元信息.
			std::vector<std::vector<std::vector<double> > > basis_gradient_v
				= the_element_v->basis_function_gradient(q_point);
			std::vector<std::vector<double> >  basis_value_v
				= the_element_v->basis_function_value(q_point);
			const std::vector<int>& element_dof_v = the_element_v->dof();
			std::vector<double> fx_value = source_v[0].value(q_point, *the_element_v);
			std::vector<double> fy_value = source_v[1].value(q_point, *the_element_v);
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
					double rhs_cont = fx_value[l] * basis_value_v[i][l];
					rhs_cont -= (vx_value[l] * vx_gradient[l][0] +
						     vy_value[l] * vx_gradient[l][1]) * basis_value_v[i][l];
					rhs_cont -= viscosity * innerProduct(basis_gradient_v[i][l], vx_gradient[l]);
					rhs_cont += p_value[l] * basis_gradient_v[i][l][0];
					rhs_cont *= Jxw;
					rhs(element_dof_v[i]) += rhs_cont;

					rhs_cont = fy_value[l] * basis_value_v[i][l];
					rhs_cont -= (vx_value[l] * vy_gradient[l][0] +
						     vy_value[l] * vy_gradient[l][1]) * basis_value_v[i][l];
					rhs_cont -= viscosity * innerProduct(basis_gradient_v[i][l], vy_gradient[l]);
					rhs_cont += p_value[l] * basis_gradient_v[i][l][1];
					rhs_cont *= Jxw;
					rhs(n_dof_v + element_dof_v[i]) += rhs_cont;
				}
			}
		}

		/// 遍历压力单元. 拼装矩阵和右端项.
		for (the_element_p = fem_space_p.beginElement();
		     the_element_p != end_element_p; ++the_element_p)
		{
			const std::vector<int>& element_dof_p = the_element_p->dof();
			int n_element_dof_p = the_element_p->n_dof();
			for (int i = 0; i < n_element_dof_p; ++i)
			{
				int idx_p = the_element_p->index();
				int n_chi = index_ele_p2v[idx_p].size();
				for (int k = 0; k < n_chi; k++)
				{
					/// 速度单元信息.
					Element<double, DIM> &v_element = fem_space_v.element(index_ele_p2v[idx_p][k]);
					/// 几何信息.
					double volume = v_element.templateElement().volume();
					const QuadratureInfo<DIM>& quad_info = v_element.findQuadratureInfo(4);
					std::vector<double> jacobian
						= v_element.local_to_global_jacobian(quad_info.quadraturePoint());
					int n_quadrature_point = quad_info.n_quadraturePoint();
					std::vector<Point<DIM> > q_point
						= v_element.local_to_global(quad_info.quadraturePoint());
					std::vector<std::vector<double> > vx_gradient
						= v_h[0].gradient(q_point, v_element);
					std::vector<std::vector<double> > vy_gradient
						= v_h[1].gradient(q_point, v_element);
					std::vector<double> vx_value = v_h[0].value(q_point, v_element);
					std::vector<double> vy_value = v_h[1].value(q_point, v_element);
					/// 压力单元信息.
					std::vector<std::vector<double> >  basis_value_p
						= the_element_p->basis_function_value(q_point);
					/// 具体拼装.
					for (int l = 0; l < n_quadrature_point; ++l)
					{
						double Jxw = quad_info.weight(l) * jacobian[l] * volume;

						/// 右端项还是零. 源项和 Neumann 条件.
						double rhs_cont = Jxw * basis_value_p[i][l]
							* (vx_gradient[l][0] + vy_gradient[l][1]);
						rhs(2 * n_dof_v + element_dof_p[i]) += rhs_cont;
					}
				}
			}

		}

		/// 初始化未知量.
		Vector<double> x(n_total_dof);

		/// 边界条件处理.
		boundaryValueNS(x);

		std::cout << "nonlinear res:" << std::endl;
		double revx = 0.0;
		for (int i = 0; i < n_dof_v; ++i)
			revx += rhs(i) * rhs(i);
		std::cout << "vx re: " << sqrt(revx) << std::endl;
		double revy = 0.0;
		for (int i = 0; i < n_dof_v; ++i)
			revy += rhs(i + n_dof_v) * rhs(i + n_dof_v);
		std::cout << "vy re: " << sqrt(revy) << std::endl;
		double rep = 0.0;
		for (int i = 0; i < n_dof_p; ++i)
			rep += rhs(i + 2 * n_dof_v) * rhs(i + 2 * n_dof_v);
		std::cout << "p re: " << sqrt(rep) << std::endl;

		double re = revx + revy +rep;
		std::cout << "total re: " << sqrt(re) << std::endl;
		std::cout << "pause ..." << std::endl;
		getchar();
		if (sqrt(re) < n_tol)
		{
			std::cout << "Covergence with residual: " << sqrt(re)
		    		  << " in step " << iteration_times << std::endl;
			break;
		}

		std::cout << "Building precondition matrix ..." << std::endl;

		/// 矩阵求解.
		SparseMatrix<double> mat_Axx(sp_vxvx);
		SparseMatrix<double> mat_Ayy(sp_vyvy);
		SparseMatrix<double> mat_Wxy(sp_vyvx);
		SparseMatrix<double> mat_Wyx(sp_vxvy);
		SparseMatrix<double> mat_BTx(sp_pvx);
		SparseMatrix<double> mat_BTy(sp_pvy);

		for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
			mat_Axx.global_entry(i) = matrix.global_entry(index_vxvx[i]);
		for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
			mat_Ayy.global_entry(i) = matrix.global_entry(index_vyvy[i]);
		for (int i = 0; i < sp_vyvx.n_nonzero_elements(); ++i)
			mat_Wxy.global_entry(i) = matrix.global_entry(index_vyvx[i]);
		for (int i = 0; i < sp_vxvy.n_nonzero_elements(); ++i)
			mat_Wyx.global_entry(i) = matrix.global_entry(index_vxvy[i]);
		for (int i = 0; i < sp_pvx.n_nonzero_elements(); ++i)
			mat_BTx.global_entry(i) = matrix.global_entry(index_pvx[i]);
		for (int i = 0; i < sp_pvy.n_nonzero_elements(); ++i)
			mat_BTy.global_entry(i) = matrix.global_entry(index_pvy[i]);
		std::cout << "Precondition matrix builded!" << std::endl;

		/// 矩阵求解.
		dealii::SolverControl solver_control (4000000, l_tol);

		SolverGMRES<Vector<double> >::AdditionalData para(2000, false, true);
		SolverGMRES<Vector<double> > gmres (solver_control, para);
		std::cout << "Begin to solve linear system ..." << std::endl;
		//	gmres.solve (matrix, x, rhs, navierstokes_preconditioner);
		gmres.solve (matrix, x, rhs, PreconditionIdentity());

		/// 调试块: 直接观测真实残量.
		//	Vector<double> tmp(n_total_dof);
		//	matrix.vmult(tmp, x);
		//	tmp -= rhs;
		//	std::cout << "linear residual: " << tmp.l2_norm() << std::endl;
		//	getchar();

		FEMFunction<double, DIM> res_vx(fem_space_v);
		FEMFunction<double, DIM> res_vy(fem_space_v);
		FEMFunction<double, DIM> res_p(fem_space_p);

		/// 更新数值解.
		for (int i = 0; i < n_dof_v; ++i)
		{
			v_h[0](i) += x(i);
			res_vx(i) = x(i);
			v_h[1](i) += x(i + n_dof_v);
			res_vy(i) = x(i+ n_dof_v);
		}
		for (int i = 0; i < n_dof_p; ++i)
		{
			p_h(i) += x(i + 2 * n_dof_v);
			res_p(i) = x(i + 2 * n_dof_v);
		}

		double r_vx = Functional::L2Norm(res_vx, 1);
		double r_vy = Functional::L2Norm(res_vy, 1);
		double r_p = Functional::L2Norm(res_p, 1);
		/// 这个其实是更新...
		error_N = r_vx + r_vy + r_p;

		std::cout.setf(std::ios::fixed);
		std::cout.precision(20);
		std::cout << "updated vx: " << r_vx << std::endl;
		std::cout << "updated vy: " << r_vy << std::endl;
		std::cout << "updated p: " << r_p << std::endl;
		std::cout << "total updated: " << error_N << std::endl;
		std::cout << "step " << iteration_times << ", total updated: " << error_N
			  << ", GMRES stpes: " << solver_control.last_step() << std::endl;

		iteration_times++;
		if (iteration_times > 10)
		{
			std::cout << "Disconvergence at step 10." << std::endl;
			break;
		}
	}

};

void ISOP2P1::buildStokesSys()
{
	int n_dof_v = fem_space_v.n_dof();
	int n_dof_p = fem_space_p.n_dof();
	int n_total_dof = 2 * n_dof_v + n_dof_p;

	matrix.reinit(sp_stokes);

	/// (0, 0)
	for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vxvx[i]) = viscosity * mat_v_stiff.global_entry(i);
	/// (1, 1)
	for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vyvy[i]) = viscosity * mat_v_stiff.global_entry(i);
	/// (0, 2)
	for (int i = 0; i < sp_pvx.n_nonzero_elements(); ++i)
		matrix.global_entry(index_pvx[i]) = mat_pvx_divT.global_entry(i);

	/// (1, 2)
	for (int i = 0; i < sp_pvy.n_nonzero_elements(); ++i)
		matrix.global_entry(index_pvy[i]) = mat_pvy_divT.global_entry(i);

	/// (2, 0)
	for (int i = 0; i < sp_vxp.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vxp[i]) =  mat_vxp_div.global_entry(i);

	/// (2, 1)
	for (int i = 0; i < sp_vyp.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vyp[i]) =  mat_vyp_div.global_entry(i);

};


void ISOP2P1::buildNewtonSys4NS()
{
	int n_dof_v = fem_space_v.n_dof();
	int n_dof_p = fem_space_p.n_dof();
	int n_total_dof = 2 * n_dof_v + n_dof_p;

	/// (0, 0)
	for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vxvx[i]) = viscosity * mat_v_stiff.global_entry(i)
			+ mat_v_convection.global_entry(i)
			+ mat_v_Jacobi_xx.global_entry(i);
	/// (0, 1)
	for (int i = 0; i < sp_vyvx.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vyvx[i]) = mat_v_Jacobi_xy.global_entry(i);
	/// (1, 0)
	for (int i = 0; i < sp_vxvy.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vxvy[i]) = mat_v_Jacobi_yx.global_entry(i);
	/// (1, 1)
	for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vyvy[i]) = viscosity * mat_v_stiff.global_entry(i)
			+ mat_v_convection.global_entry(i)
			+ mat_v_Jacobi_yy.global_entry(i);
};

void ISOP2P1::buildPicardSys4NS()
{
	int n_dof_v = fem_space_v.n_dof();
	int n_dof_p = fem_space_p.n_dof();
	int n_total_dof = 2 * n_dof_v + n_dof_p;

	/// (0, 0)
	for (int i = 0; i < sp_vxvx.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vxvx[i]) = viscosity * mat_v_stiff.global_entry(i)
			+ mat_v_convection.global_entry(i);
	/// (1, 1)
	for (int i = 0; i < sp_vyvy.n_nonzero_elements(); ++i)
		matrix.global_entry(index_vyvy[i]) = viscosity * mat_v_stiff.global_entry(i)
			+ mat_v_convection.global_entry(i);
};

#undef DIM
