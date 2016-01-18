#include "ISOP2P1.h"
#include "functions.h"
#define DIM 2

void ISOP2P1::outputMat(const std::string &prefix, SparseMatrix<double> &mat)
{
	std::stringstream result;
	result.setf(std::ios::fixed);
	result.precision(4);
	result << prefix << "1.m";
	std::ofstream matrix(result.str().c_str());
	matrix.setf(std::ios::fixed);
	matrix.precision(20);

	matrix << prefix << " = [" << std::endl;
	int n = mat.n();
	int m = mat.m();
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
			matrix << mat.el(i, j) << "\t";
		matrix << std::endl;
	}
	matrix << "];" << std::endl;
};

void ISOP2P1::outputVec(const std::string &prefix, Vector<double> &vec)
{
	std::stringstream result;
	result.setf(std::ios::fixed);
	result.precision(4);
	result << prefix << "1.m";
	std::ofstream vector(result.str().c_str());
	vector.setf(std::ios::fixed);
	vector.precision(20);

	vector << prefix << " = [" << std::endl;
	int m = vec.size();
	for (int i = 0; i < m; ++i)
		vector << vec(i) << "\t";
	vector << std::endl;
	vector << "];" << std::endl;
};
void ISOP2P1::computeError(double _t)
{
	int n_dof_v = fem_space_v.n_dof();
	int n_dof_p = fem_space_p.n_dof();

	/// 往下的程序是计算误差.
	// RealVx divx;
	// RealVy divy;
	// RealP dip(-Functional::meanValue(p_h, 3));
	DiVx divx(viscosity, _t);
	DiVy divy(viscosity, _t);
	DiP dip(viscosity, _t, 0.0);
	/// 实际上是计算的||p - p_h - mean(p - p_h)||_L2.
	double mean_p = Functional::meanValue(dip, fem_space_p, 3);
	double mean_p_h = Functional::meanValue(p_h, 3);

	/// 压力L0误差.
	double p_hL0err = 0.0;
	/// 速度L2误差。
	double v_hL2err = 0.0;
	/// 速度H1误差.
	double v_hH1err = 0.0;

	/// 计算||p - p_h - mean(p - p_h)||_L2。
	FEMSpace<double, DIM>::ElementIterator the_element_p = fem_space_p.beginElement();
	FEMSpace<double, DIM>::ElementIterator end_element_p = fem_space_p.endElement();

	the_element_p = fem_space_p.beginElement();
	for (; the_element_p != end_element_p; ++the_element_p)
	{
		/// 当前单元信息.
		double volume = the_element_p->templateElement().volume();
		/// 积分精度, u 和 p 都是 1 次, 梯度和散度 u 都是常数.
		const QuadratureInfo<DIM>& quad_info = the_element_p->findQuadratureInfo(3);
		std::vector<double> jacobian = the_element_p->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element_p->local_to_global(quad_info.quadraturePoint());
		std::vector<double> p_h_value = p_h.value(q_point, *the_element_p);
		double cont = 0.0, area = 0.0;
		for (int l = 0; l < n_quadrature_point; ++l)
		{
			double Jxw = quad_info.weight(l) * jacobian[l] * volume;
			double df_value = dip.value(q_point[l]) - mean_p - p_h_value[l] + mean_p_h;
			p_hL0err += Jxw * df_value * df_value;
		}
	}
	p_hL0err = sqrt(fabs(p_hL0err));

	/// 计算速度单元上的||\vec{u} - \vec{u_h}||_L2 误差.
	RegularMesh<DIM> &mesh_v = irregular_mesh_v->regularMesh();
	err_ele.resize(mesh_v.n_geometry(DIM), 0.0);

	FEMSpace<double, DIM>::ElementIterator the_element_v = fem_space_v.beginElement();
	FEMSpace<double, DIM>::ElementIterator end_element_v = fem_space_v.endElement();
	for (; the_element_v != end_element_v; ++the_element_v)
	{
		/// 当前单元信息.
		double volume = the_element_v->templateElement().volume();
		/// 积分精度, u 和 p 都是 1 次, 梯度和散度 u 都是常数.
		const QuadratureInfo<DIM>& quad_info = the_element_v->findQuadratureInfo(3);
		std::vector<double> jacobian = the_element_v->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element_v->local_to_global(quad_info.quadraturePoint());
		std::vector<double> u_h_value = v_h[0].value(q_point, *the_element_v);
		std::vector<double> v_h_value = v_h[1].value(q_point, *the_element_v);
		std::vector<std::vector<double> > u_h_gradient = v_h[0].gradient(q_point, *the_element_v);
		std::vector<std::vector<double> > v_h_gradient = v_h[1].gradient(q_point, *the_element_v);

		std::vector<int> element_dof_v = the_element_v->dof();
		int n_ele_dof = the_element_v->n_dof();

		double cont = 0.0;
		for (int l = 0; l < n_quadrature_point; ++l)
		{
			double Jxw = quad_info.weight(l) * jacobian[l] * volume;
			double df_uh = divx.value(q_point[l]) - u_h_value[l];
			double df_vh = divy.value(q_point[l]) - v_h_value[l];
			std::vector<double> df_uh_gradient(2, 0.0);
			df_uh_gradient[0] = (divx.gradient(q_point[l]))[0] - u_h_gradient[l][0];
			df_uh_gradient[1] = (divx.gradient(q_point[l]))[1] - u_h_gradient[l][1];
			std::vector<double> df_vh_gradient(2, 0.0);
			df_vh_gradient[0] = (divy.gradient(q_point[l]))[0] - v_h_gradient[l][0];
			df_vh_gradient[1] = (divy.gradient(q_point[l]))[1] - v_h_gradient[l][1];
			v_hH1err += Jxw * (innerProduct(df_uh_gradient, df_uh_gradient) + innerProduct(df_vh_gradient, df_vh_gradient));
			cont += Jxw * (df_uh * df_uh + df_vh * df_vh);
		}
		err_ele[the_element_v->index()] = sqrt(fabs(cont));
		/// 计算误差.
		v_hL2err += cont;
		/// 将L2误差加到u_hH1err上.
		v_hH1err += cont;
	}
	/// 得到总误差.
	v_hL2err = sqrt(fabs(v_hL2err));
	v_hH1err = sqrt(fabs(v_hH1err));

	/// 记录计算结果参数.
	std::ofstream output;
	output.open("record", std::ofstream::out | std::ofstream::app);
	output.setf(std::ios::fixed);
	output.precision(20);
	output << "nu = " << fem_space_v.n_dof() << std::endl;
	output << "np = " << fem_space_p.n_dof() << std::endl;
	output << "t = " << _t << std::endl;
	output << "dt = " << dt << std::endl;

	double error;
	error = Functional::L2Error(v_h[0], divx, 3);
	std::cout << "|| u - u_h ||_L2 = " << error << std::endl;
	output << "|| u - u_h ||_L2 = " << error << std::endl;
	error = Functional::H1SemiError(v_h[0], divx, 3);
	std::cout << "|| u - u_h ||_H1 = " << error << std::endl;
	output << "|| u - u_h ||_H1 = " << error << std::endl;
	error = Functional::L2Error(v_h[1], divy, 3);
	std::cout << "|| v - v_h ||_L2 = " << error << std::endl;
	output << "|| v - v_h ||_L2 = " << error << std::endl;
	error = Functional::H1SemiError(v_h[1], divy, 3);
	std::cout << "|| v - v_h ||_H1 = " << error << std::endl;
	output << "|| v - v_h ||_H1 = " << error << std::endl;
	std::cout << "||vec{u} - vec{u_h}||_L2 = " << v_hL2err << std::endl;
	output << "||vec{u} - vec{u_h}||_L2 = " << v_hL2err << std::endl;
	std::cout << "||vec{u} - vec{u_h}||_H1 = " << v_hH1err << std::endl;
	output << "||vec{u} - vec{u_h}||_H1 = " << v_hH1err << std::endl;
	std::cout << "|| p - p_h ||_0 = " << p_hL0err << std::endl;
	output << "|| p - p_h ||_0 = " << p_hL0err << std::endl;
	error = Functional::H1SemiError(p_h, dip, 3);
	std::cout << "|| p - p_h ||_H1 = " << error << std::endl;
	output << "|| p - p_h ||_H1 = " << error << std::endl;
}

void ISOP2P1::computeDiv_and_Vor()
{
	int n_ele = fem_space_v.n_element();
        divergence.resize(n_ele, 0.0);
        vorticity.reinit(fem_space_v);
        MassMatrix<DIM, double> mat;
        mat.reinit(fem_space_v);
        mat.algebricAccuracy() = 3;
        mat.build();
        Vector<double> rhs_vor(fem_space_v.n_dof());

        /// 用于计算整体散度.
        double div = 0.0;
        double vor = 0.0;
        /// 真实涡量.
        RealVorticity real_vor(viscosity, t);

        FEMSpace<double, DIM>::ElementIterator the_element_v = fem_space_v.beginElement();
	FEMSpace<double, DIM>::ElementIterator end_element_v = fem_space_v.endElement();
        for (; the_element_v != end_element_v; ++the_element_v)
        {
		/// 当前单元信息.
		double volume = the_element_v->templateElement().volume();
		/// 积分精度, u 和 p 都是 1 次, 梯度和散度 u 都是常数.
		const QuadratureInfo<DIM>& quad_info = the_element_v->findQuadratureInfo(3);
		std::vector<double> jacobian = the_element_v->local_to_global_jacobian(quad_info.quadraturePoint());
		int n_quadrature_point = quad_info.n_quadraturePoint();
		std::vector<Point<DIM> > q_point = the_element_v->local_to_global(quad_info.quadraturePoint());

		/// 速度单元信息.
		std::vector<std::vector<double> > basis_value = the_element_v->basis_function_value(q_point);
		std::vector<std::vector<double> > u_h_gradient = v_h[0].gradient(q_point, *the_element_v);
		std::vector<std::vector<double> > v_h_gradient = v_h[1].gradient(q_point, *the_element_v);

		std::vector<int> element_dof_v = the_element_v->dof();
		int n_ele_dof = the_element_v->n_dof();
		double cont_div = 0.0, cont_vor = 0.0;
		for (int l = 0; l < n_quadrature_point; ++l)
		{
			double Jxw = jacobian[l] * quad_info.weight(l) * volume;
			double real_vorticity = real_vor.value(q_point[l]);
			cont_div += Jxw * pow(fabs(u_h_gradient[l][0] + v_h_gradient[l][1]), 2);
                        cont_vor += Jxw * pow(fabs(v_h_gradient[l][0] - u_h_gradient[l][1] - real_vorticity), 2);
                        for (int i = 0; i < n_ele_dof; ++i)
                        	rhs_vor(element_dof_v[i]) += Jxw * basis_value[i][l] * (v_h_gradient[l][0] - u_h_gradient[l][1]);
		}
		div += cont_div;
		vor += cont_vor;
                divergence[the_element_v->index()] = cont_div;
        }
        div = sqrt(fabs(div));
        vor = sqrt(fabs(vor));
        AMGSolver solver(mat);
        solver.solve(vorticity, rhs_vor);
        double error = 0.0;
//       Operator::L2Project(real_vor, vorticity, Operator::LOCAL_LEAST_SQUARE, 5);
        error = Functional::L2Error(vorticity, real_vor, 3);
        std::cout << "vorticity L2norm = " << error << std::endl;
        std::cout << "divergence L2norm = " << div << std::endl;
        std::cout << "vorticity L2norm = " << vor << std::endl;
};
