#include "ISOP2P1.h"
#include "functions.h"
#define DIM 2

void ISOP2P1::outputTecplot(const std::string &prefix)
{
	/// 计算速度单元上的涡量和散度.
        computeDiv_and_Vor();
	RegularMesh<DIM> &mesh_p = irregular_mesh_p->regularMesh();
	RegularMesh<DIM> &mesh_v = irregular_mesh_v->regularMesh();
	int n_node = mesh_v.n_geometry(0);
	int n_ele = mesh_v.n_geometry(2);

	FEMFunction <double, DIM> p_h_refine(fem_space_v);
	Operator::L2Interpolate(p_h, p_h_refine);

	/// 减掉均值.
	p_h_refine.add(-Functional::meanValue(p_h_refine, 3));

	DiVx divx(viscosity, t);
	DiVy divy(viscosity, t);
	// RealVx divx;
	// RealVy divy;
	std::stringstream result;
	result.setf(std::ios::fixed);
	result.precision(4);
	result << prefix << ".dat";
	std::ofstream tecplot(result.str().c_str());
	tecplot.setf(std::ios::fixed);
	tecplot.precision(20);
	tecplot << "VARIABLES = \"X\", \"Y\", \"P\", \"U\", \"V\", \"L2Error\", \"divergence\", \"vorticity\"";
	tecplot << std::endl;
	tecplot << "ZONE NODES=" << n_node << ", ELEMENTS=" << n_ele << ", DATAPACKING=BLOCK," << std::endl;
	tecplot << "VARLOCATION=([6, 7]=CELLCENTERED)," << "ZONETYPE=FETRIANGLE" << std::endl;
	for (int i = 0; i < n_node; ++i)
		tecplot << mesh_v.point(i)[0] << "\n";
	tecplot << std::endl;
	for (int i = 0; i < n_node; ++i)
		tecplot << mesh_v.point(i)[1] << "\n";
	tecplot << std::endl;
	for (int i = 0; i < n_node; ++i)
		tecplot << p_h_refine(i) << "\n";
	tecplot << std::endl;
	for (int i = 0; i < n_node; ++i)
		tecplot << v_h[0](i) << "\n";
	tecplot << std::endl;
	for (int i = 0; i < n_node; ++i)
		tecplot << v_h[1](i) << "\n";
	tecplot << std::endl;
//	for (int i = 0; i < n_node; ++i)
//	{
//		double error = fabs(divx.value(mesh_v.point(i)) - v_h[0](i));
//		tecplot << error << "\n";
//	}
//	for (int i = 0; i < n_node; ++i)
//	{
//		double error = fabs(divy.value(mesh_v.point(i)) - v_h[1](i));
//		tecplot << error << "\n";
//	}
	for (int i = 0; i < n_ele; ++i)
		tecplot << err_ele[i] << "\n";
	tecplot << std::endl;
	for (int i = 0; i < n_ele; ++i)
		tecplot << divergence[i] << "\n";
	tecplot << std::endl;
	for (int i = 0; i < n_node; ++i)
		tecplot << vorticity(i) << "\n";
	tecplot << std::endl;

	for (int i = 0; i < n_ele; ++i)
	{
		std::vector<int> &vtx =  fem_space_v.element(i).geometry().vertex();
		tecplot << vtx[0] + 1 << "\n" << vtx[1] + 1 << "\n" << vtx[2] + 1 << std::endl;
	}
	tecplot.close();
};

void ISOP2P1::outputTecplotP(const std::string &prefix)
{
	computeMonitor();
	RegularMesh<DIM> &mesh_p = irregular_mesh_p->regularMesh();
	int n_node = mesh_p.n_geometry(0);
	int n_ele = mesh_p.n_geometry(2);
	/// 直接在压力数值解p上减掉均值.
	//    FEMFunction<double, DIM> _p_h(p_h);
	p_h.add(-Functional::meanValue(p_h, 3));
	//    /// 用于burgurs方程的输出，如果不是需要注释掉.
	//    Operator::L2Interpolate(v_h[0], p_h);


	std::stringstream result;
	result.setf(std::ios::fixed);
	result.precision(4);
	result << prefix << ".dat";
	std::ofstream tecplot(result.str().c_str());
	tecplot.setf(std::ios::fixed);
	tecplot.precision(20);
	tecplot << "VARIABLES = \"X\", \"Y\", \"P\", \"delta_x\", \"delta_y\", \"monitor\"";
	//    tecplot << "VARIABLES = \"X\", \"Y\", \"P\"";
	tecplot << std::endl;
	tecplot << "ZONE NODES=" << n_node << ", ELEMENTS=" << n_ele << ", DATAPACKING=BLOCK," << std::endl;
	tecplot << "VARLOCATION=([6]=CELLCENTERED)," << "ZONETYPE=FETRIANGLE" << std::endl;
	for (int i = 0; i < n_node; ++i)
		tecplot << mesh_p.point(i)[0] << "\n";
	tecplot << std::endl;
	for (int i = 0; i < n_node; ++i)
		tecplot << mesh_p.point(i)[1] << "\n";
	tecplot << std::endl;
	for (int i = 0; i < n_node; ++i)
		tecplot << p_h(i) << "\n";
	tecplot << std::endl;
	for (int i = 0; i < n_node; ++i)
		tecplot << mesh_p.point(i)[0] - mesh_bak.point(i)[0] << "\n";
	tecplot << std::endl;
	for (int i = 0; i < n_node; ++i)
		tecplot << mesh_p.point(i)[1] - mesh_bak.point(i)[1] << "\n";
	tecplot << std::endl;
	for (int i = 0; i < n_ele; ++i)
	{
		if (isMoving == 0)
			tecplot << Monitor[i]<< "\n";
		else
			tecplot << monitor(i) << "\n";
	}
	tecplot << std::endl;
	for (int i = 0; i < n_ele; ++i)
	{
		std::vector<int> &vtx =  fem_space_p.element(i).geometry().vertex();
		tecplot << vtx[0] + 1 << "\n" << vtx[1] + 1 << "\n" << vtx[2] + 1 << std::endl;
	}
	tecplot.close();
};

void ISOP2P1::computeMonitor()
{
	if (isMoving == 0)
	{
		/// 将速度有限元空间的数值解插值到压力有限元空间中，
		/// 为了方便计算压力网格单元上的monitor的值.
		FEMFunction<double, DIM> _u_h(fem_space_p);
		FEMFunction<double, DIM> _v_h(fem_space_p);
		Operator::L2Interpolate(v_h[0], _u_h);
		Operator::L2Interpolate(v_h[1], _v_h);

		RegularMesh<DIM> &mesh_p = irregular_mesh_p->regularMesh();

		int n_ele = mesh_p.n_geometry(DIM);
		Monitor.resize(n_ele, 0.0);

		/// 选取速度数值解的梯度做monitor, 跟文章中的monitor一致.
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
				d += Jxw * (innerProduct(u_h_gradient[l], u_h_gradient[l]) + innerProduct(v_h_gradient[l], v_h_gradient[l]));
			}
			Monitor[the_element->index()] = d / area;
		}
		for (int i = 0; i < n_ele; ++i)
			Monitor[i] = sqrt(1.0 + alpha * pow(Monitor[i], beta));

		/* 后验误差做monitor。
		std::vector<double> side_length(n_face);
	        std::vector<bool> flag(n_face, false);
	        std::vector<double> jump_ux(n_face);
	        std::vector<double> jump_uy(n_face);

	        FEMSpace<double,DIM>::ElementIterator the_element_p = fem_space_p.beginElement();
	        FEMSpace<double,DIM>::ElementIterator end_element_p = fem_space_p.endElement();
	        /// 遍历所有速度单元,计算每条边上的 jump / ||v||_L2.
	        for (; the_element_p != end_element_p; ++the_element_p)
	        {
	        	/// 几何信息.
	        	double volume = the_element_p->templateElement().volume();
	        	const QuadratureInfo<2>& quad_info = the_element_p->findQuadratureInfo(4);
	        	std::vector<double> jacobian = the_element_p->local_to_global_jacobian(quad_info.quadraturePoint());
	        	int n_quadrature_point = quad_info.n_quadraturePoint();
	        	std::vector<Point<2> > q_point = the_element_p->local_to_global(quad_info.quadraturePoint());

	        	const GeometryBM& geometry = the_element_p->geometry();

	        	for (int i = 0; i < geometry.n_boundary(); ++i)
	        	{
	        		int j = geometry.boundary(i);
	        		const GeometryBM& side = mesh_p.geometry(1, j);
	        		const Point<DIM>& p0 = mesh_p.point(mesh_p.geometry(0, side.boundary(0)).vertex(0));
	        		const Point<DIM>& p1 = mesh_p.point(mesh_p.geometry(0, side.boundary(1)).vertex(0));
	        		std::vector<double> vx_gradient = _u_h.gradient(midpoint(p0, p1), *the_element_p);
	        		std::vector<double> vy_gradient = _v_h.gradient(midpoint(p0, p1), *the_element_p);
	        		double vx_value = _u_h.value(midpoint(p0, p1), *the_element_p);
	        		double vy_value = _v_h.value(midpoint(p0, p1), *the_element_p);

	        		double v_L2norm = sqrt(vx_value * vx_value + vy_value * vy_value + eps);
	        		side_length[j] = distance(p0, p1);
	        		if (flag[j])
	        		{
	        			jump_ux[j] -= (vx_gradient[0] * (p0[1] - p1[1]) + vx_gradient[1] * (p1[0] - p0[0]));
	        			jump_uy[j] -= (vy_gradient[0] * (p0[1] - p1[1]) + vy_gradient[1] * (p1[0] - p0[0]));
	        			flag[j] = false;
	        		}
	        		else
	        		{
	        			jump_ux[j] = (vx_gradient[0] * (p0[1] - p1[1]) + vx_gradient[1] * (p1[0] - p0[0]));
	        			jump_uy[j] = (vy_gradient[0] * (p0[1] - p1[1]) + vy_gradient[1] * (p1[0] - p0[0]));
	        			flag[j] = true;
	        		}
	        	}
	        }
	        the_element_p = fem_space_p.beginElement();
	        for(int i = 0; the_element_p != end_element_p; ++the_element_p, ++i)
	        {
	        	const GeometryBM &geometry = the_element_p->geometry();
	        	double cont = 0.0, total_ele_side_length = 0.0;
	        	for (int l = 0; l < geometry.n_boundary(); ++l)
	        	{
	        		int j = geometry.boundary(l);
	        		/// 如果编号是j的边是边界.
	        		if (flag[j])
	        			continue;
	        		/// 将DIM个梯度沿法向方向上的跳跃平方值累加.
	        		cont += (jump_ux[j] * jump_ux[j] + jump_uy[j] * jump_uy[j]) * side_length[j];
	        		total_ele_side_length += side_length[j];
	        	}
	        	Monitor[i] = sqrt(cont / total_ele_side_length);
	        }
//		for (int i = 0; i < n_geometry(DIM); ++i)
//		        monitor(i) = 1.0 / sqrt(1.0 + alpha * pow(monitor(i) / max_monitor, 2));
		 */
	}
}
#undef DIM
