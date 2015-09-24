#include "ISOP2P1.h"
#define DIM 2

void ISOP2P1::outputTecplot(const std::string &prefix)
{
    RegularMesh<DIM> &mesh_p = irregular_mesh_p->regularMesh();
    RegularMesh<DIM> &mesh_v = irregular_mesh_v->regularMesh();
    int n_node = mesh_v.n_geometry(0);
    int n_ele = mesh_v.n_geometry(2);

    FEMFunction <double, DIM> p_h_refine(fem_space_v);

    Operator::L2Interpolate(p_h, p_h_refine);
    std::stringstream result;
    result.setf(std::ios::fixed);
    result.precision(4);
    result << prefix << ".dat";
    std::ofstream tecplot(result.str().c_str()); 
    tecplot.setf(std::ios::fixed);
    tecplot.precision(20);
    tecplot << "VARIABLES = \"X\", \"Y\", \"P\", \"U\", \"V\"";
	tecplot << std::endl;
    tecplot << "ZONE NODES=" << n_node << ", ELEMENTS=" << n_ele << ", DATAPACKING=BLOCK," << std::endl;
    tecplot << "ZONETYPE=FETRIANGLE" << std::endl;
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
    for (int i = 0; i < n_ele; ++i)
    {
	std::vector<int> &vtx =  fem_space_v.element(i).geometry().vertex();
        tecplot << vtx[0] + 1 << "\n" << vtx[1] + 1 << "\n" << vtx[2] + 1 << std::endl;
    }
    tecplot.close();
};

void ISOP2P1::outputTecplotP(const std::string &prefix)
{
    RegularMesh<DIM> &mesh_v = irregular_mesh_v->regularMesh();
    RegularMesh<DIM> &mesh_p = irregular_mesh_p->regularMesh();
    int n_node = mesh_p.n_geometry(0);
    int n_ele = mesh_p.n_geometry(2);

    std::stringstream result;
    result.setf(std::ios::fixed);
    result.precision(4);
    result << prefix << ".dat";
    std::ofstream tecplot(result.str().c_str()); 
    tecplot.setf(std::ios::fixed);
    tecplot.precision(20);
    tecplot << "VARIABLES = \"X\", \"Y\", \"P\"";
    tecplot << std::endl;
    tecplot << "ZONE NODES=" << n_node << ", ELEMENTS=" << n_ele << ", DATAPACKING=BLOCK," << std::endl;
    tecplot << "ZONETYPE=FETRIANGLE" << std::endl;
    for (int i = 0; i < n_node; ++i)
        tecplot << mesh_p.point(i)[0] << "\n";
    tecplot << std::endl;
    for (int i = 0; i < n_node; ++i)
        tecplot << mesh_p.point(i)[1] << "\n";
    tecplot << std::endl;
    for (int i = 0; i < n_node; ++i)
        tecplot << p_h(i) << "\n";
    tecplot << std::endl;
    for (int i = 0; i < n_ele; ++i)
    {
	std::vector<int> &vtx =  fem_space_p.element(i).geometry().vertex();
        tecplot << vtx[0] + 1 << "\n" << vtx[1] + 1 << "\n" << vtx[2] + 1 << std::endl;
    }
    tecplot.close();
};

#undef DIM
