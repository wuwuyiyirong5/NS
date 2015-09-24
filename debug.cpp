#include "ISOP2P1.h"
#define DIM 2

void ISOP2P1::outputMat(const std::string &prefix, SparseMatrix<double> &mat)
{
    std::stringstream result;
    result.setf(std::ios::fixed);
    result.precision(4);
    result << prefix << ".m";
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
    result << prefix << ".m";
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
