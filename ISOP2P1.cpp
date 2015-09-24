/**
 * @file   ISOP2P1.cpp
 * @author Heyu Wang <scshw@cslin107.csunix.comp.leeds.ac.uk>
 * @date   Tue Nov  4 14:33:17 2014
 * 
 * @brief ISO P2P1 有限元的网格, 空间, 矩阵信息集中封装类的实现.
 * 
 * 
 */

#include "ISOP2P1.h"
#include "functions.h"
#define DIM 2

void ISOP2P1::initialize()
{
    /// 读取配置文件.
    config("config");
    /// 构建网格.
    buildMesh();
    std::cout << "mesh builded." << std::endl;
    /// 构建混合有限元 ISOP1P2 空间.
    buildFEMSpace();
    std::cout << "space builded." << std::endl;
    /// 构建矩阵结构.
    buildMatrixStruct();
    std::cout << "matrix structure builded." << std::endl;
    buildMatrix();
    std::cout << "matrix builded." << std::endl;
};

void ISOP2P1::run()
{
	initialize();
	std::cout << "Stokes_init = " << Stokes_init << std::endl;
	if (Stokes_init == true)
	{
		solveStokes();
		outputTecplot("V0");
	}
//	else
//		stepForwardEuler();
//
//	time_step();
//
//	std::cout << "Simulation start at t = " << t << ", dt = " << dt << std::endl;
//	getchar();
//	while (t < t1)
//	{
//		if (scheme == 1)
//			stepForwardEuler();
//		else if (scheme == 2)
//			stepForwardLinearizedEuler();
//		else
//		{
//			std::cerr << "No such scheme." << std::endl;
//			exit(-1);
//		}
//
//		t += dt;
//		std::cout << "t = " << t << ", the next dt = " << dt << std::endl;
//
//		static int k = 1;
//	    std::stringstream result;
//	    result.setf(std::ios::fixed);
//	    result.precision(4);
//	    result << "V" << k++;
//		outputTecplot(result.str().c_str());
//	}
};

#undef DIM
