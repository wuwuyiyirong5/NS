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
	if (isMoving == 1)
	{
		maxStep() = 30;
		double scale_step = 0.2;
		scale = scale_step;
		do {
			DiVx real_vx(viscosity, t);
			DiVy real_vy(viscosity, t);
			Operator::L2Project(real_vx, v_h[0], Operator::LOCAL_LEAST_SQUARE, 3);
			Operator::L2Project(real_vy, v_h[1], Operator::LOCAL_LEAST_SQUARE, 3);
			outputTecplot("initial");
			v_h[0].scale(scale);
			v_h[1].scale(scale);
			movingMesh();
			v_h[0].scale(1.0 / scale);
			v_h[1].scale(1.0 / scale);
			computeError(t);
			outputTecplot("Stokes");
			scale += scale_step;
			getchar();
		} while (scale <= 1.0);
		/// 重新设置scale = 1.0。
		scale = 1.0;
	}
	DiVx real_vx(viscosity, t);
	DiVy real_vy(viscosity, t);
	Operator::L2Project(real_vx, v_h[0], Operator::LOCAL_LEAST_SQUARE, 3);
	Operator::L2Project(real_vy, v_h[1], Operator::LOCAL_LEAST_SQUARE, 3);
	
	// /// 输出一下初始网格上的初始误差.
	// computeError(t);
	std::cout << "Stokes_init = " << Stokes_init << std::endl;
	if (Stokes_init == true)
	{
		buildMatrix();
		solveStokes();
		outputTecplotP("initial_valueP");
		outputTecplot("initial_valueV");
	}
	time_step();
	std::cout << "Simulation start at t = " << t << ", dt = " << dt << std::endl;
	getchar();
	while (t < t1)
	{
		if (isMoving == 1)
		{
//			buildMatrixStruct();
			buildMatrix();
		}
		if (scheme == 1)
			stepForwardEuler();
		else if (scheme == 2)
			stepForwardLinearizedEuler(); 

		static int k = 1;
		std::stringstream result;
		result.setf(std::ios::fixed);
		result.precision(4);
		result << "V" << k++;
		outputTecplot(result.str().c_str());

		if (isMoving == 1)
		{
			movingMesh();
			outputTecplotP("monitor");
		}
		time_step();
		t += dt;
		std::cout << "t = " << t << ", the next dt = " << dt << std::endl;
		getchar();
	}
};

#undef DIM
