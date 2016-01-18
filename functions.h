#include <AFEPack/Functional.h>
#include <vector>

#include "ISOP2P1.h"

#define DIM 2

class BurgV : public Function<double>
{
private:
	double a;
	double t;
public:
BurgV(double _a, double _t) : a(_a), t(_t) 
	{};
	~BurgV()
	{};
public:
	double value(const double *p) const
	{
		return 1.0 / (1.0 + exp((p[0] + p[1] - t) / (2.0 * a)));
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		return result;
	};
};


class DiVx : public Function<double>
{
private:
	double nu;
	double t;
public:
DiVx(double _nu, double _t) : nu(_nu), t(_t) 
	{};
	~DiVx()
	{};
public:
	double value(const double *p) const
	{
//		return PI * cos(PI *p[0]) * sin(PI * p[1]) * exp(-2 * PI * PI * nu * t);
		return -cos(p[0]) * sin(p[1]) * exp(-2.0 * nu * t);
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		/* result[0] = -PI * PI * sin(PI * p[0]) * sin(PI * p[1]) * exp(-2 * PI * PI * nu * t); */
		/* result[1] = PI * PI * cos(PI * p[0]) * cos(PI * p[1]) * exp(-2 * PI * PI * nu * t); */
		result[0] = sin(p[0]) * sin(p[1]) * exp(-2.0 * nu * t);
		result[1] = -cos(p[0]) * cos(p[1]) * exp(-2.0 * nu * t);
		return result;
	};
};

class DiVy : public Function<double>
{
private:
	double nu;
	double t;
public:
DiVy(double _nu, double _t) : nu(_nu), t(_t)
	{};
	~DiVy()
	{};
public:
	double value(const double *p) const
	{
//		return -PI * sin(PI * p[0]) * cos(PI * p[1]) * exp(-2.0 * PI * PI * nu * t);

		return sin(p[0]) * cos(p[1]) * exp(-2.0 * nu * t);
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		/* result[0] = -PI * PI * cos(PI * p[0]) * cos(PI * p[1]) * exp(-2 * PI * PI * nu * t); */
		/* result[1] = PI * PI * sin(PI * p[0]) * sin(PI * p[1]) * exp(-2 * PI * PI * nu * t); */
		result[0] = cos(p[0]) * cos(p[1]) * exp(-2.0 * nu * t);
		result[1] = -sin(p[0]) * sin(p[1]) * exp(-2.0 * nu * t);
		return result;
	};
};

class DiP : public Function<double>
{
private:
	double nu;
	double t;
	double p0;
public:
        DiP(double _nu, double _t, double _p0) : nu(_nu), t(_t), p0(_p0)
	{};
	~DiP()
	{};
public:
	double value(const double *p) const
	{
		return 0.5 * (sin(p[0]) * sin(p[0]) + sin(p[1]) * sin(p[1])) * exp(-4.0 * nu * t) + p0;
//		return 0.5 * PI * PI * (sin(PI * p[0]) * sin(PI * p[0]) + sin(PI * p[1]) * sin(PI * p[1])) * exp(-4.0 * PI * PI * nu * t);

	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		/* result[0] = PI * PI * PI * cos(PI * p[0]) * sin(PI * p[0]) * exp(-4.0 * PI * PI * nu * t); */
		/* result[1] = PI * PI * PI * cos(PI * p[1]) * sin(PI * p[1]) * exp(-4.0 * PI * PI * nu * t); */
		result[0] = cos(p[0]) * sin(p[0]) * exp(-4.0 * nu * t);
		result[1] = cos(p[1]) * sin(p[1]) * exp(-4.0 * nu * t);
		return result;
	};
};
class RealVx : public Function<double>
{
public:
	RealVx()
	{};
	~RealVx()
	{};
public:
	double value(const double *p) const
	{
		///  修改一下系数,使得解的奇异性变大.
//		return 4.0 * p[0] * p[1] * p[1] * p[1];
		return 20.0 * p[0] * p[1] * p[1] * p[1];
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
//		result[0] = 4.0 * p[1] * p[1] * p[1];
//		result[1] = 12.0 * p[0] * p[1] * p[1];
		result[0] = 20.0 * p[1] * p[1] * p[1];
		result[1] = 60.0 * p[0] * p[1] * p[1];
		return result;
	};
};

class RealVy : public Function<double>
{
public:
	RealVy()
	{};
	~RealVy()
	{};
public:
	double value(const double *p) const
	{
//		return 1.0 * p[0] * p[0] * p[0] * p[0] - 1.0 * p[1] * p[1] * p[1] * p[1];
		return 5.0 * p[0] * p[0] * p[0] * p[0] - 5.0 * p[1] * p[1] * p[1] * p[1];
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
//		result[0] = 4.0 * p[0] * p[0] * p[0];
//		result[1] = -4.0 * p[1] * p[1] * p[1];
		result[0] = 20.0 * p[0] * p[0] * p[0];
		result[1] = -20.0 * p[1] * p[1] * p[1];
		return result;
	};
};

class RealP : public Function<double>
{
private:
	double average;
public:
RealP(double _a) : average(_a)
	{};
	~RealP()
	{};
public:
	double value(const double *p) const
	{
//		return 12.0 * p[0] * p[0] * p[1] - 4.0 * p[1] * p[1] * p[1] + average;
		return 60.0 * p[0] * p[0] * p[1] - 20.0 * p[1] * p[1] * p[1] + average;
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
//		result[0] = 24.0 * p[0] * p[1];
//		result[1] = 12.0 * p[0] * p[0] - 240.0 * p[1] * p[1];
		result[0] = 120.0 * p[0] * p[1];
		result[1] = 60.0 * p[0] * p[0] - 60.0 * p[1] * p[1];
		return result;
	};
};

class RealVorticity : public Function<double>
{
private:
	double nu;
	double t;
public:
RealVorticity(double _nu, double _t): nu(_nu), t(_t)
	{};
	~RealVorticity()
	{};
public:
	double value(const double *p) const
	{
//		return 20 * p[0] * p[0] * p[0] - 60 * p[0] * p[1] * p[1];
		return -2 * PI * PI * cos(PI * p[0]) * cos(PI * p[1]) * exp(-2.0 * PI * PI * nu * t);
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		result[0] = 60 * p[0] * p[0] - 60 * p[1] * p[1];
		result[1] = -120 * p[0] * p[1];
		return result;
	};
};

class PoiseuilleVx : public Function<double>
{
private:
	double y0;
	double y1;
	double c;
	double v;
public:
PoiseuilleVx(double _y0, double _y1) : y0(_y0), y1(_y1)
	{
		c = (y0 + y1) * 0.5;
		v = (y0 - y1) * 0.5;
	};
	~PoiseuilleVx()
	{};
public:
	double value(const double *p) const
	{
		return   (1.0 - (p[1] - c) * (p[1] - c)  / (v * v));
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		result[0] = 0.0;
		result[1] = -2.0 * p[1];
		return result;
	};
};

class PoiseuilleVy : public Function<double>
{
public:
	PoiseuilleVy()
	{};
	~PoiseuilleVy() 
	{};
public:
	double value(const double *p) const
	{
		return 0.0;
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		result[0] = 0.0;
		result[1] = 0.0;
		return result;
	};
};

class PoiseuilleP : public Function<double>
{
private:
	double average;
public:
PoiseuilleP(double _a) : average(_a)
	{};
	~PoiseuilleP() 
	{};
public:
	double value(const double *p) const
	{
		return -2.0 * p[0] + average;
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		result[0] = -2.0;
		result[1] = 0.0;
		return result;
	};
};

class Regularized : public Function<double>
{
public:
	Regularized() 
	{};
	~Regularized() 
	{};
public:
	double value(const double *p) const
	{
		return 1.0 - p[0] * p[0] * p[0] * p[0];
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		result[0] = -4.0 * p[0] * p[0] * p[0];
		result[1] = 0.0;
		return result;
	};
};

class ForceX : public Function<double>
{
private:
	double body_force;
	double angle;
public:
ForceX(double _bf, double _an) : body_force(_bf),
		angle(_an)
		{};
	~ForceX()
	{};
public:
	double value(const double *p) const
	{
		double deg = angle / 180 * PI;
		return (body_force * cos(deg));
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		result[0] = 0.0;
		result[1] = 0.0;
		return result;
	};
};

class ForceY : public Function<double>
{
private:
	double body_force;
	double angle;
public:
ForceY(double _bf, double _an) : body_force(_bf),
		angle(_an)
		{};
	~ForceY()
	{};
public:
	double value(const double *p) const
	{
		double deg = angle / 180 * PI;
		return (body_force * sin(deg));
	};
	std::vector<double> gradient(const double *p) const
	{
		std::vector<double> result(2);
		result[0] = 0.0;
		result[1] = 0.0;
		return result;
	};
};

#undef DIM
