#include "preconditioner.h"
#define DIM 2

LSCPreconditioner::LSCPreconditioner(SparseMatrix<double> &_BTx,
				     SparseMatrix<double> &_BTy,
				     SparseMatrix<double> &_Bx,
				     SparseMatrix<double> &_By,
				     SparseMatrix<double> &_Axx,
				     SparseMatrix<double> &_Ayy,
				     SparseMatrix<double> &_Q,
				     SchurComplement &_schur_complement,		    
				     ApproxSchurComplement _approx_schur_complement,
				     InverseMatrix &_QInv,
				     InverseMatrix &_AInvx,
				     InverseMatrix &_AInvy)
	:
	BTx(&_BTx),
	BTy(&_BTy),
	Bx(&_Bx),
	By(&_By),
	Axx(&_Axx),
	Ayy(&_Ayy),
	schur_complement(&_schur_complement),
	approx_schur_complement(_approx_schur_complement),
	AInvx(&_AInvx),
	AInvy(&_AInvy),
	QInv(&_QInv),
	Q(&_Q),
	tmp11(BTx->m()),
	tmp12(BTy->m()),
	tmp13(BTx->m()),
	tmp14(BTx->m()),
	tmp21(BTx->m()),
	tmp22(BTy->m()),
	tmp23(BTy->m()),
	tmp24(BTy->m())
{
};
void LSCPreconditioner::vmult (Vector<double>       &dst,
		const Vector<double> &src) const
{
	int n_dof_v = BTx->m();
	int n_dof_p = BTx->n();
	Vector<double> d0(n_dof_v);
	Vector<double> d1(n_dof_v);
	Vector<double> d2(n_dof_p);
	Vector<double> s0(n_dof_v);
	Vector<double> s1(n_dof_v);
	Vector<double> s2(n_dof_p);
	Vector<double> y0(n_dof_p);
	Vector<double> y1(n_dof_p);
	for (int i = 0; i < n_dof_v; ++i)
	{
		s0(i) = src(i);
		s1(i) = src(i + n_dof_v);
	}
	for (int i = 0; i < n_dof_p; ++i)
		s2(i) = src(2 * n_dof_v + i);
	
	/// 求解B T^{-1}B^T y0 = s2. 这里tol很大,不用求解精确.
	SolverControl solver_control(n_dof_p, 1e-3 * s2.l2_norm(), 1);
	SolverBicgstab<> cg(solver_control);
//	cg.solve(*schur_complement, d2, s2, approx_schur_complement);

	cg.solve(*schur_complement, y0, s2, approx_schur_complement);
	
	/// 完成(y1 = Bx T^{-1} F T^{-1} Bx^T y0).
	BTx->vmult(tmp11, y0);
        QInv->vmult(tmp12, tmp11);
 	// for (int i = 0; i < n_dof_v; ++i)
	// 	tmp12(i) = tmp11(i) / Q->diag_element(i); 
	Axx->vmult(tmp13, tmp12);
        QInv->vmult(tmp14, tmp13);
	// for (int i = 0; i < n_dof_v; ++i)
	// 	tmp14(i) = tmp13(i) / Q->diag_element(i); 
	Bx->vmult(y1, tmp14);
	

	/// 完成(y1 += By T^{-1} F T^{-1} By^T y0).
	BTy->vmult(tmp21, y0);
	QInv->vmult(tmp22, tmp21);
	// for (int i = 0; i < n_dof_v; ++i)
	// 	tmp22(i) = tmp21(i) / Q->diag_element(i);
	Ayy->vmult(tmp23, tmp22);
	QInv->vmult(tmp24, tmp23);
	// for (int i = 0; i < n_dof_v; ++i)
	// 	tmp24(i) = tmp23(i) / Q->diag_element(i);
	By->vmult_add(y1, tmp24);
	
	/// 求解B T^{-1} B^T d2 = y1. 
	cg.solve(*schur_complement, d2, y1, approx_schur_complement);
		
	tmp11.reinit(n_dof_v);
	tmp21.reinit(n_dof_v);
	BTx->vmult(tmp11, d2);
	s0.add(-1.0, tmp11);
	AInvx->vmult(d0, s0);
	
	BTy->vmult(tmp21, d2);
	s1.add(-1.0, tmp21);
	AInvy->vmult(d1, s1);

	for (int i = 0; i < n_dof_v; ++i)
	{
		dst(i) = d0(i);
		dst(i + n_dof_v) = d1(i);
	}
	for (int i = 0; i < n_dof_p; ++i)
		dst(i + 2 * n_dof_v) = d2(i);
};

SchurComplement::SchurComplement(SparseMatrix<double> &_BTx,
				 SparseMatrix<double> &_BTy,
				 SparseMatrix<double> &_Bx,
				 SparseMatrix<double> &_By,
				 SparseMatrix<double> &_Q,
				 InverseMatrix &_AInvx,
				 InverseMatrix &_AInvy)
		:
		BTx(&_BTx),
		BTy(&_BTy),
		Bx(&_Bx),
		By(&_By),
		Q(&_Q),
		AInvx(&_AInvx),
		AInvy(&_AInvy),
		tmp11(BTx->m()),
		tmp12(BTy->m()),
		tmp21(BTx->m()),
		tmp22(BTy->m())
{
};


void SchurComplement::vmult (Vector<double>       &dst,
		const Vector<double> &src) const
{
	BTx->vmult (tmp11, src);
	BTy->vmult (tmp12, src);
	// amg->solve(tmp21, tmp11, 1e-8 * tmp11.l2_norm(), 0);
	// amg->solve(tmp22, tmp12, 1e-8 * tmp11.l2_norm(), 0);
		
        AInvx->vmult(tmp21, tmp11);
        AInvy->vmult(tmp22, tmp12);
	// int n_dof_v = BTx->m();
	// for (int i = 0; i < n_dof_v; ++i)
	// {
	// 	tmp21(i) = tmp11(i) / Q->diag_element(i);
	// 	tmp22(i) = tmp12(i) / Q->diag_element(i);
	// }
	Bx->vmult(dst, tmp21);
	By->vmult_add(dst,tmp22);
};

StokesPreconditioner::StokesPreconditioner()
{
	Ax = NULL;
	Ay = NULL;
	Q = NULL;
	delta_t = 0.0;
};

StokesPreconditioner::~StokesPreconditioner()
{};

void StokesPreconditioner::initialize (const SparseMatrix<double> &_stiff_vx, 
				       const SparseMatrix<double> &_stiff_vy,
				       const SparseMatrix<double> &_mass_p_diag,
				       double _delta_t)
{
	Ax = &_stiff_vx;
	Ay = &_stiff_vy;
	Q = &_mass_p_diag;
	AMGx.reinit(*Ax);
	AMGy.reinit(*Ay);
	delta_t = _delta_t;
};

void StokesPreconditioner::vmult (Vector<double> &dst,
		const Vector<double> &src) const
{
	int n_dof_v = Ax->n();
	int n_dof_p = Q->n();
	Vector<double> d0(n_dof_v);
	Vector<double> d1(n_dof_v);
	Vector<double> d2(n_dof_p);
	Vector<double> s0(n_dof_v);
	Vector<double> s1(n_dof_v);
	Vector<double> s2(n_dof_p);

	for (int i = 0; i < n_dof_v; ++i)
		s0(i) = src(i);
	for (int i = 0; i < n_dof_v; ++i)
		s1(i) = src(n_dof_v + i);
	for (int i = 0; i < n_dof_p; ++i)
		s2(i) = src(2 * n_dof_v + i);
	for (int i = 0; i < n_dof_p; ++i)
		dst(2 * n_dof_v + i) = delta_t * src(2 * n_dof_v + i) / (*Q).diag_element(i);

	AMGx.solve(d0, s0, 1e-8 * s0.l2_norm(), 1);

	AMGy.solve(d1, s1, 1e-8 * s1.l2_norm(), 1);
	
	for (int i = 0; i < n_dof_v; ++i)
		dst(i) = d0(i);
	for (int i = 0; i < n_dof_v; ++i)
		dst(n_dof_v + i) = d1(i);
};

InverseMatrix::InverseMatrix (const SparseMatrix<double> &m,  const AMGSolver &a)
:matrix (&m), amg(&a)
{
};

InverseMatrix::~InverseMatrix()
{
};

void InverseMatrix::vmult(Vector<double> &dst,
		const Vector<double> &src) const
{
	dst = 0;
	amg->solve(dst, src, 1e-8 * src.l2_norm(), 0);
};

ApproxSchurComplement::ApproxSchurComplement(SparseMatrix<double> &_Q,
					     const AMGSolver &_amg)
	:
	Q(&_Q),
	amg(&_amg)
{
};

void ApproxSchurComplement::vmult (Vector<double> &dst,
				   const Vector<double> &src) const
{
	dst.reinit(src.size());
	/// 预处理模式比求解模式，计算时间要快，迭代步数相差无几.
	amg->solve(dst, src, 1e-8 * src.l2_norm(), 0);
	// for (int i = 0; i < src.size(); ++i)
	// 	dst(i) = src(i) / Q->diag_element(i);
};

updateSolutionPreconditioner::updateSolutionPreconditioner(SparseMatrix<double> &_BTx,
							   SparseMatrix<double> &_BTy,
							   InverseMatrix &_AInvx,
							   InverseMatrix &_AInvy,
							   SchurComplement &_schur_complement,
							   ApproxSchurComplement _approx_schur_complement)
        :
	BTx(&_BTx),
	BTy(&_BTy),
	AInvx(&_AInvx),
	AInvy(&_AInvy),
	schur_complement(&_schur_complement),
	approx_schur_complement(_approx_schur_complement),
	tmp11(BTx->m()),
	tmp12(BTy->m())
	
{};
void updateSolutionPreconditioner::vmult(Vector<double> &dst,
					 const Vector<double> &src) const
{
	int n_dof_v = BTx->m();
	int n_dof_p = BTx->n();
	Vector<double> d0(n_dof_v);
	Vector<double> d1(n_dof_v);
	Vector<double> d2(n_dof_v);
	Vector<double> s0(n_dof_v);
	Vector<double> s1(n_dof_v);
	Vector<double> s2(n_dof_v);
	for (int i = 0; i < n_dof_v; ++i)
	{
		s0(i) = src(i);
		s1(i) = src(i + n_dof_v);
	}
	for (int i = 0; i < n_dof_p; ++i)
		s2(i) = src(i + 2 * n_dof_v);

	SolverControl solver_control(n_dof_p, 1e-3 * s2.l2_norm(), 1);
	SolverBicgstab<> bicg(solver_control);
	bicg.solve(*schur_complement, d2, s2, approx_schur_complement);
	
	BTx->vmult(tmp11, d2);
	BTy->vmult(tmp12, d2);
	
	s0.add(-1.0, tmp11);
	s1.add(-1.0, tmp12);
	
	AInvx->vmult(d0, s0);
	AInvy->vmult(d1, s1);

	for (int i = 0; i < n_dof_v; ++i)
	{
		dst(i) = d0(i);
		dst(i + n_dof_v) = d1(i);
	}
	for (int i = 0; i < n_dof_p; ++i)
		dst(i + 2 * n_dof_v) = d2(i);
}; 

updateSolutionPreconditioner::~updateSolutionPreconditioner()
{};

#undef DIM

