#include <AFEPack/AMGSolver.h>
#include <lac/sparse_matrix.h>
#include <lac/sparsity_pattern.h>
#include <lac/block_vector.h>
#include <lac/full_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/solver_bicgstab.h>
#include <lac/precondition.h>

#include <vector>

#define DIM 2

/**
 * 对应的预处理.
 *
 */
class StokesPreconditioner
{
private:
	const SparseMatrix<double> *Ax; /**< 预处理矩阵各分块. */
	const SparseMatrix<double> *Ay;
	AMGSolver AMGx; /**< 预处理矩阵各分块. */
	AMGSolver AMGy;
	const SparseMatrix<double> *Q;
	double delta_t;
public:
	StokesPreconditioner();

	~StokesPreconditioner();

	/**
	 * 预处理子初始化.
	 *
	 * @param _stiff_vx vx 空间的刚度矩阵.
	 * @param _stiff_vy vy 空间的刚度矩阵.
	 * @param _mass_p_diag p 空间的质量矩阵的对角元.
	 */
	void initialize (const SparseMatrix<double> &_stiff_vx,
			 const SparseMatrix<double> &_stiff_vy,
			 const SparseMatrix<double> &_mass_p_diag,
			 double _delta_t);

	/**
	 * 实际估值 dst = M^{-1}src.
	 *
	 * @param dst
	 * @param src
	 */
	void vmult (Vector<double> &dst,
			const Vector<double> &src) const;
};

class InverseMatrix : public Subscriptor
{
public:
	InverseMatrix (const SparseMatrix<double> &m, const AMGSolver &a);
	~InverseMatrix ();

	void vmult (Vector<double>       &dst,
			const Vector<double> &src) const;

private:
	const SmartPointer<const SparseMatrix<double> > matrix;
	const AMGSolver *amg;
};



class SchurComplement : public Subscriptor
{
public:
	SchurComplement(SparseMatrix<double> &_BTx,
			SparseMatrix<double> &_BTy,
			SparseMatrix<double> &_Bx,
			SparseMatrix<double> &_By,
			SparseMatrix<double> &_Q,
			InverseMatrix &_AInvx,
			InverseMatrix &_AInvy);

	void vmult (Vector<double>       &dst,
		    const Vector<double> &src) const;

private:
	const SmartPointer<const SparseMatrix<double> > BTx;
	const SmartPointer<const SparseMatrix<double> > BTy;
	const SmartPointer<const SparseMatrix<double> > Bx;
	const SmartPointer<const SparseMatrix<double> > By;
	const SmartPointer<const SparseMatrix<double> > Q;
	const SmartPointer<const InverseMatrix> AInvx;
	const SmartPointer<const InverseMatrix> AInvy;
	
	mutable Vector<double> tmp11, tmp12, tmp21, tmp22;
};

class ApproxSchurComplement
{
private:
	const SmartPointer<const SparseMatrix<double> > Q;
	const AMGSolver *amg;

public:
	ApproxSchurComplement(SparseMatrix<double> &_Q,
			      const AMGSolver &_amg);

	/**
	 * 实际估值 dst = M^{-1}src.
	 *
	 * @param dst
	 * @param src
	 */
	void vmult (Vector<double> &dst,
		    const Vector<double> &src) const;
};

class LSCPreconditioner : public Subscriptor
{
public:
	LSCPreconditioner(SparseMatrix<double> &_BTx,
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
			  InverseMatrix &_AInvy);

	void vmult (Vector<double>       &dst,
		    const Vector<double> &src) const;

private:
        const SmartPointer<const SparseMatrix<double> > BTx;
        const SmartPointer<const SparseMatrix<double> > BTy;
        const SmartPointer<const SparseMatrix<double> > Bx;
        const SmartPointer<const SparseMatrix<double> > By;
        const SmartPointer<const SparseMatrix<double> > Axx;
	const SmartPointer<const SparseMatrix<double> > Ayy;
	const SmartPointer<const SparseMatrix<double> > Q;
	const SmartPointer<const SchurComplement> schur_complement;

	ApproxSchurComplement &approx_schur_complement;
	const SmartPointer<const InverseMatrix> QInv;
        const SmartPointer<const InverseMatrix> AInvx;
        const SmartPointer<const InverseMatrix> AInvy;

	mutable Vector<double> tmp11, tmp12, tmp13, tmp14, tmp21, tmp22, tmp23, tmp24;
	
};

class updateSolutionPreconditioner : public Subscriptor
{
private:
        const SmartPointer<const SparseMatrix<double> > BTx;
        const SmartPointer<const SparseMatrix<double> > BTy;
	const SmartPointer<const SchurComplement> schur_complement;

	ApproxSchurComplement &approx_schur_complement;
        const SmartPointer<const InverseMatrix> AInvx;
        const SmartPointer<const InverseMatrix> AInvy;
	
	mutable Vector<double> tmp11, tmp12;

public:
	updateSolutionPreconditioner(SparseMatrix<double> &_BTx, 
				     SparseMatrix<double> &_BTy,
				     InverseMatrix &_AInvx,
				     InverseMatrix &_AInvy,
				     SchurComplement &_schur_complement,
				     ApproxSchurComplement _approx_schur_complement);

	~updateSolutionPreconditioner();
public:
	void vmult (Vector<double> &dst,
		    const Vector<double> &src) const;
		
};
#undef DIM
