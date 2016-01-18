#include "ISOP2P1.h"
#define DIM 2

void ISOP2P1::buildMesh()
{
	std::cout << "The mesh tree is " << mesh_file << std::endl;
	/// 读入网格.
	h_tree.readEasyMesh(mesh_file);

	irregular_mesh_p = new IrregularMesh<DIM>;
	/// 产生宏单元网格.
	irregular_mesh_p->reinit(h_tree);
	irregular_mesh_p->globalRefine(G_refine);
	irregular_mesh_p->semiregularize();
	irregular_mesh_p->regularize(false);
	RegularMesh<DIM> &mesh_p = irregular_mesh_p->regularMesh();
	mesh_bak = mesh_p;
	mesh_p.writeEasyMesh("new");
	/// 移动网格读入. 注意：如果更换网格，则要把mesh_file中的.d文件复制给 new.d文件.
	readDomain("new");
	
	irregular_mesh_v = new IrregularMesh<DIM>(*irregular_mesh_p);

	/// 产生单元网格.
	irregular_mesh_v->globalRefine(1);
	irregular_mesh_v->semiregularize();
	irregular_mesh_v->regularize(false);
	RegularMesh<DIM> &mesh_v = irregular_mesh_v->regularMesh();

	std::cout << "v mesh refined." << std::endl;
	unsigned int n_ele_v = mesh_v.n_geometry(DIM);
	unsigned int n_ele_p = mesh_p.n_geometry(DIM);
	index_ele_v2p.resize(n_ele_v);
	index_ele_p2v.resize(n_ele_p);

	IrregularMesh<DIM, DIM>::ActiveIterator v_iterator = irregular_mesh_v->beginActiveElement();
	IrregularMesh<DIM, DIM>::ActiveIterator v_end = irregular_mesh_v->endActiveElement();
	for (; v_iterator != v_end; ++v_iterator)
	{
		HElement<DIM, DIM> *pat = v_iterator->parent;
		int ele_p_idx = pat->index;
		int n_chi = pat->n_child;
		index_ele_p2v[ele_p_idx].resize(n_chi);
		for (unsigned int i = 0; i < n_chi; ++i)
		{
			HElement<DIM, DIM> *chi = pat->child[i];
			int ele_v_idx = chi->index;
			index_ele_p2v[ele_p_idx][i] = ele_v_idx;
			index_ele_v2p[ele_v_idx] = ele_p_idx;
		}
	}

	//  测试代码, 检查如何通过树上的关系找不同网格间编号的对应关系. 但这个不能代替全局的对应关系. 因为树本身是一个链表结构.
//	v_iterator = irregular_mesh_v->beginActiveElement();
//	for (; v_iterator != v_end; ++v_iterator)
//	{
//		std::cout << "index v = " << v_iterator->index << std::endl;
//		std::cout << "parents = " << v_iterator->parent->index << std::endl;
//		std::cout << "compare: " << v_iterator->parent->index - index_ele_v2p[v_iterator->index] << std::endl;
//	}
//	IrregularMesh<DIM, DIM>::ActiveIterator p_iterator = irregular_mesh_p->beginActiveElement();
//	IrregularMesh<DIM, DIM>::ActiveIterator p_end = irregular_mesh_p->endActiveElement();
//	for (; p_iterator != p_end; ++p_iterator)
//	{
//		std::cout << "index p = " << p_iterator->index << std::endl;
//		for (int i = 0; i < p_iterator->n_child; ++i)
//		{
//			std::cout << "child[" << i << "] = " << p_iterator->h_element->child[i]->index << std::endl;
//			std::cout << "compare: " << p_iterator->h_element->child[i]->index - index_ele_p2v[p_iterator->index][i] << std::endl;
//		}
//	}
};

#undef DIM
