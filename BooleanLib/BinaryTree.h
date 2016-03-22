#pragma once
#include <map>
#include <list>
#include <vector>

#include "macroutil.h"
#include "csg.h"

namespace CSG
{
    enum BiNodeType
    {
        TYPE_UNKNOWN = 0,
        TYPE_UNION,
        TYPE_INTERSECT,
        TYPE_DIFF,
		TYPE_LEAF
    };

    template <class MPMesh>
    struct CSGTreeNode
    {
        BiNodeType Type;
		Relation relation;
        CSGTreeNode *pLeft, *pRight, *Parent;
		unsigned long long mark;

        MPMesh* pMesh;
		bool	 bInverse;

		CSGTreeNode();
		~CSGTreeNode();
    };

    template <class MPMesh>
    struct CSGTreeOld
    {
        CSGTreeNode<MPMesh>* pRoot;
		std::map<unsigned, CSGTreeNode<MPMesh>*> Leaves;

        CSGTreeOld();
        ~CSGTreeOld();
    };

    template <class MPMesh>
    struct Branch
	{
		int targetRelation;
		CSGTreeNode<MPMesh>* testTree;

        ~Branch() {SAFE_DELETE(testTree);}
	};

	typedef std::list<Branch> TestTree;
    void GetLeafList(CSGTreeNode* root, std::vector<int>& list);
    CSGTreeOld* ConvertCSGTree(CSGTree<MPMesh>* root, MPMesh*** arrMesh, int *nMes); // convert nodes.
    CSGTreeOld* ConvertToPositiveTree(const CSGTreeOld* tree);
    Relation CompressCSGTree(CSGTreeOld* tree, unsigned Id, Relation rel);
	Relation ParsingCSGTree(MPMesh* pMesh, Relation* tab, unsigned nMesh, CSGTreeNode* curTree, CSGTreeNode** leaves, TestTree& output);
    CSGTreeNode* GetNextNode(CSGTreeNode* curNode, Relation rel, Relation &output);
    CSGTreeNode* GetFirstNode(CSGTreeNode* root);
    inline bool IsLeaf(CSGTreeNode* node) {return !(node->pLeft && node->pRight);}
    CSGTreeOld* copy(const CSGTreeOld* thiz);
    CSGTreeNode* copy2(const CSGTreeNode* thiz, CSGTreeNode** leafList);
	Relation CompressCSGNodeIteration(CSGTreeNode*& root);

	/** 
	if it is a left child, return negative
	if it is a right child, return positive
	if it is a root, return 0
	*/
	inline int LeftOrRight(CSGTreeNode* node)
	{
		assert(node);

		if (!node->Parent) return 0;
		if (node->Parent->pLeft == node) return -1;
		if (node->Parent->pRight == node) return 1;

		assert(0);
		return 0;
	}

} // namespace CSG

