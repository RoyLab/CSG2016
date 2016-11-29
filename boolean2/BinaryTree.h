#pragma once
#include <map>
#include <list>
#include <vector>

#include <macros.h>
#include "global.h"
#include "RegularMesh.h"

namespace Boolean
{

	typedef RegularMesh MPMesh;

	typedef CSGNodeOp BiNodeType;

	struct CSGTreeNode
	{
		BiNodeType Type;
		Relation relation;
		CSGTreeNode *pLeft, *pRight, *Parent;
		uint64_t mark;

		MPMesh* pMesh;
		bool	 bInverse;

		CSGTreeNode();
		~CSGTreeNode();
	};

	struct CSGTreeOld
	{
		CSGTreeNode* pRoot;
		std::map<uint32_t, CSGTreeNode*> Leaves;

		CSGTreeOld();
		~CSGTreeOld();
	};

	struct Branch
	{
		int targetRelation;
		CSGTreeNode* testTree;

		~Branch() { SAFE_DELETE(testTree); }
	};

	typedef std::list<Branch> TestTree;
	void GetLeafList(CSGTreeNode* root, std::vector<uint32_t>& list);
	CSGTreeOld* ConvertToPositiveTree(const CSGTreeOld* tree);
	Relation CompressCSGTree(CSGTreeOld* tree, uint32_t Id, Relation rel);
	Relation ParsingCSGTree(MPMesh::Index meshId, Relation* tab, uint32_t nMesh, 
        CSGTreeNode* curTree, CSGTreeNode** leaves, TestTree& output);
	CSGTreeNode* GetNextNode(CSGTreeNode* curNode, Relation rel, Relation &output);
	CSGTreeNode* GetFirstNode(CSGTreeNode* root);
	inline bool IsLeaf(CSGTreeNode* node) { return !(node->pLeft && node->pRight); }
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

