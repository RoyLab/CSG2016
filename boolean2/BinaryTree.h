#pragma once
#include <map>
#include <list>
#include <vector>

#include <macros.h>
#include "global.h"

namespace Boolean
{
    typedef CSGNodeOp BiNodeType;

	template <class MeshT>
    struct CSGTreeNode
    {
        BiNodeType Type;
		Relation relation;
        CSGTreeNode *pLeft, *pRight, *Parent;
		unsigned long long mark;

		MeshT*	 pMesh;
		bool	 bInverse;

		CSGTreeNode() :
			relation(REL_UNKNOWN),
			Type(TYPE_UNKNOWN), Parent(0),
			pLeft(0), pRight(0), pMesh(0),
			bInverse(false) {}

		~CSGTreeNode()
		{
			SAFE_DELETE(pLeft);
			SAFE_DELETE(pRight);
		}
    };

	template <class MeshT>
    struct CSGTreeOld
    {
		typedef std::map<unsigned, CSGTreeNode*> LeafList;

        CSGTreeNode* pRoot;
		LeafList Leaves;

		CSGTreeOld():pRoot(nullptr) {}
		~CSGTreeOld() { SAFE_DELETE(pRoot); }
    };

	template <class MeshT>
    struct Branch
	{
		int targetRelation;
		CSGTreeNode* testTree;

        ~Branch() {SAFE_DELETE(testTree);}
	};

	template <class MeshT>
    void GetLeafList(CSGTreeNode<MeshT>* root, std::vector<int>& list)
	{
		if (root->Type == TYPE_LEAF)
		{
			list.push_back(root->pMesh->Id);
		}
		else
		{
			GetLeafList(root->pLeft, list);
			GetLeafList(root->pRight, list);
		}
	}

	template <class MeshT>
	void _GetLeafList(CSGTreeNode<MeshT>* root, typename CSGTreeOld<MeshT>::LeafList& leaves)
	{
		if (root->Type == TYPE_LEAF)
		{
			leaves[root->pMesh->Id] = root;
		}
		else
		{
			GetLeafList(root->pLeft, leaves);
			GetLeafList(root->pRight, leaves);
		}
	}

	template <class MeshT>
	CSGTreeNode<MeshT>* _ConvertToPositiveTree(const CSGTreeNode<MeshT>* root, 
		bool inverse, unsigned level, unsigned& maxLvl)
	{
		CSGTreeNode<MeshT>* res = new CSGTreeNode<MeshT>;
		if (root->Type == TYPE_LEAF) // it is a leaf node
		{
			res->bInverse = inverse;
			res->Type = TYPE_LEAF;
			res->pMesh = root->pMesh;
			res->pMesh->bInverse = res->bInverse;
			maxLvl = (maxLvl < level) ? level : maxLvl;
		}
		else
		{
			unsigned Ldepth(0), Rdepth(0);
			if (root->Type == TYPE_DIFF)
			{
				res->Type = TYPE_INTERSECT;
				res->pLeft = _ConvertToPositiveTree(root->pLeft, inverse, level + 1, Ldepth);
				res->pRight = _ConvertToPositiveTree(root->pRight, !inverse, level + 1, Rdepth);
			}
			else
			{
				res->Type = root->Type;
				res->pLeft = _ConvertToPositiveTree(root->pLeft, inverse, level + 1, Ldepth);
				res->pRight = _ConvertToPositiveTree(root->pRight, inverse, level + 1, Rdepth);
			}
			maxLvl = (Ldepth > Rdepth) ? Ldepth : Rdepth;
			if (Ldepth < Rdepth)
			{
				auto tmpNode = res->pLeft;
				res->pLeft = res->pRight;
				res->pRight = tmpNode;
			}

			if (res->pLeft) res->pLeft->Parent = res;
			if (res->pRight) res->pRight->Parent = res;

			if (inverse)
				res->Type = (res->Type == TYPE_INTERSECT) ? TYPE_UNION : TYPE_INTERSECT;
		}
		return res;
	}

	template <class MeshT>
	CSGTreeOld<MeshT>* ConvertToPositiveTree(const CSGTreeOld<MeshT>* tree)
	{
		CSGTreeOld<MeshT>* result = new CSGTreeOld<MeshT>;
		unsigned maxLvl = 0;
		result->pRoot = _ConvertToPositiveTree(myTree->pRoot, false, 0, maxLvl);
		result->Leaves.clear();
		if (!result->pRoot) return;
		_GetLeafList(result->pRoot, result->Leaves);
		return result;
	}

	template <class MeshT>
	Relation CompressCSGTreeWithInside(CSGTreeOld* tree, unsigned Id)
	{
		auto leaf = tree->Leaves[Id];

		CSGTreeNode *curPtr = leaf, *parent = leaf->Parent;
		while (parent && parent->Type == TYPE_UNION)
		{
			curPtr = parent;
			parent = parent->Parent;
		}

		if (!parent)
		{
			delete tree->pRoot;
			tree->Leaves.clear();
			tree->pRoot = nullptr;
			return REL_INSIDE;
		}

		assert(parent->Type == TYPE_INTERSECT);
		int resCur = LeftOrRight(curPtr);
		int resPar = LeftOrRight(parent);

		if (resPar == 0) // root
		{
			if (resCur < 0)
			{
				tree->pRoot = parent->pRight;
				parent->pRight = nullptr;
				delete parent;
			}
			else
			{
				assert(resCur);
				tree->pRoot = parent->pLeft;
				parent->pLeft = nullptr;
				delete parent;
			}
			tree->pRoot->Parent = nullptr;
		}
		else if (resPar < 0)
		{
			if (resCur < 0)
			{
				parent->pRight->Parent = parent->Parent;
				parent->Parent->pLeft = parent->pRight;
				parent->pRight = nullptr;
				delete parent;
			}
			else
			{
				assert(resCur);
				parent->pLeft->Parent = parent->Parent;
				parent->Parent->pLeft = parent->pLeft;
				parent->pLeft = nullptr;
				delete parent;
			}
		}
		else
		{
			if (resCur < 0)
			{
				parent->pRight->Parent = parent->Parent;
				parent->Parent->pRight = parent->pRight;
				parent->pRight = nullptr;
				delete parent;
			}
			else
			{
				assert(resCur);
				parent->pLeft->Parent = parent->Parent;
				parent->Parent->pRight = parent->pLeft;
				parent->pLeft = nullptr;
				delete parent;
			}
		}
		GetLeafList(tree);
		return REL_UNKNOWN;
	}

	template <class MeshT>
	Relation CompressCSGTreeWithOutside(CSGTreeOld* tree, unsigned Id)
	{
		auto leaf = tree->Leaves[Id];
		CSGTreeNode *curPtr = leaf, *parent = leaf->Parent;
		while (parent && parent->Type == TYPE_INTERSECT)
		{
			curPtr = parent;
			parent = parent->Parent;
		}

		if (!parent)
		{
			delete tree->pRoot;
			tree->pRoot = nullptr;
			tree->Leaves.clear();
			return REL_OUTSIDE;
		}

		assert(parent->Type == TYPE_UNION);
		int resCur = LeftOrRight(curPtr);
		int resPar = LeftOrRight(parent);

		if (resPar == 0) // root
		{
			if (resCur < 0)
			{
				tree->pRoot = parent->pRight;
				parent->pRight = nullptr;
				delete parent;
			}
			else
			{
				assert(resCur);
				tree->pRoot = parent->pLeft;
				parent->pLeft = nullptr;
				delete parent;
			}
			tree->pRoot->Parent = nullptr;
		}
		else if (resPar < 0)
		{
			if (resCur < 0)
			{
				parent->pRight->Parent = parent->Parent;
				parent->Parent->pLeft = parent->pRight;
				parent->pRight = nullptr;
				delete parent;
			}
			else
			{
				assert(resCur);
				parent->pLeft->Parent = parent->Parent;
				parent->Parent->pLeft = parent->pLeft;
				parent->pLeft = nullptr;
				delete parent;
			}
		}
		else
		{
			if (resCur < 0)
			{
				parent->pRight->Parent = parent->Parent;
				parent->Parent->pRight = parent->pRight;
				parent->pRight = nullptr;
				delete parent;
			}
			else
			{
				assert(resCur);
				parent->pLeft->Parent = parent->Parent;
				parent->Parent->pRight = parent->pLeft;
				parent->pLeft = nullptr;
				delete parent;
			}
		}
		GetLeafList(tree);
		return REL_UNKNOWN;
	}

	template <class MeshT>
	Relation CompressCSGTreeWithSame(CSGTreeOld* tree, unsigned Id)
	{
		auto leaf = tree->Leaves[Id];
		CSGTreeNode *curPtr = leaf, *parent = leaf->Parent, *neib;

		while (true)
		{
			if (!parent)
			{
				delete tree->pRoot;
				tree->pRoot = nullptr;
				tree->Leaves.clear();
				return REL_SAME;
			}
			if (LeftOrRight(curPtr) < 0.0)
				neib = parent->pRight;
			else neib = parent->pLeft;

			if (neib->relation == REL_SAME)
			{
				curPtr = parent;
				parent = parent->Parent;
			}
			else if (neib->relation == REL_OPPOSITE)
			{
				if (parent->Type == TYPE_UNION)
					return CompressCSGTreeWithInside(tree, Id);
				else if (parent->Type == TYPE_INTERSECT)
					return CompressCSGTreeWithOutside(tree, Id);
				else assert(0);
			}
			else break;
		}
		curPtr->relation = REL_SAME;
		SAFE_DELETE(curPtr->pLeft);
		SAFE_DELETE(curPtr->pRight);
		GetLeafList(tree);
		return REL_NOT_AVAILABLE;
	}

	template <class MeshT>
	Relation CompressCSGTreeWithOppo(CSGTreeOld* tree, unsigned Id)
	{
		auto leaf = tree->Leaves[Id];
		CSGTreeNode *curPtr = leaf, *parent = leaf->Parent, *neib;

		while (true)
		{
			if (!parent)
			{
				delete tree->pRoot;
				tree->pRoot = nullptr;
				tree->Leaves.clear();
				return REL_OPPOSITE;
			}
			if (LeftOrRight(curPtr) < 0.0)
				neib = parent->pRight;
			else neib = parent->pLeft;

			if (neib->relation == REL_OPPOSITE)
			{
				curPtr = parent;
				parent = parent->Parent;
			}
			else if (neib->relation == REL_SAME)
			{
				if (parent->Type == TYPE_UNION)
					return CompressCSGTreeWithInside(tree, Id);
				else if (parent->Type == TYPE_INTERSECT)
					return CompressCSGTreeWithOutside(tree, Id);
				else assert(0);
			}
			else break;
		}
		curPtr->relation = REL_OPPOSITE;
		SAFE_DELETE(curPtr->pLeft);
		SAFE_DELETE(curPtr->pRight);
		GetLeafList(tree);
		return REL_NOT_AVAILABLE;
	}

	template <class MeshT>
	Relation CompressCSGTree(CSGTreeOld<MeshT>* tree, unsigned Id, Relation rel)
	{
		// T∩N⇒N,  F∪N⇒N, F∩N⇒F, and T∪N⇒T
		if (tree->Leaves.find(Id) == tree->Leaves.end())
			return REL_UNKNOWN;

		switch (rel)
		{
		case REL_INSIDE:
			return CompressCSGTreeWithInside(tree, Id);
			break;
		case REL_OUTSIDE:
			return CompressCSGTreeWithOutside(tree, Id);
			break;
		case REL_SAME:
			return CompressCSGTreeWithSame(tree, Id);
			break;
		case REL_OPPOSITE:
			return CompressCSGTreeWithOppo(tree, Id);
			break;
		default:
			assert(0);
			break;
		}

		return REL_UNKNOWN;
	}

	template <class MeshT>
	Relation CompressCSGNodeIteration(CSGTreeNode<MeshT>*& root)
	{
		if (IsLeaf(root)) return root->relation;

		Relation rRight, rLeft;
		rRight = CompressCSGNodeIteration(root->pRight);

		if (root->Type == TYPE_UNION)
		{
			if (rRight == REL_INSIDE)
			{
				delete root;
				root = nullptr;
				return REL_INSIDE;
			}
			else if (rRight == REL_OUTSIDE)
			{
				delete root->pRight;
				auto tmp = root;
				root = root->pLeft;
				root->Parent = tmp->Parent;

				tmp->pLeft = nullptr;
				tmp->pRight = nullptr;
				delete tmp;

				return CompressCSGNodeIteration(root);
			}
		}
#ifdef _DEBUG
		else if (root->Type == TYPE_INTERSECT)
#else
		else
#endif
		{
			if (rRight == REL_OUTSIDE)
			{
				delete root;
				root = nullptr;
				return REL_OUTSIDE;
			}
			else if (rRight == REL_INSIDE)
			{
				delete root->pRight;
				auto tmp = root;
				root = root->pLeft;
				root->Parent = tmp->Parent;

				tmp->pLeft = nullptr;
				tmp->pRight = nullptr;
				delete tmp;

				return CompressCSGNodeIteration(root);
			}
		}
#ifdef _DEBUG
		else assert(0);
#endif

		rLeft = CompressCSGNodeIteration(root->pLeft);
		if (root->Type == TYPE_UNION)
		{
			if (rLeft == REL_INSIDE)
			{
				delete root;
				root = nullptr;
				return REL_INSIDE;
			}
			else if (rLeft == REL_OUTSIDE)
			{
				delete root->pLeft;
				auto tmp = root;
				root = root->pRight;
				root->Parent = tmp->Parent;

				tmp->pLeft = nullptr;
				tmp->pRight = nullptr;
				delete tmp;
				return rRight;
			}
		}
#ifdef _DEBUG
		else if (root->Type == TYPE_INTERSECT)
#else
		else
#endif
		{
			if (rLeft == REL_OUTSIDE)
			{
				delete root;
				root = nullptr;
				return REL_OUTSIDE;
			}
			else if (rLeft == REL_INSIDE)
			{
				delete root->pLeft;
				auto tmp = root;
				root = root->pRight;
				root->Parent = tmp->Parent;

				tmp->pLeft = nullptr;
				tmp->pRight = nullptr;
				delete tmp;
				return rRight;
			}
		}
#ifdef _DEBUG
		else assert(0);
#endif

		// TO DO: if no-one is outside or inside
		// Suppose we do not have on and oppo relation
		return REL_NOT_AVAILABLE;
	}

	template <class MeshT>
	Relation ParsingCSGTree(MeshT* pMesh, Relation* tab, unsigned nMesh, CSGTreeNode<MeshT>* curTree, 
		CSGTreeNode<MeshT>** leaves, std::list<Branch<MeshT>>& output)
	{
		for (unsigned i = 0; i < nMesh; i++)
		{
			if (leaves[i])
				leaves[i]->relation = tab[i];
		}

		CSGTreeNode *seed = leaves[pMesh->Id], *comp;
		int checkRel;
		bool pass = true;
		bool simple = true;
		while (seed->Parent)
		{
			// T∩N⇒N,  F∪N⇒N, F∩N⇒F, and T∪N⇒T
			if (seed->Parent->Type == TYPE_UNION)
				checkRel = REL_OUTSIDE;
			else
#ifdef _DEBUG
				if (seed->Parent->Type == TYPE_INTERSECT)
#endif
					checkRel = REL_INSIDE;
#ifdef _DEBUG
				else assert(0);
#endif
				if (LeftOrRight(seed) < 0)
				{
					comp = seed->Parent->pRight;
					checkRel ^= REL_SAME;
					seed->Parent->pRight = nullptr;
				}
				else
				{
					comp = seed->Parent->pLeft;
					seed->Parent->pLeft = nullptr;
				}
				comp->Parent = nullptr;
				Relation resRel = CompressCSGNodeIteration(comp);
				if (resRel == REL_NOT_AVAILABLE)
				{
					output.emplace_back();
					output.back().targetRelation = checkRel;
					output.back().testTree = comp;
					simple = false;
				}
				else
				{
					SAFE_DELETE(comp);
					if (!(checkRel & resRel))
					{
						pass = false;
						break;
					}
				}

				seed = seed->Parent;
		}

		SAFE_DELETE(curTree);
		if (pass)
		{
			if (!simple) return REL_NOT_AVAILABLE;
			else return REL_SAME;
		}
		output.clear();
		return REL_INSIDE; //也有可能是OutSide
	}

	template <class MeshT>
	CSGTreeNode<MeshT>* GetNextNode(CSGTreeNode<MeshT>* curNode, Relation rel, Relation &output);

	template <class MeshT>
	CSGTreeNode<MeshT>* GetFirstNode(CSGTreeNode<MeshT>* root)
	{
		assert(root);
		if (IsLeaf(root)) return root;
		else	 return GetFirstNode(root->pLeft);
	}

	template <class MeshT>
	inline bool IsLeaf(CSGTreeNode<MeshT>* node) {return !(node->pLeft && node->pRight);}

	template <class MeshT>
	CSGTreeNode<MeshT>* copy(const CSGTreeNode<MeshT>* thiz)
	{
		if (!thiz) return nullptr;

		CSGTreeNode<MeshT>* pRes = new CSGTreeNode<MeshT>(*thiz);
		pRes->pLeft = copy(thiz->pLeft);
		pRes->pRight = copy(thiz->pRight);

		if (pRes->pLeft) pRes->pLeft->Parent = pRes;
		if (pRes->pRight) pRes->pRight->Parent = pRes;

		return pRes;
	}

	template <class MeshT>
	CSGTreeOld<MeshT>* copy(const CSGTreeOld<MeshT>* thiz)
	{
		if (!thiz) return nullptr;

		CSGTreeOld<MeshT>* pCopy = new CSGTreeOld<MeshT>;
		pCopy->pRoot = copy(thiz->pRoot);
		GetLeafList(pCopy);
		return pCopy;
	}

	template <class MeshT>
	CSGTreeNode<MeshT>* copy2(const CSGTreeNode<MeshT>* thiz, CSGTreeNode<MeshT>** leafList
	{
		if (!thiz) return nullptr;

		CSGTreeNode<MeshT>* pRes = new CSGTreeNode<MeshT>(*thiz);

		if (thiz->pMesh) leafList[thiz->pMesh->Id] = pRes;

		pRes->pLeft = copy2(thiz->pLeft, leafList);
		pRes->pRight = copy2(thiz->pRight, leafList);

		if (pRes->pLeft) pRes->pLeft->Parent = pRes;
		if (pRes->pRight) pRes->pRight->Parent = pRes;

		return pRes;
	}


	/** 
	if it is a left child, return negative
	if it is a right child, return positive
	if it is a root, return 0
	*/
	template <class MeshT>
	inline int LeftOrRight(CSGTreeNode<MeshT>* node)
	{
		assert(node);

		if (!node->Parent) return 0;
		if (node->Parent->pLeft == node) return -1;
		if (node->Parent->pRight == node) return 1;

		assert(0);
		return 0;
	}

} // namespace CSG

