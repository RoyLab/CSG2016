#pragma once
#include <vector>
#include <map>
#include <cstdint>
#include <boost\smart_ptr.hpp>

#include "macroutil.h"
#include "csg.h"

namespace CSG
{
    template <class Mesh> class TrimCSGTree;

    typedef int8_t                  Indicator;
    typedef std::vector<int>        IndMap;
    typedef std::map<int, int>      IndInvMap;



    // free for copy
    class IndicatorVector
    {
    public:
        IndicatorVector(size_t sz) :codeMap(nullptr), size(sz) { inds.reset(new Indicator[size]); }
        IndicatorVector(IndMap* cm) :codeMap(cm), size(cm->size()) { inds.reset(new Indicator[size]); }
        ~IndicatorVector(){}

        Indicator& operator[](size_t meshId) {
            if (!codeInvMap)
                return inds[meshId];

            auto res = codeInvMap->find(meshId);
#ifdef _DEBUG
            if (res == codeInvMap->end())
                assert(0);
#endif
            return inds[res->second];
        }

        void copy(IndicatorVector& v)
        {
            v.codeInvMap = codeInvMap;
            v.codeMap = codeMap;
            v.size = size;
            v.inds.reset(new Indicator[size]);
            memcpy(v.inds.get(), inds.get(), sizeof(Indicator)*size);
        }

    private:

        boost::shared_array<Indicator> inds;
        int                     size = 0;
        const IndMap *     codeMap; // from meshId to index
        const IndInvMap *  codeInvMap; // inverse of IndMap 
    };

    template <class Mesh>
    class TrimCSGTree
    {
        typedef typename Mesh::Face_handle FH;
    public:
        TrimCSGTree(CSGTree<Mesh>& csg, IndicatorVector& indvec, std::vector<int>& itstPrims);
        TrimCSGTree(TrimCSGTree<Mesh>& csg, IndicatorVector& indvec, FH fh);
        ~TrimCSGTree(){}

        IndicatorVector* createIndicatorVector() const { return new IndicatorVector(2); }
        IndicatorVector* downcast(IndicatorVector& full) const
        {
            IndicatorVector* inds = createIndicatorVector();
            for (int i = 0; i < codeMap.size(); i++)
                inds[i] = full[codeMap[i]];
            return inds;
        }

        bool isSimple() const
        {
            assert(0);
            return true;
        }

    private:
        IndMap      codeMap;
        IndInvMap   codeInvMap;
    };

}
