#ifndef GRAPH_SYMMETRY
#define GRAPH_SYMMETRY

#include "vert_field_graph.h"


template<class MeshType>
class VertexFieldSymmetry
{
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename vcg::Point4<ScalarType> PlaneType;

public:

    static ScalarType PairSymmetryEval(VertexFieldGraph<MeshType>& VFGraph,
                                       CoordType& P0, 
                                       CoordType& P1)
    {
        PlaneType SymmPlane = VFGraph.SymmPlane;

        CoordType SP0 = SymmetryManager<MeshType>::SymmetryPoint(P0, SymmPlane);
        
        return Distance(SP0, P1);
    }

    static ScalarType PathSymmetryEval(VertexFieldGraph<MeshType>& VFGraph,
                                                   const std::vector<size_t>& Path0, 
                                                   const std::vector<size_t>& Path1)
    {
        std::vector<size_t> VPath[2];

        // 0 always greater than 1
        int small=1, large=0;
        if (Path0.size() >= Path1.size())
        {
            VFGraph.NodeVertI(Path0, VPath[0]);
            VFGraph.NodeVertI(Path1, VPath[1]);
        }
        else
        {
            VFGraph.NodeVertI(Path0, VPath[1]);
            VFGraph.NodeVertI(Path1, VPath[0]);
        }
        
        // init distance for samller one
        std::vector<ScalarType> DistPath;
        DistPath.resize(VPath[small].size(), 0);
        
        for (size_t i = 1; i < DistPath.size(); i++)
        {
            CoordType P0 = VFGraph.Mesh().vert[VPath[1][i-1]].P();
            CoordType P1 = VFGraph.Mesh().vert[VPath[1][i]].P();
            DistPath[i] = DistPath[i-1] + Distance(P0, P1);
        }

        // partion number
        size_t PosNum = VPath[large].size();
        size_t PartNum = PosNum - 1;

        std::vector<CoordType> PosPath[2];
        PosPath[large].resize(PosNum);
        PosPath[small].resize(PosNum);

        // part for greater one 
        for (size_t i = 0; i < PosNum; i++)
            PosPath[large][i] = VFGraph.Mesh().vert[VPath[large][i]].P();
        
        // part for smaller one
        PosPath[small][0] = VFGraph.Mesh().vert[VPath[small][0]].P();
        PosPath[small].back() = VFGraph.Mesh().vert[VPath[small].back()].P();

        int Ptr = 1;
        for (size_t i = 1; i < PosNum-1; i++)
        {
            ScalarType TargetDist = DistPath.back() * (i / (ScalarType)PartNum);
            for (; Ptr < DistPath.size(); Ptr++)
            {
                ScalarType PreDist = DistPath[Ptr-1];
                ScalarType CurDist = DistPath[Ptr];
                if (PreDist < TargetDist && CurDist >= TargetDist)
                {
                    ScalarType alpha = (TargetDist-PreDist) / (CurDist-PreDist);
                    CoordType PrePos = VFGraph.Mesh().vert[VPath[small][Ptr-1]].P();
                    CoordType CurPos = VFGraph.Mesh().vert[VPath[small][Ptr]].P();
                    PosPath[small][i] = alpha*PrePos + (1-alpha)*CurPos;
                    break;
                }
            }
            assert(Ptr < DistPath.size());
        }

        // check symmetry
        ScalarType SymmVal = 0;
        for (size_t i = 0; i < PosNum; i++)
            SymmVal += PairSymmetryEval(VFGraph, PosPath[small][i], PosPath[large][i]);
        SymmVal /= PosNum;
        
        return SymmVal;
    }


    // smaller when more symmetry
    static ScalarType PathSymmetryEval(VertexFieldGraph<MeshType>& VFGraph,
                                       const std::vector<int>& SymmVert,
                                       const std::vector<size_t>& Path,
                                       bool IsLoop)
    {
        std::vector<size_t> VPath;
        VFGraph.NodeVertI(Path, VPath);

        std::vector<size_t> SymmPtr;
        for (size_t i = 0; i < VPath.size(); i++)
        {
            if (SymmVert[VPath[i]] == VPath[(i+1)%VPath.size()])
                SymmPtr.push_back((i+1)%VPath.size());
        }
        
        // not symmetry type
        if ((SymmPtr.size() != 2 && IsLoop) ||
            (SymmPtr.size() != 1 && !IsLoop)) return std::numeric_limits<ScalarType>::max();
            
        // split path
        std::vector<size_t> SubPath[2];
        if (IsLoop)
        {
            for (size_t i = 0; i < 2; i++)
                for (size_t j = SymmPtr[i]; j != SymmPtr[(i+1)%2]; j=(j+1)%Path.size())
                    SubPath[i].push_back(Path[j]);   
        }
        else
        {
            SubPath[0] = std::vector<size_t>(Path.begin(), Path.begin()+SymmPtr[0]);
            SubPath[1] = std::vector<size_t>(Path.begin()+SymmPtr[0], Path.end());
        }
        
        std::reverse(SubPath[0].begin(), SubPath[0].end());
        return PathSymmetryEval(VFGraph, SubPath[0], SubPath[1]);   
    }

};








#endif