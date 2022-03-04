#ifndef ADD_VERTEX_H
#define ADD_VERTEX_H

#include "ray_intersection.h"
#include "bvh.h"
#include <vector>
#include <set>

// typedef FieldTriMesh MeshType; 
template<class MeshType>
class ExtraVertexProcess
{
public:
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FacePointer FacePointer;
    typedef typename MeshType::FaceIterator FaceIterator;
    typedef typename MeshType::VertexPointer VertexPointer;
    typedef typename MeshType::VertexIterator VertexIterator;

    typedef Intersection<ScalarType> InterType;
    typedef std::vector<InterType> InterListType;
    typedef vcg::Point3<ScalarType> CoordType;
    typedef vcg::Point2<ScalarType> TexCoordType;

private:

    static bool ConfirmFaceTopology(MeshType& mesh, size_t topo[][3], VertexPointer* vp, size_t faces)
    {
        return true;
    }

public:

    static size_t OnEdgeOrInternal(MeshType& mesh, InterType& inter)
    {
        // judge on edge / on vertex / internal
        const double EPSILON = 0.1;
        size_t zeroc = 0, zeroid = 0;
        for (size_t i = 0; i < 3; i++)
        {
            if (inter.bary[i] < EPSILON)
            {
                zeroid = i;
                zeroc++;
            }
        }

        // on edge
        if (zeroc == 1)
        {
            inter.edge = (zeroid + 1) % 3;
            // reset bary
            CoordType& bary = inter.bary;
            bary[zeroid] = 0;
            ScalarType sum = bary[0] + bary[1] + bary[2];

            for (size_t m = 0; m < 3; m++)
                bary[m] /= sum;

            FacePointer fp = &mesh.face[inter.idx];
            inter.pos = fp->P(0)*bary[0] + fp->P(1)*bary[1] + fp->P(2)*bary[2];
        }

        return zeroc;
    }

    static bool SelectFeatureEdge(MeshType& mesh, InterType& inter)
    {
        FacePointer fp = &mesh.face[inter.idx];
        if (fp->IsD())
            return false;

        // judge if close to edge
        size_t zeroc = OnEdgeOrInternal(mesh, inter);

        if (zeroc != 1)
            return false;

        size_t edge = inter.edge;
        size_t edgep = fp->FFi(edge);
        FacePointer fpp = fp->FFp(edge);

        if (fp->IsFaceEdgeS(edge))
        {
            fp->ClearFaceEdgeS(edge);
            fpp->ClearFaceEdgeS(edgep);
        }
        else
        {
            fp->SetFaceEdgeS(edge);
            fpp->SetFaceEdgeS(edgep);
        }

        return true;
    }

    static size_t AddExtraVertex(MeshType& mesh, BVHT<MeshType>& bvh, const InterListType& inters)
    {   
        size_t count = 0;

        // topology
        const size_t new_idx[3][3] = {{1,3,0},{1,2,3},{3,2,0}};
        const size_t edge_idx[3][2][3] = {{{0,1,3},{0,3,2}},
                                            {{1,2,3},{1,3,0}},
                                            {{2,0,3},{2,3,1}}};
                                            
        const double EPSILON = 0.1;

        // face history
        std::set<size_t> fp_set;

        // traverse and add vertices
        for (size_t i = 0; i < inters.size(); i++)
        {
            // face...
            size_t idx = inters[i].idx;

            if (fp_set.find(idx) != fp_set.end())
            {
                std::cout << "Failed! Same face...\n";
                continue;
            }

            // not found, simply split face 
            FacePointer fp = &mesh.face[idx];

            // judge if close to edge
            CoordType bary = inters[i].bary;
            size_t zeroc = 0, zeroid = 0;
            for (size_t m = 0; m < 3; m++)
            {
                if (bary[m] < EPSILON)
                {
                    zeroid = m;
                    zeroc++;
                }
            }

            // close to original vertex
            if (zeroc >= 2)
            {
                std::cout << "Collide on vertex!\n";
                continue;
            }

            // insert into history
            fp_set.insert(idx);

            // set coord
            VertexIterator vi = vcg::tri::Allocator<MeshType>::AddVertex(mesh, inters[i].pos);

            // origin vertex
            VertexPointer vp[4] = {fp->V(0),fp->V(1),fp->V(2),&(*vi)};
            
            // close to edge
            if (zeroc == 1)
            {
                std::cout << "Collide on edge!\n";
                
                size_t edge = (zeroid + 1) % 3;

                // reset bary
                bary[zeroid] = 0;
                double sum = bary[0] + bary[1] + bary[2];
                for (size_t m = 0; m < 3; m++)
                    bary[m] = bary[m] / sum;
                vi->P() = (fp->P(0))*bary[0] + (fp->P(1))*bary[1] + (fp->P(2))*bary[2];               

                // face adj
                FacePointer fpp = fp->FFp(edge);
                size_t zeroidp = (fp->FFi(edge) + 2) % 3;
                size_t fpidx = vcg::tri::Index(mesh, fpp);
                VertexPointer vpp[4] = {fpp->V(0),fpp->V(1),fpp->V(2),&(*vi)};

                // what if find again?
                if (fp_set.find(fpidx) != fp_set.end())
                {
                    //todo
                }
                else
                {
                    // inset into history
                    fp_set.insert(fpidx);

                    // cur face
                    FaceIterator fi = vcg::tri::Allocator<MeshType>::AddFaces(mesh, 2);
                    for (size_t m = 0; m < 2; m++)
                        for (size_t n = 0; n < 3; n++)
                            (fi+m)->V(n) = vp[edge_idx[zeroid][m][n]];

                    // opp face
                    FaceIterator fip = vcg::tri::Allocator<MeshType>::AddFaces(mesh, 2);
                    for (size_t m = 0; m < 2; m++)
                        for (size_t n = 0; n < 3; n++)
                            (fip+m)->V(n) = vpp[edge_idx[zeroidp][m][n]];
                    
                    // reset fp
                    fp = &mesh.face[idx];
                    fpp = &mesh.face[fpidx];
                    
                    if (fp->HasWedgeTexCoord())
                    {
                        TexCoordType texs[4] = { fp->WT(0).P(), fp->WT(1).P(), fp->WT(2).P() };
                        texs[3] = bary[0]*texs[0] + bary[1]*texs[1] + bary[2]*texs[2];
                        for (size_t m = 0; m < 2; m++)
                            for (size_t n = 0; n < 3; n++)
                                (fi+m)->WT(n).P() = texs[edge_idx[zeroid][m][n]];

                        TexCoordType texsp[4] = { fpp->WT(0).P(), fpp->WT(1).P(), fpp->WT(2).P(), texs[3] };
                        for (size_t m = 0; m < 2; m++)
                            for (size_t n = 0; n < 3; n++)
                                (fip+m)->WT(n).P() = texsp[edge_idx[zeroidp][m][n]];
                    }
                    
                }
            }
            
            // internal
            else
            {
                // add 3 faces, set vertices
                FaceIterator fi = vcg::tri::Allocator<MeshType>::AddFaces(mesh, 3);
                for (size_t m = 0; m < 3; m++)
                    for (size_t n = 0; n < 3; n++)
                        (fi+m)->V(n) = vp[new_idx[m][n]];
                
                // reset fp
                fp = &mesh.face[idx];

                // set texcoord
                if (fp->HasWedgeTexCoord())
                {
                    TexCoordType texs[4] = { fp->WT(0).P(), fp->WT(1).P(), fp->WT(2).P() };
                    texs[3] = bary[0]*texs[0] + bary[1]*texs[1] + bary[2]*texs[2];
                        
                    for (size_t m = 0; m < 3; m++)
                        for (size_t n = 0; n < 3; n++)
                            (fi+m)->WT(n).P() = texs[new_idx[m][n]];
                }
            }

            // how many vertices added
            count++;
        }

        // update mesh
        for (auto idx: fp_set)
            vcg::tri::Allocator<MeshType>::DeleteFace(mesh, mesh.face[idx]);
        mesh.UpdateDataStructures();
        mesh.InitEdgeType();

        // update bvh
        bvh.BuildTree(mesh);

        return count;
    }
};












#endif
