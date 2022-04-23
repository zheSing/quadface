#ifndef ADAPTIVE_EVAL_H
#define ADAPTIVE_EVAL_H

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/simplex/face/topology.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <opencv2/opencv.hpp>
#include <wrap/io_trimesh/export_field.h>
#include <iostream>
#include <fstream>
#include <vcg/complex/algorithms/attribute_seam.h>
#include <vcg/complex/algorithms/crease_cut.h>
#include <vcg/complex/algorithms/polygonal_algorithms.h>
#include <vcg/complex/algorithms/point_sampling.h>
#include <vcg/complex/algorithms/voronoi_processing.h>
#include "triangle_mesh_type.h"

typedef FieldTriMesh MeshType;
// template<typename MeshType>
class AdaptProcess
{
public:
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::CoordType  CoordType;

    typedef typename MeshType::VertexType    VertexType;
    typedef typename MeshType::VertexPointer VertexPointer;

    typedef typename MeshType::FaceType    FaceType;
    typedef typename MeshType::FacePointer FacePointer;

    typedef std::pair<CoordType,CoordType> CoordPair;
    typedef vcg::face::Pos<FaceType> PosType;
    typedef vcg::face::JumpingPos<FaceType> JpPosType; // for the border

    // average edge length
    static ScalarType AvgEdge(MeshType& mesh)
    {
        ScalarType ret = 0;
        size_t Num = 0;
        std::set<CoordPair> EdgeSet;

        for (size_t i = 0; i < mesh.face.size(); i++)
        {
            FacePointer fp = &mesh.face[i];
            for (size_t j = 0; j < fp->VN(); j++)
            {
                CoordType P0 = fp->P0(j);
                CoordType P1 = fp->P1(j);
                CoordPair CurEdge(std::min(P0,P1), std::max(P0,P1));

                if (EdgeSet.find(CurEdge) != EdgeSet.end()) continue;

                ret += (P0-P1).Norm();
                Num++;
            }
        }
            
        return ret / Num;
    }

    // get adj vertices of a vertex
    static void GetVVStar(VertexPointer vp, std::vector<VertexPointer>& VV)
    {
        if (vp == nullptr) return;

        FacePointer fp = vp->VFp();
        if (fp == nullptr) return;

        VV.clear();
        JpPosType startPos(fp, vp);
        JpPosType curPos = startPos;

        // 1-ring
        do
        {
            VV.push_back(curPos.VFlip());
            curPos.NextFE();
        } while (startPos != curPos);
    }

    // get adj faces of a vertex
    static void GetVFStar(VertexPointer vp, std::set<FacePointer>& VF)
    {
        if (vp == nullptr) return;

        FacePointer fp = vp->VFp();
        if (fp == nullptr) return;

        VF.clear();
        JpPosType startPos(fp, vp);
        JpPosType curPos = startPos;

        // 1-ring
        do
        {
            VF.insert(curPos.F());
            curPos.NextFE();
        } while (startPos != curPos);
    }

    static void EvalAdaptiveness(MeshType& mesh, ScalarType scale=3)
    {
        // update normal
        mesh.UpdateDataStructures();

        // compute tolerance
        ScalarType Tau = AvgEdge(mesh)*scale;

        // traverse every vertex
        for (size_t i = 0; i < mesh.vert.size(); i++)
        {
            // preparation
            vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);
            VertexPointer vp = &mesh.vert[i];
            vp->SetV();

            // quque
            std::queue<VertexPointer> VertexQ;
            VertexQ.push(vp);

            // normal list
            std::vector<CoordType> Normals;

            while (VertexQ.size())
            {
                VertexPointer curV = VertexQ.front();
                Normals.push_back(curV->N());
                VertexQ.pop();

                // get 1-ring
                std::vector<VertexPointer> vvStar;
                GetVVStar(curV, vvStar);

                // add vertex
                for (size_t j = 0; j < vvStar.size(); j++)
                {
                    if (vvStar[j]->IsV()) continue;
                    if ((vvStar[j]->P()-vp->P()).Norm() > Tau) continue;

                    vvStar[j]->SetV();
                    VertexQ.push(vvStar[j]);
                }
            }

            // normal distribution by variance
            
        }
    }

    static void EvalAdaptiveness(MeshType& mesh,
                                 ScalarType refScale,
                                 ScalarType minScale,
                                 ScalarType maxScale,
                                 size_t step=8,
                                 ScalarType maxScaleDiag=0.2,
                                 bool colorDbg=false)
    {
        // update normal
        mesh.UpdateDataStructures();
        // compute face roughness
        std::vector<ScalarType> FaceRoughness;
        EvalFaceRoughness(mesh, FaceRoughness);
        // sample vertices
        ScalarType avgE = AvgEdge(mesh);
        ScalarType radius = avgE * refScale;
        SampleVerticesPoisson(mesh, radius);

        // compute tolerance
        size_t count = 0;
        ScalarType Tau = 0;
        for (size_t i = 0; i < mesh.vert.size(); i++)
        {
            VertexPointer vp = &mesh.vert[i];
            if (!vp->IsS()) continue;

            count++;
            Tau += EvalVertexRoughness(mesh, vp, radius, FaceRoughness);
        }
        assert(count);
        Tau /= (ScalarType)count;
        vcg::tri::UpdateFlags<MeshType>::VertexClearS(mesh);

        // set radius sequence
        assert(step > 1);
        assert(minScale < maxScale);
        assert(refScale < maxScale);
        std::vector<ScalarType> sampleSequence(step);
        ScalarType minRadius = avgE * minScale;
        ScalarType maxRadius = std::min(avgE*maxScale, mesh.bbox.Diag()*maxScaleDiag);
        ScalarType areaStep = (maxRadius*maxRadius - minRadius*minRadius) / (step-1);
        for (size_t i = 0; i < step; i++)
        {
            ScalarType curArea = minRadius*minRadius + i*areaStep;
            sampleSequence[i] = sqrt(curArea);
        }

        // compute max tolerance area for each vertex
        assert(vcg::tri::HasPerVertexQuality(mesh));
        for (size_t i = 0; i < mesh.vert.size(); i++)
        {
            VertexPointer vp = &mesh.vert[i];
            vp->Q() = EvalAdaptiveness(mesh, vp, Tau, sampleSequence);
        }

        // set adaptiveness value as color for debugging
        if (colorDbg)
        {
            ScalarType minAdpt = std::numeric_limits<ScalarType>::max();
            ScalarType maxAdpt = 0;
            for (size_t i = 0; i < mesh.vert.size(); i++)
            {
                ScalarType curAdpt = mesh.vert[i].Q();
                if (curAdpt < minAdpt) minAdpt = curAdpt;
                if (curAdpt > maxAdpt) maxAdpt = curAdpt;
            }
            ScalarType factor = (ScalarType)std::numeric_limits<u_char>::max() / (maxAdpt - minAdpt);
            for (size_t i = 0; i < mesh.vert.size(); i++)
            {
                mesh.vert[i].C() = (u_char)((mesh.vert[i].Q()-minAdpt) * factor);
            }
            
        }
    }

    // for given radius list, return max area of roughness that less than Tau
    static ScalarType EvalAdaptiveness(MeshType& mesh,
                                       VertexPointer vp, 
                                       ScalarType Tau,
                                       std::vector<ScalarType>& sampleSeq)
    {
        
    }

    static ScalarType EvalVertexRoughness(MeshType& mesh,
                                          VertexPointer vp,
                                          ScalarType radius,
                                          const std::vector<ScalarType>& FaceRoughness)
    {
        std::set<FacePointer> Faces1ring;
        GetVFStar(vp, Faces1ring);

        // BFS
        std::set<FacePointer> SphereCoverFaces;
        VertexBFS(mesh, vp, radius, SphereCoverFaces);

        // 1-ring must be included
        for (auto fp: Faces1ring) 
            SphereCoverFaces.insert(fp);
        
        // compute vertex roughness
        ScalarType area2 = 0;
        ScalarType roughness = 0;
        for (auto fp: SphereCoverFaces)
        {
            ScalarType triArea2;
            size_t indexF = vcg::tri::Index(mesh, fp);

            if (Faces1ring.find(fp) != Faces1ring.end())
            {
                size_t indexFV = FaceIndexV(fp, vp);
                assert(indexFV >=0 && indexFV <= 2);
                CoordType P0 = vp->P();
                CoordType E1 = fp->P1(indexFV) - P0;
                CoordType E2 = fp->P2(indexFV) - P0;

                if (E1.Norm() > radius) E1 = E1.normalized() * radius;
                if (E2.Norm() > radius) E2 = E2.normalized() * radius;

                triArea2 = (E1^E2).Norm();
            }
            else
            {
                triArea2 += vcg::DoubleArea(*fp);
            }

            roughness += triArea2*FaceRoughness[indexF];
        }
        assert(area2 > 0);

        return roughness / area2;
    }

    // vertex index inside a face
    static int FaceIndexV(FacePointer fp, VertexPointer vp)
    {
        for (size_t i = 0; i < 3; i++)
        {
            if (fp->V(i) == vp) return i;
        }
        return -1;
    }


    // bfs start from a vertex, inside a sphere of radius
    static ScalarType VertexBFS(MeshType& mesh, VertexPointer vp, ScalarType radius, std::set<FacePointer>& CoverFaces)
    {
        // preparation
        CoverFaces.clear();
        std::set<FacePointer> CandidateFaces;
        vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);

        // quque
        std::queue<VertexPointer> VertexQ;
        VertexQ.push(vp);
        vp->SetV();

        // BFS
        while (VertexQ.size())
        {
            VertexPointer curV = VertexQ.front();
            VertexQ.pop();

            // get 1-ring
            FacePointer fp = curV->VFp();
            assert(fp != nullptr);
            JpPosType startPos(fp, curV);
            JpPosType curPos = startPos;
            do
            {
                VertexPointer adjV = curPos.VFlip();
                if (adjV->IsV()) continue;
                if ((adjV->P()-vp->P()).Norm() > radius) continue;

                // for covering faces
                adjV->SetV();
                VertexQ.push(adjV);
                CandidateFaces.insert(curPos.F());
                CandidateFaces.insert(curPos.FFlip());

                // next junmp pos
                curPos.NextFE();
            } while (startPos != curPos);   
        }

        // cover faces
        for (auto fp: CandidateFaces)
        {
            if (!fp->V(0)->IsV() || !fp->V(1)->IsV() || !fp->V(2)->IsV()) continue;
            CoverFaces.insert(fp);
        }
    }

    /*
    // give a point, calculate normal variance inside circle of radi
    static void EvalVariance(MeshType& mesh,
                             VertexPointer vp,
                             ScalarType Radi,
                             ScalarType& Var,
                             ScalarType& Area)
    {
        // preparation
        vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);
        vp->SetV();
        // quque
        std::queue<VertexPointer> VertexQ;
        VertexQ.push(vp);
        // normal list
        std::vector<CoordType> Normals;
        while (VertexQ.size())
        {
            VertexPointer curV = VertexQ.front();
            Normals.push_back(curV->N());
            VertexQ.pop();
            // get 1-ring
            std::vector<VertexPointer> vvStar;
            GetVVStar(curV, vvStar);
            // add vertex
            for (size_t j = 0; j < vvStar.size(); j++)
            {
                if (vvStar[j]->IsV()) continue;
                if ((vvStar[j]->P()-vp->P()).Norm() > Radi) continue;

                vvStar[j]->SetV();
                VertexQ.push(vvStar[j]);
            }
        }
        // calculate area
        Area = 0;
        for (size_t i = 0; i < mesh.face.size(); i++)
        {
            FacePointer fp = &mesh.face[i];
            if (fp->V(0)->IsV() && fp->V(1)->IsV() && fp->V(2)->IsV())
                Area += vcg::DoubleArea(*fp);
        }
        // calculate variance
        Var = 0;
        for (size_t i = 0; i < Normals.size(); i++)
        {   
        }
    }*/

    static void EvalFaceRoughness(MeshType& mesh, std::vector<ScalarType>& Roughness)
    {
        typedef std::pair<ScalarType,ScalarType> AvgVarType;

        mesh.UpdateDataStructures();

        Roughness.clear();
        Roughness.resize(mesh.face.size());
        std::vector<AvgVarType> Vertex1ring(mesh.vert.size());
        
        // vertex roughness in 1-ring
        for (size_t i = 0; i < mesh.vert.size(); i++)
        {
            VertexPointer vp = &mesh.vert[i];
            FacePointer fp = vp->VFp();
            assert(fp != nullptr);
            std::vector<ScalarType> Dials;

            JpPosType startPos(fp, vp);
            JpPosType curPos = startPos;

            // 1-ring
            do
            {
                FacePointer curFp = curPos.F();
                FacePointer nextFp = curPos.FFlip();
                // not border, push dial angle
                if (curFp != nextFp) Dials.push_back(acos(curFp->N()*nextFp->N()));
                // next pos (jump border)
                curPos.NextFE();
            } while (startPos != curPos);

            // 
            assert(Dials.size());

            // calculate average and variance of dial angle
            ScalarType avg = 0, var = 0;
            for (size_t j = 0; j < Dials.size(); j++)
                avg += Dials[j];
            avg /= Dials.size();
            for (size_t j = 0; j < Dials.size(); j++)
                var += (Dials[j]-avg) * (Dials[j]-avg);
            var /= Dials.size();
            Vertex1ring[i] = AvgVarType(avg, var);
        }

        // face roughness
        for (size_t i = 0; i < mesh.face.size(); i++)
        {
            FacePointer fp = &mesh.face[i];
            ScalarType nume = 0;
            ScalarType deno = 0;
            for (size_t j = 0; j < 3; j++)
            {
                size_t indexV = vcg::tri::Index(mesh, fp->V(j));
                nume += Vertex1ring[indexV].first * Vertex1ring[indexV].second;
                deno += Vertex1ring[indexV].second;
            }
            if (deno < 1e-5) Roughness[i] = 0;
            else Roughness[i] = nume / deno;
        }
    }


    static void SampleVerticesPoisson(MeshType& mesh, ScalarType radius)
    {   
        std::vector<CoordType> points;
        
        vcg::tri::PoissonSampling<MeshType>(mesh, points, 0, radius, 1, 0.04f, 276519752);

        std::vector<VertexPointer> vps;
        vcg::tri::VoronoiProcessing<MeshType>::SeedToVertexConversion(mesh, points, vps);

        vcg::tri::UpdateFlags<MeshType>::VertexClearS(mesh);
        for (size_t i = 0; i < vps.size(); i++)
        {
            size_t indexV = vcg::tri::Index(mesh, vps[i]);
            mesh.vert[indexV].SetS();
        }
    }
    

};









#endif //!ADAPTIVE_EVAL_H