#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "ray_intersection.h"
#include "triangle_mesh_type.h"
#include <Eigen/Core>
#include <Eigen/SVD>
#include <vector>
#include <queue>
#include <set>

// typedef FieldTriMesh MeshType;
template<class MeshType>
class SymmetryManager
{
private:
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::CoordType  CoordType;

    typedef typename MeshType::VertexType  VertexType;
    typedef typename MeshType::VertexPointer VertexPointer;

    typedef typename MeshType::FaceType    FaceType;
    typedef typename MeshType::FacePointer FacePointer;

    typedef std::pair<CoordType,CoordType> CoordPair;
    typedef vcg::Point4<ScalarType> PlaneType;
    typedef Intersection<ScalarType> InterType;
    typedef vcg::face::Pos<FaceType> PosType;

    static ScalarType CosOfDir(const CoordType& a, const CoordType& b)
    {
        return (a * b) / (a.Norm() * b.Norm());
    }

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

    static ScalarType PointToPlane(const CoordType& p, const PlaneType& plane, ScalarType deno=-1)
    {
        ScalarType nume = abs(p[0]*plane[0] + p[1]*plane[1] + p[2]*plane[2] + plane[3]);
        if (deno < 0)
            deno = sqrt(plane[0]*plane[0] + plane[1]*plane[1] + plane[2]*plane[2]);
        return nume / deno;
    }

    static ScalarType AvgOffset(const std::set<CoordPair>& symmAxis, const PlaneType& plane)
    {
        ScalarType avg = 0;
        ScalarType deno = sqrt(plane[0]*plane[0] + plane[1]*plane[1] + plane[2]*plane[2]);
        for (auto edge: symmAxis)
        {
            CoordType mid = (edge.first + edge.second) / 2.;
            avg += PointToPlane(mid, plane);
        }
        return avg / symmAxis.size();
    }

    static int EdgeIndex(FacePointer fp, VertexPointer vp0, VertexPointer vp1)
    {
        for (size_t i = 0; i < 3; i++)
        {
            if (fp->V0(i) == vp0 && fp->V1(i) == vp1)
                return i;
            if (fp->V1(i) == vp0 && fp->V0(i) == vp1)
                return i;
        }
        return -1;
    }

    static void SetSymmAxis(FacePointer fp, size_t edge, int* symmbit)
    {
        fp->SetFaceEdgeS(edge);
        fp->SetUserBit(symmbit[edge]);

        if (fp->IsB(edge))
        {
            std::cout << "Warning: marked a border edge as symmetry axis!\n";
            return;
        }

        size_t edgep = fp->FFi(edge);   
        FacePointer fpp = fp->FFp(edge);  
        fpp->SetFaceEdgeS(edgep);
        fpp->SetUserBit(symmbit[edgep]);
    }


public:
    static PlaneType PlaneSVD(const std::set<CoordType>& pointSet)
    {
        assert(pointSet.size() > 3);

        CoordType avg(0, 0, 0);
        for (auto point: pointSet)
            avg += point;
        avg /= pointSet.size();

        Eigen::MatrixXd A(pointSet.size(), 3);
        int cnt = 0;
        for (auto point: pointSet)
        {
            for (size_t i = 0; i < 3; i++)
                A(cnt, i) = point[i] - avg[i];
            cnt++;
        }

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinV);

        PlaneType plane(0, 0, 0, 0);
        for (size_t i = 0; i < 3; i++)
        {
            plane[i] = svd.matrixV()(i, 2);
            plane[3] -= plane[i] * avg[i];
        }

        return plane;
    }

    static bool GetPlaneQuad(MeshType& mesh, std::vector<CoordType>& quad)
    {
        quad.clear();

        CoordType center;
        std::set<CoordType> pointSet;
        
        for (size_t i = 0; i < mesh.face.size(); i++)
        {
            FacePointer fp = &mesh.face[i];
            for (size_t j = 0; j < 3; j++)
            {
                if (fp->IsUserBit(mesh.symmbit[j]))
                {
                    pointSet.insert(fp->P0(j));
                    pointSet.insert(fp->P1(j));
                }
            }
        }

        if (pointSet.size() < 3)
            return false;
        
        PlaneType plane = PlaneSVD(pointSet);

        center.SetZero();
        for (auto point: pointSet)
            center += point;
        center /= pointSet.size();
        
        
        ScalarType Ext = mesh.bbox.Diag();
        CoordType vec1 = mesh.bbox.Dim();
        
        vec1[2] = (-plane[3] - plane[0]*vec1[0] - plane[1]*vec1[1]) / plane[2];
        vec1.normalize();

        CoordType vec2 = CoordType(plane[0], plane[1], plane[2]).normalized() ^ vec1;

        quad.push_back(center + vec1 * Ext);
        quad.push_back(center + vec2 * Ext);
        quad.push_back(center - vec1 * Ext);
        quad.push_back(center - vec2 * Ext);
        return true;
    }

    static void FindSymmetryAxis(MeshType& mesh, PlaneType& plane, bool trans=false)
    {
        // int *mesh.symmbit = mesh.mesh.symmbit;

        // init two sets
        std::set<CoordType> pointSet;
        std::set<CoordPair> pairSet;

        for (size_t i = 0; i < mesh.face.size(); i++)
        {
            FacePointer fp = &mesh.face[i];
            for (size_t j = 0; j < 3; j++)
            {
                if (fp->IsUserBit(mesh.symmbit[j]))
                {
                    pointSet.insert(fp->P0(j));
                    pointSet.insert(fp->P1(j));
                    
                    int idx0 = vcg::tri::Index(mesh, fp->V0(j));
                    int idx1 = vcg::tri::Index(mesh, fp->V1(j));

                    assert(idx0 != idx1); 
                    if (idx0 < idx1)
                        pairSet.insert(CoordPair(fp->P0(j), fp->P1(j)));
                    else
                        pairSet.insert(CoordPair(fp->P1(j), fp->P0(j)));
                }
            }
        }

        // init plane and error
        PlaneType plane0 = PlaneSVD(pointSet);
        ScalarType offset = AvgOffset(pairSet, plane0) * 2;

        std::queue<VertexPointer> vertexQ;
        vcg::tri::UpdateFlags<MeshType>::VertexClearV(mesh);

        for (size_t i = 0; i < mesh.face.size(); i++)
        {
            FacePointer fp = &mesh.face[i];
            for (size_t j = 0; j < 3; j++)
            {
                if (fp->IsUserBit(mesh.symmbit[j]))
                {
                    if (!fp->V0(j)->IsV())
                    {
                        vertexQ.push(fp->V0(j));
                        fp->V0(j)->SetV();
                    }
                    if (!fp->V1(j)->IsV())
                    {
                        vertexQ.push(fp->V1(j));
                        fp->V1(j)->SetV();
                    }
                }
            }
        }
        
        while (vertexQ.size())
        {
            VertexPointer curV = vertexQ.front();
            vertexQ.pop();

            FacePointer fp = curV->VFp();
            PosType startPos(fp, curV);
            PosType curPos = startPos;

            // traverse 1-ring
            do
            {
                CoordType mid = (curPos.V()->P() + curPos.VFlip()->P()) * 0.5;
                // check if less than error
                if (PointToPlane(mid, plane0) < offset)
                {
                    FacePointer fp = curPos.F();
                    VertexPointer vp = curPos.V(), vpp = curPos.VFlip();
                    
                    // get current edge index
                    int edgeIdx = EdgeIndex(fp, vp, vpp);
                    assert(edgeIdx >= 0 && edgeIdx < 3);

                    // set symm axis
                    fp->SetUserBit(mesh.symmbit[edgeIdx]);
                    fp->SetFaceEdgeS(edgeIdx);
                    
                    // insert into set
                    pointSet.insert(vpp->P());
                    int idx0 = vcg::tri::Index(mesh, vp);
                    int idx1 = vcg::tri::Index(mesh, vpp);
                    
                    // push into queue
                    if (!vpp->IsV())
                    {
                        vertexQ.push(vpp);
                        vpp->SetV();
                    }
                }
                
                // next pos
                curPos.FlipE();
                curPos.FlipF();
            
            } while (startPos != curPos);
        }


        // set final plane
        plane = PlaneSVD(pointSet);

        // todo rotate the plane to fit z
        if (trans)
        {

        }

    }

    static void SetUpForSymmetry(MeshType& mesh)
    {
        // int *mesh.symmbit = mesh.mesh.symmbit;
        for (size_t i = 0; i < 3; i++)
        {
            mesh.symmbit[i] = FaceType::NewBitFlag();
            vcg::tri::UpdateFlags<MeshType>::FaceClear(mesh, mesh.symmbit[i]);
        }
    }

    static bool SelectSymmAxis(MeshType& mesh, InterType& inter)
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

        // int *mesh.symmbit = mesh.mesh.symmbit;

        if (fp->IsUserBit(mesh.symmbit[edge]))
        {
            fp->ClearFaceEdgeS(edge);
            fp->ClearUserBit(mesh.symmbit[edge]);
            fpp->ClearFaceEdgeS(edgep);
            fpp->ClearUserBit(mesh.symmbit[edgep]);
        }
        else
        {
            fp->SetFaceEdgeS(edge);
            fp->SetUserBit(mesh.symmbit[edge]);
            fpp->SetFaceEdgeS(edgep);
            fpp->SetUserBit(mesh.symmbit[edgep]);
            std::cout << "DBG\n";
        }

        return true;
    }

    static size_t TraceAxisByDirection(MeshType& mesh, InterType& inter)
    {
        FacePointer fp = &mesh.face[inter.idx];
        if (fp->IsD())
            return 0;
        
        size_t zeroc = OnEdgeOrInternal(mesh, inter);

        if (zeroc != 1)
            return 0;
        
        size_t edge  = inter.edge;
        SetSymmAxis(fp, edge, mesh.symmbit);

        VertexPointer vp0 = fp->V0(edge);
        VertexPointer vp1 = fp->V1(edge);
        CoordType dir0 =  vp0->P() - vp1->P();
        CoordType dir1 = -dir0;

        return 1 + TraceAxisByDirection(vp0, dir0, mesh.symmbit) 
                 + TraceAxisByDirection(vp1, dir1, mesh.symmbit);
    }

    static size_t TraceAxisByDirection(VertexPointer vp, CoordType dir, int* symmbit, int smooth=5)
    {
        if (vp->IsB())
            return 0;

        size_t cnt = 0;
        VertexPointer startVp = vp;
        CoordType* dirQ = new CoordType[smooth];
        for (size_t i = 0; i < smooth; i++)
            dirQ[i].SetZero();
        dirQ[0] = dir;
        int avgCnt = 1, ptr = 1;
        
        // start tracing
        do
        {
            FacePointer fp = vp->VFp();
            PosType startPos(fp, vp);
            PosType curPos = startPos;

            ScalarType maxCos = -1;
            VertexPointer nextVp = nullptr;
            PosType posToSet;

            do
            {
                CoordType curDir = curPos.VFlip()->P() - vp->P();
                ScalarType curCos = CosOfDir(dir, curDir);
                
                if (maxCos <= curCos && curCos > 0)
                {
                    maxCos = curCos;
                    nextVp = curPos.VFlip();
                    posToSet = curPos;
                }
                
                curPos.FlipE();
                curPos.FlipF();
            } while (curPos != startPos);
            

            if (nextVp == nullptr)
            {
                std::cout << "DBG1\n";
                break;
            }

            // set symm axis
            SetSymmAxis(posToSet.F(), posToSet.E(), symmbit);
            dirQ[(ptr++)%smooth] = nextVp->P() - vp->P();
            dir.SetZero();
            for (size_t i = 0; i < smooth; i++)
                dir += dirQ[i];
            dir /= avgCnt;
            avgCnt = std::min(avgCnt+1, smooth);

            vp = nextVp;
            cnt++;

            if (vp == startVp)
                std::cout << "DBG2\n";

            if (vp->IsB())
                std::cout << "DBG3\n"; 

        } while (vp != startVp && !vp->IsB());
        
        delete[] dirQ;
        return cnt;
    }


};

#endif