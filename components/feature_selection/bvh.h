#ifndef BVH_H
#define BVH_H

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include "triangle_mesh_type.h"
#include "ray_intersection.h"

// typedef FieldTriMesh MeshType;
template<class MeshType>
class BVHT
{
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FacePointer FacePointer;
    typedef typename vcg::Box3<ScalarType> BoxType;
    typedef typename vcg::Point3<ScalarType> PointType;

    typedef  Ray<ScalarType> RayType;
    typedef  Intersection<ScalarType> InterType;

public:
    BVHT() 
    {
        root = nullptr;
    }

    ~BVHT()
    {
        Clear(root);
    }

    void BuildTree(MeshType& mesh)
    {
        mptr = &mesh;
        Clear(root);
        std::vector<FacePointer> fplist;
        for (size_t i = 0; i < mesh.face.size(); i++)
            fplist.push_back(&mesh.face[i]);
        root = BuildTreeRecur(fplist);
    }
    
    bool Intersect(const RayType& ray, InterType& inter)
    {
        inter.happened = false;
        return IntersectRecur(root, ray, inter);
    }

    bool Replace(FacePointer fp, std::vector<FacePointer>fplist)
    {
        if (fp == nullptr)
            return false;

        BVHN** pptr;
        PointType c = GetBBox(fp).Center();
        if (FindRecur(root, fp, c, pptr))
        {
            Clear(*pptr);
            *pptr = BuildTreeRecur(fplist);
            return true;
        }
        return false;
    }

private:

    struct BVHN
    {   
        BoxType bbox;
        BVHN* left;
        BVHN* right;
        FacePointer obj;
        BVHN(): bbox(), left(nullptr), right(nullptr), obj(nullptr) {   bbox.SetNull(); }
    };

    void Clear(BVHN* node)
    {
        if (node != nullptr)
        {
            Clear(node->left);
            Clear(node->right);
            delete node;
        }
    }

    bool FindRecur(BVHN* node, FacePointer fp, PointType c, BVHN**& ret)
    {
        if (node == nullptr)
            return false;

        if (node->bbox.IsIn(c))
        {
            if (node->obj == fp)
            {
                ret = &node;
                return true;
            }
            if (FindRecur(node->left, fp, c, ret))
                return true;
            if (FindRecur(node->right, fp, c, ret))
                return true;
        }
        return false;
    }

    bool Intersect(const BoxType& bbox, const RayType& ray)
    {
        PointType tIns = bbox.min - ray.ori;
        PointType tOuts = bbox.max - ray.ori;
        
        for (size_t i = 0; i < 3; i++)
        {
            tIns[i] *= ray.dir_inv[i];
            tOuts[i] *= ray.dir_inv[i];
        }

        for (size_t i = 0; i < 3; i++) 
        {
            if (ray.dir[i] <= 0) 
            {
                ScalarType tmp = tIns[i];
                tIns[i] = tOuts[i];
                tOuts[i] = tmp;
            }
        }
        ScalarType tIn = std::max(tIns.X(), std::max(tIns.Y(), tIns.Z()));
        ScalarType tOut = std::min(tOuts.X(), std::min(tOuts.Y(), tOuts.Z()));

        if (tIn <= tOut && tOut >= 0)
            return true;
        
        return false;
    }

    bool Intersect(FacePointer fp, const RayType& ray, InterType& inter)
    {
        inter.happened = false;

        PointType a = fp->P(0);
        PointType b = fp->P(1);
        PointType c = fp->P(2);

        PointType e1 = b - a;
        PointType e2 = c - a;

        ScalarType u, v, t_tmp = 0;
        PointType pvec = ray.dir ^ e2;

        ScalarType det = e1 * pvec;

        if (fabs(det) < 1e-10)
            return false;

        ScalarType det_inv = 1. / det;
        PointType tvec = ray.ori - a;
        u = (tvec * pvec) * det_inv;
        if (u < 0 || u > 1)
            return false;

        PointType qvec = tvec ^ e1;
        v = (ray.dir * qvec) * det_inv;
        if (v < 0 || u + v > 1)
            return false;

        t_tmp = (e2 * qvec) * det_inv;
        if (t_tmp < 0)
            return false;

        inter.happened = true;
        inter.pos = ray(t_tmp);
        inter.dist = t_tmp;
        inter.idx = vcg::tri::Index(*mptr, fp);
        inter.bary = PointType(1-u-v, u, v);
        return true;
    }

    BoxType GetBBox(FacePointer fp)
    {
        BoxType box;
        fp->GetBBox(box);
        return box;
    }

    bool IntersectRecur(BVHN* node, const RayType& ray, InterType& inter)
    {
        inter = InterType();

        if (node == nullptr)
            return false;
        
        if (Intersect(node->bbox, ray))
        {
            if (node->obj != nullptr)
                return Intersect(node->obj, ray, inter);

            InterType interl, interr;
            bool happenl = IntersectRecur(node->left, ray, interl);
            bool happenr = IntersectRecur(node->right, ray, interr);
            
            if (happenl && happenr)
                inter = interl.dist < interr.dist ? interl : interr;
            else if (happenl)
                inter = interl;
            else if (happenr)
                inter = interr;
            
            return inter.happened;
        }
        return false;
    }

    BVHN* BuildTreeRecur(std::vector<FacePointer> fplist)
    {
        if( fplist.size() == 0)
            return nullptr;

        BVHN* node = new BVHN;

        if (fplist.size() == 1)
        {
            node->obj = fplist[0];
            node->left= node->right = nullptr;
            node->bbox = GetBBox(fplist[0]);
        }

        else if (fplist.size() == 2){
            node->left = BuildTreeRecur(std::vector<FacePointer>{fplist[0]});
            node->right = BuildTreeRecur(std::vector<FacePointer>{fplist[1]});
            node->bbox.Add(node->left->bbox);
            node->bbox.Add(node->right->bbox);
        }

        else {
            BoxType centerBox;
            for (size_t i = 0; i < fplist.size(); i++)
                centerBox.Add(GetBBox(fplist[i]).Center());
            
            int longestDim = centerBox.DimX() > centerBox.DimY() ? (centerBox.DimX() > centerBox.DimZ() ? 0 : 2) : (centerBox.DimY() > centerBox.DimZ() ? 1 : 2);
            switch (longestDim)
            {
            case 0:
                std::sort(fplist.begin(), fplist.end(), [&](auto f1, auto f2){
                    return GetBBox(f1).Center().X() < GetBBox(f2).Center().X();
                });
                break;
            case 1:
                std::sort(fplist.begin(), fplist.end(), [&](auto f1, auto f2){
                    return GetBBox(f1).Center().Y() < GetBBox(f2).Center().Y();
                });
                break;
            case 2:
                std::sort(fplist.begin(), fplist.end(), [&](auto f1, auto f2){
                    return GetBBox(f1).Center().Z() < GetBBox(f2).Center().Z();
                });
                break;
            }

            auto start = fplist.begin();
            auto mid = fplist.begin() + fplist.size() / 2;
            auto end = fplist.end();

            node->left = BuildTreeRecur(std::vector<FacePointer>(start, mid));
            node->right = BuildTreeRecur(std::vector<FacePointer>(mid, end));

            node->bbox.Add(node->left->bbox);
            node->bbox.Add(node->right->bbox);
        }
        return node;
    }

    BVHN* root;
    MeshType* mptr;
};


#endif