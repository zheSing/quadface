#include "triangle_mesh_type.h"
#include <opencv2/opencv.hpp>

typedef FieldTriMesh MeshType;

class TextureProcess
{
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::CoordType CoordType;

    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::VertexPointer VertexPointer;

    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::FacePointer FacePointer;
    typedef typename FaceType::TexCoordType TexCoordType;
    typedef typename FaceType::ColorType ColorType;

    static void BottomFlatTriangle(TexCoordType *t,
                                   vcg::Point3i &sum,
                                   size_t &count,
                                   cv::Mat &img,
                                   bool debugmsg = false)
    {
        if (debugmsg)
        {
            std::cout << "In bottom flat triangle:\n";
            std::cout << "Texture Coordinate after sorting: \n";
            for (size_t i = 0; i < 3; i++)
            {
                std::cout << t[i].U() << " " << t[i].V() << std::endl;
            }
        }

        size_t width = img.size().width;
        size_t height = img.size().height;

        size_t start = height * t[1].V() + 0.5;
        size_t end = height * t[2].V() + 0.5;
        ScalarType ratio = width / (ScalarType)height;

        // inverse slope
        ScalarType s0 = ratio * (t[2].U() - t[0].U()) / (t[2].V() - t[0].V());
        ScalarType s1 = ratio * (t[2].U() - t[1].U()) / (t[2].V() - t[1].V());
        if (s0 < s1)
            std::swap(s0, s1);

        ScalarType x0 = width * t[2].U() + 0.5;
        ScalarType x1 = width * t[2].U() + 0.5;

        if (debugmsg)
        {
            std::cout << "Begin drawing...\n";
            std::cout << "start: " << start << ", end: " << end << "\n";
            std::cout << "slope0: " << s0 << ", slope1: " << s1 << "\n";
        }

        //           t2     --------end
        //          /  \
        //         /    \ 
        //       t0-----t1  -------start
        for (size_t i = end; i >= start; i--)
        {
            if (debugmsg)
            {
                std::cout << "Iter" << end - i << ": y=" << i
                          << ", from=" << (int)x0 << " to=" << (int)x1 << "\n";
            }
            count += (int)x1 - (int)x0 + 1;
            for (size_t j = (int)x0; j <= (int)x1; j++)
            {
                if (i < 0 || j < 0 || i > height || j > width)
                {
                    std::cout << "out of image bound\n";
                    continue;
                }
                cv::Vec3b tmp = img.at<cv::Vec3b>(i, j);
                for (size_t k = 0; k < 3; k++)
                {
                    sum[k] += (int)tmp[k];
                }
                count++;
            }
            x0 -= s0;
            x1 -= s1;
        }
    }

    static void TopFlatTriangle(TexCoordType *t,
                                vcg::Point3i &sum,
                                size_t &count,
                                cv::Mat &img,
                                bool debugmsg = false)
    {
        if (debugmsg)
        {
            std::cout << "In top flat triangle:\n";
            std::cout << "Texture Coordinate after sorting: \n";
            for (size_t i = 0; i < 3; i++)
            {
                std::cout << t[i].U() << " " << t[i].V() << std::endl;
            }
        }
        size_t width = img.size().width;
        size_t height = img.size().height;

        size_t start = height * t[0].V() + 0.5;
        size_t end = height * t[2].V() + 0.5;
        ScalarType ratio = width / (ScalarType)height;

        // inverse slope
        ScalarType s1 = ratio * (t[0].U() - t[1].U()) / (t[0].V() - t[1].V());
        ScalarType s2 = ratio * (t[0].U() - t[2].U()) / (t[0].V() - t[2].V());
        if (s1 > s2)
            std::swap(s1, s2);

        ScalarType x1 = width * t[0].U() + 0.5;
        ScalarType x2 = width * t[0].U() + 0.5;

        if (debugmsg)
        {
            std::cout << "Begin drawing...\n";
            std::cout << "start: " << start << ", end: " << end << "\n";
            std::cout << "slope1: " << s1 << ", slope2: " << s2 << "\n";
        }

        //    t1-----t2   --------end
        //     \    /
        //      \  /
        //       t0       -------start
        for (size_t i = start; i <= end; i++)
        {
            if (debugmsg)
            {
                std::cout << "Iter" << i - start << ": y=" << i
                          << ", from=" << (int)x1 << " to=" << (int)x2 << "\n";
            }
            count += (int)x2 - (int)x1 + 1;
            for (size_t j = (int)x1; j <= (int)x2; j++)
            {
                if (i < 0 || j < 0 || i > height || j > width)
                {
                    std::cout << "out of image bound\n";
                    continue;
                }
                cv::Vec3b tmp = img.at<cv::Vec3b>(i, j);
                for (size_t k = 0; k < 3; k++)
                {
                    sum[k] += (int)tmp[k];
                }
                count++;
            }
            x1 += s1;
            x2 += s2;
        }
    }

    static ColorType AverageFaceColor(FacePointer F,
                                      cv::Mat &img,
                                      bool debugmsg = false)
    {
        TexCoordType t[4];
        for (size_t i = 0; i < 3; i++)
        {
            t[i] = F->WT(i);
        }

        // sort
        if (t[0].V() > t[1].V())
            std::swap(t[0], t[1]);
        if (t[0].V() > t[2].V())
            std::swap(t[0], t[2]);
        if (t[1].V() > t[2].V())
            std::swap(t[1], t[2]);

        if (debugmsg)
        {
            std::cout << "Texture Coordinate after sorting: \n";
            for (size_t i = 0; i < 3; i++)
            {
                std::cout << t[i].U() << " " << t[i].V() << std::endl;
            }
        }

        // discretization
        size_t height = img.size().height;
        size_t y[3];
        for (size_t i = 0; i < 3; i++)
        {
            y[i] = height * t[i].V() + 0.5;
        }

        // return
        vcg::Point3i sum = vcg::Point3i::Zero();
        size_t count = 0;

        // bottom flat
        if (y[0] == y[1])
        {
            BottomFlatTriangle(t, sum, count, img, debugmsg);
        }

        // top flat
        else if (y[1] == y[2])
        {
            TopFlatTriangle(t, sum, count, img, debugmsg);
        }

        // universal
        else
        {
            ScalarType mid = (t[0].V() - t[1].V()) / (t[0].V() - t[2].V()) * (t[2].U() - t[0].U()) + t[0].U();
            t[3] = TexCoordType(mid, t[1].V());
            std::swap(t[2], t[3]);
            // bottom flat
            BottomFlatTriangle(t + 1, sum, count, img, debugmsg);
            // top flat
            TopFlatTriangle(t, sum, count, img, debugmsg);
        }

        if (count == 0)
        {
            std::cout << "div by 0\n";
        }

        ColorType ret = ColorType::Black;
        for (size_t i = 0; i < 3; i++)
        {
            ret[i] = sum[i] / count;
        }

        return ret;
    }

    // area of triangle
    static ScalarType FaceArea(FacePointer F)
    {
        CoordType A = F->V(0)->P();
        CoordType B = F->V(1)->P();
        CoordType C = F->V(2)->P();
        CoordType AB = B - A;
        CoordType AC = C - A;
        return 0.5 * (AB ^ AC).Norm();
    }

    struct WeightedEdge
    {
        size_t fIdx, fOpIdx, eIdx;
        ScalarType diff;
        bool deleted;

        WeightedEdge(size_t f, size_t fop, size_t e, ScalarType d) : fIdx(f), fOpIdx(fop), eIdx(e), diff(d), deleted(false) {}
    };

    struct FindSet
    {
    private:
        std::vector<size_t> parent;
        std::vector<size_t> size;
        std::vector<size_t> actsize;
        std::vector<size_t> rank;
        std::vector<ScalarType> thr;
        size_t count;

        size_t Find(size_t i)
        {
            return i == parent[i] ? i : parent[i] = Find(parent[i]);
        }

        ScalarType tau(size_t s, ScalarType k)
        {
            return k / s;
        }

        void Union(size_t i, size_t j, ScalarType diff, ScalarType k)
        {
            ScalarType noDiff = 1.5;
            size_t ri = Find(i);
            size_t rj = Find(j);
            if (rank[ri] < rank[rj])
                std::swap(ri, rj);
            if (rank[ri] == rank[rj])
                rank[ri]++;
            parent[rj] = ri;
            actsize[ri] += actsize[rj];
            if (diff > noDiff)
            {
                size[ri] += size[rj];
                thr[ri] = diff + tau(size[ri], k);
            }
        }

        void Union(size_t i, size_t j)
        {
            size_t ri = Find(i);
            size_t rj = Find(j);
            if (rank[ri] < rank[rj])
                std::swap(ri, rj);
            if (rank[ri] == rank[rj])
                rank[ri]++;
            parent[rj] = ri;
            actsize[ri] += actsize[rj];
        }

    public:
        FindSet(size_t c) : count(c)
        {
            parent = std::vector<size_t>(count);
            rank = std::vector<size_t>(count, 0);
            size = std::vector<size_t>(count, 1);
            actsize = std::vector<size_t>(count, 1);
            std::iota(parent.begin(), parent.end(), 0);
        }

        void Segment(std::vector<WeightedEdge> &edges, ScalarType k)
        {
            thr = std::vector<ScalarType>(count, tau(1, k));
            for (auto &edge : edges)
            {
                size_t fIdx = edge.fIdx;
                size_t fOpIdx = edge.fOpIdx;
                size_t fRoot = Find(fIdx);
                size_t fOpRoot = Find(fOpIdx);
                if (fRoot == fOpRoot)
                    edge.deleted = true;
                else if (edge.diff <= std::min(thr[fRoot], thr[fOpRoot]))
                {
                    Union(fRoot, fOpRoot, edge.diff, k);
                    edge.deleted = true;
                }
            }
        }

        bool MergeSmallSet(std::vector<WeightedEdge> &edges, size_t minSize)
        {
            bool flag = false;
            for (auto &edge : edges)
            {
                size_t fIdx = edge.fIdx;
                size_t fOpIdx = edge.fOpIdx;
                size_t fRoot = Find(fIdx);
                size_t fOpRoot = Find(fOpIdx);
                if (fRoot == fOpRoot)
                    edge.deleted = true;
                else if (std::min(actsize[fRoot], actsize[fOpRoot]) < minSize)
                {
                    Union(fRoot, fOpRoot);
                    edge.deleted = true;
                    flag = true;
                }
            }
            return flag;
        }
    };

    static void CreateGraph(MeshType &mesh, std::vector<WeightedEdge> &edges, bool debugmsg = false)
    {
        for (size_t i = 0; i < mesh.face.size(); i++)
        {
            FacePointer fp = &mesh.face[i];
            for (size_t j = 0; j < 3; j++)
            {
                if (vcg::face::IsBorder(*fp, j))
                {
                    fp->SetFaceEdgeS(j);
                    continue;
                }
                size_t indexV0 = vcg::tri::Index(mesh, fp->V0(j));
                size_t indexV1 = vcg::tri::Index(mesh, fp->V1(j));
                if (indexV0 >= indexV1)
                    continue;
                // calculate difference
                FacePointer fopp = fp->FFp(j);
                vcg::Point3<ScalarType> diff;
                for (size_t k = 0; k < 3; k++)
                    diff[k] = (int)fp->C()[k] - (int)fopp->C()[k];
                // append a edge
                edges.push_back(WeightedEdge(i, vcg::tri::Index(mesh, fopp), j, diff.Norm()));
            }
        }
    }

public:
    static void SetUpAvgColor(MeshType &mesh, cv::Mat &img, ScalarType &thereshold)
    {
        for (size_t i = 0; i < mesh.face.size(); i++)
        {
            // todo: validation judge
            ColorType color;
            color = AverageFaceColor(&mesh.face[i], img);
            mesh.face[i].C() = color;
        }
        thereshold = 10;
    }

    static bool SegmentTexture(MeshType &mesh,
                               cv::Mat &img,
                               ScalarType k = 100,
                               bool debugmsg = false)
    {
        if (img.empty())
            return false;

        ScalarType thereshold;
        SetUpAvgColor(mesh, img, thereshold);

        mesh.UpdateDataStructures();
        vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(mesh);

        std::vector<WeightedEdge> edges;
        CreateGraph(mesh, edges);
        std::sort(edges.begin(), edges.end(), [](const WeightedEdge &a, const WeightedEdge &b)
                  { return a.diff < b.diff; });

        FindSet forest(mesh.face.size());
        forest.Segment(edges, k);
        while (forest.MergeSmallSet(edges, 20))
        {
        };

        for (auto edge : edges)
        {
            if (!edge.deleted)
            {
                FacePointer fp = &mesh.face[edge.fIdx];
                FacePointer fopp = &mesh.face[edge.fOpIdx];
                size_t eOpIdx = fp->FFi(edge.eIdx);
                fp->SetFaceEdgeS(edge.eIdx);
                fopp->SetFaceEdgeS(eOpIdx);
            }
        }
        mesh.InitEdgeType();
        mesh.SaveSharpFeatures("segment.sharp");
        return true;
    }

    static bool ExtractTexFeature(MeshType &mesh,
                                  cv::Mat &img,
                                  size_t erodilastep = 0,
                                  ScalarType alpha = 0.5,
                                  bool debugmsg = false)
    {
        if (img.empty())
            return false;

        ScalarType thereshold;
        SetUpAvgColor(mesh, img, thereshold);

        mesh.UpdateDataStructures();
        vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(mesh);

        // loop over edges
        for (size_t i = 0; i < mesh.face.size(); i++)
        {
            FacePointer fp = &mesh.face[i];
            for (size_t j = 0; j < 3; j++)
            {
                if (fp->IsFaceEdgeS(j))
                    continue;
                if (vcg::face::IsBorder(*fp, j))
                    continue;
                if (fp->IsFaceEdgeS((j + 1) % 3) && fp->IsFaceEdgeS((j + 2) % 3))
                    continue;
                size_t indexV0 = vcg::tri::Index(mesh, fp->V0(j));
                size_t indexV1 = vcg::tri::Index(mesh, fp->V1(j));
                if (indexV0 >= indexV1)
                    continue;
                // calculate difference between face
                FacePointer fopp = fp->FFp(j);
                vcg::Point3i diff;
                for (size_t k = 0; k < 3; k++)
                {
                    diff[k] = (int)fp->C()[k] - (int)fopp->C()[k];
                }
                if (diff.Norm() > thereshold)
                {
                    fp->SetFaceEdgeS(j);
                    size_t jj = fp->FFi(j);
                    fopp->SetFaceEdgeS(jj);
                }
            }
        }
        mesh.ErodeDilate(erodilastep);
        mesh.InitEdgeType();
        mesh.SaveSharpFeatures("test.sharp");
        return true;
    }
};
