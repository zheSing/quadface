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
    typedef std::vector<Intersection<ScalarType>> InterListType;
    typedef vcg::Point3<ScalarType> CoordType;
    typedef vcg::Point2<ScalarType> TexCoordType;

    static size_t AddExtraVertex(MeshType& mesh, BVHT<MeshType>& bvh, const InterListType& inters)
    {   
        size_t count = 0;

        // topology
        const size_t new_idx[3][3] = {{1,3,0},{1,2,3},{3,2,0}};

        // face history
        std::set<size_t> fp_set;

        std::cout << "DEBUG: 1"  << std::endl;

        // traverse and add vertices
        for (size_t i = 0; i < inters.size(); i++)
        {
            // set coord
            VertexIterator vi = vcg::tri::Allocator<MeshType>::AddVertex(mesh, inters[i].pos);
            
            // face...
            size_t idx = inters[i].idx;

            std::cout << "DEBUG: 2"  << std::endl;

            // not found, simply split face 
            if (fp_set.find(idx) == fp_set.end())
            {
                std::cout << "DEBUG: 3"  << std::endl;

                // insert into history
                fp_set.insert(idx);
                FacePointer fp = &mesh.face[idx];

                // origin vertex
                VertexPointer vp[4] = {fp->V(0),fp->V(1),fp->V(2),&(*vi)};

                // add 3 faces, set vertices
                FaceIterator fi = vcg::tri::Allocator<MeshType>::AddFaces(mesh, 3);
                for (size_t m = 0; m < 3; m++)
                    for (size_t n = 0; n < 3; n++)
                        (fi+m)->V(n) = vp[new_idx[m][n]];
                
                // reset fp
                fp = &mesh.face[idx];
                
                std::cout << "DEBUG: 4"  << std::endl;

                // set texcoord
                if (fp->HasWedgeTexCoord())
                {
                    std::cout << "DEBUG: 5"  << std::endl;
                 
                    CoordType bary = inters[i].bary;
                    std::cout << "DEBUG: 5.1"  << std::endl;

                    TexCoordType texs[4] = {};
                    for (size_t m = 0; m < 3; m++)
                    {
                        typename MeshType::FaceType::TexCoordType t = fp->WT(m);
                        std::cout << fp->WT(m).P()[0] << std::endl;
                        texs[m] = fp->WT(m).P();
                    }
                    
                    std::cout << "DEBUG: 5.2"  << std::endl;
                    texs[3] = bary[0]*texs[0] + bary[1]*texs[1] + bary[2]*texs[2];
                    std::cout << "DEBUG: 5.3"  << std::endl;
                     
                    for (size_t m = 0; m < 3; m++)
                        for (size_t n = 0; n < 3; n++)
                            (fi+m)->WT(n).P() = texs[new_idx[m][n]];
                }

                std::cout << "DEBUG: 6"  << std::endl;

                std::cout << "DEBUG: 7"  << std::endl;

                count++;
            }
            
            else
            {
                std::cout << "Failed! Same face...\n";
            }
        }

        // update mesh
        for (auto idx: fp_set)
            mesh.face[idx].SetD();
        mesh.UpdateDataStructures();
        vcg::tri::Allocator<MeshType>::CompactEveryVector(mesh);


        // update bvh
        bvh.BuildTree(mesh);

        std::cout << "DEBUG: 9"  << std::endl;

        return count;
    }
};












#endif
