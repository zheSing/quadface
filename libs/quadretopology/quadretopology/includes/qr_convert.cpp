/***************************************************************************/
/* Copyright(C) 2021


The authors of

Reliable Feature-Line Driven Quad-Remeshing
Siggraph 2021


 All rights reserved.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
****************************************************************************/

#include "qr_convert.h"

#include <vcg/complex/complex.h>

namespace QuadRetopology {
namespace internal {

template<class PolyMeshType>
void VCGToEigen(
        PolyMeshType& vcgMesh,
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F,
        std::vector<int>& vMap,
        std::vector<int>& fMap,
        bool selectedOnly,
        int numVerticesPerFace,
        int dim)
{
    assert(dim >= 2);
    assert(numVerticesPerFace > 2);

    int nSelectedVertices = 0;
    for (size_t i = 0; i < vcgMesh.vert.size(); i++){
        if ((!selectedOnly || vcgMesh.vert[i].IsS()) && !vcgMesh.vert[i].IsD()) {
            nSelectedVertices++;
        }
    }
    int nSelectedFaces = 0;
    for (size_t i = 0; i < vcgMesh.face.size(); i++){
        if ((!selectedOnly || vcgMesh.face[i].IsS()) && !vcgMesh.face[i].IsD()) {
            nSelectedFaces++;
        }
    }

    V.resize(nSelectedVertices, dim);
    F.resize(nSelectedFaces, numVerticesPerFace);

    vMap.resize(vcgMesh.vert.size(), -1);
    int vId = 0;
    for (size_t i = 0; i < vcgMesh.vert.size(); i++){
        if ((!selectedOnly || vcgMesh.vert[i].IsS()) && !vcgMesh.vert[i].IsD()) {
            vMap[i] = vId;
            for (int j = 0; j < dim; j++) {
                V(vId, j) = vcgMesh.vert[i].P()[j];
            }
            vId++;
        }
    }

    fMap.resize(vcgMesh.face.size(), -1);
    int fId = 0;
    for (size_t i = 0; i < vcgMesh.face.size(); i++){
        if ((!selectedOnly || vcgMesh.face[i].IsS()) && !vcgMesh.face[i].IsD()) {
            fMap[i] = fId;
            for (int j = 0; j < vcgMesh.face[i].VN(); j++) {
                F(fId, j) = vMap[vcg::tri::Index(vcgMesh, vcgMesh.face[i].V(j))];
            }
            fId++;
        }
    }
}

template<class PolyMeshType>
void VCGToEigenUV(
        PolyMeshType& vcgMesh,
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F,
        std::vector<std::vector<double>>& UV,
        bool selectedOnly)
{
    int dbg = 0;
    assert (vcg::tri::HasPerWedgeTexCoord(vcgMesh));
    UV.resize(F.rows());
    int fId = 0;
    for (size_t i = 0; i < vcgMesh.face.size(); i++){
        if ((!selectedOnly || vcgMesh.face[i].IsS()) && !vcgMesh.face[i].IsD()) {
            for (int j = 0; j < vcgMesh.face[i].VN(); j++) {
                for (int k = 0; k < 2; k++) {
                    UV[fId].push_back(vcgMesh.face[i].WT(j).P()[k]);
                }
            }
            fId++;
        }
    }
}

template<class PolyMeshType>
void eigenToVCG(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        PolyMeshType& vcgMesh,
        int numVertices,
        int dim)
{
    assert(dim >= 2);
    assert(numVertices > 2);

    vcgMesh.Clear();

    vcg::tri::Allocator<PolyMeshType>::AddVertices(vcgMesh, V.rows());
    for (int i = 0; i < V.rows(); i++) {
        typename PolyMeshType::CoordType vv(V(i,0), V(i,1), V(i,2));
        vcgMesh.vert[static_cast<size_t>(i)].P() = vv;
    }

    vcg::tri::Allocator<PolyMeshType>::AddFaces(vcgMesh, static_cast<size_t>(F.rows()));
    for (int i = 0; i < F.rows(); i++) {
        vcgMesh.face[static_cast<size_t>(i)].Alloc(numVertices);
        for (int j = 0; j < numVertices; j++) {
            vcgMesh.face[static_cast<size_t>(i)].V(j) = &(vcgMesh.vert[static_cast<size_t>(F(i,j))]);
        }
    }
}

template<class PolyMeshType>
void eigenUVToVCG(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::MatrixXd& UV,
        PolyMeshType& vcgMesh,
        int numVertices)
{
    assert(vcg::tri::HasPerWedgeTexCoord(vcgMesh));

    for (int i = 0; i < F.rows(); i++) {
        for (int j = 0; j < numVertices; j++) {
            size_t vidx = static_cast<size_t>(F(i,j));
            vcgMesh.face[static_cast<size_t>(i)].WT(j).P() = vcg::Point2f((float)UV(vidx, 0), (float)UV(vidx, 1));
        }
    }
}

}
}
