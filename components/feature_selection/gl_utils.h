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

#ifndef GL_UTILS
#define GL_UTILS


#include <wrap/gl/trimesh.h>
#include "ray_intersection.h"
#include "symmetry.h"

template <class MeshType>
void GLDrawSharpEdges(MeshType &mesh)
{
    typedef typename MeshType::CoordType  CoordType;

    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDisable(GL_LIGHTING);
    glDepthRange(0,0.9999);
    glLineWidth(5);
    glBegin(GL_LINES);
    for (size_t i=0;i<mesh.face.size();i++)
        for (size_t j=0;j<3;j++)
        {
            if (!mesh.face[i].IsFaceEdgeS(j))continue;

            if (mesh.face[i].FKind[j]==ETConcave)
                vcg::glColor(vcg::Color4b(0,0,255,255));
            else
                vcg::glColor(vcg::Color4b(255,0,0,255));

            CoordType Pos0=mesh.face[i].P0(j);
            CoordType Pos1=mesh.face[i].P1(j);
            vcg::glVertex(Pos0);
            vcg::glVertex(Pos1);
        }
    glEnd();
    glPopAttrib();
    
    std::vector<CoordType> quad;
    if (SymmetryManager<MeshType>::GetPlaneQuad(mesh, quad))
    {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_LIGHTING);
        glDepthRange(0,0.9999);
        glBegin(GL_QUADS);
        vcg::glColor(vcg::Color4b(125, 50, 199, 128));
        for (size_t i = 0; i < 4; i++)
            vcg::glVertex(quad[i]);
        glEnd();
        glPopAttrib();
    }

}

template<typename ScalarType>
void GLDrawAddedVertices(const std::vector<Intersection<ScalarType>>& interList)
{    
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glEnable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
    glDepthRange(0,0.9999);
    glPointSize(10);
    glBegin(GL_POINTS);

    for (auto v: interList)
    {
        if (v.edge == -1)
            vcg::glColor(vcg::Color4b(0,255,0,255));
        else
            vcg::glColor(vcg::Color4b(255,0,255,255));

        vcg::glVertex(v.pos);
    }
    
    glEnd();
    glPopAttrib();
}


#endif
