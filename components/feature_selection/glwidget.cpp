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

/* LZ change
- class QGLWidget -> QOpenGLWidget
- function updateGL() -> update()
- function delta() -> -angleDelta().y()
- OpenGL draw mode to DMFlat

- add variable "has_texture", "pathT"
- add header "feature_texture.h"
*/ 

#include <GL/glew.h>
#include <QMouseEvent>

#include <chrono>

#include <math.h>
#include "glwidget.h"
#include "mesh_manager.h"
#include "ray_intersection.h"
#include "bvh.h"
#include "fields/field_smoother.h"
#include "feature_texture.h"
#include "history_queue.h"
#include "add_vertex.h"
#include "symmetry.h"
#include "triangle_mesh_type.h"
// #include "poly_mesh_type.h"
#include <wrap/qt/trackball.h>
#include <wrap/gl/picking.h>
#include <wrap/qt/anttweakbarMapper.h>
#include <wrap/io_trimesh/import_field.h>
#include <wrap/io_trimesh/export_field.h>
#include <wrap/io_trimesh/export.h>
#include <wrap/gl/trimesh.h>
//#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>
#include <wrap/gl/gl_field.h>
#include <AutoRemesher.h>
#include <vcg/complex/algorithms/hole.h>
#include "mesh_field_smoother.h"
#include "gl_utils.h"

#ifdef MIQ_QUADRANGULATE
#include <wrap/igl/miq_parametrization.h>
#include <vcg/complex/algorithms/quadrangulator.h>
#endif

std::string pathM="";
std::string pathS="";
std::string pathT="";

std::string projM="";
std::string projS="";
std::string projSymm="";

vcg::Trackball track;//the active manipulator

bool drawfield=false;

FieldTriMesh tri_mesh;

#ifdef MIQ_QUADRANGULATE
PMesh quad_mesh;
bool quadrangulated=false;
vcg::tri::MiQParametrizer<FieldTriMesh>::MIQParameters MiqP;
#endif

TwBar *barQuad;

vcg::GlTrimesh<FieldTriMesh> glWrap;
vcg::GLField<FieldTriMesh> glField;

typedef typename FieldTriMesh::ScalarType ScalarType;
typedef typename FieldTriMesh::CoordType CoordType;

int Iterations;
ScalarType EdgeStep;
ScalarType Multiplier=2;
ScalarType SharpFactor=6;
ScalarType alpha=0.3;
//ScalarType sharp_feature_thr=45;


int xMouse,yMouse;

//vcg::GridStaticPtr<FieldTriMesh::FaceType,FieldTriMesh::ScalarType> Gr;

typedef vcg::tri::FieldSmoother<FieldTriMesh> FieldSmootherType;
FieldSmootherType::SmoothParam FieldParam;

AutoRemesher<FieldTriMesh>::Params RemPar;
bool do_batch=false;
bool has_features=false;
bool has_features_fl=false;
bool has_texture=false;
//size_t ErodeDilateSteps=4;

bool do_remesh=true;
bool do_init_feature=false;
int remesher_iterations=15;
ScalarType remesher_aspect_ratio=0.3;
int remesher_termination_delta=10000;
bool surf_dist_check=true;
ScalarType sharp_feature_thr=35;
int feature_erode_dilate=4;

//MeshPrepocess<FieldTriMesh> MP(tri_mesh);

u_int tID;
cv::Mat texture;

enum gui_mode{ None, Vertex, Edge, Symm } uimode, uimode_;

typedef MyIntersection<ScalarType> InterType;
typedef Ray<ScalarType> RayType;

BVHT<FieldTriMesh> bvh_tree;
HistoryQueue<InterType> vertices_added;

void InitSharp()
{
    assert(!(has_features || has_features_fl));
    tri_mesh.UpdateDataStructures();
    tri_mesh.InitSharpFeatures(sharp_feature_thr);
    //MP.InitSharpFeatures(sharp_feature_thr);
}



void SaveAllData()
{
    MeshPrepocess<FieldTriMesh>::SaveAllData(tri_mesh,pathM);
}

void DoBatchProcess ()
{
    typename MeshPrepocess<FieldTriMesh>::BatchParam BPar;
    BPar.DoRemesh=do_remesh;
    BPar.feature_erode_dilate=feature_erode_dilate;
    BPar.remesher_aspect_ratio=remesher_aspect_ratio;
    BPar.remesher_iterations=remesher_iterations;
    BPar.remesher_termination_delta=remesher_termination_delta;
    BPar.SharpFactor=SharpFactor;
    BPar.sharp_feature_thr=sharp_feature_thr;
    BPar.surf_dist_check=surf_dist_check;
    BPar.UpdateSharp=do_init_feature;
    MeshPrepocess<FieldTriMesh>::BatchProcess(tri_mesh,BPar,FieldParam);
    drawfield=true;
}


void TW_CALL VertexUndo(void*)
{
    vertices_added.Undo();
}

void TW_CALL VertexRedo(void*)
{
    vertices_added.Redo();
}

void TW_CALL VertexConfirm(void*)
{
    uimode_ = None;
    ExtraVertexProcess<FieldTriMesh>::AddExtraVertex(tri_mesh, bvh_tree, vertices_added.HistoryList());
    vertices_added.Clear();
}

void TW_CALL VertexCancel(void*)
{
    uimode_ = None;
    vertices_added.Clear();
}

void TW_CALL AddVertices(void*)
{
    uimode_ = Vertex;
}

void TW_CALL EdgeConfirm(void*)
{
    uimode_ = None;
}

void TW_CALL SelectFeatures(void*)
{
    uimode_ = Edge;
}

void TW_CALL SelectSymmetry(void *)
{
    uimode_ = Symm;
}

void TW_CALL SaveFeatures(void*)
{
    tri_mesh.SaveSharpFeatures(projS);
    std::cout<<"Saving Sharp TO:" << projS.c_str() << std::endl;
}

void TW_CALL SaveMesh(void*)
{
    // size_t indexExt=pathM.find_last_of(".");       
    // std::string projM = pathM.substr(0,indexExt) + "_ref.obj";
    tri_mesh.SaveTriMesh(projM);
    std::cout<<"Saving Mesh TO:" << projM.c_str() << std::endl;
}

void TW_CALL SaveSymmetryAxis(void*)
{
    tri_mesh.SaveSymmetryAxis(projSymm);
    std::cout<<"Saving Symmetry Axis TO:" << projSymm.c_str() << std::endl;
}

void DoAutoRemesh()
{
    RemPar.iterations   = remesher_iterations;
    RemPar.targetAspect = remesher_aspect_ratio;
    RemPar.targetDeltaFN= remesher_termination_delta;
    RemPar.creaseAngle  = sharp_feature_thr;
    RemPar.erodeDilate  = feature_erode_dilate;

    //RemPar.userSelectedCreases = true;
    RemPar.surfDistCheck = surf_dist_check;

    AutoRemesher<FieldTriMesh>::RemeshAdapt(tri_mesh,RemPar);
    //AutoRemesher<FieldTriMesh>::RemeshAdapt(tri_mesh,RemPar);

    tri_mesh.InitEdgeType();
    tri_mesh.InitFeatureCoordsTable();
//    tri_mesh.InitSharpFeatures(sharp_feature_thr);
//    tri_mesh.ErodeDilate(feature_erode_dilate);

    //MeshPrepocess<FieldTriMesh>::AutoRemesh(tri_mesh,RemPar);
//    Gr.Set(tri_mesh.face.begin(),tri_mesh.face.end());
}

void TW_CALL CleanMesh(void *)
 {
   MeshPrepocess<FieldTriMesh>::SolveGeometricArtifacts(tri_mesh);
 }

void TW_CALL AutoRemesh(void *)
{
    DoAutoRemesh();
}

void TW_CALL InitSharpFeatures(void *)
{
   InitSharp();
}

void TW_CALL RefineIfNeeded(void *)
{
   // tri_mesh.RefineIfNeeded();
    MeshPrepocess<FieldTriMesh>::RefineIfNeeded(tri_mesh);
    //MP.RefineIfNeeded();
    //Gr.Set(tri_mesh.face.begin(),tri_mesh.face.end());
}

void TW_CALL BatchProcess(void *)
{
    DoBatchProcess();
}

void TW_CALL SmoothField(void *)
{
//    //tri_mesh.SplitFolds();
//    MeshPrepocess<FieldTriMesh>::SplitFolds(tri_mesh);
//    //tri_mesh.RemoveFolds();
//    MeshPrepocess<FieldTriMesh>::RemoveFolds(tri_mesh);
//    //tri_mesh.RemoveSmallComponents();
//    MeshPrepocess<FieldTriMesh>::RemoveSmallComponents(tri_mesh);
//    MP.SplitFolds();
//    MP.RemoveFolds();
//    MP.RemoveSmallComponents();

    //vcg::tri::io::ExporterPLY<FieldTriMesh>::Save(tri_mesh,"test0.ply");
    //tri_mesh.SmoothField(FieldParam);
    MeshFieldSmoother<FieldTriMesh>::SmoothField(tri_mesh,FieldParam);
    //vcg::tri::io::ExporterPLY<FieldTriMesh>::Save(tri_mesh,"test1.ply");
    drawfield=true;

    //Gr.Set(tri_mesh.face.begin(),tri_mesh.face.end());
}



#ifdef MIQ_QUADRANGULATE
void TW_CALL MiqQuadrangulate(void *)
{

    vcg::tri::MiQParametrizer<FieldTriMesh>::MIQParametrize(tri_mesh,MiqP);//,MaxFeature);
    vcg::tri::Quadrangulator<FieldTriMesh,PMesh> Quadr;

    std::cout<<"Quadrangulating"<<std::endl;
    FieldTriMesh splitted_mesh;
    vcg::tri::Append<FieldTriMesh,FieldTriMesh>::Mesh(splitted_mesh,tri_mesh);
    std::vector< std::vector< short int> > AxisUV;

    Quadr.Quadrangulate(splitted_mesh,quad_mesh,AxisUV,false);
    quadrangulated=true;
    quad_mesh.UpdateAttributes();
    vcg::tri::Allocator<PMesh>::CompactEveryVector(quad_mesh);
    std::cout<<"Num Faces: "<<quad_mesh.face.size()<<std::endl;
    vcg::tri::io::ExporterOBJ<PMesh>::Save(quad_mesh,"quadrangulated.obj",vcg::tri::io::Mask::IOM_BITPOLYGONAL);
}
#endif



void TW_CALL SaveData(void *)
{
    SaveAllData();
}

void TW_CALL AutoSetupField(void *)
{
    MeshFieldSmoother<FieldTriMesh>::AutoSetupParam(tri_mesh,FieldParam);
}

void TW_CALL ErodeDilateFeatureStep(void *)
{
    tri_mesh.ErodeDilate(feature_erode_dilate);
    //MP.ErodeDilate(feature_erode_dilate);
}

void TW_CALL InitTexFeature(void *)
{
    TextureProcess<FieldTriMesh>::TexCoordFeature(tri_mesh);
}

void TW_CALL InitBorderFeature(void *)
{
    std::vector<FieldTriFace>& face = tri_mesh.face;
    for (size_t i = 0; i < face.size(); i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            if (face[i].IsD()) continue;
            if (vcg::face::IsBorder(face[i],j))
            {
                face[i].SetFaceEdgeS(j);
                face[i].FKind[j]=ETConvex;
            }
            else if (!vcg::face::IsManifold(face[i],j))
            {
                face[i].SetFaceEdgeS(j);
                face[i].FKind[j]=ETConvex;
            }
        }
    }
}

void TW_CALL SymmConfirm(void *)
{
    // vcg::Point4<ScalarType> plane;
    // SymmetryManager<FieldTriMesh>::FindSymmetryAxis(tri_mesh, plane);
    uimode_ = None;
}

void TW_CALL InitThrFeature(void *)
{
    // backup
    std::vector<FieldTriFace>& face = tri_mesh.face;
    std::vector<int> back(face.size(), 0);
    for (size_t i = 0; i < face.size(); i++)
    {
        if (face[i].IsD()) continue;
        for (size_t j = 0; j < 3; j++)
        {
            if (face[i].IsFaceEdgeS(j))
                back[i] |= (1<<j); 
        }
    }
    
    if (sharp_feature_thr > 0)
    {
        vcg::tri::UpdateFlags<FieldTriMesh>::FaceEdgeSelCrease(tri_mesh, vcg::math::ToRad(sharp_feature_thr));
    }

    for (size_t i = 0; i < face.size(); i++)
    {
        if (face[i].IsD()) continue;
        for (size_t j = 0; j < 3; j++)
        {
            if (back[i] & (1<<j))
                face[i].SetFaceEdgeS(j);; 
        }
    }
}


void SetFieldBarSizePosition(QWidget *w)
{
    int params[2];
    params[0] = QTDeviceWidth(w) / 3;
    params[1] = QTDeviceHeight(w) / 1.8;
    TwSetParam(barQuad, NULL, "size", TW_PARAM_INT32, 2, params);
    params[0] = QTLogicalToDevice(w, 10);
    params[1] = 30;//QTDeviceHeight(w) - params[1] - QTLogicalToDevice(w, 10);
    TwSetParam(barQuad, NULL, "position", TW_PARAM_INT32, 2, params);
}

void InitFieldBar(QWidget *w)
{
    if (uimode == None)
    {
        TwDeleteAllBars();

        barQuad = TwNewBar("QuadWild");

        // SetFieldBarSizePosition(w);

        TwAddButton(barQuad, "AddVertices", AddVertices, 0, "label='AddVertices'");
        TwAddButton(barQuad, "SelectFeatures", SelectFeatures, 0, "label='SelectFeatures'");
        TwAddButton(barQuad, "SelectSymmetry", SelectSymmetry, 0, "label='SelectSymmetry'");

        TwAddSeparator(barQuad, "sep1", "label=''");
        TwAddButton(barQuad, "InitTexFeature", InitTexFeature, 0, "label='InitTexFeature'");
        TwAddButton(barQuad, "InitBorderFeature", InitBorderFeature, 0, "label='InitBorderFeature'");
        TwAddVarRW(barQuad, "sharp_feature_thr", TW_TYPE_DOUBLE, &sharp_feature_thr," label='Sharp Degree'");
        TwAddButton(barQuad, "InitThrFeature", InitThrFeature, 0, "label='InitThrFeature'");

        TwAddSeparator(barQuad, "sep2", "label=''");

        // TwAddVarRW(barQuad,"LimitConcave",TW_TYPE_DOUBLE, &tri_mesh.LimitConcave," label='Limit Concave'");
        //TwAddVarRW(barQuad,"LimitConcave",TW_TYPE_DOUBLE, &MP.LimitConcave," label='Limit Concave'");

        // TwAddButton(barQuad,"CleanMesh",CleanMesh,0,"label='CleanMesh'");

        // TwAddButton(barQuad,"SetSharp",InitSharpFeatures,0,"label='InitSharp'");

        TwAddVarRW(barQuad,"ErodeDilSteps",TW_TYPE_INT32,&feature_erode_dilate,"label='ErodeDilateSteps'");

        TwAddButton(barQuad,"Erode Dilate",ErodeDilateFeatureStep,0,"label='ErodeDilateSharp'");

        TwAddSeparator(barQuad, "sep3", "label=''");

        TwAddVarRW(barQuad, "SharpName", TW_TYPE_STDSTRING, &projS, "label='SharpName'");
        TwAddButton(barQuad, "SaveFeatures", SaveFeatures, 0, "label='SaveFeatures'");
        TwAddVarRW(barQuad, "MeshName", TW_TYPE_STDSTRING, &projM, "label='MeshName'");
        TwAddButton(barQuad, "SaveMesh", SaveMesh, 0, "label='SaveMesh'");
        TwAddVarRW(barQuad, "SymmName", TW_TYPE_STDSTRING, &projSymm, "label='SymmName'");
        TwAddButton(barQuad, "SaveSymmetryAxis", SaveSymmetryAxis, 0, "label='SaveSymmetryAxis'");

        TwAddSeparator(barQuad, "sep4", "label=''");
        TwAddVarRW(barQuad, "Init_Feature", TW_TYPE_BOOLCPP, &do_init_feature, "label='Init_Feature'");
        TwAddVarRW(barQuad, "Do_Remesh", TW_TYPE_BOOLCPP, &do_remesh, "label='Do_Remesh'");
        TwAddButton(barQuad, "BatchProcess", BatchProcess, 0, "label='Batch Process'");
        TwAddButton(barQuad, "SaveData", SaveData, 0, "label='SaveData'");

        // TwAddButton(barQuad,"AutoRemesh",AutoRemesh,0,"label='AutoRemesh'");

        // TwAddButton(barQuad,"Refine",RefineIfNeeded,0,"label='Refine if needed'");

        // TwAddVarRW(barQuad,"Alpha",TW_TYPE_DOUBLE, &FieldParam.alpha_curv," label='Alpha Curvature'");
        // TwAddVarRW(barQuad,"HardCT",TW_TYPE_DOUBLE, &FieldParam.curv_thr," label='Hard Curv Thr'");
        // TwAddVarRW(barQuad,"CurvRing",TW_TYPE_INT32,&FieldParam.curvRing,"label='Curvature Ring'");

        // TwEnumVal smoothmodes[3] = {
        //     {vcg::tri::SMMiq,"MIQ"},
        //     {vcg::tri::SMNPoly,"NPoly"},
        //     {vcg::tri::SMIterative,"Ite"}
        // };
        // TwType smoothMode = TwDefineEnum("SmoothMode", smoothmodes, 3);
        // TwAddVarRW(barQuad, "Smooth Mode", smoothMode, &FieldParam.SmoothM," label='Smooth Mode' ");


        // TwAddButton(barQuad,"AutoSetup",AutoSetupField,0,"label='Auto Setup Field'");

        // TwAddButton(barQuad,"ComputeField",SmoothField,0,"label='Compute Field'");

        // TwAddButton(barQuad,"BatchProcess",BatchProcess,0,"label='Batch Process'");

        // TwAddButton(barQuad,"SaveData",SaveData,0,"label='Save Data'");

#ifdef MIQ_QUADRANGULATE
    TwAddSeparator(barQuad,"","");
    TwAddVarRW(barQuad,"Gradient",TW_TYPE_DOUBLE, &MiqP.gradient," label='Gradient'");
    TwAddVarRW(barQuad,"Direct Round",TW_TYPE_BOOLCPP, &MiqP.directRound," label='Direct Round'");
    TwAddVarRW(barQuad,"Round Singularities",TW_TYPE_BOOLCPP, &MiqP.round_singularities," label='Round Singularities'");
    TwAddVarRW(barQuad,"IsotropyVsAlign",TW_TYPE_DOUBLE, &MiqP.miqAnisotropy," label='Isotropy Vs Align'");
    TwAddVarRW(barQuad,"Align Sharp",TW_TYPE_BOOLCPP, & MiqP.crease_as_feature," label='Align Sharp'");
    TwAddButton(barQuad,"Quadrangulate",MiqQuadrangulate,0,"label='Miq Quadrangulate'");
#endif
    }

    else if (uimode == Vertex)
    {
        TwDeleteAllBars();
        barQuad = TwNewBar("Press right mouse button to add vertex...");
        TwAddButton(barQuad, "undo", VertexUndo, 0, "label='VertexUndo'");
        TwAddButton(barQuad, "redo", VertexRedo, 0, "label='VertexRedo'");
        TwAddButton(barQuad, "cancel", VertexCancel, 0, "label='VertexCancel'");
        TwAddButton(barQuad, "confirm", VertexConfirm, 0, "label='VertexConfirm'");
    }

    else if (uimode == Edge)
    {
        TwDeleteAllBars();
        barQuad = TwNewBar("Press right mouse to select/unselect sharp edge...");
        TwAddButton(barQuad, "confirm", EdgeConfirm, 0, "label='EdgeConfirm'");
    }

    else if (uimode == Symm)
    {
        TwDeleteAllBars();
        barQuad = TwNewBar("Symmetry Edge Selection");
        TwAddButton(barQuad, "SymmConfirm", SymmConfirm, 0, "label='Generate'");
    }
}


GLWidget::GLWidget(QWidget *parent)
    : QOpenGLWidget(parent)
{
    uimode = uimode_ = None;
    hasToPick=false;
    bool AllQuad=false;
    bool Loaded=tri_mesh.LoadTriMesh(pathM,AllQuad);
    FieldParam.alpha_curv=alpha;
    FieldParam.curv_thr=0.8;

    if (!Loaded)
    {
        std::cout<<"Error Loading Mesh"<<std::endl;
        exit(0);
    }
    if (AllQuad)
    {
        drawfield=true;
    }
    for (auto &f: tri_mesh.face)
    {
        f.C() = FieldTriMesh::FaceType::ColorType::White;
    }
    
    size_t indexExt=pathM.find_last_of(".");       
    projM = pathM.substr(0, indexExt) + "_ref.obj";
    projS = pathM.substr(0, indexExt) + "_ref.sharp";
    projSymm = pathM.substr(0, indexExt) + "_ref.symm";

    std::cout<<"Loaded "<<tri_mesh.face.size()<<" faces "<<std::endl;
    std::cout<<"Loaded "<<tri_mesh.vert.size()<<" vertices "<<std::endl;

    std::cout<<"Initializing BVH tree..."<<std::endl;
    bvh_tree.BuildTree(tri_mesh); 
    std::cout<<"BVH tree established."<<std::endl;


    glWrap.m=&tri_mesh;

    tri_mesh.UpdateDataStructures();
    tri_mesh.InitEdgeType();

    SymmetryManager<FieldTriMesh>::SetUpForSymmetry(tri_mesh);

    if (has_texture)
    {
        cv::flip(cv::imread(pathT, cv::ImreadModes::IMREAD_COLOR), texture, 0);

        // glGenTextures(1, &tID);
        // glBindTexture(GL_TEXTURE_2D, tID);
        // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        // glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texture.cols, texture.rows, 0, GL_BGR_EXT, GL_UNSIGNED_BYTE, texture.ptr());
        // std::cout << texture.cols << ", " << texture.rows << std::endl;
        // FieldTriMesh::ScalarType thereshold;
        // TextureProcess::SetUpAvgColor(tri_mesh, imgrgb, thereshold);
        // TextureProcess::ExtractTexFeature(tri_mesh, imgrgb);
        // tri_mesh.UpdateDataStructures();
        // TextureProcess::TexCoordFeature(tri_mesh);
        // TextureProcess::SegmentTexture(tri_mesh, imgrgb, 250, true);
        // std::cout << "Color setup finished.\n";
    }

    tri_mesh.LimitConcave=0;
    //MP.LimitConcave=0;
    if (has_features)
    {
        std::cout<<"*** Loading SHARP FEATURES ***"<<std::endl;
        bool HasRead=tri_mesh.LoadSharpFeatures(pathS);
        //bool HasRead=MP.LoadSharpFeatures(pathS);
        assert(HasRead);
    }

    if (has_features_fl)
    {
        std::cout<<"*** Loading SHARP FEATURES FL***"<<std::endl;
        bool HasRead=tri_mesh.LoadSharpFeaturesFL(pathS);
        //bool HasRead=MP.LoadSharpFeaturesFL(pathS);
        assert(HasRead);
    }
    tri_mesh.InitEdgeType();
    
//    if ((do_batch)&&(!has_features))
//    {
////        bool SufficientFeatures=mesh.SufficientFeatures(SharpFactor);
////        if (SufficientFeatures)

//        DoBatchProcess();
//        SaveAllData();
//        exit(0);
//    }
//    if ((do_batch)&&(has_features))
//    {
//        //tri_mesh.PrintSharpInfo();
//        DoBatchProcess();
//        SaveAllData();
//        exit(0);
//    }
    if (do_batch)
    {

        DoBatchProcess();
        SaveAllData();
        exit(0);
    }

    //remeshed_mesh.UpdateDataStructures();
    //Gr.Set(tri_mesh.face.begin(),tri_mesh.face.end());
}



void GLWidget::initializeGL ()
{
    //initialize Glew
    glewInit();
    //CaptInt.GLInit( GLWidget::width(),GLWidget::height());
    glClearColor(0, 0, 0, 0);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);

    if (has_texture)
    {
        glGenTextures(1, &tID);
        glBindTexture(GL_TEXTURE_2D, tID);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texture.cols, texture.rows, 0, GL_BGR_EXT, GL_UNSIGNED_BYTE, texture.ptr());
    }
}


void GLWidget::resizeGL (int w, int h)
{
    glViewport (0, 0, (GLsizei) w, (GLsizei) h);
    TwWindowSize(w, h);
    InitFieldBar(this);
    initializeGL();
}

void GLWidget::paintGL ()
{

    //    if (RType!=OldRType)
    //    {
    //        PatchDeco.ColorPatches(RType);
    //        OldRType=RType;
    //    }
    if (uimode != uimode_)
    {
        uimode = uimode_;
        InitFieldBar(this);
    }

    glClearColor(255,255,255,255);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40, GLWidget::width()/(float)GLWidget::height(), 0.1, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,3.5f,   0,0,0,   0,1,0);
    track.center=vcg::Point3f(0, 0, 0);
    track.radius= 1;
    track.GetView();


    glPushMatrix();
    track.Apply();
    glPushMatrix();

    glDisable(GL_CULL_FACE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

    if(tri_mesh.vert.size()>0)
    {
        vcg::glScale(2.0f/tri_mesh.bbox.Diag());
        glTranslate(-tri_mesh.bbox.Center());
        //tri_mesh.GLDrawSharpEdges();
        GLDrawSharpEdges(tri_mesh);
        GLDrawAddedVertices(vertices_added.HistoryList());
        //MP.GLDrawSharpEdges();
        if (has_texture)
        {
            glEnable(GL_TEXTURE_2D);
            glBindTexture(GL_TEXTURE_2D, tID);
        }
        glWrap.Draw(vcg::GLW::DMFlatWire,vcg::GLW::CMPerFace,vcg::GLW::TMPerWedge);
        //glWrap.Draw(vcg::GLW::DMSmooth,vcg::GLW::CMNone,vcg::GLW::TMNone);
    }

#ifdef MIQ_QUADRANGULATE
    if (quadrangulated)
    {
        quad_mesh.GLDraw();
    }
#endif

    if (drawfield)
    {
        vcg::GLField<FieldTriMesh>::GLDrawFaceField(tri_mesh,false,false,0.007);
        vcg::GLField<FieldTriMesh>::GLDrawSingularity(tri_mesh);
    }

//    if(hasToPick)
//    {
//        hasToPick=false;
//        typename FieldTriMesh::CoordType pp;
//        if(vcg::Pick<typename FieldTriMesh::CoordType>(xMouse,yMouse,pp))
//        {
//            typename FieldTriMesh::CoordType closPt,bary;
//            typename FieldTriMesh::ScalarType minD;
//            typename FieldTriMesh::FaceType *f=vcg::tri::GetClosestFaceBase(tri_mesh,Gr,pp,tri_mesh.bbox.Diag(),minD,closPt);
//            vcg::InterpolationParameters(*f,closPt,bary);
//            size_t EdgeI=1;
//            if ((bary.Y()<bary.X())&&(bary.Y()<bary.Z()))EdgeI=2;
//            if ((bary.Z()<bary.X())&&(bary.Z()<bary.Y()))EdgeI=0;

////            FieldTriMesh::FaceType *fOpp=f->FFp(EdgeI);
////            int eOpp=f->FFi(EdgeI);

//            if (f->IsFaceEdgeS(EdgeI))
//            {
//                tri_mesh.ClearSharp((*f),EdgeI);
//                //MP.ClearSharp((*f),EdgeI);
////                {
////                f->ClearFaceEdgeS(EdgeI);
////                if (fOpp!=f)
////                    fOpp->ClearFaceEdgeS(eOpp);
//            }else
//            {

//                tri_mesh.SetSharp((*f),EdgeI);
//                //MP.SetSharp((*f),EdgeI);
////                f->SetFaceEdgeS(EdgeI);
////                if (fOpp!=f)
////                    fOpp->SetFaceEdgeS(eOpp);
//            }
//        }
//    }

    glPopMatrix();
    glPopMatrix();


    TwDraw();

}

void GLWidget::keyReleaseEvent (QKeyEvent * e)
{
    e->ignore ();
    if (e->key () == Qt::Key_Control)  track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ControlModifier));
    if (e->key () == Qt::Key_Shift)  track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
    if (e->key () == Qt::Key_Alt) track.ButtonUp (QT2VCG (Qt::NoButton, Qt::AltModifier));
    update ();
}


void GLWidget::keyPressEvent (QKeyEvent * e)
{
    e->ignore ();
    if (e->key () == Qt::Key_Control) track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ControlModifier));
    if (e->key () == Qt::Key_Shift)  track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
    if (e->key () == Qt::Key_Alt)  track.ButtonDown (QT2VCG (Qt::NoButton, Qt::AltModifier));

    TwKeyPressQt(e);
    update ();
}

void GLWidget::mousePressEvent (QMouseEvent * e)
{
    if(!TwMousePressQt(this,e))
    {
        e->accept ();
        setFocus ();
        track.MouseDown(QT2VCG_X(this, e), QT2VCG_Y(this, e), QT2VCG (e->button (), e->modifiers ()));

        if (uimode != None && e->button() == Qt::RightButton)
        {
            int x = QT2VCG_X(this, e);
            int y = QT2VCG_Y(this, e);

            double coordx = (2 * (x + 0.5f) / (double)width() - 1) * tan(M_PI / 9) * 0.1 * ((double)width() / height());
            double coordy = (2 * (y + 0.5f) / (double)height() - 1) * tan(M_PI / 9) * 0.1;

            // from world to model        
            vcg::Matrix44f srt_inv = track.track.InverseMatrix();
            vcg::Matrix44f s_inv, t_inv;
            float inv_s = tri_mesh.bbox.Diag()/2.0f;
            s_inv.SetScale(inv_s, inv_s, inv_s);
            t_inv.SetTranslate(tri_mesh.bbox.Center());
            vcg::Matrix44f m_inv = t_inv * s_inv *srt_inv;

            // ray construction
            vcg::Point3f pixel(coordx, coordy, 3.4f);
            vcg::Point3f eye(0, 0, 3.5f);
            pixel = m_inv * pixel;
            eye = m_inv * eye;
            RayType ray(eye, pixel-eye);            

            // collide detection
            if (uimode == Vertex)
            {
                InterType inter;
                if (bvh_tree.Intersect(ray, inter))
                {
                    std::cout << "Intersection happended!\n";
                    if (ExtraVertexProcess<FieldTriMesh>::OnEdgeOrInternal(tri_mesh, inter) < 2)
                        vertices_added.Insert(inter);
                }
            }

            // edge select detection
            else if (uimode == Symm)
            {
                InterType inter;
                if (bvh_tree.Intersect(ray, inter))
                {
                    size_t cnt = SymmetryManager<FieldTriMesh>::TraceAxisByDirection(tri_mesh, inter);
                    std::cout << "Marked " << cnt << " symmetry axis.\n";
                }
                else 
                    std::cout << "Didn't select any edge!\n";
            }
            else if (uimode == Edge)
            {
                InterType inter;
                if (bvh_tree.Intersect(ray, inter) && ExtraVertexProcess<FieldTriMesh>::SelectFeatureEdge(tri_mesh, inter));
                else 
                    std::cout << "Didn't select any edge!\n";
            }
        }
    }
    update ();
}

void GLWidget::mouseMoveEvent (QMouseEvent * e)
{
    if (e->buttons ()) {
        track.MouseMove(QT2VCG_X(this, e), QT2VCG_Y(this, e));
        update ();
    }
    TwMouseMotion(QTLogicalToDevice(this, e->x()), QTLogicalToDevice(this, e->y()));
}

void GLWidget::mouseDoubleClickEvent (QMouseEvent * e)
{
    if (e->buttons ())
    {
        xMouse=QT2VCG_X(this, e);
        yMouse=QT2VCG_Y(this, e);
        //pointToPick=Point2i(e->x(),height()-e->y());
        //pointToPick=Point2i(xMouse,yMouse);
        hasToPick=true;
        update ();
    }
    update();
}

void GLWidget::mouseReleaseEvent (QMouseEvent * e)
{
    track.MouseUp(QT2VCG_X(this, e), QT2VCG_Y(this, e), QT2VCG(e->button (), e->modifiers ()));
    TwMouseReleaseQt(this,e);
    update ();
}

void GLWidget::wheelEvent (QWheelEvent * e)
{
    const int WHEEL_STEP = 120;
    track.MouseWheel (-e->angleDelta ().y() / float (WHEEL_STEP), QTWheel2VCG (e->modifiers ()));
    update ();
}
