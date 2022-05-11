//***************************************************************************/
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


//**************************************************************************/
/* Modified by LinZhe 2022

Receive a triangular mesh along with 
- sharp feature
- symmetry feature
- texture
and convert it to quad mesh.

****************************************************************************/

#ifndef QUADRANGULATION_H
#define QUADRANGULATION_H

#include <iomanip>

#include <triangle_mesh_type.h>
#include <mesh_manager.h>
#include <vcg/space/box3.h>
#include <tracing/mesh_type.h>
#include <tracing/tracer_interface.h>

#include <load_save.h>
#include <mesh_types.h>
#include <smooth_mesh.h>
#include <quad_from_patches.h>
#include <quad_mesh_tracer.h>

#include <clocale>

template<typename TriMesh>
class Quadrangulation
{   
public:

    typedef typename TriMesh::ScalarType ScalarType;

    struct QuadPara
    {
        float Alpha=0.02; 
        float ScaleFact=1;
        float minFact=0.5;
        float maxFact=2;
        std::string PathM;
    };

    static ScalarType avgEdge(const TriangleMesh& trimesh)
    {
        ScalarType AvgVal=0;
        size_t Num=0;
        for (size_t i=0;i<trimesh.face.size();i++)
            for (size_t j=0;j<3;j++)
            {
                AvgVal+=(trimesh.face[i].cP0(j)-trimesh.face[i].cP1(j)).Norm();
                Num++;
            }
        return (AvgVal/Num);    
    }

    static void TracingPath(std::string PathM)
    {
        std::string pathProject=PathM;
        pathProject.erase(pathProject.find_last_of("."));

        PathM=pathProject;
        PathM.append("_rem.obj");
        std::cout<<"Loading Remeshed M:"<<PathM.c_str()<<std::endl;

        std::string PathF=pathProject;
        PathF.append("_rem.rosy");
        std::cout<<"Loading Rosy Field:"<<PathF.c_str()<<std::endl;

        std::string PathS=pathProject;
        PathS.append("_rem.sharp");
        std::cout<<"Loading Sharp F:"<<PathS.c_str()<<std::endl;

        std::string PathSymm=pathProject;
        PathSymm.append("_rem.symm");
        std::cout<<"Loading Symmetry Axis:"<<PathSymm.c_str()<<std::endl;

        //MESH LOAD
        TraceMesh trace_mesh;
        printf("Loading the mesh \n");
        bool loadedMesh=trace_mesh.LoadMesh(PathM);
        assert(loadedMesh);
        trace_mesh.UpdateAttributes();

        //FIELD LOAD
        bool loadedField=trace_mesh.LoadField(PathF);
        assert(loadedField);
        trace_mesh.UpdateAttributes();

        //SHARP LOAD
        bool loadedFeatures=trace_mesh.LoadSharpFeatures(PathS);
        assert(loadedFeatures);

        //SYMM LOAD
        bool loadedSymmetry=trace_mesh.LoadSymmetryAxis(PathSymm);
        if (!loadedSymmetry) std::cout << "Symmetry Axis loading failed!\n";

        // bool loadedAdpt=trace_mesh.LoadVertexAdapt(PathA);
        // assert(loadedAdpt);

        trace_mesh.SolveGeometricIssues();
        trace_mesh.UpdateSharpFeaturesFromSelection();

        //preprocessing mesh
        PreProcessMesh(trace_mesh);

        //initializing graph
        VertexFieldGraph<TraceMesh> VGraph(trace_mesh);
        VGraph.InitGraph(false);

        //INIT TRACER
        typedef PatchTracer<TraceMesh> TracerType;
        TracerType PTr(VGraph);
        TraceMesh::ScalarType Drift=100;
        bool add_only_needed=true;
        bool final_removal=true;
        bool meta_mesh_collapse=true;
        bool force_split=false;
        bool force_symmetry=true;
        PTr.sample_ratio=0.01;
        PTr.CClarkability=1;
        PTr.split_on_removal=true;
        PTr.away_from_singular=true;
        PTr.match_valence=true;
        PTr.check_quality_functor=false;
        PTr.MinVal=3;
        PTr.MaxVal=5;
        PTr.Concave_Need=1;

        //TRACING
        PTr.InitTracer(Drift,false);
        RecursiveProcess<TracerType>(PTr,Drift, add_only_needed,final_removal,force_symmetry,true,meta_mesh_collapse,force_split,true,false);
        // recover
        trace_mesh.RecoverAdapt();
        PTr.SmoothPatches();
        SaveAllData(PTr,pathProject,0,false,false);
    }

    static void Patch2Quad(PolyMesh& quadmesh, const QuadPara& QuadP)
    {
        ScalarType ScaleFact=QuadP.ScaleFact;
        ScalarType Alpha=QuadP.Alpha;
        std::string PathM=QuadP.PathM;

        TriangleMesh to_quad_trimesh;
        std::vector<std::vector<size_t>> trimeshPartitions;
        std::vector<std::vector<size_t>> trimeshCorners;
        std::vector<std::pair<size_t,size_t> > trimeshFeatures;
        std::vector<size_t> trimeshFeaturesC;

        // PolyMesh quadmesh;
        std::vector<std::vector<size_t>> quadmeshPartitions;
        std::vector<std::vector<size_t>> quadmeshCorners;
        std::vector<int> ilpResult;

        PathM.erase(PathM.find_last_of("."));
        std::string pathProject = PathM;
        PathM.append("_p0.obj");
        int mask=vcg::tri::io::Mask::IOM_WEDGTEXCOORD;
        vcg::tri::io::ImporterOBJ<TriMesh>::LoadMask(PathM.c_str(), mask);
        std::cout << "Mask: " << mask << std::endl;
        int err = vcg::tri::io::ImporterOBJ<TriangleMesh>::Open(to_quad_trimesh, PathM.c_str(), mask);
        if ((err!=0)&&(err!=5))
            assert(0);

        //FACE PARTITIONS
        std::string partitionFilename = pathProject;
        partitionFilename.append("_p0.patch");
        trimeshPartitions = loadPatches(partitionFilename);
        std::cout<<"Loaded "<<trimeshPartitions.size()<<" patches"<<std::endl;

        //PATCH CORNERS
        std::string cornerFilename = pathProject;
        cornerFilename.append("_p0.corners");
        trimeshCorners = loadCorners(cornerFilename);
        std::cout<<"Loaded "<<trimeshCorners.size()<<" corners set"<<std::endl;

        //FEATURES
        std::string featureFilename = pathProject;
        featureFilename.append("_p0.feature");
        trimeshFeatures = LoadFeatures(featureFilename);
        std::cout<<"Loaded "<<trimeshFeatures.size()<<" features"<<std::endl;

        //FEATURE CORNERS
        std::string featureCFilename = pathProject;
        featureCFilename.append("_p0.c_feature");
        trimeshFeaturesC = loadFeatureCorners(featureCFilename);
        std::cout<<"Loaded "<<featureCFilename.size()<<" corner features"<<std::endl;
        loadFeatureCorners(featureCFilename);

        // adpt
        std::string adaptiveness = pathProject;
        adaptiveness.append("_p0.adpt");
        bool adpt = LoadVertexAdapt(to_quad_trimesh, adaptiveness);
        assert(adpt);

        OrientIfNeeded(to_quad_trimesh,trimeshPartitions,trimeshCorners,trimeshFeatures,trimeshFeaturesC);

        //COMPUTE QUADRANGULATION
        QuadRetopology::internal::updateAllMeshAttributes(to_quad_trimesh);

        QuadRetopology::Parameters parameters;
        float scaleFactor;
        int fixedChartClusters;

        parameters.alpha=Alpha;
        parameters.ilpMethod=QuadRetopology::ILPMethod::LEASTSQUARES;
        parameters.timeLimit=200;
        parameters.gapLimit=0.0;
        parameters.callbackTimeLimit.push_back(3.0);
        parameters.callbackTimeLimit.push_back(5.0);
        parameters.callbackTimeLimit.push_back(10.0);
        parameters.callbackTimeLimit.push_back(20.0);
        parameters.callbackTimeLimit.push_back(30.0);
        parameters.callbackTimeLimit.push_back(60.0);
        parameters.callbackTimeLimit.push_back(90.0);
        parameters.callbackTimeLimit.push_back(120.0);

        parameters.callbackGapLimit.push_back(0.005);
        parameters.callbackGapLimit.push_back(0.02);
        parameters.callbackGapLimit.push_back(0.05);
        parameters.callbackGapLimit.push_back(0.1);
        parameters.callbackGapLimit.push_back(0.15);
        parameters.callbackGapLimit.push_back(0.20);
        parameters.callbackGapLimit.push_back(0.25);
        parameters.callbackGapLimit.push_back(0.3);

        parameters.minimumGap=0.4;

        parameters.isometry=true;

        parameters.regularityQuadrilaterals=true;

        parameters.regularityNonQuadrilaterals=true;

        parameters.regularityNonQuadrilateralsWeight=0.9;

        parameters.alignSingularities=true;

        parameters.alignSingularitiesWeight=0.1;

        parameters.repeatLosingConstraintsIterations=true;

        parameters.repeatLosingConstraintsQuads=false;

        parameters.repeatLosingConstraintsNonQuads=false;

        parameters.repeatLosingConstraintsAlign=true;

        parameters.hardParityConstraint=true;
        
        scaleFactor=ScaleFact;

        fixedChartClusters=300;

        parameters.chartSmoothingIterations = 0; //Chart smoothing
        parameters.quadrangulationFixedSmoothingIterations = 0; //Smoothing with fixed borders of the patches
        parameters.quadrangulationNonFixedSmoothingIterations = 0; //Smoothing with fixed borders of the quadrangulation
        parameters.feasibilityFix = false;

        double EdgeSize=avgEdge(to_quad_trimesh)*scaleFactor;
        std::cout<<"Edge Size "<<EdgeSize<<std::endl;
        const std::vector<double> edgeFactor(trimeshPartitions.size(), EdgeSize); //adaptiveness
        
        float scaleFactorMin=EdgeSize*QuadP.minFact;
        float scaleFactorMax=EdgeSize*QuadP.maxFact;

        qfp::AdaptParam aPar(true, scaleFactorMin, scaleFactorMax);

        qfp::quadrangulationFromPatches(to_quad_trimesh, trimeshPartitions, trimeshCorners, edgeFactor, aPar, parameters, fixedChartClusters, quadmesh, quadmeshPartitions, quadmeshCorners, ilpResult);

        //SAVE OUTPUT
        std::string outputFilename = pathProject;
        outputFilename+=std::string("_quadrangulation")+std::string(".obj");
        mask = vcg::tri::io::Mask::IOM_WEDGTEXCOORD;
        vcg::tri::io::ExporterOBJ<PolyMesh>::Save(quadmesh, outputFilename.c_str(),mask);


        //SMOOTH
        std::vector<size_t> QuadPart(quadmesh.face.size(),0);
        for (size_t i=0;i<quadmeshPartitions.size();i++)
            for (size_t j=0;j<quadmeshPartitions[i].size();j++)
                QuadPart[quadmeshPartitions[i][j]]=i;

        std::vector<size_t> TriPart(to_quad_trimesh.face.size(),0);
        for (size_t i=0;i<trimeshPartitions.size();i++)
            for (size_t j=0;j<trimeshPartitions[i].size();j++)
                TriPart[trimeshPartitions[i][j]]=i;

        std::vector<size_t> QuadCornersVect;
        for (size_t i=0;i<quadmeshCorners.size();i++)
            for (size_t j=0;j<quadmeshCorners[i].size();j++)
                QuadCornersVect.push_back(quadmeshCorners[i][j]);

        std::sort(QuadCornersVect.begin(),QuadCornersVect.end());
        auto last=std::unique(QuadCornersVect.begin(),QuadCornersVect.end());
        QuadCornersVect.erase(last, QuadCornersVect.end());

        std::cout<<"** SMOOTHING **"<<std::endl;
        MultiCostraintSmooth(quadmesh,to_quad_trimesh,trimeshFeatures,trimeshFeaturesC,TriPart,QuadCornersVect,QuadPart,0.5,EdgeSize,30,1);

        //SAVE OUTPUT
        outputFilename = pathProject;
        outputFilename+=std::string("_quadrangulation_smooth")+std::string(".obj");

        vcg::tri::io::ExporterOBJ<PolyMesh>::Save(quadmesh, outputFilename.c_str(),mask);
    }

    static void QuadrangulateFromTri(TriMesh& Mesh, const QuadPara& Para, PolyMesh& quadmesh)
    {
        std::cout << "1- Saving Tri Data to \"$Path_rem...\"\n";
        MeshPrepocess<TriMesh>::SaveAllData(Mesh, Para.PathM);
        
        std::cout << "2- Tracing...\n";
        TracingPath(Para.PathM);

        std::cout << "3- Quadding...\n";
        Patch2Quad(quadmesh, Para);
    }
};



















#endif //!QUADRANGULATION_H