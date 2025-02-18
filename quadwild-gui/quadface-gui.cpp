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
- remove header QDesktopWidget
- add variable extern "has_texture", "pathT"
- add args "*.png", "*.jpg"
*/ 

#include <QApplication>
// #include <QDesktopWidget>
#include <GL/glew.h>
#include "glwidget-gui.h"
#include <wrap/qt/anttweakbarMapper.h>
#include <QWindow>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QTextStream>
//#include <triangle_mesh_type.h>
//#include <wrap/io_trimesh/import_field.h>
//#include <wrap/io_trimesh/import.h>

#include<vcg/complex/algorithms/hole.h>
#include <cstdio>
#include <clocale>
//#include <locale>


extern bool do_batch;
extern bool has_features;
extern bool has_features_fl;
extern bool has_texture;
extern bool has_symm;
extern std::string pathM,pathS,pathT,pathSymm;
extern std::string pathProj, MeshName;

extern bool do_remesh;
extern int remesher_iterations;
extern double remesher_aspect_ratio;
extern double alpha;
extern int remesher_termination_delta;

extern double sharp_feature_thr;
extern int feature_erode_dilate;

extern double AlphaSing;
extern double ScaleFact;

void correctOrAbort(const bool ok, const std::string & line, const int linenumber)
{
    if(!ok)
    {
        std::cerr << "[fieldComputation] Malformerd config file |-> " << line << " (" << linenumber << ")" << std::endl;
        std::cerr << "[fieldComputation] ...aborting..." << std::endl;
        abort();
    }
}

void parsePath()
{
    std::string name;
    int pos = pathM.find_last_of('/');
    if (pos == std::string::npos) name = pathM;
    else
    {
        pathProj = pathM.substr(0, pos);
        name = pathM.substr(pos+1);
    }
    
    pos = name.find_last_of('.');
    MeshName = name.substr(0, pos);
}

bool loadConfigFile(const std::string & filename)
{

    // FILE *f=fopen(filename.c_str(),"rt");
    // if (f==NULL)return false;
    // std::cout<<"READ CONFIG FILE"<<std::endl;
    // int IntVar;
    // fscanf(f,"do_remesh %d\n",&IntVar);
    // if (IntVar==0)
    //     do_remesh=false;
    // std::cout<<"do_remesh "<<do_remesh<<std::endl;
    // fscanf(f,"remesh_iterations %d\n",&remesher_iterations);
    // std::cout<<"remesh_iterations "<<remesher_iterations<<std::endl;
    // float remesher_aspect_ratiof;
    // fscanf(f,"remesh_target_aspect_ratio %f\n",&remesher_aspect_ratiof);
    // remesher_aspect_ratio=remesher_aspect_ratiof;
    // std::cout<<"remesher_aspect_ratio "<<remesher_aspect_ratio<<std::endl;
    // fscanf(f,"remesh_termination_delta %d\n",&remesher_termination_delta);
    // std::cout<<"remesh_termination_delta "<<remesher_termination_delta<<std::endl;
    // float sharp_feature_thrf;
    // fscanf(f,"sharp_feature_thr %f\n",&sharp_feature_thrf);
    // sharp_feature_thr=sharp_feature_thrf;
    // std::cout<<"sharp_feature_thr "<<sharp_feature_thr<<std::endl;
    // fscanf(f,"sharp_feature_erode_dilate %d\n",&feature_erode_dilate);
    // std::cout<<"sharp_feature_erode_dilate "<<feature_erode_dilate<<std::endl;
    // float alphaf;
    // fscanf(f,"alpha %f\n",&alphaf);
    // alpha=(double)alphaf;
    // std::cout<<"alpha "<<alpha<<std::endl;
    // fclose(f);
    // std::cout << "[fieldComputation] Successful config import" << std::endl;

    FILE *f=fopen(filename.c_str(),"rt");

    if (f==NULL)return false;

    std::cout<<"READ CONFIG FILE"<<std::endl;

    int IntVar;
    fscanf(f,"do_remesh %d\n",&IntVar);
    if (IntVar==0)
        do_remesh=false;
    else
        do_remesh=true;

    fscanf(f,"sharp_feature_thr %f\n",&sharp_feature_thr);

    fscanf(f,"alpha %lf\n",&AlphaSing);

    fscanf(f,"scaleFact %lf\n",&ScaleFact);

    fclose(f);

    std::cout << "[fieldComputation] Successful config import" << std::endl;

    return true;

}

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    //Use "." as decimal separator
    std::setlocale(LC_NUMERIC, "en_US.UTF-8");

    QWindow dummy;
    QString def_string = QString("GLOBAL fontscaling=%1").arg((int)dummy.devicePixelRatio());
    TwDefine(def_string.toStdString().c_str());
    printf("%s\n",qPrintable(def_string));
    fflush(stdout);

    // Set functions to handle string copy
    TwCopyCDStringToClientFunc(CopyCDStringToClient);
    TwCopyStdStringToClientFunc(CopyStdStringToClient);

    if( !TwInit(TW_OPENGL, NULL) )
    {
        fprintf(stderr, "AntTweakBar initialization failed: %s\n", TwGetLastError());
        return 1;
    }

    //load default config

    loadConfigFile("basic_setup.txt");
    assert(argc>1);
    pathM=std::string(argv[1]);
    parsePath();
    for (size_t i=2;i<argc;i++)
    {
        // if(std::string(argv[i])==std::string("batch"))
        // {
        //     do_batch=true;
        //     continue;
        // }

        std::string pathTest=std::string(argv[i]);
        int position=pathTest.find(".sharp");
        if (position!=-1)
        {
            pathS=pathTest;
            has_features=true;
            continue;
        }
        position=pathTest.find(".fl");
        if (position!=-1)
        {
            pathS=pathTest;
            has_features_fl=true;
            continue;
        }
        position=pathTest.find(".txt");
        if (position!=-1)
        {
           loadConfigFile(pathTest.c_str());
           continue;
        }
        position=pathTest.find(".png");
        if(position!=-1)
        {
            pathT=pathTest;
            has_texture=true;
            continue;
        }
        position=pathTest.find(".jpg");
        if(position!=-1)
        {
            pathT=pathTest;
            has_texture=true;
            continue;
        }
        position=pathTest.find("symm");
        if(position!=-1)
        {
            pathSymm=pathTest;
            has_symm=true;
            continue;   
        }
    }
    //    loadConfigFile("basic_setup.txt");

    //    std::cout << pathM << std::endl;

    //    if ((argc>2)&&(std::string(argv[2])==std::string("batch")))
    //        do_batch=true;
    //    else
    //    {
    //        if (argc>2)
    //        {
    //            std::string pathTest=std::string(argv[2]);
    //            int position=pathTest.find(".sharp");

    //            if (position!=-1)
    //            {
    //                pathS=pathTest;
    //                has_features=true;
    //            }
    //            position=pathTest.find(".fl");
    //            if (position!=-1)
    //            {
    //                pathS=pathTest;
    //                has_features_fl=true;
    //            }
    //        }
    //    }

    //    if ((has_features || has_features_fl)&&(argc>3))
    //    {
    //        if (std::string(argv[3])==std::string("batch"))
    //            do_batch=true;
    //    }
    GLWidget window;

    window.show();
    return app.exec();
}
