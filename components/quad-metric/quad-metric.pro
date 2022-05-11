CONFIG += c++11

TARGET = quad-metric

CONFIG += console

CONFIG(release, debug|release): DEFINES += NDEBUG

include(../../libs/libs.pri)

#vcglib
INCLUDEPATH += $$VCGLIB_PATH

#eigen
INCLUDEPATH += $$EIGEN_PATH

SOURCES += main.cpp

HEADERS += mesh_def.h

#Parallel computation (just in release)
unix:!mac {
    QMAKE_CXXFLAGS += -fopenmp
    LIBS += -fopenmp
}