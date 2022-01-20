############################ CONFIGURATION ############################

DEFINES += COMISO_FIELD


############################ LIBRARY PATHS ############################

#Libraries
QUADRETOPOLOGY_PATH = $$PWD/quadretopology/
LIBIGL_PATH         = $$PWD/libigl/
VCGLIB_PATH         = $$PWD/vcglib/
GLEW_PATH           = $$PWD/glew/
COMISO_PATH         = $$PWD/CoMISo/
GMM_PATH            = $$PWD/CoMISo/gmm/
EIGEN_PATH          = $$PWD/eigen/

#GUI external libraries (needed only for field_computation and field_tracing projects)
ANTTWEAKBAR_PATH    = ~/AntTweakBar/

#External libraries
BOOST_PATH          = /usr/include/boost/
GUROBI_PATH         = ~/gurobi950/linux64/
GUROBI_COMPILER     = gurobi_g++5.2
GUROBI_LIB          = gurobi95
OPENCV_PATH         = /usr/local/include/opencv4
