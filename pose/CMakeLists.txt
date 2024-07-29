cmake_minimum_required( VERSION 2.8 )
project( g2o_test )

set( CMAKE_BUILD_TYPE "Debug" )
set( CMAKE_CXX_FLAGS "-std=c++11 -O3" )

list( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules )

find_package( G2O REQUIRED )
find_package( Sophus REQUIRED )
# find_package( Cholmod REQUIRED )

include_directories( 
    ${G2O_INCLUDE_DIRS}
    ${Sophus_INCLUDE_DIRS}
    # ${Cholmod_INCLUDE_DIRS}
    "/usr/include/eigen3"
    "include/"
    # "/usr/include/suitesparse"
)

SET(G2O_LIBS g2o_cli g2o_ext_freeglut_minimal g2o_simulator g2o_solver_slam2d_linear 
    g2o_types_icp g2o_types_slam2d g2o_types_sba g2o_types_slam3d g2o_core g2o_interface 
    g2o_solver_csparse g2o_solver_structure_only g2o_csparse_extension g2o_opengl_helper g2o_solver_dense 
    g2o_stuff g2o_types_sclam2d g2o_parser g2o_solver_pcg g2o_types_data g2o_types_sim3 cxsparse )

add_executable(pose_graph_left src/pose_graph_left.cpp)
target_link_libraries( pose_graph_left 
    ${G2O_LIBS}
    ${Sophus_LIBS}
)

add_executable(pose_graph_right src/pose_graph_right.cpp)
target_link_libraries( pose_graph_right 
    ${G2O_LIBS}
    ${Sophus_LIBS}
)

add_executable(pose_graph_so3_t src/pose_graph_so3_t.cpp)
target_link_libraries( pose_graph_so3_t 
    ${G2O_LIBS}
    ${Sophus_LIBS}
    # ${Cholmod_LIBRARIES}
)

add_executable(pose_graph_so2_t src/pose_graph_so2_t.cpp)
target_link_libraries( pose_graph_so2_t 
    ${G2O_LIBS}
    ${Sophus_LIBS}
)
# message("cmake_module_path: " ${CMAKE_MODULE_PATH})
# message("sophus include path: ${Sophus_INCLUDE_DIRS}")
# message("sophus library path: ${Sophus_LIBRARIES}")

add_executable(pose_graph_SE3 src/pose_graph_SE3.cpp)
target_link_libraries( pose_graph_SE3 
    ${G2O_LIBS}
    # ${Sophus_LIBS}
)