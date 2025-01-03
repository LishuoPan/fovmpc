cmake_minimum_required(VERSION 3.16)
project(mpc_cbf)
set_property(GLOBAL PROPERTY mpc_cbf_ADDED "ON")

set(CMAKE_CXX_STANDARD 17)

# add local library
get_property(is_mpc_ADDED GLOBAL PROPERTY mpc_ADDED)
if(NOT is_mpc_ADDED)
    add_subdirectory(../mpc ${CMAKE_CURRENT_BINARY_DIR}/mpc)
endif()

get_property(is_cbf_ADDED GLOBAL PROPERTY cbf_ADDED)
if(NOT is_cbf_ADDED)
    add_subdirectory(../cbf ${CMAKE_CURRENT_BINARY_DIR}/cbf)
endif()

get_property(is_pf_ADDED GLOBAL PROPERTY pf_ADDED)
if(NOT is_pf_ADDED)
    add_subdirectory(../particle_filter ${CMAKE_CURRENT_BINARY_DIR}/particle_filter)
endif()

# add libraries
add_library(mpc_cbf_optimization_PiecewiseBezierMPCCBFQPOperations SHARED src/optimization/PiecewiseBezierMPCCBFQPOperations.cpp)
target_include_directories(mpc_cbf_optimization_PiecewiseBezierMPCCBFQPOperations PUBLIC include)
target_link_libraries(mpc_cbf_optimization_PiecewiseBezierMPCCBFQPOperations PUBLIC
        mpc_optimization_PiecewiseBeizerMPCQPOperations
        cbf_detail_FovCBF)

add_library(mpc_cbf_optimization_PiecewiseBezierMPCCBFQPGenerator SHARED src/optimization/PiecewiseBezierMPCCBFQPGenerator.cpp)
target_include_directories(mpc_cbf_optimization_PiecewiseBezierMPCCBFQPGenerator PUBLIC include)
target_link_libraries(mpc_cbf_optimization_PiecewiseBezierMPCCBFQPGenerator PUBLIC
        mpc_cbf_optimization_PiecewiseBezierMPCCBFQPOperations
        mpc_optimization_PiecewiseBezierMPCQPGenerator)

add_library(mpc_cbf_controller_BezierMPCCBF SHARED src/controller/BezierMPCCBF.cpp)
target_include_directories(mpc_cbf_controller_BezierMPCCBF PUBLIC include)
target_link_libraries(mpc_cbf_controller_BezierMPCCBF PUBLIC
        mpc_cbf_optimization_PiecewiseBezierMPCCBFQPGenerator
        qpcpp_solvers_CPLEX)

add_library(mpc_cbf_controller_BezierIMPCCBF SHARED src/controller/BezierIMPCCBF.cpp)
target_include_directories(mpc_cbf_controller_BezierIMPCCBF PUBLIC include)
target_link_libraries(mpc_cbf_controller_BezierIMPCCBF PUBLIC
        mpc_cbf_optimization_PiecewiseBezierMPCCBFQPGenerator
        qpcpp_solvers_CPLEX
        separating_hyperplanes_Voronoi)

# add examples
add_executable(mpc_cbf_examples_BezierMPCCBFXYYaw_example examples/BezierMPCCBFXYYaw_example.cpp)
target_include_directories(mpc_cbf_examples_BezierMPCCBFXYYaw_example PUBLIC include ../../third_party/json/single_include)
target_link_libraries(mpc_cbf_examples_BezierMPCCBFXYYaw_example PUBLIC
        mpc_cbf_controller_BezierMPCCBF
        model_DoubleIntegratorXYYaw
        math_collision_shapes_AlignedBoxCollisionShape)

add_executable(mpc_cbf_examples_BezierIMPCCBFXYYaw_example examples/BezierIMPCCBFXYYaw_example.cpp)
target_include_directories(mpc_cbf_examples_BezierIMPCCBFXYYaw_example PUBLIC include ../../third_party/json/single_include)
target_link_libraries(mpc_cbf_examples_BezierIMPCCBFXYYaw_example PUBLIC
        mpc_cbf_controller_BezierIMPCCBF
        model_DoubleIntegratorXYYaw
        math_collision_shapes_AlignedBoxCollisionShape)

add_executable(mpc_cbf_pf_example examples/pf_ex.cpp)
target_include_directories(mpc_cbf_pf_example PUBLIC include ../../third_party/json/single_include)
target_link_libraries(mpc_cbf_pf_example PUBLIC
        mpc_cbf_controller_BezierIMPCCBF
        model_DoubleIntegratorXYYaw
        math_collision_shapes_AlignedBoxCollisionShape
        ParticleFilter)