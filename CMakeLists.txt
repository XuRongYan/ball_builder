cmake_minimum_required(VERSION 3.16)
project(ball_builder)

set(CMAKE_CXX_STANDARD 17)

#------------Eigen--------------
find_package(Eigen3 CONFIG REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

#----------libigl--------------------
find_package(libigl CONFIG REQUIRED)

#----------spdlog-----------------
find_package(spdlog CONFIG REQUIRED)

#---------SURFACE MESH----------
include_directories("$ENV{Surface_Mesh_DIR}\\include")

aux_source_directory(. SRC)

add_executable(ball_builder ${SRC})

target_link_libraries(${PROJECT_NAME}
        Eigen3::Eigen
        spdlog::spdlog
        igl::core
        igl::common)

target_link_libraries(${PROJECT_NAME} "$ENV{Surface_Mesh_DIR}\\lib\\surface_mesh_release_double.lib")