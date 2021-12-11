project("collective_migration")
cmake_minimum_required(VERSION 3.1)

find_package(AMDIS REQUIRED SEQUENTIAL)
find_package(mpi14 REQUIRED)

add_executable(random_model 
  src/AdaptInstationarySeq.cc 
  include/ContourWriter.cc 
  src/Coupling_phase.h
  src/box_migration.cc)
target_include_directories(random_model PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(random_model PRIVATE AMDiS mpi14)
target_compile_definitions(random_model PRIVATE HAVE_MPI)