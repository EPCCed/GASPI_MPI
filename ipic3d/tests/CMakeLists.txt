#This release was prepared by Dana Akhmetova <danaak@kth.se>/<danieka@gmail.com> on behalf of the INTERTWinE European Exascale Project <http://www.intertwine-project.eu>
cmake_minimum_required(VERSION 2.8.8)
project (ipic_tests)
find_package(MPI REQUIRED)
include_directories(
    ../include
    ${MPI_INCLUDE_PATH}
)
add_executable(
  testspawn
  testspawn.cpp)
target_link_libraries(
    testspawn
    ${MPI_LIBRARIES}
    #${CMAKE_DL_LIBS}
)

