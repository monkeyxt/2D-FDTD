cmake_minimum_required(VERSION 3.10)
project(EMSolver CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Add source files
set(SOURCES
    src/1DBareBone.h
    src/snapshot.h
    src/constants.h
)

# Create header-only library target
add_library(EMSolver INTERFACE)
target_include_directories(EMSolver INTERFACE 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include>
)

# Add examples directory and compile FDTD1D examples
add_executable(FDTD1D_1 examples/FDTD1D_1.cpp)
add_executable(FDTD1D_2 examples/FDTD1D_2.cpp) 
add_executable(FDTD1D_3 examples/FDTD1D_3.cpp)
add_executable(FDTD2D_1 examples/FDTD2D_1.cpp)
add_executable(FDTD2D_2 examples/FDTD2D_2.cpp)
add_executable(FDTD2D_3 examples/FDTD2D_3.cpp)

target_link_libraries(FDTD1D_1 PRIVATE EMSolver)
target_link_libraries(FDTD1D_2 PRIVATE EMSolver)
target_link_libraries(FDTD1D_3 PRIVATE EMSolver)
target_link_libraries(FDTD2D_1 PRIVATE EMSolver)
target_link_libraries(FDTD2D_2 PRIVATE EMSolver)
target_link_libraries(FDTD2D_3 PRIVATE EMSolver)

# Installation rules
install(TARGETS EMSolver FDTD1D_1 FDTD1D_2 FDTD1D_3 FDTD2D_1 FDTD2D_2 FDTD2D_3
    EXPORT EMSolverTargets
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

install(FILES ${SOURCES}
    DESTINATION include/EMSolver
)

install(EXPORT EMSolverTargets
    FILE EMSolverTargets.cmake
    NAMESPACE EMSolver::
    DESTINATION lib/cmake/EMSolver
)

# Add tests directory if it exists
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    enable_testing()
    add_subdirectory(tests)
endif()
