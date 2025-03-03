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
    src/2DBareBone.h
    src/snapshot.h
    src/constants.h
    src/pulses.h
)

# Create header-only library target
add_library(EMSolver INTERFACE)
target_include_directories(EMSolver INTERFACE 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include>
)

# Add examples directory and compile FDTD1D examples
add_executable(2DPointSource examples/2DPointSource.cpp) 
add_executable(waveguideTE1 examples/waveguideTE1.cpp)
add_executable(waveguideTE5 examples/waveguideTE5.cpp)
add_executable(waveguideAttenuation examples/waveguideAttenuation.cpp)
add_executable(2DPointSourcePML examples/2DPointSourcePML.cpp)

#target_link_libraries(1DPointSource PRIVATE EMSolver)
#target_link_libraries(2DPointSource PRIVATE EMSolver)
#target_link_libraries(2DPlaneWave PRIVATE EMSolver)
#target_link_libraries(2DPlaneWaveWLoss PRIVATE EMSolver)
target_link_libraries(waveguideTE1 PRIVATE EMSolver)
target_link_libraries(waveguideTE5 PRIVATE EMSolver)
target_link_libraries(waveguideAttenuation PRIVATE EMSolver)
target_link_libraries(2DPointSource PRIVATE EMSolver)
target_link_libraries(2DPointSourcePML PRIVATE EMSolver)

# Plane Wave Examples
add_executable(2DPlaneWave examples/2DPlaneWave.cpp)
add_executable(2DPlaneWaveWLoss examples/2DPlaneWaveWLoss.cpp)
add_executable(2DPlaneWaveIncident examples/2DPlaneWaveIncident.cpp)
add_executable(2DPlaneWaveIncidentTotal examples/2DPlaneWaveIncidentTotal.cpp)

target_link_libraries(2DPlaneWave PRIVATE EMSolver)
target_link_libraries(2DPlaneWaveWLoss PRIVATE EMSolver)
target_link_libraries(2DPlaneWaveIncident PRIVATE EMSolver)
target_link_libraries(2DPlaneWaveIncidentTotal PRIVATE EMSolver)

# Installation rules
install(TARGETS EMSolver 2DPointSource waveguideTE1 waveguideTE5 waveguideAttenuation 2DPointSourcePML 2DPlaneWaveIncident 2DPlaneWaveIncidentTotal 2DPlaneWave 2DPlaneWaveWLoss
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
