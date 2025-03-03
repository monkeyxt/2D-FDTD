#include <iostream>
#include <string>

#include "1DBareBone.h"
#include "snapshot.h"

///============================================================================
/// Simple 1D FDTD simulation with a point source, but with a PML boundary, 
/// to absorb the wave at the end of the domain.
///============================================================================
int main(int argc, char* argv[]) {
    std::size_t Nx = 1000;
    std::size_t Npml = 20;

    double Lx = 0.10;
    double CourantFactor = 0.99;
    double tMax = 400e-12;
    double fSrc = 30e9;
    double sourcePosition = 0.05;

    std::string dataDir = "../data/1DPointSourcePML"; // Default directory
    if (argc > 1) {
        dataDir = argv[1];
    }
    Snapshot<double> snapshot(dataDir);

    PointSourcePulse<double> source(fSrc, tMax, {{sourcePosition, 0.0}});
    BareBone1D<double> simulation(Nx, 
                                  Npml, 
                                  Lx, 
                                  CourantFactor, 
                                  tMax, 
                                  source, 
                                  snapshot);
    
    // PML parameters
    const double R0 = 1e-8;  // Theoretical reflection coefficient
    const int m = 3;         // PML order
    const double sigma_max 
        = -(m + 1) * eps0 * c0 * std::log(R0) / (2 * Npml * Lx / Nx);
    std::cout << "sigma_max: " << sigma_max << std::endl;
    simulation.setupPMLBoundary(sigma_max, m);
    simulation.setupMaterialEps(std::vector<double>{0.04, 4.0, 1.0});
    simulation.computeUpdateCoefficients();

    std::cout << "Running simulation..." << std::endl;
    simulation.runSimulation();
    std::cout << "Simulation completed." << std::endl;
    return 0;
}