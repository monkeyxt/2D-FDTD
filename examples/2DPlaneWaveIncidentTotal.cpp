#include "2DBareBone.h"
#include "snapshot.h"
#include "pulses.h"
///============================================================================
/// Simple 2D FDTD simulation with a plane wave source and incident on a 
/// dielectric interface. 
///============================================================================
int main(int argc, char* argv[]) {
    constexpr std::size_t Nx = 2000;
    constexpr std::size_t Ny = 2000;
    constexpr std::size_t Npmlx = 50;
    constexpr std::size_t Npmly = 50;

    const double Lx = 0.008;
    const double Ly = 0.008;
    const double CourantFactor = 0.99;
    const double tMax = 10e-12;
    const double fSrc = 500e9;
    const double sourcePositionx = 0.002;
    const double sourcePositiony = 0.002;
    const PositionType sourcePosition = {sourcePositionx, sourcePositiony};

    const double theta = 30 * std::numbers::pi / 180;

    std::string dataDir1 = "../data/2DPlaneWaveIncidentTotal/1";
    std::string dataDir2 = "../data/2DPlaneWaveIncidentTotal/2";
    std::string dataDir3 = "../data/2DPlaneWaveIncidentTotal/3";
    std::string dataDir4 = "../data/2DPlaneWaveIncidentTotal/4";
    if (argc > 1) {
        dataDir1 = argv[1];
        dataDir2 = argv[2];
        dataDir3 = argv[3];
        dataDir4 = argv[4];
    }
    Snapshot<double> snapshot1(dataDir1);
    Snapshot<double> snapshot2(dataDir2);
    Snapshot<double> snapshot3(dataDir3);
    Snapshot<double> snapshot4(dataDir4);

    /// Initialize a plane wave source
    PlaneWavePulse<double> source(fSrc, tMax, {sourcePosition});

    BareBone2D<double> simulation1(
        Nx, Ny, Npmlx, Npmly, Lx, Ly, CourantFactor, tMax, source, snapshot1
    );
    // Setup PML boundary
    const double m = 10.0;  // PML order
    const double R0 = 1e-8; // Theoretical reflection coefficient
    const double sigmaMax 
        = -(m + 1) * eps0 * c0 * std::log(R0) / (2.0 * std::min(Lx / Nx, Ly / Ny));
    std::cout << "sigmaMax: " << sigmaMax << std::endl;
    simulation1.setupPMLBoundary(sigmaMax, m);
    simulation1.computeUpdateCoefficients();
    std::cout << "Running simulation 1..." << std::endl;
    simulation1.runSimulation();
    std::cout << "Simulation 1 completed." << std::endl;


    BareBone2D<double> simulation2(
        Nx, Ny, Npmlx, Npmly, Lx, Ly, CourantFactor, tMax, source, snapshot2
    );

    std::vector<std::vector<double>> epsR2(Nx, std::vector<double>(Ny, 1.0));
    for (std::size_t i = 0; i < Nx; i++) {
        for (std::size_t j = 0; j < Ny; j++) {
            // Create a diagonal interface at 30 degrees
            // Using y = tan(30°) * (x - x0) + y0 to define the interface line
            double x = i * Lx / Nx;
            double y = j * Ly / Ny;
            double x0 = 0.002;  // x-position of interface at y=0
            double y0 = 0.000;  // y-position offset
            if (y > std::tan(50.0 * std::numbers::pi / 180.0) * (x - x0) + y0) {
                epsR2[i][j] = 2.4;
            }
        }
    }
    simulation2.setupMaterialEps(std::move(epsR2));
    simulation2.setupPMLBoundary(sigmaMax, m);
    simulation2.computeUpdateCoefficients();
    std::cout << "Running simulation 2..." << std::endl;
    simulation2.runSimulation();
    std::cout << "Simulation 2 completed." << std::endl;

    
    BareBone2D<double> simulation3(
        Nx, Ny, Npmlx, Npmly, Lx, Ly, CourantFactor, tMax, source, snapshot3
    );

    std::vector<std::vector<double>> epsR3(Nx, std::vector<double>(Ny, 1.0));
    for (std::size_t i = 0; i < Nx; i++) {
        for (std::size_t j = 0; j < Ny; j++) {
            // Create a diagonal interface at 30 degrees
            // Using y = tan(30°) * (x - x0) + y0 to define the interface line
            double x = i * Lx / Nx;
            double y = j * Ly / Ny;
            double x0 = 0.002;  // x-position of interface at y=0
            double y0 = 0.000;  // y-position offset
            if (y > std::tan(60.0 * std::numbers::pi / 180.0) * (x - x0) + y0) {
                epsR3[i][j] = 2.4;
            }
        }
    }
    simulation3.setupMaterialEps(std::move(epsR3));
    simulation3.setupPMLBoundary(sigmaMax, m);
    simulation3.computeUpdateCoefficients();
    std::cout << "Running simulation 3..." << std::endl;
    simulation3.runSimulation();
    std::cout << "Simulation 3 completed." << std::endl;

    
    BareBone2D<double> simulation4(
        Nx, Ny, Npmlx, Npmly, Lx, Ly, CourantFactor, tMax, source, snapshot4
    );

    std::vector<std::vector<double>> epsR4(Nx, std::vector<double>(Ny, 1.0));
    for (std::size_t i = 0; i < Nx; i++) {
        for (std::size_t j = 0; j < Ny; j++) {
            // Create a diagonal interface at 30 degrees
            // Using y = tan(30°) * (x - x0) + y0 to define the interface line
            double x = i * Lx / Nx;
            double y = j * Ly / Ny;
            double x0 = 0.002;  // x-position of interface at y=0
            double y0 = 0.000;  // y-position offset
            if (y > std::tan(40.0 * std::numbers::pi / 180.0) * (x - x0) + y0) {
                epsR4[i][j] = 2.4;
            }
        }
    }
    simulation4.setupMaterialEps(std::move(epsR4));
    simulation4.setupPMLBoundary(sigmaMax, m);
    simulation4.computeUpdateCoefficients();
    std::cout << "Running simulation 4..." << std::endl;
    simulation4.runSimulation();
    std::cout << "Simulation 4 completed." << std::endl;

    
    
    return 0;
}