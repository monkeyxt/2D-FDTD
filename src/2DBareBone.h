#ifndef EMSOLVER_2DBAREBONE_H
#define EMSOLVER_2DBAREBONE_H

#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <functional>
#include <numbers>
#include <cassert>
#include <omp.h>

#include "pulses.h"
#include "constants.h"
#include "snapshot.h"

///============================================================================
/// 2D Barebone Solver for the 2D Wave Equation
///============================================================================
template<FloatingPoint T>
class BareBone2D {
public:
    /// Constructor
    explicit BareBone2D(const std::size_t Nx_,
                        const std::size_t Ny_,
                        const std::size_t Npmlx_,
                        const std::size_t Npmly_,
                        const T Lx_,
                        const T Ly_,
                        const T courantFactor_,
                        const T tMax_,
                        Pulse<T>& source_,
                        Snapshot<T>& snapshot_);

    ~BareBone2D() = default;

    BareBone2D() = delete;
    BareBone2D(const BareBone2D&) = delete;
    BareBone2D(BareBone2D&&) = delete;
    BareBone2D& operator=(const BareBone2D&) = delete;
    BareBone2D& operator=(BareBone2D&&) = delete;

    /// Setup material dielectric parameters
    void setupMaterialEps(std::vector<std::vector<T>>&& epsR_);

    /// Setup material magnetic parameters
    void setupMaterialMu(std::vector<std::vector<T>>&& muR_);

    /// Setup PML boundary
    void setupPMLBoundary(const T sigmaMax,
                          const T factor);

    /// Setup material conductivity parameters
    void setupMaterialSigmaE(std::vector<std::vector<T>>&& sigmaEx_,
                             std::vector<std::vector<T>>&& sigmaEy_);

    /// Setup material sigma for magnetic field
    void setupMaterialSigmaM(std::vector<std::vector<T>>&& sigmaMx_,
                             std::vector<std::vector<T>>&& sigmaMy_);

    /// Compute update coefficients
    void computeUpdateCoefficients();

    /// Run simulation
    void runSimulation();

    /// Write snapshots of the fields
    void writeFields();

private:
    /// Simulation parameters
    const std::size_t Nx;      /// number of grid points in x-direction
    const std::size_t Ny;      /// number of grid points in y-direction
    const std::size_t Npmlx;   /// number of PML grid points in x-direction
    const std::size_t Npmly;   /// number of PML grid points in y-direction
    const T Lx;                /// total x-dimension in meters
    const T Ly;                /// total y-dimension in meters
    const T courantFactor;     /// Courant factor
    T dt;                      /// time step in seconds
    const T tMax;              /// maximum time in seconds

    T dx;                      /// grid spacing in meters
    T dy;                      /// grid spacing in meters
    std::size_t Nt;            /// number of time steps

    /// Source parameters
    Pulse<T>& source;
    std::vector<PositionTypeIdx> srcIndex;  /// computed source index

    /// Field arrays
    std::vector<std::vector<T>> Ez;
    std::vector<std::vector<T>> Hy;
    std::vector<std::vector<T>> Hx;

    /// Material parameters
    std::vector<std::vector<T>> epsR;
    std::vector<std::vector<T>> muR;
    std::vector<std::vector<T>> sigmaEx;
    std::vector<std::vector<T>> sigmaEy;
    std::vector<std::vector<T>> sigmaMx;
    std::vector<std::vector<T>> sigmaMy;

    /// Update coefficients
    std::vector<std::vector<T>> chxh;
    std::vector<std::vector<T>> chxe;
    std::vector<std::vector<T>> chyh;
    std::vector<std::vector<T>> chye;
    std::vector<std::vector<T>> ceze;
    std::vector<std::vector<T>> cezh;

    /// Snapshot handler
    Snapshot<T>& snapshot;

    /// Initialization
    void initialize();

    /// Update electric field
    void updateElectricField();

    /// Update magnetic field
    void updateMagneticField();

    /// Update source
    void updateSource(const std::size_t t);
};

template<FloatingPoint T>
BareBone2D<T>::BareBone2D(const std::size_t Nx_,
                          const std::size_t Ny_,
                          const std::size_t Npmlx_,
                          const std::size_t Npmly_,
                          const T Lx_,
                          const T Ly_,
                          const T courantFactor_,
                          const T tMax_,
                          Pulse<T>& source_,
                          Snapshot<T>& snapshot_)
    : Nx(Nx_), 
      Ny(Ny_), 
      Npmlx(Npmlx_), 
      Npmly(Npmly_), 
      Lx(Lx_), 
      Ly(Ly_), 
      courantFactor(courantFactor_), 
      tMax(tMax_), 
      source(source_), 
      snapshot(snapshot_) {
    if (Nx <= 2) {
        throw std::runtime_error("Nx must be greater than 2");
    }
    if (Ny <= 2) {
        throw std::runtime_error("Ny must be greater than 2"); 
    }
    if (Npmlx <= 0) {
        throw std::runtime_error("Npmlx must be greater than 0");
    }
    if (Npmly <= 0) {
        throw std::runtime_error("Npmly must be greater than 0");
    }
    if (Nx <= Npmlx) {
        throw std::runtime_error("Nx must be greater than Npmlx");
    }
    if (Ny <= Npmly) {
        throw std::runtime_error("Ny must be greater than Npmly");
    }

    /// Setup the grid for simulation
    initialize();
}

template<FloatingPoint T>
void BareBone2D<T>::initialize() {
    dx = Lx / Nx;
    dy = Ly / Ny;

    T cfl2d = 1.0 / (c0 * std::sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy)));
    dt = courantFactor * cfl2d;
    Nt = static_cast<std::size_t>(std::floor(tMax / dt));

    /// Zero initial fields
    Ez.resize(Nx, std::vector<T>(Ny, 0.0));
    Hy.resize(Nx, std::vector<T>(Ny, 0.0));
    Hx.resize(Nx, std::vector<T>(Ny, 0.0));
    
    /// Setup material parameters
    epsR.resize(Nx, std::vector<T>(Ny, 1.0));
    muR.resize(Nx, std::vector<T>(Ny, 1.0));
    
    /// Setup conductivity parameters
    sigmaEx.resize(Nx, std::vector<T>(Ny, 0.0));
    sigmaEy.resize(Nx, std::vector<T>(Ny, 0.0));
    sigmaMx.resize(Nx, std::vector<T>(Ny, 0.0));
    sigmaMy.resize(Nx, std::vector<T>(Ny, 0.0));

    /// Setup update coefficients
    chxh.resize(Nx, std::vector<T>(Ny, 0.0));
    chxe.resize(Nx, std::vector<T>(Ny, 0.0));
    chyh.resize(Nx, std::vector<T>(Ny, 0.0));
    chye.resize(Nx, std::vector<T>(Ny, 0.0));
    ceze.resize(Nx, std::vector<T>(Ny, 0.0));
    cezh.resize(Nx, std::vector<T>(Ny, 0.0));

    /// Calculate location of the source index
    srcIndex.clear();
    if (auto* planeWave = dynamic_cast<PlaneWavePulse<T>*>(&source)) {
        std::size_t xIndex = static_cast<std::size_t>(planeWave->getSourcePositions()[0].first / dx);
        if (xIndex >= Nx) {
            throw std::runtime_error("Source position is out of bounds");
        }
        // For plane wave, create source points along entire y-axis at x position
        for (std::size_t j = 0; j < Ny; j++) {
            srcIndex.push_back({xIndex, j});
        }
    }
    else if (auto* pointSource = dynamic_cast<PointSourcePulse<T>*>(&source)) {
        // For point source, use the exact position
        for (const auto& pos : pointSource->getSourcePositions()) {
            std::size_t xIdx = static_cast<std::size_t>(pos.first / dx);
            std::size_t yIdx = static_cast<std::size_t>(pos.second / dy);
            if (xIdx >= Nx || yIdx >= Ny) {
                throw std::runtime_error("Source position is out of bounds");
            }
            srcIndex.push_back({xIdx, yIdx});
        }
    }

    std::cout << "Simulation parameters:" << std::endl;
    std::cout << "--------------------------------" << std::endl;
    std::cout << "Lx = " << Lx << " m" << std::endl;
    std::cout << "Ly = " << Ly << " m" << std::endl;
    std::cout << "dx = " << dx << " m" << std::endl;
    std::cout << "dy = " << dy << " m" << std::endl;
    std::cout << "dt = " << dt << " s" << std::endl;
    std::cout << "steps = " << Nt << std::endl;
    std::cout << "--------------------------------" << std::endl;
}

template<FloatingPoint T>
void BareBone2D<T>::setupMaterialEps(std::vector<std::vector<T>>&& epsR_) {
    epsR = std::move(epsR_);
}

template<FloatingPoint T>
void BareBone2D<T>::setupMaterialMu(std::vector<std::vector<T>>&& muR_) {
    muR = std::move(muR_);
}

template<FloatingPoint T>
void BareBone2D<T>::setupMaterialSigmaE(std::vector<std::vector<T>>&& sigmaEx_,
                                        std::vector<std::vector<T>>&& sigmaEy_) {
    sigmaEx = std::move(sigmaEx_);
    sigmaEy = std::move(sigmaEy_);
}

template<FloatingPoint T>
void BareBone2D<T>::setupMaterialSigmaM(std::vector<std::vector<T>>&& sigmaMx_,
                                        std::vector<std::vector<T>>&& sigmaMy_) {
    sigmaMx = std::move(sigmaMx_);
    sigmaMy = std::move(sigmaMy_);
}

template<FloatingPoint T>
void BareBone2D<T>::setupPMLBoundary(const T sigmaMax,
                                     const T factor) {
    // Left PML region - x-direction conductivity only
    for(std::size_t i = 0; i < Npmlx; i++) {
        T x = static_cast<T>(Npmlx - i) / Npmlx;
        T sigma = sigmaMax * std::pow(x, factor);
        for(std::size_t j = 0; j < Ny; j++) {
            sigmaEx[i][j] = sigma;
            sigmaMx[i][j] = sigma * mu0/eps0;
        }
    }

    // Right PML region - x-direction conductivity only
    for(std::size_t i = Nx - Npmlx; i < Nx; i++) {
        T x = static_cast<T>(i - (Nx - Npmlx)) / Npmlx;
        T sigma = sigmaMax * std::pow(x, factor);
        for(std::size_t j = 0; j < Ny; j++) {
            sigmaEx[i][j] = sigma;
            sigmaMx[i][j] = sigma * mu0/eps0;
        }
    }

    // Bottom PML region - y-direction conductivity only
    for(std::size_t j = 0; j < Npmly; j++) {
        T y = static_cast<T>(Npmly - j) / Npmly;
        T sigma = sigmaMax * std::pow(y, factor);
        for(std::size_t i = 0; i < Nx; i++) {
            sigmaEy[i][j] = sigma;
            sigmaMy[i][j] = sigma * mu0/eps0;
        }
    }

    // Top PML region - y-direction conductivity only
    for(std::size_t j = Ny - Npmly; j < Ny; j++) {
        T y = static_cast<T>(j - (Ny - Npmly)) / Npmly;
        T sigma = sigmaMax * std::pow(y, factor);
        for(std::size_t i = 0; i < Nx; i++) {
            sigmaEy[i][j] = sigma;
            sigmaMy[i][j] = sigma * mu0/eps0;
        }
    }

    // Corner regions - combine x and y conductivities
    // This creates a smoother transition in the corners
    for(std::size_t i = 0; i < Npmlx; i++) {
        for(std::size_t j = 0; j < Npmly; j++) {
            T x = static_cast<T>(Npmlx - i) / Npmlx;
            T y = static_cast<T>(Npmly - j) / Npmly;
            T sigmaX = sigmaMax * std::pow(x, factor);
            T sigmaY = sigmaMax * std::pow(y, factor);
            sigmaEx[i][j] = sigmaX * std::sqrt(x*x + y*y) / std::sqrt(2.0);
            sigmaEy[i][j] = sigmaY * std::sqrt(x*x + y*y) / std::sqrt(2.0);
            sigmaMx[i][j] = sigmaEx[i][j] * mu0/eps0;
            sigmaMy[i][j] = sigmaEy[i][j] * mu0/eps0;
        }
    }
}

template<FloatingPoint T>
void BareBone2D<T>::computeUpdateCoefficients() {
    /// Conduction coefficients for electric field
    #pragma omp simd
    for(std::size_t i = 0; i < Nx; i++) {
        for(std::size_t j = 0; j < Ny; j++) {
            T denomE = 1.0 + (sigmaEx[i][j] + sigmaEy[i][j]) * dt / (2.0 * epsR[i][j] * eps0);
            T numerE = 1.0 - (sigmaEx[i][j] + sigmaEy[i][j]) * dt / (2.0 * epsR[i][j] * eps0);
            ceze[i][j] = numerE / denomE;
            cezh[i][j] = (dt / (epsR[i][j] * eps0 * dx)) / denomE;
        } 
    }

    /// Conduction coefficients for magnetic field
    #pragma omp simd
    for(std::size_t i = 0; i < Nx; i++) {
        for(std::size_t j = 0; j < Ny - 1; j++) {
            T denomH = 1.0 + (sigmaMx[i][j] + sigmaMy[i][j]) * dt / (2.0 * muR[i][j] * mu0);
            T numerH = 1.0 - (sigmaMx[i][j] + sigmaMy[i][j]) * dt / (2.0 * muR[i][j] * mu0);
            chxh[i][j] = numerH / denomH;
            chxe[i][j] = (dt / (muR[i][j] * mu0 * dx)) / denomH;
        } 
    }

    #pragma omp simd
    for(std::size_t i = 0; i < Nx - 1; i++) {
        for(std::size_t j = 0; j < Ny; j++) {
            T denomH = 1.0 + (sigmaMx[i][j] + sigmaMy[i][j]) * dt / (2.0 * muR[i][j] * mu0);
            T numerH = 1.0 - (sigmaMx[i][j] + sigmaMy[i][j]) * dt / (2.0 * muR[i][j] * mu0);
            chyh[i][j] = numerH / denomH;
            chye[i][j] = (dt / (muR[i][j] * mu0 * dx)) / denomH;
        }
    }
}

template<FloatingPoint T>
void BareBone2D<T>::updateElectricField() {
    #pragma omp simd
    for(std::size_t i = 1; i < Nx - 1; i++) {
        for(std::size_t j = 1; j < Ny - 1; j++) {
            Ez[i][j] = ceze[i][j] * Ez[i][j] 
                + cezh[i][j] * ((Hy[i][j] - Hy[i-1][j]) - (Hx[i][j] - Hx[i][j-1]));
        }
    }
}

template<FloatingPoint T>
void BareBone2D<T>::updateMagneticField() {
    #pragma omp simd
    for(std::size_t i = 0; i < Nx; i++) {
        for(std::size_t j = 0; j < Ny - 1; j++) {
            Hx[i][j] = chxh[i][j] * Hx[i][j] - chxe[i][j] * (Ez[i][j+1] - Ez[i][j]);
        }
    }
    #pragma omp simd
    for(std::size_t i = 0; i < Nx - 1; i++) {
        for(std::size_t j = 0; j < Ny; j++) {
            Hy[i][j] = chyh[i][j] * Hy[i][j] + chye[i][j] * (Ez[i+1][j] - Ez[i][j]);
        }
    }
}

template<FloatingPoint T>
void BareBone2D<T>::updateSource(const std::size_t t) {
    for(const auto& src : srcIndex) {
        auto pulseInc = source.computePulse(t * dt, src);
        Ez[src.first][src.second] += pulseInc;
    }
}

template<FloatingPoint T>
void BareBone2D<T>::runSimulation() {
    /// Run simulation for Nt time steps
    for (std::size_t t = 0; t < Nt; t++) {
        updateMagneticField();
        updateElectricField();
        updateSource(t);
    }
    writeFields();
}

template<FloatingPoint T>
void BareBone2D<T>::writeFields() {
    snapshot.write({"Ez", "Hy", "Hx"}, Ez, Hy, Hx);
}

#endif // EMSOLVER_2DBAREBONE_H