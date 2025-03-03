#ifndef EMSOLVER_CONSTANTS_H
#define EMSOLVER_CONSTANTS_H

#include <cmath>
#include <numbers>

///============================================================================
/// Physical constants in SI Units
///============================================================================
inline constexpr double c0   = 299792458.0;         // Speed of light in vacuum
inline constexpr double eps0 = 8.8541878128e-12;    // Free space permittivity
inline constexpr double mu0  = 1.25663706212e-6;    // Free space permeability
inline constexpr double eta0 = 376.730313412;       // Free space impedance

///============================================================================
/// Concept for floating point types
///============================================================================
template<typename T>
concept FloatingPoint = std::is_floating_point_v<T>;

///============================================================================
/// Type alias for position type
///============================================================================
template<FloatingPoint T>
using PositionType = std::pair<T, T>;
using PositionTypeIdx = std::pair<std::size_t, std::size_t>;

#endif // EMSOLVER_CONSTANTS_H
