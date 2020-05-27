#ifndef SAMPLING_TYPES_HPP
#define SAMPLING_TYPES_HPP

#include <cassert>
#include <complex>
#include <ostream>

#include <valarray>
#include <vector>

using Real = float;
using Complex = std::complex<Real>;

template <typename T> using Array = std::vector<T>;
using RealArray = Array<Real>;
using ComplexArray = Array<Complex>;

inline Real operator""_r(long double v) { return static_cast<Real>(v); }

#endif // SAMPLING_TYPES_HPP
