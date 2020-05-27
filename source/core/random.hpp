#ifndef SAMPLING_RANDOM_HPP
#define SAMPLING_RANDOM_HPP

#include <chrono>
#include <random>

#include "types.hpp"

class RandomNumberGenerator
{
public:
    RandomNumberGenerator()
        : m_engine{compute_seed()}
    {
    }

    explicit RandomNumberGenerator(size_t seed)
        : m_engine{seed}
    {
    }

    Real uniform_real() { return m_dist(m_engine); }

private:
    size_t compute_seed() const
    {
        auto now = std::chrono::system_clock::now();
        return std::chrono::system_clock::to_time_t(now);
    }

    std::mt19937 m_engine;
    std::uniform_real_distribution<Real> m_dist{0.0_r, 1.0_r};
};

#endif // SAMPLING_RANDOM_HPP
