#ifndef SAMPLING_SAMPLER_HPP
#define SAMPLING_SAMPLER_HPP

#include <core/types.hpp>

struct Sample
{
    Real x;
    Real y;
};
using SampleArray = Array<Sample>;

///
/// Sampler interface
///
/// Generates samples within [0.0 - 1.0) domain.
///
template <typename Impl> class Sampler
{
public:
    Sampler(const Sampler<Impl>&) = delete;
    Sampler(Sampler<Impl>&&) = delete;

    Sampler<Impl>& operator=(const Sampler<Impl>&) = delete;
    Sampler<Impl>& operator=(Sampler<Impl>&&) = delete;

    SampleArray request_samples(const size_t num_samples) { return impl()->request_samples(num_samples); }

protected:
    Sampler() = default;

private:
    // cast to implementation
    Impl* impl() { return static_cast<Impl*>(this); }
};

///
/// Helper to write Sample to out stream
///
inline std::ostream& operator<<(std::ostream& os, const Sample& sample)
{
    os << sample.x << " " << sample.y;
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const SampleArray& sample_array)
{
    for (const Sample& sample : sample_array)
    {
        os << sample << '\n';
    }
    return os;
}

#endif // SAMPLING_SAMPLER_HPP
