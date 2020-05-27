#ifndef SAMPLING_ELEMENTARY_INTERVAL_H
#define SAMPLING_ELEMENTARY_INTERVAL_H

#include "sampler.hpp"

#include <core/random.hpp>

class ElementaryIntervalSampler : public Sampler<ElementaryIntervalSampler>
{
public:
    ElementaryIntervalSampler()
        : Sampler<ElementaryIntervalSampler>{}
    {
    }


    SampleArray request_samples(const size_t num_samples)
    {
        // we need to be sure we work in power of two strata, otherwise fractions appear in elementary intervals.
        const double lx_ly = std::ceil(std::log2(num_samples));
        const size_t lx = std::ceil(lx_ly / 2);
        const size_t ly = lx; // enforce equal strata in x and y

        const size_t nx = std::pow(2, lx);
        const size_t ny = std::pow(2, ly);
        const size_t n = nx * ny;

        // allocate output samples storage
        SampleArray samples;
        samples.resize(n);

        // number of all elementary intervals is n + n*log2(n), but it easier to split them into two separate occupancy
        // intervals, n*log2(nx*ny), can be broken down
        Array<bool> x_occupancy(n * (lx_ly + 1), false);
        Array<bool> y_occupancy(n * (lx_ly + 1), false);

        std::vector<NodeOffset> x_offsets;
        std::vector<NodeOffset> y_offsets;
        x_offsets.reserve(nx);
        y_offsets.reserve(ny);

        const Real step = 1.0_r / n;

        size_t index{};
        for (size_t y{}; y < ny; ++y)
        {
            for (size_t x{}; x < nx; ++x)
            {
                // 2) initialize binary trees
                initialize_x_tree(x_occupancy, 0, x, y, nx, ny);
                initialize_y_tree(y_occupancy, 0, x, y, nx, ny);

                // get valid offsets
                x_offsets.clear();
                y_offsets.clear();

                get_valid_x_offsets(x_occupancy, 0, x, y, nx, ny, x_offsets);
                get_valid_y_offsets(y_occupancy, 0, x, y, nx, ny, y_offsets);

                assert(x_offsets.size() > 0 && x_offsets.size() <= nx);
                assert(y_offsets.size() > 0 && y_offsets.size() <= ny);

                // pick a random strata
                const size_t x_random_index = std::floor(m_rng.uniform_real() * x_offsets.size());
                const size_t y_random_index = std::floor(m_rng.uniform_real() * y_offsets.size());

                const NodeOffset& x_node_offset = x_offsets[x_random_index];
                const NodeOffset& y_node_offset = y_offsets[y_random_index];

                samples[index].x = x_node_offset.offset * step + m_rng.uniform_real()*step;
                samples[index].y = y_node_offset.offset * step+ m_rng.uniform_real()*step;

                // update tree, mark visited nodes as occupied
                x_occupancy[x_node_offset.node] = true;
                y_occupancy[y_node_offset.node] = true;

                ++index;
            }
        }
        return samples;
    }

private:
    struct NodeOffset
    {
        size_t node;
        size_t offset;
    };

    void initialize_x_tree(Array<bool>& x_occupancy, size_t offset, //
                           const size_t x, const size_t y,          //
                           const size_t nx, const size_t ny) const
    {
        assert(offset < x_occupancy.size());
        assert(x < nx);
        assert(y < ny);

        // check if node has been already filled
        const size_t node = offset + nx * y + x;
        assert(node < x_occupancy.size());

        if (x_occupancy[node])
        {
            return;
        }

        // if not continue drilling the hierarchy, to mark nodes occupied
        bool is_leaf = ny == 1;
        if (!is_leaf)
        {
            // compute new intervals
            const size_t new_x{2 * x};
            const size_t new_y = std::floor(y / 2);

            // compute new stratification
            const size_t new_nx{nx * 2};
            const size_t new_ny{ny / 2};
            assert(new_nx * new_ny == nx * ny);

            // process children
            const size_t new_offset{nx * ny};
            initialize_x_tree(x_occupancy, new_offset, new_x, new_y, new_nx, new_ny);
            initialize_x_tree(x_occupancy, new_offset, new_x + 1, new_y, new_nx, new_ny);

            // compute offset for children in global array, if children are occupied then this node has to be marked as
            // occupied
            const size_t l_child{new_offset + new_nx * y + new_x};
            const size_t r_child{new_offset + new_nx * y + new_x + 1};

            assert(l_child < x_occupancy.size());
            assert(r_child < x_occupancy.size());

            x_occupancy[offset] = x_occupancy[l_child] && x_occupancy[r_child];
        }
    }

    void get_valid_x_offsets(Array<bool>& x_occupancy, size_t offset, //
                            const size_t x, const size_t y,          //
                            const size_t nx, const size_t ny,
                             std::vector<NodeOffset>& x_offsets) const
    {
        assert(offset < x_occupancy.size());
        assert(x < nx);
        assert(y < ny);

        // check if node has been already filled
        const size_t node = offset + nx * y + x;
        assert(node < x_occupancy.size());

        if (!x_occupancy[node])
        {
            bool is_leaf = ny == 1;
            if(is_leaf)
            {
                x_offsets.push_back(NodeOffset{node, x});
            }
            else
            {
                // compute new intervals
                const size_t new_x{2 * x};
                const size_t new_y = std::floor(y / 2);

                // compute new stratification
                const size_t new_nx{nx * 2};
                const size_t new_ny{ny / 2};
                assert(new_nx * new_ny == nx * ny);

                // process children
                const size_t new_offset{offset + (nx * ny)};
                get_valid_x_offsets(x_occupancy, new_offset, new_x, new_y, new_nx, new_ny, x_offsets);
                get_valid_x_offsets(x_occupancy, new_offset, new_x + 1, new_y, new_nx, new_ny, x_offsets);
            }
        }
    }

    void initialize_y_tree(Array<bool>& y_occupancy, size_t offset, //
                           const size_t x, const size_t y,          //
                           const size_t nx, const size_t ny) const
    {
        assert(offset < y_occupancy.size());
        assert(x < nx);
        assert(y < ny);

        // check if node has been already filled
        const size_t node = offset + nx * y + x;
        assert(node < y_occupancy.size());

        if (y_occupancy[node])
        {
            return;
        }

        // if not continue drilling the hierarchy, to mark nodes occupied
        bool is_leaf = nx == 1;
        if (!is_leaf)
        {
            // compute new intervals
            const size_t new_x = std::floor(x / 2);;
            const size_t new_y{2 * y};

            // compute new stratification
            const size_t new_nx{nx / 2};
            const size_t new_ny{ny * 2};
            assert(new_nx * new_ny == nx * ny);

            // process children
            const size_t new_offset{nx * ny};
            initialize_y_tree(y_occupancy, new_offset, new_x, new_y, new_nx, new_ny);
            initialize_y_tree(y_occupancy, new_offset, new_x, new_y + 1, new_nx, new_ny);

            // compute offset for children in global array, if children are occupied then this node has to be marked as
            // occupied
            const size_t l_child{new_offset + new_nx * y + new_x};
            const size_t r_child{new_offset + new_nx * (y + 1) + new_x};

            assert(l_child < y_occupancy.size());
            assert(r_child < y_occupancy.size());

            y_occupancy[offset] = y_occupancy[l_child] && y_occupancy[r_child];
        }
    }

    void get_valid_y_offsets(Array<bool>& y_occupancy, size_t offset, //
                             const size_t x, const size_t y,          //
                             const size_t nx, const size_t ny,
                             std::vector<NodeOffset>& y_offsets) const
    {
        assert(offset < y_occupancy.size());
        assert(x < nx);
        assert(y < ny);

        // check if node has been already filled
        const size_t node = offset + nx * y + x;
        assert(node < y_occupancy.size());

        if (!y_occupancy[node])
        {
            bool is_leaf = nx == 1;
            if(is_leaf)
            {
                y_offsets.push_back(NodeOffset{node, y});
            }
            else
            {
                // compute new intervals
                const size_t new_x = std::floor(x / 2);;
                const size_t new_y{2 * y};

                // compute new stratification
                const size_t new_nx{nx / 2};
                const size_t new_ny{ny * 2};
                assert(new_nx * new_ny == nx * ny);

                // process children
                const size_t new_offset{offset + (nx * ny)};
                get_valid_y_offsets(y_occupancy, new_offset, new_x, new_y, new_nx, new_ny, y_offsets);
                get_valid_y_offsets(y_occupancy, new_offset, new_x, new_y+1, new_nx, new_ny, y_offsets);
            }
        }
    }

    RandomNumberGenerator m_rng;
};

#endif // SAMPLING_ELEMENTARY_INTERVAL_H
