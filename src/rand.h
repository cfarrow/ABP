#include <functional>
#include <random>

#ifndef ABP_RAND_H
#define ABP_RAND_H


using bounded_rng_type = const std::function<size_t(void)>;
using rng_type = const std::function<size_t(size_t, size_t)>;


/* Make an RNG with fixed bounds */
bounded_rng_type makeRNG(int seed, const size_t& min_val, const size_t& max_val);

/* Make an RNG with passable bounds */
rng_type makeRNG(int seed);


#endif