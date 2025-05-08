#include <functional>
#include <random>

#include "rand.h"

/* Make an RNG with fixed bounds */
bounded_rng_type makeRNGRange(int seed, const size_t& min_val, const size_t& max_val)
{
	static std::mt19937_64 gen(seed);
	static std::uniform_int_distribution<size_t> dist(min_val, max_val);
	auto rng = [&]() { return dist(gen); };
	return rng;
}


/* Make an RNG with variable bounds */
rng_type makeRNG(int seed)
{
	static std::mt19937_64 gen(seed);
    static std::uniform_int_distribution<size_t> dist(0);
	auto rng = [&](const size_t& min_val, const size_t& max_val)
    {
        std::uniform_int_distribution<size_t>::param_type params{min_val, max_val};
        dist.param(params);
        return dist(gen); 
    };
	return rng;
}


bounded_rng_type makeRNG(int seed, const size_t& min_val, const size_t& max_val)
{
	static std::mt19937_64 gen(seed);
	static std::uniform_int_distribution<size_t> dist(min_val, max_val);
	auto rng = [&]() { return dist(gen); };
	return rng;
}