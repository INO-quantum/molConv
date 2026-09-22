// random.hpp
// implements random number generator (RNG) either from random C++11 library or custom implementation here.
// created 4/11/2023 by Andi
// last change 19/9/2026 by Andi

#ifndef RANDOM_NUMBER_GENERATOR
#define RANDOM_NUMBER_GENERATOR

#include <stdint.h>
#include <random>
#include <iomanip>

#include "threads.h"

// convert 64/32bit unsigned integer to double/float
// - the individual manipulation of mantissa and exponent is the standard way but it generates only half of the possible random numbers!
// - the multiplication used to be slower but on modern computers should be ok. I use this since v1.6.
// - the simple division by 2^32 (MODULUS in Lehmer32) used in old codes does not distribute the floats equally and should not be used!

#define TO_DOUBLE    MULT // MANT_EXP (until v1.5), MULT (since v1.6)

#if ( TO_DOUBLE == MANT_EXP )

// IEEE 754 double precision 64bit: {1bit sign, 11bit exponent, 52bit mantissa}
// see https://prng.di.unimi.it/
// note: this generates only half of the possible random numbers!
static inline double to_double(uint64_t x) {
    union { uint64_t i; double d; } u;
    u.i = (UINT64_C(0x3FF) << 52) | (x >> 12);
    return u.d - 1.0;
}

// IEEE 754 single precision 32bit: {1bit sign, 8bit exponent, 23bit mantissa}
// https://blog.bithole.dev/blogposts/random-float/
// note: this generates only half of the possible random numbers!
static inline float to_double(uint32_t x) {
    union { uint32_t i; float f; } u;
    u.i = (0x7F << 23) | (x >> 9);
    return u.f - 1.0f;
}

// inverse functions double/float to uint64_t/uint32_t
// gives only MSB 52/23 bits!
static inline uint64_t from_double(double x) {
    union { uint64_t i; double d; } u;
    u.d = x + 1.0;
    return (u.i & (((uint64_t)1<<52)-1)) << 12;
}

static inline uint32_t from_double(float x) {
    union { uint32_t i; float d; } u;
    u.d = x + 1.0f;
    return (u.i & (((uint32_t)1<<23)-1)) << 9;
}

#elif ( TO_DOUBLE == MULT )

// IEEE 754 double precision 64bit: {1bit sign, 11bit exponent, 52bit mantissa}
// this gives full resolution of random numbers but might be slightly slower than MANT_EXP
// https://prng.di.unimi.it/
#define TO_DOUBLE_MULT_CONST    0x1.0p-53
static inline double to_double(uint64_t x) {
    return (x >> 11) * TO_DOUBLE_MULT_CONST;
}

// IEEE 754 single precision 32bit: {1bit sign, 8bit exponent, 23bit mantissa}
// https://stackoverflow.com/questions/79581403/how-do-i-convert-a-random-uint32-t-into-a-random-float-in-the-interval-0-1-us
#define TO_FLOAT_MULT_CONST     0x1.0p-24f
static inline float to_double(uint32_t x) {
    return (x >> 8) * 0x1.0p-24f;
}

// inverse functions double/float to uint64_t/uint32_t
// gives only MSB 53/24 bits! but 1 bit more than with MANT_EXP
static inline uint64_t from_double(double x) {
    return (((uint64_t)(x/TO_DOUBLE_MULT_CONST))<<11);
}

static inline uint32_t from_double(float x) {
    return (((uint64_t)(x/TO_FLOAT_MULT_CONST))<<8);
}

#endif

// TODO: somehow the definition of CLASS_UINT128 in threads.h is not regognized here???
//       probably since it is a circular import: threads.h imports random.hpp and vice versa.
#if defined(CLASS_UINT128) || defined(_WIN32) || defined (_WIN64)

// MSVC does not support 128bit integers, so need to implement it here
// best is using the gcc implementation here:
// https://github.com/gcc-mirror/gcc/blob/master/libstdc%2B%2B-v3/src/c%2B%2B17/uint128_t.h
// TODO: newer versions of Visual Studio can use clang which implements uint128_t

#ifndef CLASS_UINT128
#define CLASS_UINT128
#endif

#include "uint128_t.h"

#define uint128       my_uint128_t

#else

// GCC knows __uint128_t
#define uint128     __uint128_t

#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
// virtual basic template class which all custom random generators need to implement
// implements operator(), seed(), min() and max() as given for std::UniformRandomBitGenerator
// which can be passed to any instance of std::RandomNumberDistribution
// T must be an integer type. this is the result_type of the operator().
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename result_type, class=std::enable_if_t<std::is_integral<result_type>::value>>
class random_generator_basic {
public:
    // generate a single random number
    virtual result_type operator()() = 0;
    
    // generate a bunch of random numbers
    void operator()(result_type *vrnd, long num) { 
        for (long i = 0; i < num; ++i) *vrnd++ = this(); 
    }

    // reset to seed values
    virtual void seed(std::seed_seq &seq) = 0;

#if defined(_WIN32) || defined (_WIN64) // windows
    static result_type min(void) { return 0; };
    static result_type max(void) { return 0; };
#else
    // note: each new compiler version breaks this code? cannot use static since expects constexpr.
    virtual result_type min(void) const = 0;
    virtual result_type max(void) const = 0;
#endif
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// this is a interface class to access standard and custom random number generators with a single pointer
// T can be float or integer type.
// TODO: test if using this vs. direct call to RNG might be slower? 
//       even if this would be slower, calling operator(ptr, num) to get several random numbers is better,
//       and additionally, atom creation is not the bottleneck of Monte Carlo but molecule search.
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
class random_generator {
public:
    random_generator() {};
    virtual ~random_generator() {};
    
    virtual T operator()() = 0;
    virtual void operator()(T *vrnd, long num) = 0;
    virtual void generate(std::vector<T> &vec) = 0;
    virtual void seed(std::vector<uint32_t> &vec) = 0;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// direct implementation of random number generator to generate uniform distributed data.
// notes: 
// - these are specific implementations for random number generators returning result_type (real or integer) data.
// - this should be faster than RNG_generic but gives only uniform distributed data.
// - give int_type appropriate for result_type, i.e. uint64_t for double or uint32_t for float.
// - the RNG must be seeded with std::seed_seq& which is the only way to give variable-size seed to standard RNGs.
//   in case needed, add template parameter for seed_type of constructor.
////////////////////////////////////////////////////////////////////////////////////////////////////

template <
    typename result_type = double, 
    typename int_type    = uint64_t, 
    class Generator      = std::mt19937_64, 
    class = std::enable_if_t<std::is_floating_point<result_type>::value>
    >
class RNG_uni_direct_real : public random_generator<result_type> {
private:
    Generator engine;
public:
    explicit RNG_uni_direct_real(std::seed_seq & seed) : random_generator<result_type>(), engine(seed) {};
    ~RNG_uni_direct_real() {};
    
    result_type operator()() { 
        // note: PractRand test fails when casting float from double! this obviously correlates some bits.
        //       but the test passes with to_double called with uint32_t casted (truncated) from uint64_t!
        //       therefore, we explicitly cast to int_type before calling to_double.
        //       this ensures that proper to_double function is called for the result_type!
        return to_double(static_cast<int_type>(engine())); 
    }

    void operator()(result_type *vrnd, long num) { 
        for (long i = 0; i < num; ++i) *vrnd++ = to_double(static_cast<int_type>(engine()));
    }

    // fill vector with random numbers
    void generate(std::vector<result_type> &vec) {
        for (auto &value : vec) {
            value = to_double(static_cast<int_type>(engine()));
        }
    }

    // reset to seed values
    void seed(std::vector<uint32_t> &vec) {
        std::seed_seq seq(vec.begin(), vec.end());
        engine.seed(seq);
    }
};

template <
    typename result_type = uint64_t,
    class    Generator   = std::mt19937_64, 
    class = std::enable_if_t<std::is_integral<result_type>::value>
    >
class RNG_uni_direct_int : public random_generator<result_type> {
private:
    Generator engine;
public:
    explicit RNG_uni_direct_int(std::seed_seq & s) : random_generator<result_type>(), engine(s) {};
    ~RNG_uni_direct_int() {};

    result_type operator()() {
        // note: a larger integer type (uint64_t) might be cast (truncated) here to a smaller one (uint32_t)
        return static_cast<result_type>(engine()); 
    }

    void operator()(result_type *vrnd, long num) { 
        for (long i = 0; i < num; ++i) *vrnd++ = static_cast<result_type>(engine()); 
    }

    // fill vector with random numbers
    void generate(std::vector<result_type> &vec) {
        for (auto &value : vec) {
            value = static_cast<result_type>(engine());
        }
    }

    // reset to seed values
    void seed(std::vector<uint32_t> &vec) {
        std::seed_seq seq(vec.begin(), vec.end());
        engine.seed(seq);
    }
    
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// wrapper class for C++11 random library and other classes to generate
// arbitrary distributed data from any kind of random number generator.
// this is not the most efficient but allows to easily generate non-uniform distributed data.
// note: result_type must be the appropriate type for the given distribution, i.e.
//       floating point for std::uniform_real_distribution and std_normal,
//       integer type for std::uniform_integer_distribution (for integer use RNG_uni_direct_int)
////////////////////////////////////////////////////////////////////////////////////////////////////

template<
    typename result_type = double,
    class Generator      = std::mt19937_64,
    class D              = std::uniform_real_distribution<result_type>
    > 
class RNG_generic : public random_generator<result_type> {
private:
    Generator engine;
    D dist;
public:
    explicit RNG_generic(std::seed_seq & seed) : random_generator<result_type>(), engine(seed), dist() {};
    
    // generate next random value from distribution
    result_type operator()() { return dist(engine); }
    void operator()(result_type *vrnd, long num) { 
        for (long i = 0; i < num; ++i) *vrnd++ = dist(engine); 
    }
    // fill vector with random numbers
    void generate(std::vector<result_type> &vec) {
        for (auto &value : vec) {
            value = dist(engine);
        }
    }

    // reset to seed values
    void seed(std::vector<uint32_t> &vec) {
        std::seed_seq seq(vec.begin(), vec.end());
        engine.seed(seq);
    }

};

////////////////////////////////////////////////////////////////////////////////////////////////////
// Lehmer32 generator
// DO NOT USE!!! I keep it here only for reference.
// this is the original generator used by Florian. 
// it is fast but its period of 2^32 is not sufficient for a Monte Carlo simulation! 
// especially in combination with inefficient rejection sampling, using 2 floats per random number, and signed integer,
// already for 2x15k atoms the produced data has correlations leading to increased 'noise' in the distribution on the y-axis!
// switching between streams of Lehmer32 but seeded differently does not improve anything (I guess the correlations persist).
// this implementation is similar (or even the same?) as for LCG 'minstd_rand' in the C++ random library.
// see Schrage's Method: https://en.wikipedia.org/wiki/Lehmer_random_number_generator
// note: the original implementation generates the float by directly dividing by MODULUS which is not recommended.
//       the min/max functions are needed for RNG_generic to generate non-uniform random numbers.
////////////////////////////////////////////////////////////////////////////////////////////////////
 
#define MODULUS    0x7fffffff
#define MULTIPLIER 48271

class Lehmer32 : public random_generator_basic<int32_t> {
private:
    const int32_t Q = MODULUS / MULTIPLIER;
    const int32_t R = MODULUS % MULTIPLIER;
    int32_t state;
public:
    Lehmer32(int32_t seed)       : random_generator_basic<int32_t>() { 
        state = seed;
        //if (state < 0) state += MODULUS;
        if (state <= 0) state += MODULUS; // update v1.6 to ensure seed != 0
    };
    Lehmer32(std::seed_seq &seq) : random_generator_basic<int32_t>() { 
        seed(seq);
    };
    ~Lehmer32() {};

    int32_t operator()() {
        int32_t t = MULTIPLIER * (state % Q) - R * (state / Q);
        state = (t > 0) ? t : t + MODULUS;
        // note: next line would convert the state to float/double which we not need here.
        //return ((T) state / MODULUS);
        return state;
    };
    
    // reset to seed values
    void seed(std::seed_seq & seq) {
        std::vector<int32_t> seed_values(1);
        seq.generate(seed_values.begin(), seed_values.end());
        state = seed_values[0];
        //if (state < 0) state += MODULUS;
        if (state <= 0) state += MODULUS; // update v1.6 to ensure seed != 0
    };

#if defined(_WIN32) || defined (_WIN64) // windows
    static int32_t min(void) { return 0; };
    static int32_t max(void) { return MODULUS; };
#else
    constexpr int32_t min(void) { return 0; };
    constexpr int32_t max(void) { return MODULUS; };
#endif
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// Lehmer64 generator
// this is very fast. I have not found what is the period?
// this is the default generator. 
// however, it fails in PractRand 0.94 statistical test
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
// D. H. Lehmer, Mathematical methods in large-scale computing units.
// Proceedings of a Second Symposium on Large Scale Digital Calculating
// Machinery;
// Annals of the Computation Laboratory, Harvard Univ. 26 (1951), pp. 141-146.
////////////////////////////////////////////////////////////////////////////////////////////////////

// this number is from the blog below but have not found a published paper or its periodiciy
// this is 64bit multiplication with a 128bit number. this should be slightly faster than Lehmer128.
// https://lemire.me/blog/2019/03/19/the-fastest-conventional-random-number-generator-that-can-pass-big-crush/
// https://github.com/lemire/testingRNG/blob/master/source/lehmer64.h
#define LEHMER64_MUL UINT64_C(0xda942042e4dd58b5)

// TODO: seed cannot be 0! and maybe has to be odd always?
//       therefore, |1 is used which ensures both.

class Lehmer64 : public random_generator_basic<uint64_t> {
private:
    uint128 state;
public:
    Lehmer64(uint128 seed)       : random_generator_basic<uint64_t>() { 
        state = seed | 1;
    };
    Lehmer64(std::seed_seq& seq) : random_generator_basic<uint64_t>() { 
        seed(seq);
    };
    ~Lehmer64() {};

    uint64_t operator()() { 
        state *= LEHMER64_MUL;
        return (uint64_t)(state >> 64);
    };

    // reset to seed value
    void seed(std::seed_seq &seq) {
        std::vector<uint32_t> seed_values(4);
        seq.generate(seed_values.begin(), seed_values.end());
        state = (((uint128)  seed_values[3]     ) << 96) | 
                (((uint128)  seed_values[2]     ) << 64) | 
                (((uint128)  seed_values[1]     ) << 32) | 
                 ((uint128) (seed_values[0] | 1))        ; 
    };

#if defined(_WIN32) || defined (_WIN64) // windows
    static uint64_t min(void) { return 0; };
    static uint64_t max(void) { return UINT64_C(0xffffffffffffffff); };
#else
    constexpr virtual uint64_t min(void) const { return 0; };
    constexpr virtual uint64_t max(void) const { return UINT64_C(0xffffffffffffffff); };
#endif
};

// this is similar to Lehmer64 but it uses a 128x128-bit multiplication
// this might be slower than the 64x128-bit multiplication but it still very fast.
// according to Wikipedia this has a period of 2^126.
// it fails PraceRand test but which is expected for a modulo-2 PRNG of this kind.
// https://en.wikipedia.org/wiki/Lehmer_random_number_generator
// with python calculated hex -> decimal conversion:
// (0x12e15e35b500f16e << 64) | 0x2e714eb2b37916a5
// decimal 25096281518912105342191851917838718629
// this is from the paper below for 2^128 with c=0

////////////////////////////////////////////////////////////////////////////////////////////////////
// P L'Ecuyer,  Tables of linear congruential generators of different sizes and good lattice structure. 
// Mathematics of Computation of the American Mathematical Society 68.225 (1999): 249-260.
////////////////////////////////////////////////////////////////////////////////////////////////////

// TODO: seed cannot be 0! here I could not find anything that it must be odd?
//       to be sure, |1 is used which ensures both.

#ifdef CLASS_UINT128
const uint128 Lehmer128_mult(UINT64_C(0x2e714eb2b37916a5), UINT64_C(0x12e15e35b500f16e));
#else
const uint128 Lehmer128_mult = (((uint128)UINT64_C(0x12e15e35b500f16e)) << 64) | UINT64_C(0x2e714eb2b37916a5);
#endif

class Lehmer128 : public random_generator_basic<uint64_t> {
private:
    uint128 state;
public:
    Lehmer128(uint128 seed)                : random_generator_basic<uint64_t>() { 
        state = seed | 1;
    };
    Lehmer128(std::seed_seq& seq)          : random_generator_basic<uint64_t>() { 
        seed(seq);
    };
    ~Lehmer128() {};

    uint64_t operator()() { 
        state *= Lehmer128_mult;
        return (state >> 64);
    };

    // reset to seed value
    void seed(std::seed_seq &seq) {
        std::vector<uint32_t> seed_values(4);
        seq.generate(seed_values.begin(), seed_values.end());
#ifdef _DEBUG
    std::cout << std::hex << "create Lehmer128 seed 0x" << std::setw(8) << seed_values[0] 
                                              << ", 0x" << std::setw(8) << seed_values[1] 
                                              << ", 0x" << std::setw(8) << seed_values[2] 
                                              << ", 0x" << std::setw(8) << seed_values[3] << std::dec << std::endl; 
#endif          
        state = (((uint128)  seed_values[3]     ) << 96) | 
                (((uint128)  seed_values[2]     ) << 64) | 
                (((uint128)  seed_values[1]     ) << 32) | 
                 ((uint128) (seed_values[0] | 1))        ; 
        if (state == ((uint128)0)) state += 1; // ensure seed != 0 
    };

#if defined(_WIN32) || defined (_WIN64) // windows
    static uint64_t min(void) { return 0; };
    static uint64_t max(void) { return UINT64_C(0xffffffffffffffff); };
#else
    constexpr virtual uint64_t min(void) const { return 0; };
    constexpr virtual uint64_t max(void) const { return UINT64_C(0xffffffffffffffff); };
#endif

};

// code from here:
// https://lemire.me/blog/2019/03/19/the-fastest-conventional-random-number-generator-that-can-pass-big-crush/
// https://github.com/lemire/testingRNG/blob/master/source/lehmer64.h
// adapted to this project by D. Lemire, from https://github.com/wangyi-fudan/wyhash/blob/master/wyhash.h
// This uses mum hashing.
// The state can be seeded with any value!

class Wyhash64 : public random_generator_basic<uint64_t> {
private:
    uint_fast64_t state;
public:
    Wyhash64(uint64_t seed)      : random_generator_basic<uint64_t>() { 
        state = seed;
    };
    Wyhash64(std::seed_seq& seq) : random_generator_basic<uint64_t>() { 
        seed(seq);
    };
    ~Wyhash64() {};

    uint64_t operator()() { 
        state += UINT64_C(0x60bee2bee120fc15);
        uint128 tmp = (uint128)state; 
        tmp *= UINT64_C(0xa3b195354a39b70d);
        uint64_t m1 = (tmp >> 64) ^ tmp;
        tmp = (uint128)m1;
        tmp *= UINT64_C(0x1b03738712fad5c9);
        uint64_t m2 = (tmp >> 64) ^ tmp;
        return m2;
    };

    // reset to seed values
    void seed(std::seed_seq &seq) {
        std::vector<uint_fast32_t> seed_values(2);
        seq.generate(seed_values.begin(), seed_values.end());
        state = (((uint_fast64_t) seed_values[1]) << 32) | 
                 ((uint_fast64_t) seed_values[0])        ; 
    };

#if defined(_WIN32) || defined (_WIN64) // windows
    static uint64_t min(void) { return 0; };
    static uint64_t max(void) { return UINT64_C(0xffffffffffffffff); };
#else
    constexpr virtual uint64_t min(void) const { return 0; };
    constexpr virtual uint64_t max(void) const { return UINT64_C(0xffffffffffffffff); };
#endif

};

// implementation of WELL1024a by FRANCOIS PANNETON and PIERRE L'ECUYER
// paper "Generators Based on Linear Recurrences Modulo 2"
// http://www.iro.umontreal.ca/~lecuyer/myftp/papers/lfsr04.pdf
// init WELL with WELL1024a_R random uint32_t integers!
// Andi: copied from one of my old projectes

#define WELL1024a_R 32
#define M1 3
#define M2 24
#define M3 10
#define MAT3POS(t,v) (v^(v>>t))
#define MAT3NEG(t,v) (v^(v<<(t)))
#define Identity(v) (v)
#define V0 STATE[ state_n ]
#define VM1 STATE[ (state_n+M1) & 0x0000001fUL ]
#define VM2 STATE[ (state_n+M2) & 0x0000001fUL ]
#define VM3 STATE[ (state_n+M3) & 0x0000001fUL ]
#define VRm1 STATE[ (state_n+31) & 0x0000001fUL ]
#define newV0 STATE[ (state_n+31) & 0x0000001fUL ]
#define newV1 STATE[ state_n ]

class WELL1024 : public random_generator_basic<uint32_t> {
private:
    uint32_t z0, z1, z2, state_n;
    uint32_t STATE[WELL1024a_R];
public:
    WELL1024(uint32_t seed[WELL1024a_R]) : random_generator_basic<uint32_t>() { 
        state_n = 0;
        for (int j = 0; j < WELL1024a_R; j++) STATE[j] = seed[j];
    };
    WELL1024(std::seed_seq& seq)         : random_generator_basic<uint32_t>() { 
        seed(seq);
    };
    ~WELL1024() {};

    uint32_t operator()() { 
        z0 = VRm1;
        z1 = Identity(V0) ^ MAT3POS (8, VM1);
        z2 = MAT3NEG (19, VM2) ^ MAT3NEG(14,VM3);
        newV1 = z1 ^ z2;
        newV0 = MAT3NEG (11,z0) ^ MAT3NEG(7,z1) ^ MAT3NEG(13,z2) ;
        state_n = (state_n + 31) & 0x0000001fUL; // note: dependency on WELL1024a_R is hard-coded here!
        // note: next line would convert the state to double which we not need here
        //return ((double) STATE[state_n] * 2.32830643653869628906e-10);
        return STATE[state_n];
    };

    // reset to seed values
    void seed(std::seed_seq &seq) {
        std::vector<uint_fast32_t> seed_values(WELL1024a_R);
        seq.generate(seed_values.begin(), seed_values.end());
        state_n = 0;
        for (int j = 0; j < WELL1024a_R; j++) STATE[j] = seed_values[j];
    };

#if defined(_WIN32) || defined (_WIN64) // windows
    static uint32_t min(void) { return 0; };
    static uint32_t max(void) { return (uint32_t)0xffffffff; };
#else
    constexpr virtual uint32_t min(void) const { return 0; };
    constexpr virtual uint32_t max(void) const { return (uint32_t)0xffffffff; };
#endif

};

#endif // RANDOM_NUMBER_GENERATOR

////////////////////////////////////////////////////////////////////////////////////////////////////

