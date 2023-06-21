#ifndef BMUTIL__H__INCLUDED__
#define BMUTIL__H__INCLUDED__
/*
Copyright(c) 2002-2017 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

For more information please visit:  http://bitmagic.io
*/

/*! \file bmutil.h
    \brief Bit manipulation primitives (internal)
*/

#include "bmdef.h"
#include "bmconst.h"

#if defined(_M_AMD64) || defined(_M_X64) 
#include <intrin.h>
#elif defined(BMSSE2OPT) || defined(BMSSE42OPT)
#include <emmintrin.h>
#elif defined(BMAVX2OPT)
#include <emmintrin.h>
#include <avx2intrin.h>
#endif

#ifdef __GNUG__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#endif

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4146)
#endif


namespace bm
{

        /**
        bit-block array wrapped into union for correct interpretation of
        32-bit vs 64-bit access vs SIMD
        @internal
        */
        struct bit_block_t
        {
            union bunion_t
            {
                bm::word_t BM_VECT_ALIGN w32[bm::set_block_size] BM_VECT_ALIGN_ATTR;
                bm::id64_t BM_VECT_ALIGN w64[bm::set_block_size / 2] BM_VECT_ALIGN_ATTR;
                
#if defined(BMAVX512OPT)
                __m512i  BM_VECT_ALIGN w512[bm::set_block_size / 16] BM_VECT_ALIGN_ATTR;
#endif
#if defined(BMAVX2OPT)
                __m256i  BM_VECT_ALIGN w256[bm::set_block_size / 8] BM_VECT_ALIGN_ATTR;
#endif
#if defined(BMSSE2OPT) || defined(BMSSE42OPT)
                __m128i  BM_VECT_ALIGN w128[bm::set_block_size / 4] BM_VECT_ALIGN_ATTR;
#endif
            } b_;

            operator bm::word_t*() { return &(b_.w32[0]); }
            operator const bm::word_t*() const { return &(b_.w32[0]); }
            explicit operator bm::id64_t*() { return &b_.w64[0]; }
            explicit operator const bm::id64_t*() const { return &b_.w64[0]; }
#ifdef BMAVX512OPT
            explicit operator __m512i*() { return &b_.w512[0]; }
            explicit operator const __m512i*() const { return &b_.w512[0]; }
#endif
#ifdef BMAVX2OPT
            explicit operator __m256i*() { return &b_.w256[0]; }
            explicit operator const __m256i*() const { return &b_.w256[0]; }
#endif
#if defined(BMSSE2OPT) || defined(BMSSE42OPT)
            explicit operator __m128i*() { return &b_.w128[0]; }
            explicit operator const __m128i*() const { return &b_.w128[0]; }
#endif

            const bm::word_t* begin() const { return (b_.w32 + 0); }
            bm::word_t* begin() { return (b_.w32 + 0); }
            const bm::word_t* end() const { return (b_.w32 + bm::set_block_size); }
            bm::word_t* end() { return (b_.w32 + bm::set_block_size); }
        };
    
/**
    Get minimum of 2 values
*/
template<typename T>
T min_value(T v1, T v2) BMNOEXCEPT
{
    return v1 < v2 ? v1 : v2;
}

/**
    \brief ad-hoc conditional expressions
    \internal
*/
template <bool b> struct conditional
{
    static bool test() { return true; }
};
template <> struct conditional<false>
{
    static bool test() { return false; }
};


/**
    Fast loop-less function to find LOG2
*/
template<typename T>
T ilog2(T x) BMNOEXCEPT
{
    unsigned int l = 0;
    
    if (x >= 1<<16) { x = (T)(x >> 16); l |= 16; }
    if (x >= 1<<8)  { x = (T)(x >> 8);  l |= 8; }
    if (x >= 1<<4)  { x = (T)(x >> 4);  l |= 4; }
    if (x >= 1<<2)  { x = (T)(x >> 2);  l |= 2; }
    if (x >= 1<<1)  l |=1;
    return (T)l;
}

template<>
inline bm::gap_word_t ilog2(gap_word_t x) BMNOEXCEPT
{
    unsigned int l = 0;
    if (x >= 1<<8)  { x = (bm::gap_word_t)(x >> 8); l |= 8; }
    if (x >= 1<<4)  { x = (bm::gap_word_t)(x >> 4); l |= 4; }
    if (x >= 1<<2)  { x = (bm::gap_word_t)(x >> 2); l |= 2; }
    if (x >= 1<<1)  l |=1;
    return (bm::gap_word_t)l;
}

/**
     Mini auto-pointer for internal memory management
     @internal
*/
template<class T>
class ptr_guard
{
public:
    ptr_guard(T* p) BMNOEXCEPT : ptr_(p) {}
    ~ptr_guard() { delete ptr_; }
private:
    ptr_guard(const ptr_guard<T>& p);
    ptr_guard& operator=(const ptr_guard<T>& p);
private:
    T* ptr_;
};

/**
    Portable LZCNT with (uses minimal LUT)
    @ingroup bitfunc
    @internal
*/
inline unsigned count_leading_zeros(unsigned x) BMNOEXCEPT
{
    unsigned n =
        (x >= (1U << 16)) ?
        ((x >= (1U << 24)) ? ((x >= (1 << 28)) ? 28u : 24u) : ((x >= (1U << 20)) ? 20u : 16u))
        :
        ((x >= (1U << 8)) ? ((x >= (1U << 12)) ? 12u : 8u) : ((x >= (1U << 4)) ? 4u : 0u));
    return unsigned(bm::lzcnt_table<true>::_lut[x >> n]) - n;
}

/**
    Portable TZCNT with (uses 37-LUT)
    @ingroup bitfunc
    @internal
*/
inline
unsigned count_trailing_zeros(unsigned v) BMNOEXCEPT
{
    // (v & -v) isolates the last set bit
    return unsigned(bm::tzcnt_table<true>::_lut[(-v & v) % 37]);
}

/**
    Lookup table based integer LOG2
*/
template<typename T>
T ilog2_LUT(T x) BMNOEXCEPT
{
    unsigned l = 0;
    if (x & 0xffff0000) 
    {
        l += 16; x >>= 16;
    }
    
    if (x & 0xff00) 
    {
        l += 8; x >>= 8;
    }
    return l + T(first_bit_table<true>::_idx[x]);
}

/**
    Lookup table based short integer LOG2
*/
template<>
inline bm::gap_word_t ilog2_LUT<bm::gap_word_t>(bm::gap_word_t x) BMNOEXCEPT
{
    bm::gap_word_t l = 0;
    if (x & 0xff00) 
    {
        l = bm::gap_word_t( + 8u);
        x = bm::gap_word_t(x >> 8u);
    }
    return bm::gap_word_t(l + bm::gap_word_t(first_bit_table<true>::_idx[x]));
}


// if we are running on x86 CPU we can use inline ASM 

#ifdef BM_x86
#ifdef __GNUG__

BMFORCEINLINE
unsigned bsf_asm32(unsigned int v) BMNOEXCEPT
{
    unsigned r;
    asm volatile(" bsfl %1, %0": "=r"(r): "rm"(v) );
    return r;
}
 
BMFORCEINLINE
unsigned bsr_asm32(unsigned int v) BMNOEXCEPT
{
    unsigned r;
    asm volatile(" bsrl %1, %0": "=r"(r): "rm"(v) );
    return r;
}

#endif  // __GNUG__

#ifdef _MSC_VER

#if defined(_M_AMD64) || defined(_M_X64) // inline assembly not supported

BMFORCEINLINE
unsigned int bsr_asm32(unsigned int value) BMNOEXCEPT
{
    unsigned long r;
    _BitScanReverse(&r, value);
    return r;
}

BMFORCEINLINE
unsigned int bsf_asm32(unsigned int value) BMNOEXCEPT
{
    unsigned long r;
    _BitScanForward(&r, value);
    return r;
}

#else

BMFORCEINLINE
unsigned int bsr_asm32(unsigned int value) BMNOEXCEPT
{   
  __asm  bsr  eax, value
}

BMFORCEINLINE
unsigned int bsf_asm32(unsigned int value) BMNOEXCEPT
{   
  __asm  bsf  eax, value
}

#endif

#endif // _MSC_VER

#endif // BM_x86


// From:
// http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.37.8562
//
template<typename T>
T bit_scan_fwd(T v) BMNOEXCEPT
{
    return
        DeBruijn_bit_position<true>::_multiply[(((v & -v) * 0x077CB531U)) >> 27];
}

inline
unsigned bit_scan_reverse32(unsigned value) BMNOEXCEPT
{
    BM_ASSERT(value);
#if defined(BM_USE_GCC_BUILD)
    return (unsigned) (31 - __builtin_clz(value));
#else
# if defined(BM_x86) && (defined(__GNUG__) || defined(_MSC_VER))
    return bm::bsr_asm32(value);
# else
    return bm::ilog2_LUT<unsigned int>(value);
# endif
#endif
}

inline
unsigned bit_scan_forward32(unsigned value) BMNOEXCEPT
{
    BM_ASSERT(value);
#if defined(BM_USE_GCC_BUILD)
    return (unsigned) __builtin_ctz(value);
#else
# if defined(BM_x86) && (defined(__GNUG__) || defined(_MSC_VER))
    return bm::bsf_asm32(value);
# else
        return bit_scan_fwd(value);
# endif
#endif
}


BMFORCEINLINE
unsigned long long bmi_bslr_u64(unsigned long long w) BMNOEXCEPT
{
#if defined(BMAVX2OPT) || defined (BMAVX512OPT)
    return _blsr_u64(w);
#else
    return w & (w - 1);
#endif
}

BMFORCEINLINE
unsigned long long bmi_blsi_u64(unsigned long long w)
{
#if defined(BMAVX2OPT) || defined (BMAVX512OPT)
    return _blsi_u64(w);
#else
    return w & (-w);
#endif
}

/// 64-bit bit-scan reverse
inline
unsigned count_leading_zeros_u64(bm::id64_t w) BMNOEXCEPT
{
    BM_ASSERT(w);
#if defined(BMAVX2OPT) || defined (BMAVX512OPT)
    return (unsigned)_lzcnt_u64(w);
#else
    #if defined(BM_USE_GCC_BUILD)
        return (unsigned) __builtin_clzll(w);
    #else
        unsigned z;
        unsigned w1 = unsigned(w >> 32);
        if (!w1)
        {
            z = 32;
            w1 = unsigned(w);
            z += 31 - bm::bit_scan_reverse32(w1);
        }
        else
        {
            z = 31 - bm::bit_scan_reverse32(w1);
        }
        return z;
    #endif
#endif
}

/// 64-bit bit-scan fwd
inline
unsigned count_trailing_zeros_u64(bm::id64_t w) BMNOEXCEPT
{
    BM_ASSERT(w);

#if defined(BMAVX2OPT) || defined (BMAVX512OPT)
    return (unsigned)_tzcnt_u64(w);
#else
    #if defined(BM_USE_GCC_BUILD)
        return (unsigned) __builtin_ctzll(w);
    #else
        unsigned z;
        unsigned w1 = unsigned(w);
        if (!w1)
        {
            z = 32;
            w1 = unsigned(w >> 32);
            z += bm::bit_scan_forward32(w1);
        }
        else
        {
            z = bm::bit_scan_forward32(w1);
        }
        return z;
    #endif
#endif
}



/*!
    Returns BSR value
    @ingroup bitfunc
*/
template <class T>
unsigned bit_scan_reverse(T value) BMNOEXCEPT
{
    BM_ASSERT(value);

    if (bm::conditional<sizeof(T)==8>::test())
    {
    #if defined(BM_USE_GCC_BUILD)
        return (unsigned) (63 - __builtin_clzll(value));
    #else
        bm::id64_t v8 = value;
        v8 >>= 32;
        unsigned v = (unsigned)v8;
        if (v)
        {
            v = bm::bit_scan_reverse32(v);
            return v + 32;
        }
    #endif
    }
    return bm::bit_scan_reverse32((unsigned)value);
}

/*! \brief and functor
    \internal
 */
struct and_func
{
    static
    BMFORCEINLINE unsigned op(unsigned v1, unsigned v2) BMNOEXCEPT2
        { return v1 & v2; }
};
/*! \brief xor functor
    \internal
 */
struct xor_func
{
    static
    BMFORCEINLINE unsigned op(unsigned v1, unsigned v2) BMNOEXCEPT2
        { return v1 ^ v2; }
};
/*! \brief or functor
    \internal
 */
struct or_func
{
    static
    BMFORCEINLINE unsigned op(unsigned v1, unsigned v2) BMNOEXCEPT2
        { return v1 | v2; }
};
/*! \brief sub functor
    \internal
 */
struct sub_func
{
    static
    BMFORCEINLINE unsigned op(unsigned v1, unsigned v2) BMNOEXCEPT2
        { return v1 & ~v2; }
};



#ifdef __GNUG__
#pragma GCC diagnostic pop
#endif
#ifdef _MSC_VER
#pragma warning( pop )
#endif


} // bm

#endif
