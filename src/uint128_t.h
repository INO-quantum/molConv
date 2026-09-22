// A relatively minimal unsigned 128-bit integer class type, used by the
// floating-point std::to_chars implementation on targets that lack __int128.

// Copyright (C) 2021-2023 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

// downloaded 3/12/2023 by Andi
// https://github.com/gcc-mirror/gcc/blob/master/libstdc%2B%2B-v3/src/c%2B%2B17/uint128_t.h
// simplified for only the functions I need here
// see also test_uint128 in MolecularConversion.cpp for testing.
// TODO: try to compile with clang on Visual Studio, then this code is not needed any longer.

// last change 19/4/2024 by Andi
// comment about clang added on 20/9/2026

#ifndef _UINT128_IMPLEMENTATION
#define _UINT128_IMPLEMENTATION

#if defined(_WIN32) || defined(_WIN64)
#include <stdint.h>
#include <intrin.h>
#endif

struct my_uint128_t
{
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    uint64_t lo, hi;
#else
    uint64_t hi, lo;
#endif

    // default constructor
    my_uint128_t() = default;

    // constructor
    constexpr
    my_uint128_t(uint64_t lo, uint64_t hi = 0)
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
        : lo(lo), hi(hi)
#else
        : hi(hi), lo(lo)
#endif
    { }
    
    // comparison
    friend bool 
    operator==(const my_uint128_t &x, const my_uint128_t &y) {
        return ((x.lo==y.lo) && (x.hi==y.hi));
    }

    friend bool
    operator==(const my_uint128_t& x, const uint64_t& y) {
        return ((x.lo == y) && (x.hi == 0));
    }

    friend bool
    operator!=(const my_uint128_t &x, const my_uint128_t &y) {
        return (!(x == y));
    }
    
    friend bool
        operator!=(const my_uint128_t& x, const uint64_t& y) {
        return (!((x.lo == y) && (x.hi == 0)));
    }

    // not really needed
    friend bool 
    operator<(const my_uint128_t &l, const my_uint128_t &r) {
        return (l-r);
    }

    // cast to uint64_t
    operator uint64_t() const {
        return lo;
    }
    
    // addition
    friend my_uint128_t
    operator+(my_uint128_t x, const my_uint128_t y)
    {
#if defined (_WIN32) || defined(_WIN64)
        //unsigned char carry = 0;
        x.hi += _addcarry_u64(0, x.lo, y.lo, &x.lo);
        x.hi += y.hi;
#else
        x.hi += __builtin_add_overflow(x.lo, y.lo, &x.lo);
        x.hi += y.hi;
#endif
        return x;
    }
    
    my_uint128_t&
    operator+=(const my_uint128_t y)
    { 
        return *this = *this + y; 
    }

    // pre-increment
    my_uint128_t&
    operator++()
    { 
        return *this += 1; 
    }

    // post-increment
    my_uint128_t
    operator++(int)
    { 
        return *this += 1; 
    }
    
    // multiplication
    static my_uint128_t
    umul64_64_128(const uint64_t x, const uint64_t y)
    {
#if defined(_WIN64)
        uint64_t l, h;
        l = _umul128(x,y,&h);
#else 
        const uint64_t xl = x & 0xffffffff;
        const uint64_t xh = x >> 32;
        const uint64_t yl = y & 0xffffffff;
        const uint64_t yh = y >> 32;
        const uint64_t ll = xl * yl;
        const uint64_t lh = xl * yh;
        const uint64_t hl = xh * yl;
        const uint64_t hh = xh * yh;
        const uint64_t m = (ll >> 32) + lh + (hl & 0xffffffff);
        const uint64_t l = (ll & 0xffffffff ) | (m << 32);
        const uint64_t h = (m >> 32) + (hl >> 32) + hh;
#endif
        return {l, h};
    }

    friend my_uint128_t
    operator*(const my_uint128_t x, const my_uint128_t y)
    {
        my_uint128_t z = umul64_64_128(x.lo, y.lo);
        z.hi += x.lo * y.hi + x.hi * y.lo;
        return z;
    }
    my_uint128_t&
    operator*=(const my_uint128_t y)
    { 
        return *this = *this * y; 
    }

    // shift-right
    friend my_uint128_t
    operator>>(my_uint128_t x, const int y)
    {
        if (y >= 128) 
        {
            x.lo = 0;
            x.hi = 0;
        } 
        else if (y >= 64)
        {
            x.lo = x.hi >> (y - 64);
            x.hi = 0;
        }
        else if (y != 0)
        {
            x.lo >>= y;
            x.lo |= x.hi << (64 - y);
            x.hi >>= y;
        }
        return x;
    }  
    
    my_uint128_t&
    operator>>=(const int y)
    { 
        return *this = *this >> y; 
    }

    // shift-left
    friend my_uint128_t
    operator<<(my_uint128_t x, const int y)
    {
        if (y >= 128) 
        {
            x.lo = 0;
            x.hi = 0;
        } 
        else if (y >= 64)
        {
            x.hi = x.lo << (y - 64);
            x.lo = 0;
        }
        else if (y != 0)
        {
            x.hi <<= y;
            x.hi |= x.lo >> (64 - y);
            x.lo <<= y;
        }
        return x;
    }
 
    my_uint128_t&
    operator<<=(const int y)
    { 
        return *this = *this << y; 
    }

    // AND
    friend my_uint128_t
    operator&(my_uint128_t x, const int y)
    {
        x.lo &= (uint64_t)y;
        x.hi = 0;
        return x;
    }

    friend my_uint128_t
    operator&(my_uint128_t x, const uint64_t y)
    {
        x.lo &= y;
        x.hi = 0;
        return x;
    }

    friend my_uint128_t
    operator&(my_uint128_t x, const my_uint128_t y)
    {
        x.lo &= y.lo;
        x.hi &= y.hi;
        return x;
    }
    
    my_uint128_t&
    operator&=(const my_uint128_t y)
    { 
        return *this = *this & y; 
    }

    // OR
    friend my_uint128_t
    operator|(my_uint128_t x, const int y)
    {
        x.lo |= (uint64_t)y;
        x.hi = 0;
        return x;
    }

    friend my_uint128_t
    operator|(my_uint128_t x, const my_uint128_t y)
    {
        x.lo |= y.lo;
        x.hi |= y.hi;
        return x;
    }
    
    my_uint128_t&
    operator|=(const my_uint128_t y)
    { 
        return *this = *this | y; 
    }

    // XOR
    friend my_uint128_t
    operator^(my_uint128_t x, const my_uint128_t y)
    {
        x.lo ^= y.lo;
        x.hi ^= y.hi;
        return x;
    }
    
    my_uint128_t&
    operator^=(const my_uint128_t y)
    { 
        return *this = *this ^ y; 
    }

};

#endif // _UINT128_IMPLEMENTATION
