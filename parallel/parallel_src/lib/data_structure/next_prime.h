/******************************************************************************
 * next_prime.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Graph Partitioning 
 ******************************************************************************
 * Copyright (C) 2017 Christian Schulz 
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef NEXT_PRIME
#define NEXT_PRIME

static const std::size_t small_primes[] =
{
    2,
    3,
    5,
    7,
    11,
    13,
    17,
    19,
    23,
    29
};

static const std::size_t indices[] =
{
    1,
    7,
    11,
    13,
    17,
    19,
    23,
    29
};

bool
is_prime(std::size_t x)
{
    const size_t N = sizeof(small_primes) / sizeof(small_primes[0]);
    for (std::size_t i = 3; i < N; ++i)
    {
        const std::size_t p = small_primes[i];
        const std::size_t q = x / p;
        if (q < p)
            return true;
        if (x == q * p)
            return false;
    }
    for (std::size_t i = 31; true;)
    {
        std::size_t q = x / i;
        if (q < i)
            return true;
        if (x == q * i)
            return false;
        i += 6;

        q = x / i;
        if (q < i)
            return true;
        if (x == q * i)
            return false;
        i += 4;

        q = x / i;
        if (q < i)
            return true;
        if (x == q * i)
            return false;
        i += 2;

        q = x / i;
        if (q < i)
            return true;
        if (x == q * i)
            return false;
        i += 4;

        q = x / i;
        if (q < i)
            return true;
        if (x == q * i)
            return false;
        i += 2;

        q = x / i;
        if (q < i)
            return true;
        if (x == q * i)
            return false;
        i += 4;

        q = x / i;
        if (q < i)
            return true;
        if (x == q * i)
            return false;
        i += 6;

        q = x / i;
        if (q < i)
            return true;
        if (x == q * i)
            return false;
        i += 2;
    }
    return true;
}

std::size_t
next_prime(std::size_t n)
{
    const size_t L = 30;
    const size_t N = sizeof(small_primes) / sizeof(small_primes[0]);
    // If n is small enough, search in small_primes
    if (n <= small_primes[N-1])
        return *std::lower_bound(small_primes, small_primes + N, n);
    // Else n > largest small_primes
    // Start searching list of potential primes: L * k0 + indices[in]
    const size_t M = sizeof(indices) / sizeof(indices[0]);
    // Select first potential prime >= n
    //   Known a-priori n >= L
    size_t k0 = n / L;
    size_t in = std::lower_bound(indices, indices + M, n - k0 * L) - indices;
    n = L * k0 + indices[in];
    while (!is_prime(n))
    {
        if (++in == M)
        {
            ++k0;
            in = 0;
        }
        n = L * k0 + indices[in];
    }
    return n;
}


#endif /* end of include guard: */
