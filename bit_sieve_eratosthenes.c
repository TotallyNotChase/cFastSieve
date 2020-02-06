/*
* @file = bit_sieve_eratosthenes.c
* @author = https://github.com/TotallyNotChase
* @licence = public domain
*/

/*
* The following sieve of eratosthenes is a C implementation of an
* improved sieve of eratosthenes algorithm written by Kim Wilsch in C++
* View the documentation of the algo here - [https://github.com/kimwalisch/primesieve/wiki/Segmented-sieve-of-Eratosthenes]
* ///////////////////////////////////////////////////////////////////////
* COPYRIGHT NOTICE:-
* BSD 2-Clause License
*
* Copyright (c) 2010 - 2019, Kim Walisch.
* All rights reserved.
* ///////////////////////////////////////////////////////////////////////
*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdbool.h>
#include<stdint.h>
#include<inttypes.h>
#include<math.h>
#include<time.h>

// Defining custom structs to imitate vectors in C++

typedef struct uint64_vector
{
    // For integer values
    uint64_t* data;
    uint64_t capacity, size;
} uint64vec_t;

typedef struct uint8_vector
{
    // For boolean values
    uint8_t* data;
    uint64_t size, count;
} uint8vec_t;

// The Level 1 Data cache for the user's CPU (must be per core)
uint64_t L1D_CACHE;

const int64_t unset_bit[16] = {
                                ~(1 << 0), ~(1 << 0),
                                ~(1 << 1), ~(1 << 1),
                                ~(1 << 2), ~(1 << 2),
                                ~(1 << 3), ~(1 << 3),
                                ~(1 << 4), ~(1 << 4),
                                ~(1 << 5), ~(1 << 5),
                                ~(1 << 6), ~(1 << 6),
                                ~(1 << 7), ~(1 << 7)
                              };

const int64_t popcnt[256] = {
                              0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
                              1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
                              1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
                              2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                              1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
                              2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                              2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                              3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
                              1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
                              2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                              2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                              3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
                              2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
                              3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
                              3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
                              4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
                            };

uint64vec_t create_uint64_vector(uint64_t size)
{
    // Creates a uint64_t vector and returns it
    uint64vec_t vector;
    vector.data = malloc(size * sizeof(uint64_t));
    if (vector.data == NULL)
    {
        printf("\nAn error occured while allocating memory for uint64_vector\n");
        exit(1);
    }
    vector.capacity = size;
    vector.size = 0;
    return vector;
}

uint8vec_t create_uint8_vector(uint64_t size)
{
    // Creates a bool vector, assigns all slots to True, and returns it
    uint8vec_t vector;
    vector.data = malloc(size * sizeof(uint8_t));
    if (vector.data == NULL)
    {
        printf("\nAn error occured while allocating memory for uint8_vector\n");
        exit(1);
    }
    memset(vector.data, 0xff, sizeof(uint8_t) * size);
    vector.size = size;
    vector.count = 0;
    return vector;
}

void uint64_vector_append(uint64vec_t* vector, uint64_t data)
{
    // Appends values to uint8 vector, reallocates if necessary
    if ((vector->size + 1) < vector->capacity)
    {
        vector->data[vector->size++] = data;
        return;
    }
    vector->data = realloc(vector->data, (vector->capacity *= 2) * sizeof(uint64_t));
    if (vector->data == NULL)
    {
        printf("\nAn error occured while re-allocating memory for uint64_vector\n");
        exit(1);
    }
    vector->data[vector->size++] = data;
}

size_t approximate_size(uint64_t limit)
{
    int i;
    float x = 1;
    for (i = log10(limit); i > 0; i--)
    {
        x *= 2.4;
    }
    return x;
}

void segmented_sieve(uint64_t limit)
{
    // A detailed explaination of this algo can be found in the wiki mentioned above
    int64_t low, high, i = 3, j, k, n, s = 3, bits, count = (limit == 1) ? -1 : 0;
    size_t i_size, approx_arr_size = approximate_size(limit);
    uint64_t sqrtval = (uint64_t)sqrt(limit);
    uint64_t sieve_size = sqrtval < L1D_CACHE ? L1D_CACHE : sqrtval;            // This is a imitation of std::max()
    uint64_t segment_size = sieve_size * 16;
    uint64vec_t prime_arr = create_uint64_vector(approx_arr_size);               // Assuming the vectors to have 10 slots, no real formula to approximate exact amount
    uint64vec_t multiples = create_uint64_vector(approx_arr_size);
    uint8vec_t sieve = create_uint8_vector(segment_size);
    uint8vec_t is_prime = create_uint8_vector(sqrtval + 1);
    for (low = 0; low <= limit; low += segment_size)
    {
        memset(sieve.data, 0xff, sizeof(uint8_t) * sieve.size);
        high = low + segment_size - 1;
        high = high < limit ? high : limit;
        sieve_size = (high - low) / 16 + 1;
        for (; i * i <= high; i += 2)
        {
            if (is_prime.data[i])
            {
                for (j = i * i; j <= sqrtval; j += i)
                {
                    is_prime.data[j] = false;
                }
            }
        }
        for (; s * s <= high; s += 2)
        {
            if (is_prime.data[s])
            {
                uint64_vector_append(&prime_arr, s);
                uint64_vector_append(&multiples, s * s - low);
            }
        }
        for (i_size = 0; i_size < prime_arr.size; i_size++)
        {
            j = multiples.data[i_size];
            for (k = prime_arr.data[i_size] * 2; j < segment_size; j += k)
            {
                sieve.data[j >> 4] &= unset_bit[j & 15];
            }
            multiples.data[i_size] = j - segment_size;
        }
        if (high == limit)
        {
            bits = 0xff << (limit % 16 + 1) / 2;
            sieve.data[sieve_size - 1] &= ~bits;
        }
        for (n = 0; n < sieve_size; n++)
        {
            count += popcnt[sieve.data[n]];
        }
    }
    printf("\nFound primes: %" SCNd64, count);
}

int main()
{
    uint64_t N;
    printf("Enter your CPU's L1D_CACHE per thread (in bytes): ");
    scanf("%" SCNu64, &L1D_CACHE);
    printf("Enter upper limit for prime check: ");
    scanf("%" SCNu64, &N);
    clock_t t0 = clock();
    segmented_sieve(N);
    clock_t t1 = clock();
    double time_taken = (double)(t1 - t0) / CLOCKS_PER_SEC;
    printf("\nDone! Time taken: %f\n", time_taken);
    return 0;
}
