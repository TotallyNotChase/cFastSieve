/*
* @file = sieveprime.c
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

typedef struct bool_vector
{
    // For boolean values
    bool* data;
    uint64_t size, count;
} boolvec_t;

// The Level 1 Data cache for the user's CPU (must be per core)
uint64_t L1D_CACHE;

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

boolvec_t create_bool_vector(uint64_t size)
{
    // Creates a bool vector, assigns all slots to True, and returns it
    boolvec_t vector;
    vector.data = malloc(size * sizeof(bool));
    if (vector.data == NULL)
    {
        printf("\nAn error occured while allocating memory for bool_vector\n");
        exit(1);
    }
    memset(vector.data, true, sizeof(bool) * size);
    vector.size = size;
    vector.count = 0;
    return vector;
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

void uint64_vector_append(uint64vec_t *vector, uint64_t data)
{
    // Appends values to uint64 vector, reallocates if necessary
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

void segmented_sieve(uint64_t limit)
{
    int count = 1;
    // A detailed explaination of this algo can be found in the wiki mentioned above
    int64_t low, high, i = 3, j, k, n = 3, s = 3;
    size_t i_size, approx_arr_size = approximate_size(limit);
    uint64_t sqrtval = (uint64_t) sqrt(limit);
    uint64_t segment_size = sqrtval < L1D_CACHE ? L1D_CACHE : sqrtval;          // This is a imitation of std::max()
    uint64vec_t prime_arr = create_uint64_vector(approx_arr_size);              // An assumption on approx size
    uint64vec_t multiples = create_uint64_vector(approx_arr_size);
    boolvec_t sieve = create_bool_vector(segment_size);
    boolvec_t is_prime = create_bool_vector(sqrtval + 1);
    printf("2 ");
    for (low = 0; low <= limit; low += segment_size)
    {
        memset(sieve.data, true, sizeof(bool) * sieve.size);
        high = low + segment_size - 1;
        high = high < limit ? high : limit;
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
                sieve.data[j] = false;
            }
            multiples.data[i_size] = j - segment_size;
        }
        for (; n <= high; n += 2)
        {
            if (sieve.data[n - low])
            {
                printf("%" SCNu64 " ", n);
                count++;
            }
        }
    }
    printf("\nFound primes: %d", count);
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
    double time_taken = (double) (t1 - t0) / CLOCKS_PER_SEC;
    printf("\nDone! Time taken: %f\n", time_taken);
    return 0;
}
