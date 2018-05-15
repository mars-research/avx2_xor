#include <immintrin.h>
#include <stdio.h>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cstdint>
#include "rdtsc_helper.h"
#include "stat.h"

using namespace std;

#define XOR_MASK	_mm256_set_epi32(0x77777777, 0x77777777, 0x77777777, 0x77777777, 0x77777777, 0x77777777, 0x77777777, 0x77777777)

#define NUM_TRANSACTIONS	10000000
#define LENGTH			64
#define REV_COMP_ASM
#define REV_COMP_LL

int64_t times[NUM_TRANSACTIONS] = { 0 };

struct Kmer {
	union {
		char data[LENGTH];	/* 512 bits */
		struct {
			union {
				char data1[LENGTH / 2];	/* 256 bits */
				__m256i d1;
			};
			union {
				char data2[LENGTH / 2];	/* 256 bits */
				__m256i d2;
			};
		};
	};
};

// for storing kmer data
struct Kmer kmer_data[NUM_TRANSACTIONS];

// uses the native xor instruction to xor 8-bit or 64-bit wide
// data according to the last parameter
void reverse_complement_naive(char *a, int len, int width)
{
	switch (width) {
	case 8:
		{
			for (int i = 0; i < len; i++)
				a[i] ^= 0x77;
			break;
		}
	case 64:
		{
			uint64_t *_a = (uint64_t *) a;
			for (int i = 0; i < len / 8; i++) {
				_a[i] ^= 0x7777777777777777;
			}
			break;
		}
	default:
		cerr << "Unknown bitwidth requested" << endl;
		break;
	}
}

// gcc intrinisic for computing 256-bit avx2 xor
inline __m256i reverse_complement_intrinsic(__m256i a, __m256i b)
{
	return _mm256_xor_si256(a, b);
}

// generate random data from std::rand
void generate_random_data(void)
{
	int i;
	std::string s(LENGTH, 0);

	std::srand(0);
	cout << "generating random data ... " << endl;
	for (i = 0; i < NUM_TRANSACTIONS; i++) {
		std::generate_n(s.begin(), LENGTH, std::rand);
		memcpy(kmer_data[i].data1, s.data(), s.length() / 2);
		memcpy(kmer_data[i].data2, s.data() + (LENGTH / 2),
		       s.length() / 2);
	}
}

void do_naive_xor8(void)
{
	int i;
	uint64_t start, end;

	start = RDTSC_START();
	for (i = 0; i < NUM_TRANSACTIONS; i++) {
		reverse_complement_naive(kmer_data[i].data, LENGTH, 8);
	}
	end = RDTSCP();
	cout << "(xor 8-bit) Time taken for " << NUM_TRANSACTIONS <<
	    " transactions: " << (float)(end -
					 start) / NUM_TRANSACTIONS << endl;
}

void do_naive_xor64(void)
{
	int i;
	uint64_t start, end;

	start = RDTSC_START();
	for (i = 0; i < NUM_TRANSACTIONS; i++) {
		reverse_complement_naive(kmer_data[i].data, LENGTH, 64);
	}
	end = RDTSCP();
	cout << "(xor 64bit) Time taken for " << NUM_TRANSACTIONS <<
	    " transactions: " << (float)(end -
					 start) / NUM_TRANSACTIONS << endl;
}

void do_simd_xor_intrinsic(void)
{
	int i;
	uint64_t start, end;

	start = RDTSC_START();
	for (i = 0; i < NUM_TRANSACTIONS; i++) {
		reverse_complement_intrinsic(kmer_data[i].d1, XOR_MASK);
		reverse_complement_intrinsic(kmer_data[i].d2, XOR_MASK);
	}
	end = RDTSCP();
	cout << "(xor avx2-gccintrinsic) Time taken for " << NUM_TRANSACTIONS <<
	    " transactions: " << (float)(end -
					 start) / NUM_TRANSACTIONS << endl;
	cout <<
	    "Warning: Do NOT trust the above number as the compiler reorders" <<
	    " the vpxor instruction after rdtscp. Barriers didn't help!" <<
	    endl;
}

void do_simd_xor_asm(void)
{
	int i;
	uint64_t start, end;
	__m256i b = XOR_MASK;
	start = RDTSC_START();

	for (i = 0; i < NUM_TRANSACTIONS - 1; i = i + 2) {
		asm volatile ("vpxor %%ymm1, %%ymm0, %%ymm0":"=x"
			      (kmer_data[i].d1):"x"(kmer_data[i].d1), "x"(b));
		asm volatile ("vpxor %%ymm2, %%ymm3, %%ymm3":"=x"
			      (kmer_data[i].d2):"x"(kmer_data[i].d2), "x"(b));
		asm volatile ("vpxor %%ymm4, %%ymm5, %%ymm5":"=x"
			      (kmer_data[i + 1].d1):"x"(kmer_data[i + 1].d1),
			      "x"(b));
		asm volatile ("vpxor %%ymm6, %%ymm7, %%ymm7":"=x"
			      (kmer_data[i + 1].d2):"x"(kmer_data[i + 1].d2),
			      "x"(b));
	}
	end = RDTSCP();
	cout << "(AVX2 xor) Time taken for " << NUM_TRANSACTIONS <<
	    " transactions: " << (float)(end -
					 start) / NUM_TRANSACTIONS << endl;

}

int main(int argc, char **argv)
{
	generate_random_data();
	do_naive_xor8();
	do_naive_xor64();
	do_simd_xor_intrinsic();
	do_simd_xor_asm();
}
