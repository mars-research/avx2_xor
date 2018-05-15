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

#ifdef REV_COMP1
__m256i reverse_complement(__m256i a, __m256i b)
{
	return _mm256_xor_si256(a, b);
}
#endif

#ifdef REV_COMP_NAIVE_BYTE
void reverse_complement_naive(char *a, int len)
{
	for (int i = 0; i < len; i++)
		a[i] ^= 0x77;
}
#elif defined(REV_COMP_LL)
void reverse_complement_naive(char *a, int len)
{
	uint64_t *_a = (uint64_t *) a;
	for (int i = 0; i < len / 8; i++) {
		_a[i] ^= 0x7777777777777777;
	}
}
#endif

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

struct Kmer kmer_data[NUM_TRANSACTIONS];

#ifdef REV_COMP1
static __attribute__ ((noinline))
void calculate_complement(void)
{
	for (int i = 0; i < NUM_TRANSACTIONS; i++) {
		reverse_complement(kmer_data[i].d1, XOR_MASK);
		reverse_complement(kmer_data[i].d2, XOR_MASK);
	}
}
#endif

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

int main(int argc, char **argv)
{
	int i;
	uint64_t start, end;
	__m256i b = XOR_MASK;

	generate_random_data();

	start = RDTSC_START();

#ifdef REV_COMP1
	for (i = 0; i < NUM_TRANSACTIONS; i++) {
		reverse_complement(kmer_data[i].d1, XOR_MASK);
		reverse_complement(kmer_data[i].d2, XOR_MASK);
	}
#elif defined(REV_COMP2)
	calculate_complement();
#elif defined(REV_COMP_ASM)
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
#endif
	end = RDTSCP();
	cout << "(AVX2 xor) Time taken for " << NUM_TRANSACTIONS <<
	    " transactions: " << (float)(end -
					 start) / NUM_TRANSACTIONS << endl;

	start = RDTSC_START();
	for (i = 0; i < NUM_TRANSACTIONS; i++) {
		reverse_complement_naive(kmer_data[i].data, LENGTH);
	}
	end = RDTSCP();
	cout << "(byte-wise xor) Time taken for " << NUM_TRANSACTIONS <<
	    " transactions: " << (float)(end -
					 start) / NUM_TRANSACTIONS << endl;
}
