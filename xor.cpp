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

#define XOR_MASK	_mm256_set_epi32(0x77777777,0x77777777, 0x77777777, 0x77777777, 0x77777777, 0x77777777, 0x77777777, 0x77777777)

#define SET_MM256(a)								\
		_mm256_set_epi8(						\
		a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], 		\
		a[8], a[9], a[10], a[11], a[12], a[13], a[14], a[15], 		\
		a[16], a[17], a[18], a[19], a[20], a[21], a[22], a[23], 	\
		a[24], a[25], a[26], a[27], a[28], a[29], a[30], a[31]	 	\
		)

#define SETR_MM256(a)								\
		_mm256_setr_epi8(						\
		a[31], a[30], a[29], a[28], a[27], a[26], a[25], a[24], a[23], 	\
		a[22], a[21], a[20], a[19], a[18], a[17], a[16], a[15], 	\
		a[14], a[13], a[12], a[11], a[10], a[9], a[8], a[7], 		\
		a[6], a[5], a[4], a[3], a[2], a[1], a[0]			\
		)

#define NUM_TRANSACTIONS	10000000
#define LENGTH			64
#define REV_COMP_ASM
#define REV_COMP_LL

int64_t times[NUM_TRANSACTIONS] = {0};

#ifdef REV_COMP1
__m256i reverse_complement(__m256i a,__m256i b) {
	return _mm256_xor_si256 (a,b);
}
#endif

#ifdef REV_COMP_NAIVE_BYTE
void reverse_complement_naive(char *a, int len) {
	for (int i = 0; i < len; i++)
		a[i] ^= 0x77;
}
#elif defined(REV_COMP_LL)
void reverse_complement_naive(char *a, int len) {
	uint64_t *_a = (uint64_t*) a;
	for (int i = 0; i < len/8; i++) {
		_a[i] ^= 0x7777777777777777;
	}
}
#endif

struct Kmer {
	union {
		char data[LENGTH]; /* 512 bits */
		struct {
			union{
				char data1[LENGTH/2]; /* 256 bits */
				__m256i d1;
			};
			union{
				char data2[LENGTH/2]; /* 256 bits */
				__m256i d2;
			};
		};
	};
};

struct Kmer kmer_data[NUM_TRANSACTIONS];

#ifdef REV_COMP1
static __attribute__((noinline)) void calculate_complement(void)
{
    for (int i = 0; i < NUM_TRANSACTIONS; i++) {
	    reverse_complement(kmer_data[i].d1, XOR_MASK);
	    reverse_complement(kmer_data[i].d2, XOR_MASK);
    }
}
#endif

int main(int argc, char **argv) {
    std::string s(LENGTH,0);
    int i;
    uint64_t start, end;
    __m256i b = XOR_MASK;

    std::srand(0);

    for (i = 0; i < NUM_TRANSACTIONS; i++) {
	    std::generate_n(s.begin(), LENGTH, std::rand);
	    memcpy(kmer_data[i].data1, s.data(), s.length()/2);
	    memcpy(kmer_data[i].data2, s.data()+(LENGTH/2), s.length()/2);
    }

    start = RDTSC_START();

#ifdef REV_COMP1
    for (i = 0; i < NUM_TRANSACTIONS; i++) {
	    reverse_complement(kmer_data[i].d1, XOR_MASK);
	    reverse_complement(kmer_data[i].d2, XOR_MASK);
    }
#elif defined(REV_COMP2)
    calculate_complement();
#elif defined(REV_COMP_ASM)
    for (i = 0; i < NUM_TRANSACTIONS; i++) {
	    asm volatile("vpxor %%ymm0, %%ymm0, %%ymm1"
				: "=x" (kmer_data[i].d1)
				: "x"(kmer_data[i].d1), "x"(b));
	    asm volatile("vpxor %%ymm0, %%ymm0, %%ymm1"
				: "=x" (kmer_data[i].d2)
				: "x"(kmer_data[i].d2), "x"(b));
    }
#endif
    end = RDTSCP();
    cout << "(AVX2 xor) Time taken for " << NUM_TRANSACTIONS << " transactions: " << (float)(end - start)/NUM_TRANSACTIONS << endl;

    start = RDTSC_START();
    for (i = 0; i < NUM_TRANSACTIONS; i++) {
	    reverse_complement_naive(kmer_data[i].data, LENGTH);
    }
    end = RDTSCP();
    cout << "(byte-wise xor) Time taken for " << NUM_TRANSACTIONS << " transactions: " << (float)(end - start)/NUM_TRANSACTIONS << endl;
}
