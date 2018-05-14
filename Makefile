all:
	g++ -o xor xor.cpp -mavx2 -O2 -std=c++11 -Wall


clean:
	rm -f xor
