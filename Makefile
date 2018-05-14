all:
	g++ -o xor xor.cpp -O2 -mavx2 -std=c++11 -Wall -g


clean:
	rm -f xor
