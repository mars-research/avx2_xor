all:
	g++ -o xor xor.cpp -O2 -mavx2 -std=c++11 -Wall -g

clean:
	rm -f xor


indent:
	indent -npro -kr -i8 -ts8 -sob -l80 -ss -ncs -cp1 ./xor.cpp
