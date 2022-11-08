CXX=clang++
CXXFLAGS=-ggdb -O2 -W -Wall -Wextra -std=c++20 # -fsanitize=address -fanalyzer
LDFLAGS=-lgsl -lcblas -lm

main: main.cpp

run : main
	./$<

clean:
	$(RM) *~ *.o main
