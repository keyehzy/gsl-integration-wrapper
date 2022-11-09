CXX=clang++
CXXFLAGS= -g -O2 -W -Wall -Wextra -std=c++20 # -fsanitize=address -fanalyzer
LDFLAGS= -lgsl -lgslcblas -lm

main: main.cpp

run : main
	./$<

clean:
	$(RM) *~ *.o main
