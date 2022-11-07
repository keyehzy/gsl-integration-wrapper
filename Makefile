CXX=clang++
CXXFLAGS=-ggdb -O0 -W -Wall -Wextra # -fsanitize=address -fanalyzer
LDFLAGS=-lgsl -lcblas -lm

main: main.cpp

run : main
	./$<

clean:
	$(RM) *~ *.o main
