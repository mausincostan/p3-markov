CXX ?= g++
CXXFLAGS ?= -Wall -Werror -pedantic -g --std=c++17 -Wno-sign-compare -Wno-comment -fsanitize=address -fsanitize=undefined -D_GLIBCXX_DEBUG

# Run regression test
test: main.exe
	./main.exe < main_test.in > main_test.out
	diff main_test.out main_test.out.correct
	echo PASS

# Compile the main executable
main.exe: main.cpp
	$(CXX) $(CXXFLAGS) main.cpp -o $@

matrix_tests.exe: matrix_tests.cpp matrix.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

# Remove automatically generated files
clean :
	rm -rvf *.exe *~ *.out *.dSYM *.stackdump
