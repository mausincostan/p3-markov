#include "matrix.hpp"
#include "unit_test_framework.hpp"

using namespace std;

TEST(test_matrix_size) {
  Matrix matrix(9);

  ASSERT_EQUAL(matrix.get_matrix_size(), 100);
  matrix.print_all_points();
  matrix.print_markov_matrix_size();
}



TEST_MAIN()