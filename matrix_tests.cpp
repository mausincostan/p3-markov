#include "matrix.hpp"
#include "unit_test_framework.hpp"

using namespace std;
/*
TEST(test_matrix_size) {
  Matrix matrix(9);

  ASSERT_EQUAL(matrix.get_matrix_size(), 100);
  matrix.print_all_points();
  matrix.print_markov_matrix_size();
}

TEST(test_print_markov_matrix) {
  Matrix matrix(9);

  matrix.print_all_points();
  matrix.sort_all_points();
  matrix.print_all_vectors();
}
*/
TEST(test_master) {
  Matrix m(9);

  m.print_all_points();
  m.sort_all_points();
  m.print_all_vectors();
  m.generate_markov_matrix();
  cout << endl << endl;
  m.print_markov_matrix();
}



TEST_MAIN()