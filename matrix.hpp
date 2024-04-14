#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>

using namespace std;
using Eigen::MatrixXd;
class Matrix {
  public:
  
  Matrix(int grid_in)
    : grid_size(grid_in), 
      matrix_size(calc_matrix_size(grid_in)) {

        for (int i = 0; i < matrix_size; ++i) {
          points.push_back(i);
        }
      
      markov_matrix.resize(matrix_size, matrix_size);
      }

  // EFFECTS: calculates matrix_size by using grid_size
  //          matrix_size = (n+1)^2
  int calc_matrix_size(int n) {
    return pow((n + 1), 2);
  }

  // EFFECTS: grid_size getter
  int get_grid_size() {
    return grid_size;
  }

  // EFFECTS: matrix_size getter
  int get_matrix_size() {
    return matrix_size;
  }

  void print_all_points() {
    for (int p : points) {
      cout << p << " ";
    }
    cout << endl;
  }

  void print_markov_matrix_size() {
    cout << "Markov matrix size is " << markov_matrix.rows() 
          << "x" << markov_matrix.cols() << endl;
  }





  private:
  int grid_size;
  int matrix_size;

  vector<int> points; // vector of all points

  // vectors for different classification of points
  vector<int> bot_left_corner_points;
  vector<int> bot_right_corner_points;
  vector<int> top_left_corner_points;
  vector<int> top_right_corner;

  vector<int> top_edge;
  vector<int> bot_edge;
  vector<int> right_edge;
  vector<int> left_edge;

  vector<int> middle;

  // alleged markov matrix 
  MatrixXd markov_matrix;
  
  
};

#endif
