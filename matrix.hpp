#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <cmath>
#include <vector>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

using namespace std;
using Eigen::MatrixXd;

class Matrix {
  public:
  
  // custom constructor, only ctor used for Matrix class
  Matrix(int grid_in)
    : grid_size(grid_in), 
      matrix_size(calc_matrix_size(grid_in)) {

        for (int i = 0; i < matrix_size; ++i) {
          points.push_back(i);
        }
      
      m.resize(matrix_size, matrix_size);
      fill_markov_matrix(0);
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

  // EFFECTS: prints all points
  // MODIFIES: cout
  void print_all_points() {
    for (int p : points) {
      cout << p << " ";
    }
    cout << endl;
  }

  // EFFECTS: prints all different vectors for corners, edges, and middle
  // MODIFIES: cout
  void print_all_vectors() {
    // corners
    cout << endl << "bot right corner: ";
    for (int i : bot_right_corner_points) { cout << i << " "; }
    cout << endl << "bot left corner: ";
    for (int i : bot_left_corner_points) { cout << i << " "; }
    cout << endl << "top right corner: ";
    for (int i : top_right_corner_points) { cout << i << " "; }
    cout << endl << "top left corner: ";
    for (int i : top_left_corner_points) { cout << i << " "; }

    // edges
    cout << endl << "bottom edge: ";
    for (int i : bot_edge_points) { cout << i << " "; }
    cout << endl << "top edge: ";
    for (int i : top_edge_points) { cout << i << " "; }
    cout << endl << "left edge: ";
    for (int i : left_edge_points) { cout << i << " "; }
    cout << endl << "right edge: ";
    for (int i : right_edge_points) { cout << i << " "; }

    // middle
    cout << endl << "middle points: ";
    for (int i : middle_points) { cout << i << " ";}
  }

  // EFFECTS: print markov matrix size
  // MODIFIES: cout
  void print_markov_matrix_size() {
    cout << "Markov matrix size is " << m.rows() 
          << "x" << m.cols() << endl;
  }

  // EFFECTS: fill markov matrix with num
  //          used in Matrix constructor to init m to fill with 0
  void fill_markov_matrix(int num) {
    for (int i = 0; i < matrix_size; ++i) {
      for (int j = 0; j < matrix_size; ++j) {
        m(i, j) = num;
      }
    }
  }

  // EFFECTS: prints markov matrix
  // MODIFIES: cout
  void print_markov_matrix() {
    cout << m << endl;
  }

  // EFFECTS: sorts all points from points vector into different categories:
  //          corners, eedges, and the rest are middle points
  // MODIFIES: corner, edges, and middle vectors
  void sort_all_points() {
    for (int i = 0; i < matrix_size; ++i) {
      // corners
      if (is_bot_left_corner(i)) { bot_left_corner_points.push_back(i); }
      else if (is_bot_right_corner(i)) { bot_right_corner_points.push_back(i); }
      else if (is_top_left_corner(i)) { top_left_corner_points.push_back(i); }
      else if (is_top_right_corner(i)) { top_right_corner_points.push_back(i); }

      // edges
      else if (is_left_edge(i)) { left_edge_points.push_back(i); }
      else if (is_right_edge(i)) { right_edge_points.push_back(i); }
      else if (is_bot_edge(i)) { bot_edge_points.push_back(i); }
      else if (is_top_edge(i)) { top_edge_points.push_back(i); }

      // middle (the rest and the majority of points are middle)
      else {
        middle_points.push_back(i);
      }
    }
  }

  // conditionals to sort points
  bool is_bot_left_corner(int num) {
    return num == 0;
  }
  bool is_bot_right_corner(int num) {
    return num == grid_size;
  }
  bool is_top_left_corner(int num) {
    int n = grid_size;
    return num == (n * (n + 1));
  }
  bool is_top_right_corner(int num) {
    int n = grid_size;
    return num == ((pow(n + 1, 2)) - 1);
  }
  bool is_left_edge(int p) {
    int n = grid_size;
    int condition = p % (n + 1);
    return condition == 0;
  }
  bool is_right_edge(int p) {
    int n = grid_size;
    int condition = p % (n + 1);
    return condition == n;
  }
  bool is_bot_edge(int p) {
    int n = grid_size;
    return (0 <= p) && (p <= n);
  }
  bool is_top_edge(int p) {
    int n = grid_size;
    int left_condition = n * (n + 1);
    int right_condition = pow((n + 1), 2) - 1;
    return (left_condition <= p) && (p <= right_condition);
  }

  bool is_corner(int num) {
    return is_bot_left_corner(num) || is_bot_right_corner(num) ||
           is_top_left_corner(num) || is_top_right_corner(num); 
  }


  Eigen::VectorXd get_coordinate_distribution() {
    Eigen::EigenSolver<Eigen::MatrixXd> solver(m);
    Eigen::VectorXd eigenvalues_real = solver.eigenvalues().real();
    double epsilon = 1e-6;
    int index = -1;
    for (int i = 0; i < eigenvalues_real.size(); ++i) {
        if (std::abs(eigenvalues_real(i) - 1.0) < epsilon) {
            index = i;
            break;
        }
    }
    return solver.eigenvectors().col(index).real();
  }


  private:
  int grid_size;
  int matrix_size;

  vector<int> points; // vector of all points

  // vectors for different classification of points
  vector<int> bot_left_corner_points;
  vector<int> bot_right_corner_points;
  vector<int> top_left_corner_points;
  vector<int> top_right_corner_points; 

  vector<int> top_edge_points;
  vector<int> bot_edge_points;
  vector<int> right_edge_points;
  vector<int> left_edge_points; 

  vector<int> middle_points; 

  // alleged markov matrix 
  MatrixXd m;
  
  
};

#endif
