#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <cmath>
#include <vector>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

using namespace std;
using namespace Eigen;

class Matrixx {
  public:
  
  // custom constructor, only ctor used for Matrix class
  Matrixx(int grid_in)
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

  void generate_markov_matrix() {
    int n = grid_size;
    double half = (1.0 / 2.0);
    //double third = 0.33; // for purposes of viewing couted version better
    double third = 1.0 / 3.0;
    double fourth = 1.0 / 4.0;
    
    // iterate for all p in points
    for (int p : points) {
      // shortcut vars for different directions points can go
      int up = p + (n + 1);
      int down = p - (n + 1);
      int left = p - 1;
      int right = p + 1;

      // corners (1/2)
      if (is_bot_left_corner(p)) {
        m(p, 1) = half;
        m(p, (n + 1)) = half;
        continue;
      }
      if (is_bot_right_corner(p)) {
        m(p, (n - 1)) = half;
        m(p, p + (n + 1)) = half;
        continue;
      }
      if (is_top_left_corner(p)) {
        int top_left_1 = (pow(n, 2) + n) + 1;
        int top_left_2 = (pow(n, 2) + n) - (n + 1);
        m(p, top_left_1) = half;
        m(p, top_left_2) = half;
        continue;
      }
      if (is_top_right_corner(p)) {
        int top_right_1 = (pow((n + 1), 2) - 1) - 1;
        int top_right_2 = (pow((n + 1), 2) - 1) - (n + 1);
        m(p, top_right_1) = half;
        m(p, top_right_2) = half;
        continue;
      }
  
      // edges (1/3)
      if (is_bot_edge(p)) {
        m(p, right) = third;
        m(p, left) = third;
        m(p, up) = third;
        continue;
      }
      if (is_top_edge(p)) {
        m(p, right) = third;
        m(p, left) = third;
        m(p, down) = third;
        continue;
      }
      if (is_left_edge(p)) {
        m(p, right) = third;
        m(p, up) = third;
        m(p, down) = third;
        continue;
      }
      if (is_right_edge(p)) {
        m(p, left) = third;
        m(p, up) = third;
        m(p, down) = third;
        continue;
      }

      // middle points (1/4)
      for (int point : middle_points) {
        if (point == p) {
          m(p, up) = fourth;
          m(p, down) = fourth;
          m(p, left) = fourth;
          m(p, right) = fourth;
        }
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

  /* WORKS BUT DOESNT ACTUALLY CALCULATE STEADY STATE VECTOR
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
  */

  /* FAILED ATTEMPT DOES NOT WORK
  VectorXd steady_state_vector() {
    int n = matrix_size; // Dimension of transition matrix P
    MatrixXd P = m;
    MatrixXd I = MatrixXd::Identity(n, n); // Identity matrix of same size as P
    MatrixXd Q = P - I; // Compute matrix Q
    VectorXd e = VectorXd::Ones(n); // Vector of all ones

    // Append vector e to Q
    Q.conservativeResize(Q.rows() + 1, NoChange);
    Q.bottomRows(1) = e.transpose();

    // Transpose Q
    MatrixXd QT = Q.transpose();
    std::cout << "Size of QT: " << QT.rows() << "x" << QT.cols() << std::endl;
    
    // Construct vector b with correct dimensions
    VectorXd b = VectorXd::Zero(QT.rows()); // Size of b should match the rows of QT
    std::cout << "Size of b: " << b.size() << std::endl;
    // Set the last element of b to 1
    b(n - 1) = 1;

    // Solve for x using least squares
    VectorXd x = QT.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);
    
    

    return x;
}
*/

  VectorXd computeSteadyState() {
    EigenSolver<MatrixXd> es(m.transpose());

    // Find the index of the eigenvalue with unit magnitude
    int index;
    double magnitude = 0.0;
    for (int i = 0; i < es.eigenvalues().rows(); ++i) {
        double currMagnitude = abs(es.eigenvalues()(i).real());
        if (abs(currMagnitude - 1.0) < 1e-9 && currMagnitude > magnitude) {
            magnitude = currMagnitude;
            index = i;
        }
    }
  
    // Extract the corresponding eigenvector
    VectorXd pi = es.eigenvectors().col(index).real();

    // Normalize the steady-state vector
    pi /= pi.sum();

    return pi;
}

  double sumVector(const VectorXd& v) {
      double sum = 0.0;
      for (int i = 0; i < v.size(); ++i) {
          sum += v(i);
      }
      return sum;
  }

  VectorXd vector_multiply(const VectorXd& v, double num) {
    return v * num;
  }


  MatrixXd get_markov_matrix() {
    return m;
  }

  void print_vector(VectorXd mat) {
    cout << mat << endl;
  }

  void test_master(int n) {
    sort_all_points();
    generate_markov_matrix();
    cout << "Markov Transition Matrix Below: \n\n";
    print_markov_matrix();
    cout << "\nSteady State Vector Below: \n\n";
    print_vector(computeSteadyState());
    cout << endl << endl;
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
