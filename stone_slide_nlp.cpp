// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "stone_slide_nlp.hpp"
//#include "Differentiation_Matrix.hpp"

#include <cassert>
#include <fstream>
#include <iostream>
#include <cmath>
#include "Differentiation_Matrix.cpp"
using namespace Ipopt;

// define no of constraints and size of decision vector

//#define N_ 10 // no of grid points

// define size of stepsize for finite difference scheme to find gradient
const Number step_size = 1e-8;
#define PI acos(-1)
// define computation of objective function to be used in finite differenc scheme
inline Number Obj_func(Number* X, Index n)
{

  Number value = X[n - 1];

  return value;
}

//Default Constructor
STONE_SLIDE_NLP::STONE_SLIDE_NLP
(
  Index N
)
{
  N_ = N;
  T = define_time_stamps<Number, Index>(N_ + 1);
  for (Index i = 0; i <= N_; i++)
  {
    std::cout << "T[" << i << "]" << " : " << T[i] << "\n";
  }
}

// constructor
STONE_SLIDE_NLP::STONE_SLIDE_NLP() {}

// destructor
STONE_SLIDE_NLP::~STONE_SLIDE_NLP() {}


// returns the size of the problem
bool STONE_SLIDE_NLP::get_nlp_info(Index & n,         // size of problem
                                   Index & m,         // no of constraintsno of constraints
                                   Index & nnz_jac_g, // no of non zero elements in jacobain
                                   Index & nnz_h_lag, // no of non zero elements in hessian
                                   IndexStyleEnum & index_style)
{
  // The problem described
  n = (3 * (N_ + 1) ) + 2;

  // No of constraints
  m = (2 * (N_ + 1)) + 3;

  // size of jacobian matrix
  nnz_jac_g = m * n;

  // size of hessian of lagrangian
  // nnz_h_lag = 10;

  // style of indexing
  index_style = TNLP::C_STYLE; // use the C style indexing (0-based)

  return true;
}

// returns the variable bounds
bool STONE_SLIDE_NLP::get_bounds_info(Index n,      // size of problem
                                      Number * x_l, // lower limits for decision variables
                                      Number * x_u, // uppe limits for decision variables
                                      Index m,      // no of constraints
                                      Number * g_l, // lower limits for constraints
                                      Number * g_u) // upper limits for constraints
{
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  assert(n == (3 * (N_ + 1)) + 2);
  assert(m == (2 * (N_ + 1)) + 3);

  // Lower bounds
  // for X
  for (Index i = 0; i <= N_; i++) {
    x_l[i] = -10.0;
  }
  // for Y
  for (Index i = (N_ + 1); i <= ((2 * N_) + 1); i++) {
    x_l[i] = -40.0;
  }
  // // for Theta
  for (Index i = ((2 * N_) + 2); i <= ((3 * N_) + 3); i++) {
    x_l[i] = -PI;
  }
  //for tow
  x_l[n - 1] = -1.0;

  // Upper Bounds
  // for X
  for (Index i = 0; i <= N_; i++) {
    x_u[i] = +10.0;
  }
  // for Y
  for (Index i = (N_ + 1); i <= ((2 * N_) + 1); i++) {
    x_u[i] = +40.0;
  }
  // for Theta
  for (Index i = ((2 * N_) + 2); i <= ((3 * N_) + 3); i++) {
    x_u[i] = -PI;
  }
  //for tow
  x_u[n - 1] = +1.0;

  // set bounds on constraints

  return true;
}

// returns the initial point for the problem
bool STONE_SLIDE_NLP::get_starting_point(Index n,           //
    bool init_x,       //
    Number * x,        //
    bool init_z,       //
    Number * z_L,      //
    Number * z_U,      //
    Index m,           //
    bool init_lambda,  //
    Number * lambda)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);
  //std::cout << "Initialization\n";
  // for distance X
  for (Index i = 0; i <= N_; i++)
  {
    x[i] = T[i];
    //std::cout << T[i] << "," << x[i] << "\n";
  }
  // for Y
  for (Index i = (N_ + 1); i <= ((2 * N_) + 1); i++) {
    x[i] = T[i - (N_ + 1)];
    //std::cout << T[i - (N_ + 1)] << "," << x[i] << "\n";
  }
  // for Theta
  for (Index i = ((2 * N_) + 2); i <= ((3 * N_) + 3); i++) {
    x[i] = atan(x[i - (N_ + 1) ] / x[i - (2 * (N_ + 1))]);
    //std::cout << T[i - 2 * (N_ + 1)] << "," << x[i] << "\n";
  }
  //for tow f
  x[n - 1] = 0.5;
  return true;
}

// returns the value of the objective function
bool STONE_SLIDE_NLP::eval_f(Index n,          //
                             const Number * x, //
                             bool new_x,       //
                             Number & obj_value)
{
  assert(n == (3 * (N_ + 1)) + 2);


  return true;
}

/*
Number grad_at_x(Number Obj_func(Number * X, Index n), //
                 Number * X, //
                 Index pos,//
                 Index n, //
                 Number h)
{
  X[pos] = X[pos] + h;
  Number f_x_p_h = Obj_func(X, n);
  X[pos] = X[pos] - h;

  X[pos] = X[pos] - h;
  Number f_x_m_h = Obj_func(X, n);
  X[pos] = X[pos] + h;

  Number grad = (f_x_p_h - f_x_m_h) / (2.0 * h);
  return grad;
}
*/

// return the gradient of the objective function grad_{x} f(x)
bool STONE_SLIDE_NLP::eval_grad_f
(
  Index n, //
  const Number * x, //
  bool new_x,       //
  Number * grad_f)
{
  assert(n == (3 * (N_ + 1)) + 2);
  // Approx Gradient using FDS
  // declare X a copy array of same size of x
  /*
  Number X[n];
  // make a copy of x in X
  for (Index k = 0; k < n; k++) {
    X[k] = x[k];
    //std::cout << X[k] << " ";
  }
  std::cout << std::endl;
  Number h = step_size;

  for (Index at = 0; at < n; at++) {
    Number val = grad_at_x(Obj_func, X, at, n, h);
    grad_f[at] = val;

    // std::cout << "grad[" << at << "]: " << grad_f[at] << "\n";
  }
  // std::cout << "--------------------" << std::endl;

  */

  return true;
}

// return the value of the constraints: g(x)
bool STONE_SLIDE_NLP::eval_g(Index n,          //
                             const Number * x, //
                             bool new_x,       //
                             Index m,          //
                             Number * g)
{
  assert(n == (3 * (N_ + 1)) + 2);
  assert(m == (2 * (N_ + 1)) + 3);

  return true;
}

// return the structure or values of the Jacobian
bool STONE_SLIDE_NLP::eval_jac_g(Index n,          //
                                 const Number * x, //
                                 bool new_x,       //
                                 Index m,          //
                                 Index nele_jac,   //
                                 Index * iRow,     //
                                 Index * jCol,     //
                                 Number * values)
{
  if (values == NULL) {
    // return the structure of the Jacobian

    // this particular Jacobian is dense
    Index nnz = 0;
    for (Index i = 0; i < m; i++) {
      for (Index j = 0; j < n; j++) {
        iRow[nnz] = i;
        jCol[nnz] = j;
        nnz += 1;
      }
    }
  }
  // else
  // {
  //    // return the values of the Jacobian of the constraintss
  // }

  return true;
}

/*
//return the structure or values of the Hessian
bool STONE_SLIDE_NLP::eval_h(
   Index         n,
   const Number* x,
   bool          new_x,
   Number        obj_factor,
   Index         m,
   const Number* lambda,
   bool          new_lambda,
   Index         nele_hess,
   Index*        iRow,
   Index*        jCol,
   Number*       values
   )
{
   if( values == NULL )
   {
      // return the structure. This is a symmetric matrix, fill the lower left
      // triangle only.

      // the hessian for this problem is actually dense
      Index idx = 0;
      for( Index row = 0; row < 4; row++ )
      {
         for( Index col = 0; col <= row; col++ )
         {
            iRow[idx] = row;
            jCol[idx] = col;
            idx++;
         }
      }

      assert(idx == nele_hess);
   }
   else
   {
      // return the values. This is a symmetric matrix, fill the lower left
      // triangle only

      // fill the objective portion

   }

   return true;
}
*/

void STONE_SLIDE_NLP::finalize_solution(SolverReturn status,       //
                                        Index n,                   //
                                        const Number * x,          //
                                        const Number * z_L,        //
                                        const Number * z_U,        //
                                        Index m,                   //
                                        const Number * g,          //
                                        const Number * lambda,     //
                                        Number obj_value,          //
                                        const IpoptData * ip_data, //
                                        IpoptCalculatedQuantities * ip_cq)
{
  // here is where we would store the solution to variables, or write to a
  // file,
  // etc
  // so we could use the solution.
  /*
  // For this example, we write the solution to the console
  std::cout << std::endl
            << std::endl
            << "Numerical Solution " << std::endl;

  // std::cout << std::endl << "Time t " << std::endl;
  // for (Index i = 0; i < n / 3; i++) {
  //   std::cout << "t[" << i << "] = " << i << std::endl;
  // }

  std::ofstream myfile;
  myfile.open("output.txt");
  // myfile << "Writing this to a file.\n";

  std::cout << std::endl << "Distance x " << std::endl;
  for (Index i = 0; i < n / 3; i++) {
    std::cout << "x[" << i << "] = " << x[i] << std::endl;
    myfile << "" << i * h_k << "," << x[i] << "\n";
  }

  std::cout << std::endl << "Velocity v " << std::endl;
  for (Index i = n / 3; i < ((2 * n) / 3); i++) {
    std::cout << "v[" << i << "] = " << x[i] << std::endl;
  }

  std::cout << std::endl << "Acc u " << std::endl;
  for (Index i = ((2 * n) / 3); i < n; i++) {
    std::cout << "u[" << i << "] = " << x[i] << std::endl;
  }

  std::cout << std::endl;

  // std::cout << std::endl << std::endl << "Solution of the bound
  // multipliers,
  // z_L and z_U" << std::endl;
  // for ( Index i = 0; i < n; i++ )
  // {
  //    std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
  // }
  // for ( Index i = 0; i < n; i++ )
  // {
  //    std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
  // }

  std::cout << std::endl << std::endl << "Objective value" << std::endl;
  std::cout << "f(x*) = " << obj_value << std::endl;

  // std::cout << std::cout << std::endl
  //          << "Final value of the constraints:" << std::endl;
  //for (Index i = 0; i < m; i++) {
  //   std::cout << "g(" << i << ") = " << g[i] << std::endl;
  // }
  //
  //

  // Analytical Solution
  std::cout << std::endl
            << std::endl
            << "Numerical Solution " << std::endl;

  // std::cout << std::endl << "Time t " << std::endl;
  // for ( Index i = 0; i < n / 3; i++ )
  // {
  //    std::cout << "t[" << i << "] = " << i << std::endl;
  // }

  std::cout << std::endl << "Distance x " << std::endl;
  for ( double i = 0.0; i < 1; i += 0.1 )
  {
    //myfile << "" << i << "," << (3 * (Square(i))) - ((2 * i) * (Square(i))) << "\n";
    std::cout << "x[" << i << "] = " << (3 * (Square(i))) - ((2 * i) * (Square(i))) << std::endl;
  }

  std::cout << std::endl << "Velocity v " << std::endl;
  for ( double i = 0.0; i < 1; i += 0.1 )
  {
    std::cout << "v[" << i << "] = " << (((-6) * (Square(i))) + (6 * i))
              <<
              std::endl;
  }

  std::cout << std::endl << "Acc u " << std::endl;
  for ( double i = 0.0; i < 1; i += 0.1 )
  {
    std::cout << "u[" << i << "] = " << (((-12)*i) + 6) << std::endl;
  }
  myfile.close();
  */
}
