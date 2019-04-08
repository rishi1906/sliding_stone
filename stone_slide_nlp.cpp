// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "stone_slide_nlp.hpp"

#include <cassert>
#include <fstream>
#include <iostream>
#include <cmath>

using namespace Ipopt;

// define no of constraints and size of decision vector
#define no_of_cons (2 * (n / 3 - 1)) + 4
#define N_ 60

// define sie of stepsize for finite difference scheme to find gradient
const Number step_size = 1e-8;

// define computation of objective function to be used in finite differenc scheme
inline Number Obj_func(Number* X, Index n)
{

  Numer value = 0.0;

  return value;
}

// constructor
STONE_SLIDE_NLP::STONE_SLIDE_NLP() {}

// destructor
STONE_SLIDE_NLP::~STONE_SLIDE_NLP() {}

// returns the size of the problem
bool STONE_SLIDE_NLP::get_nlp_info(Index& n,          // size of problem
                                   Index& m,          // no of constraintsno of constraints
                                   Index& nnz_jac_g,  // no of non zero elements in jacobain 
                                   Index& nnz_h_lag,  // no of non zero elements in hessian
                                   IndexStyleEnum& index_style)
{
  // The problem described 
  n = N_;  

  // No of constraints
  m = no_of_cons;

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
                                      Number* x_l,  // lower limits for decision variables
                                      Number* x_u,  // uppe limits for decision variables
                                      Index m,      // no of constraints
                                      Number* g_l,  // lower limits for constraints
                                      Number* g_u)  // upper limits for constraints
{
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  assert(n == N_);
  assert(m == no_of_cons);

  // Lower bounds
  // for x
  for (Index i = 0; i < n; i++) {
    x_l[i] = -5.0;
  }
  // for v
  // for (Index i = 0; i < (2 * n) / 3; i++) {
  //    x_l[i] = 0.0;
  // }
  // // for u
  // for (Index i = 0; i < (2 * n) / 3; i++) {
  //    x_l[i] = 0.0;
  // }

  // Upper Bounds
  for (Index i = 0; i < n; i++) {
    x_u[i] = +5.0;
  }

  // Set bounds for constraints


  return true;
}

// returns the initial point for the problem
bool STONE_SLIDE_NLP::get_starting_point(Index n,           //
    bool init_x,       //
    Number* x,         //
    bool init_z,       //
    Number* z_L,       //
    Number* z_U,       //
    Index m,           //
    bool init_lambda,  //
    Number* lambda)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);
  //std::cout << "Initialization\n";
  // for distance x
  for (Index i = 0; i < n / 3; i++) {
    x[i] = i * h_k;
    //std::cout << "x[" << i << "]: " << x[i] << "\n";
  }
  //std::cout << std::endl;

  // hard coded
  /*x[10] = 0;
  x[11] = 0;
  x[12] = -10;
  x[13] = -36;
  x[14] = -72;
  x[15] = -120;
  x[16] = -180;
  x[17] = -252;
  x[18] = -336;
  x[19] = -432;*/

  // for veclocity v
  for (Index i = n / 3; i < ((2 * n) / 3); i++) {
    x[i] = 1;
    //std::cout << "x[" << i << "]: " << x[i] << "\n";
  }
  //std::cout << std::endl;

  // for acceleration u
  //double t = 0.0;
  for (Index i = ((2 * n) / 3); i < n; i++) {
    //x[i] = 6 - (12 * t);
    //t = t + h_k;
    x[i] = 0.0;
    //std::cout << "x[" << i << "]: " << x[i] << "\n";
  }
  //std::cout << std::endl;

  return true;
}

// returns the value of the objective function
bool STONE_SLIDE_NLP::eval_f(Index n,          //
                             const Number* x,  //
                             bool new_x,       //
                             Number& obj_value)
{
  assert(n == N_);

  // obj_value = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
  Index shift = (2 * n) / 3;
  obj_value = ((h_k) * (Square(x[shift]) + Square(x[shift + 1]))) / 2;
  for (Index k = shift + 1; k < n - 1; k++) {
    obj_value = obj_value + ((h_k) * (Square(x[k]) + Square(x[k + 1]))) / 2;
  }

  return true;
}

Number grad_at_x(Number Obj_func(Number* X, Index n), //
                 Number* X,//
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



// return the gradient of the objective function grad_{x} f(x)
bool STONE_SLIDE_NLP::eval_grad_f(Index n,          //
                                  const Number* x,  //
                                  bool new_x,       //
                                  Number* grad_f)
{
  assert(n == N_);
  // Approx Gradient using FDS
  // declare X a copy array of same size of x
  Number X[n];
  // make a copy of x in X
  for (Index k = 0; k < n; k++) {
    X[k] = x[k];
    //std::cout << X[k] << " ";
  }
  std::cout << std::endl;
  Number h = step_size;


  //std::cout << " approx gradient of function " << std::endl;
  //std::cout << "--------------------" << std::endl;
  // set gradt upto (2*n)/3 0
  for (Index at = 0; at < ((2 * n) / 3); at++) {
    grad_f[at] = 0;

    // std::cout << "grad[" << at << "]: " << grad_f[at] << "\n";
  }
  //std::cout << " approx gradient " << std::endl;
  //std::cout << "--------------------" << std::endl;
  for (Index at = ((2 * n) / 3); at < n; at++) {
    Number val = grad_at_x(Obj_func, X, at, n, h);
    grad_f[at] = val;

    // std::cout << "grad[" << at << "]: " << grad_f[at] << "\n";
  }
  // std::cout << "--------------------" << std::endl;


  // Exact Gradient with size n
  // //std::cout << "Exact Gradient\n";
  // // /set gradt upto (2 * n) / 3 0
  // for (Index at = 0; at < ((2 * n) / 3); at++) {
  //    grad_f[at] = 0;

  //    //std::cout << "grad[" << at << "]: " << grad_f[at] << "\n";
  // }
  // grad_f[((2 * n) / 3)] = h_k * x[((2 * n) / 3)];
  // grad_f[n - 1] = h_k * x[n - 1];
  // //std::cout << "grad[" << ((2 * n) / 3) << "]: " << grad_f[((2 * n) / 3)] << "\n";
  // for (Index at = ((2 * n) / 3) + 1; at < n - 1; at++) {
  //    grad_f[at] = 2 * h_k * x[at];
  //    //std::cout << "grad[" << at << "]: " << grad_f[at] << "\n";
  // }
  //std::cout << "grad[" << (n - 1) << "]: " << grad_f[n - 1] << "\n";



  ///*
  /* Exact Gradient with size n/3
  // set gradt upto (2*n)/3 0
  std::cout << "Exact Gradient\n";
  grad_f[0] = h_k * x[((2 * n) / 3)];
  grad_f[(n / 3) - 1] = h_k * x[n - 1];
  Index cnt = 1;
  std::cout << "grad[" << "0" << "]: " << grad_f[0] << "\n";
  for (Index at = ((2 * n) / 3) + 1 ; at < n - 1; at++) {
     grad_f[cnt++] = 2 * h_k * x[at];
     std::cout << "grad[" << (cnt - 1) << "]: " << grad_f[(cnt - 1)] << "\n";
  }
  std::cout << "grad[" << ((n / 3) - 1) << "]: " << grad_f[(n / 3) - 1] << "\n";
  */
  return true;
}

// return the value of the constraints: g(x)
bool STONE_SLIDE_NLP::eval_g(Index n,          //
                             const Number* x,  //
                             bool new_x,       //
                             Index m,          //
                             Number* g)
{
  assert(n == N_);
  assert(m == no_of_cons);

  Index k_th = 0;
  Index shift = n / 3;
  for (Index at = 0; at < (n / 3) - 1; at++) {
    g[k_th] =
      x[at + 1] - x[at] - (((h_k) * (x[at + shift + 1] + x[at + shift])) / 2);
    k_th = k_th + 1;
  }

  for (Index at = n / 3; at < ((2 * n) / 3) - 1; at++) {
    g[k_th] =
      x[at + 1] - x[at] - (((h_k) * (x[at + shift + 1] + x[at + shift])) / 2);
    k_th = k_th + 1;
  }
  g[k_th] = x[0];  // x(0) = 0
  g[++k_th] = x[(n / 3) - 1] - 1; // x(1) = 1
  g[++k_th] = x[n / 3]; // v(0) = 0
  g[++k_th] = x[(2 * (n / 3)) - 1]; // v(1) = 0
  //assert(k_th == no_of_cons);
  return true;
}

// return the structure or values of the Jacobian
bool STONE_SLIDE_NLP::eval_jac_g(Index n,          //
                                 const Number* x,  //
                                 bool new_x,       //
                                 Index m,          //
                                 Index nele_jac,   //
                                 Index* iRow,      //
                                 Index* jCol,      //
                                 Number* values)
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
                                        const Number* x,           //
                                        const Number* z_L,         //
                                        const Number* z_U,         //
                                        Index m,                   //
                                        const Number* g,           //
                                        const Number* lambda,      //
                                        Number obj_value,          //
                                        const IpoptData* ip_data,  //
                                        IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a
  // file,
  // etc
  // so we could use the solution.

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

}
