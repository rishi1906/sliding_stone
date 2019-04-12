// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "stone_slide_nlp.hpp"
#include "Differentiation_Matrix.cpp"
#include <cassert>
#include <fstream>
#include <iostream>
#include <cmath>
//#include "Differentiation_Matrix.cpp"
using namespace Ipopt;

// define no of constraints and size of decision vector

//#define N_ 10 // no of grid points

// define size of stepsize for finite difference scheme to find gradient
//const Number step_size = 1e-8;

const Number grav = 1.0; // value of gravitational constant
const Number t0 = 0.00;
const Number tf = 3.0;
#define PI acos(-1)
Number RAD_to_DEG(Number r) { return r * 180.0 / PI; }

/*
// define computation of objective function to be used in finite differenc scheme
inline Number Obj_func(Number* X, Index n)
{

  Number value = X[n - 1];

  return value;
}*/

//Default Constructor
STONE_SLIDE_NLP::STONE_SLIDE_NLP
(
  Index N
)
{
  N_ = N;
  T = define_time_stamps<Number, Index>(N_ + 1);
  //T = define_time_stamps(N_ + 1);
  /*for (Index i = 0; i <= N_; i++)
  {
    std::cout << "T[" << i << "]" << " : " << T[i] << "\n";
  }*/
}

// constructor
STONE_SLIDE_NLP::STONE_SLIDE_NLP() {}

// destructor
STONE_SLIDE_NLP::~STONE_SLIDE_NLP() {}


// returns the size of the problem
bool STONE_SLIDE_NLP::get_nlp_info
(
  Index & n,         // size of problem
  Index & m,         // no of constraintsno of constraints
  Index & nnz_jac_g, // no of non zero elements in jacobain
  Index & nnz_h_lag, // no of non zero elements in hessian
  IndexStyleEnum & index_style
)
{
  // size of problem
  n = (3 * (N_ + 1) ) + 1;

  // size of constraints
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
bool STONE_SLIDE_NLP::get_bounds_info
(
  Index n,      // size of problem

  Number * x_l, // lower limits for decision variables
  Number * x_u, // uppe limits for decision variables
  Index m,      // no of constraints
  Number * g_l, // lower limits for constraints
  Number * g_u  // upper limits for constraints
) // upper limits for constraints
{
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  //assert(n == (3 * (N_ + 1)) + 1);
  //assert(m == (2 * (N_ + 1)) + 3);

  // Lower bounds
  // for X
  for (Index i = 0; i <= N_; i++) {
    x_l[i] = -5.0;
  }
  // for Y
  for (Index i = (N_ + 1); i <= ((2 * N_) + 1); i++) {
    x_l[i] = -5.0;
  }

  // // for Theta
  for (Index i = ((2 * N_) + 2); i <= ((3 * N_) + 2); i++) {
    x_l[i] = -5.0;
  }

  //for tow
  x_l[n - 1] = 0.00;

  // Upper Bounds
  // for X
  for (Index i = 0; i <= N_; i++) {
    x_u[i] = +5.0;
  }
  // for Y
  for (Index i = (N_ + 1); i <= ((2 * N_) + 1); i++) {
    x_u[i] = +5.0;
  }

  // for Theta
  for (Index i = ((2 * N_) + 2); i <= ((3 * N_) + 2); i++) {
    x_u[i] = +5.0;
  }

  //for tow
  x_u[n - 1] = +5.0;

  // set bounds on constraints for ineuality constraints
  /*g_l[m - 1] = 0.0;
  g_u[m - 1] = 0.5;*/

  return true;
}

// returns the initial point for the problem
bool STONE_SLIDE_NLP::get_starting_point
(
  Index n,           //
  bool init_x,       //
  Number * x,        //
  bool init_z,       //
  Number * z_L,      //
  Number * z_U,      //
  Index m,           //
  bool init_lambda,  //
  Number * lambda
)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);


  std::ofstream myfile;


  std::vector<Number > time(N_ + 1);

  for (Index i = 0; i <= N_; i++)
  {
    time[i] = ((tf - t0) * (T[i]) + (tf + t0)) / 2.0;
    //std::cout << "time[" << i << "]" << " : " << time[i] << std::endl;
  }
  //myfile << "Initialization\n";
  //myfile << "X\n";
  myfile.open("init_X.txt");
  //for X
  for (Index i = 0; i <= N_; i++)
  {
    //x[i] = time[i];

    //x[i] = ( ((grav * tf) / (PI)) * ((time[i]) - ((tf / PI) * sin((PI) * (1.0 - (time[i] / tf)) ))) );
    x[i] = ((tf - t0) * i) / (4 * N_);
    //x[i] = i * delta1 + 0.001;
    //std::cout << T[i] << "," << x[i] << "\n";
    myfile << time[i] << "," << x[i] << "\n";
  }
  //myfile << "\n";
  myfile.close();

  //myfile << "Y\n";
  myfile.open("init_Y.txt");
  // for Y
  //Number delta2 = 0.3 / (N_ + 1);
  for (Index i = 0; i <= N_; i++) {
    x[i + ((N_ + 1))] = time[i];
    //x[i] = ((tf-t0)*i)/(4*N_);
    //x[i + ((N_ + 1))] =  ( (2.0 * grav * tf * tf / (PI * PI)) * ((cos((PI / 2.0) * (1.0 - (time[i] / tf)))) * (cos((PI / 2.0) * (1.0 - (time[i] / tf))))) );
    //x[i + ((N_ + 1))] = 0.03 * (N_ - i);
    //std::cout << T[i - (N_ + 1)] << "," << x[i] << "\n";
    myfile << time[i] << "," << x[i + ((N_ + 1))] << "\n";
  }
  //myfile << "\n";
  myfile.close();

  //myfile << "Theta\n";
  myfile.open("init_theta.txt");
  // for Theta
  for (Index i = 0; i <= N_; i++) {
    x[i + ((2 * N_) + 2)] = PI / 4;
    //x[i + ((2 * N_) + 2)] = atan((x[(N_ + 1)]) / (x[(2 * N_ + 2) ]));
    //x[i + ((2 * N_) + 2)] = (PI / 2.0) * (1.0 - (time[i] / tf));
    //std::cout << T[i - 2 * (N_ + 1)] << "," << x[i] << "\n";
    myfile << time[i] << "," << x[i + ((2 * N_) + 2)] << "\n";
  } //myfile << "\n";
  myfile.close();
  //for tow f
  x[n - 1] = tf;
  //std::cout << x[n - 1] << "\n";
  return true;
}

// returns the value of the objective function
bool STONE_SLIDE_NLP::eval_f
(
  Index n,          //
  const Number * x, //
  bool new_x,       //
  Number & obj_value
)
{
  //assert(n == (3 * (N_ + 1)) + 1);

  obj_value = x[n - 1];
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
  //assert(n == (3 * (N_ + 1)) + 1);

  for (Index i = 0; i <= n - 2; i++)
  {
    grad_f[i] = 0;
  }

  grad_f[n - 1] = 1;
  return true;
}

// return the value of the constraints: g(x)
bool STONE_SLIDE_NLP::eval_g
(
  Index n,          //

  const Number * x, //
  bool new_x,       //
  Index m,          //
  Number * g
)
{
  //assert(n == (3 * (N_ + 1)) + 1);
  //assert(m == (2 * (N_ + 1)) + 3);

  Index nth = 0 ;
  //Number tf = 1.25;
  //Number c1 = 2.0 / x[n - 1], c2 = sqrt(2 * grav);
  Number c1 = 2.0 / (x[n - 1] - t0) , c2 = sqrt(2.0 * grav);

  std::vector<Number > P1(N_ + 1), P2(N_ + 1), C(N_ + 1), X(N_ + 1), Y(N_ + 1);
  std::vector < std::vector<Number > > D(N_ + 1 , std::vector <Number > (N_ + 1));
  C = compute_c<Number, Index>(N_ + 1);
  D = formulate_differentiation_matrix<Number, Index>(C, T, N_ + 1);

  // form X
  for (Index i = 0; i <= N_; i++) {
    X[i] = x[i];
  }
  // form Y
  for (Index i = (N_ + 1); i <= ((2 * N_) + 1); i++) {
    Y[i - (N_ + 1)] = x[i];
  }

  P1 = multiply_D_X<Number, Index>(D, X, N_ + 1);
  P2 = multiply_D_X<Number, Index>(D, Y, N_ + 1);
  //std::cout << std::endl;
  /*for (Index i = 0; i <= N_; i++) {
    std::cout << "P1[" << i << "] : " << P1[i] << "\n";
  }  //std::cout << std::endl;
  */
  Index shift = ((2 * N_) + 2);
  //std::cout << "\nConstraints\n";
  for (Index i = 0; i <= N_; i++) {
    g[nth] = (c1 * P1[i]) - (c2 * sqrt(Y[i]) * cos(x[i + shift])) ;
    /*std::cout << "\n --\n";
    std::cout << "|\n";
    std::cout << "\tX[" << nth << "] = " << X[nth] << "\n";
    std::cout << "\tY[" << nth << "] = " << Y[nth] << "\n";
    std::cout << "\tP[" << nth << "] = " << P1[nth] << "\n";
    std::cout << "\tC[" << nth << "] = " << cos(x[i + shift]) << "\n";
    std::cout << "\tg[" << nth << "] = " << g[nth] << "\n";
    std::cout << "\t\t\t  |\n";
    std::cout << "\t\t\t--\n";*/
    //std::cout << "g[" << nth << "] : " << g[nth] << "\n";
    nth = nth + 1;

  }

  for (Index i = 0; i <= N_; i++) {
    g[nth] = (c1 * P2[i]) - (c2 * sqrt(Y[i]) * sin(x[i + shift])) ;
    //std::cout << "g[" << nth << "]" << g[nth] << "\n";
    //std::cout << "g[" << nth << "] : " << g[nth] << "\n";
    nth = nth + 1;

  }
  g[nth] = x[N_]; // X[0] = 0
  //std::cout << "g[" << nth << "] : " << g[nth] << "\n";
  nth = nth + 1;

  g[nth] = x[(2 * N_) + 1]; // Y[0]=0
  //std::cout << "g[" << nth << "] : " << g[nth] << "\n";
  nth = nth + 1;

  g[nth] = x[0] - 0.5; // X(N)-0.5 = 0
  //nth = nth + 1;

  //g[nth] = x[n+1]-0.5;
  //std::cout << "g[" << (nth) << "] : " << g[nth] << "\n";
  //assert(nth == (2 * (N_ + 1)) + 3);std::ofstream myfile;
  //myfile.open("output.txt");
  //P1.clear();
  //P2.clear();
  //X.clear();
  //Y.clear();
  //C.clear();
  //D.clear();
  return true;
}

// return the structure or values of the Jacobian
bool STONE_SLIDE_NLP::eval_jac_g
(
  Index n,          //

  const Number * x, //
  bool new_x,       //
  Index m,          //
  Index nele_jac,   //
  Index * iRow,     //
  Index * jCol,     //
  Number * values   //
)
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

void STONE_SLIDE_NLP::finalize_solution
(
  SolverReturn status,       //

  Index n,                   //
  const Number * x,          //
  const Number * z_L,        //
  const Number * z_U,        //
  Index m,                   //
  const Number * g,          //
  const Number * lambda,     //
  Number obj_value,          //
  const IpoptData * ip_data, //
  IpoptCalculatedQuantities * ip_cq
)
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




  // myfile << "Writing this to a file.\n";

  std::vector<Number > time(N_ + 1);
  //Number tf = 1.25;
  for (Index i = 0; i <= N_; i++)
  {
    time[i] = (x[n - 1] / 2.0) * (T[i] + 1.0);
    std::cout << "t[" << i << "]" << " : " << time[i] << std::endl;
  }

  std::ofstream myfile;
  myfile.open("X.txt");
  // for X
  std::cout << "X \n";
  for (Index i = 0; i <= N_; i++) {
    std::cout << x[i] << std::endl;
    myfile << time[i] << "," << x[i] << "\n";
  }
  myfile.close();
  //

  myfile.open("Y.txt");
  // for Y
  std::cout << "\nY\n";
  for (Index i = (N_ + 1); i <= ((2 * N_) + 1); i++) {
    std::cout << x[i] << std::endl;
    myfile << time[i - (N_ + 1)] << "," << x[i] << "\n";
  }
  myfile.close();
  myfile.open("theta.txt");
  // // for Theta
  std::cout << "\nTheta \n";
  for (Index i = ((2 * N_) + 2); i <= ((3 * N_) + 2); i++) {
    std::cout << x[i] << std::endl;
    myfile << time[i - ((2 * N_) + 2)] << "," << x[i] << "\n";
  }
  myfile.close();
  //for tow
  std::cout << "\nt_f \n";
  std::cout << x[n - 1] << std::endl;

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

  std::cout << std::cout << std::endl
            << "Final value of the constraints:" << std::endl;
  for (Index i = 0; i < m; i++) {
    std::cout << "g(" << i << ") = " << g[i] << std::endl;
  }

  // Analytical Solution
  /*std::cout << std::endl
            << std::endl
            << "Numerical Solution " << std::endl;
  */


}
