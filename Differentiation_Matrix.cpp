//#include <iostream>
#include "Differentiation_Matrix.h"
#include <cmath>
#define PI acos(-1)
//decimal DEG_to_RAD(decimal d) { return d * PI / 180.0; }
//decimal RAD_to_DEG(decimal r) { return r * 180.0 / PI; }

template<class decimal, class integer>
std::vector<decimal > compute_c
(
  integer N
)
{
	std::vector<decimal > c(N);
	c[0] = c[N - 1] = 2;
	for (integer i = 1 ; i < N - 1 ; i++) {
		c[i] = 1;
	}
	return c;
}

template<class decimal, class integer>
std::vector<decimal > define_time_stamps
(
  integer N
)
{
	std::vector<decimal > t(N);
	for (integer k = 0 ; k < N ; k++) {
		t[k] = (decimal)cos(((PI * k) / (N - 1) ));
	
	}
	return t;
}

template<class decimal, class integer>
std::vector<decimal > multiply_D_X
(
  std::vector<std::vector<decimal > > 	D,	// matrix D
  std::vector<decimal > 								X,	// vector X
  integer 																N 	// length of vector X
)
{
	std::vector<decimal > prod(N);

	for ( integer i = 0; i < N; i++ )
	{
		prod[i] = 0.0;                        // <===== Needs initialising
		for ( integer j = 0; j < N; j++ )
		{
			prod[i] += D[j][i] * X[j];       // <===== Add terms to sum for ith element
		}
	}
	return prod;
}

template<class decimal, class integer>
std::vector<std::vector<decimal > > formulate_differentiation_matrix
(
  std::vector<decimal > c, //
  std::vector<decimal > t,
  integer 				 N //
)
{
	std::vector<std::vector<decimal > > D (N , std::vector <decimal > (N));
	for (integer k = 0 ; k < N; k++) {
		for (integer j = 0 ; j < N ; j++)
		{
			if (j == k)
			{
				if (j == 0)
				{
					D[j][k] = decimal((2 * (N - 1) * (N - 1)) + 1) / 6;
				} else if (j == (N - 1) )
				{
					D[j][k] = -1.0 * (decimal((2 * (N - 1) * (N - 1)) + 1) / 6);
				} else
				{
					D[j][k] = (-1.0 * t[k]) / (2.0 * (1 - (t[k] * t[k])));
				}
			} else if (j != k)
			{
				D[j][k] = (c[k] / c[j]) * ( (pow(-1, (j + k))) / (t[k] - t[j]) );
			}
			//D[j][k] = -D[j][k];
		}
	}
	return D;
}
