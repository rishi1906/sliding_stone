// #ifndef __Differentiation_Matrix__
// #define __Differentitaion_Matrix__

// //#include "IpTNLP.hpp"
// #include <iostream>
// #include <vector>
// //using namespace Ipopt;
// using namespace std;
std::vector<double > compute_c
(
    int N 	//
);

std::vector<double > define_time_stamps
(
    int N 	//
);

std::vector<double > multiply_D_X
(
    std::vector<std::vector<double > > D,	// matrix D
    std::vector<double > 			   X,	// vector X
    int 							   N	// length of vector X
);

std::vector<std::vector<double > > formulate_differentiation_matrix
(
    std::vector<double > c,	//
    std::vector<double > t,	//
    int 				 N	//
);

//#endif