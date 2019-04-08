// to test Differentiation Matrix over Chebyshev polynomial
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "Differentiation_Matrix.cpp"
#define PI acos(-1)
using namespace std;
int main() {
	ios_base::sync_with_stdio(false), cin.tie(0), cout.tie(0);
	freopen("test_input.txt", "r", stdin);
	freopen("test_output.txt", "w", stdout);
	int n;
	cin >> n;
	std::vector<double > X(n + 1), P(n + 1), C(n + 1), T(n + 1);
	std::vector < std::vector<double > > D(n + 1 , std::vector <double > (n + 1));

	T = define_time_stamps<double, int>(n + 1);

	//cout << endl;
	// for (int i = 0 ; i <=n ; i++) {
	//cout <<double(i)<< ","<<T[i] << "\n";
	//} //cout << endl;

	for (int i = 0 ; i <= n ; i++) {
		X[i] = sin(T[i]);
		//X[i] = 1;
		//cout << T[i] << "," << X[i] << "\n";
	}

	C = compute_c<double, int>(n + 1);
	// cout << endl;
	// for (int i = 0 ; i < n ; i++) {
	// 	cout << i << "," << C[i] << "\n";
	// } cout << endl;

	D = formulate_differentiation_matrix<double, int>(C, T, n + 1);
	std::cout.precision(4);
	/*
	for (int k = 0 ; k <=n; k++)
	{
		double sum=0;
		for (int j = 0 ; j <=n ; j++)
		{
			sum+=D[j][k];
			//cout << D[j][k] << " ";
		} cout <<sum<< endl<<endl;
	}
	*/
	P = multiply_D_X<double, int>(D, X, n + 1);
	for (int i = 0 ; i <= n ; i++) {
		cout << T[i] << "," << P[i] << "\n";
	}
	return 0;
}


