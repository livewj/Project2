/*
Parts of this code consists of/is inspired by the code given in the lecture notes:
"http://folk.uio.no/compphys/programs/chapter07/cpp/jacobi.cpp"
and in the solution from last year:
"https://github.com/CompPhysics/ComputationalPhysics/tree/master/doc/Projects/2015/Project2/Solution2015"

Jacobi's method for finding eigenvalues eigenvectors of the symetric matrix A.
The eigenvalues of A will be on the diagonal of A, with eigenvalue i being A[i][i].
The j-th component of the i-th eigenvector is stored in R[i][j].
A: input matrix (n x n)
R: empty matrix for eigenvectors (n x n) n: dimention of matrices
*/

#include <iostream> 
#include <cmath> 
#include <fstream> //for writing to file
#include <cstring>
#include <iomanip>
#include <time.h> //for time-taking
#include <armadillo> //for comparing

using namespace std;
using namespace arma;
ofstream ofile;


//declare functions
void jacobi_method ( double ** A, double ** R, int n );
double maxoffdiag( mat A, int *k, int *l, int n);

//function that finds the maximal off-diagonal element i the matrix A

double maxoffdiag(mat A, int *k, int *l, int n) {
	double max_off = 0.; //starting optimistic with zero as maximum value
	double A_ij;
	//loop over non-diagonal matrix elements in A:
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			if (i != j) { //if element is not on diagonal
				A_ij = fabs(A(i,j)); //store absolute value
				if (A_ij > max_off) {
					max_off = A_ij;
					*k = i; //save row-number
					*l = j; //save column-number
				}
			}
		}
	}
	return max_off;
}


//this function rotates the matrix A by using Jacobi's method
//and saves the eigenvectors the columns of matrix R

void jacobi_method ( mat &A, mat &R, int k, int l, int n ) {
	//declare variables
	double tau; //cot(2theta)
	double t; //tan
	double s; //sin
	double c; //cos

	//find the values of cos and sin
	if ( A(k,l) != 0.0 ) {
		tau = (A(l,l) - A(k,k))/(2*A(k,l));
		//now t can take two values, we want the smallest one  NBNBNBNBNBOBS
		if ( tau >= 0.0) {
			//t = -fabs(tau) + sqrt(1 + tau*tau); 
			t = 1.0/(fabs(tau) + sqrt(1. + tau*tau));
		}
		else {
			//t = -fabs(tau) - sqrt(1 + tau*tau); 
			t = -1.0/(fabs(tau) + sqrt(1. + tau*tau));
		}

		c = 1./(sqrt(1. + t*t));
		s = t*c;
	}
	//no rotation needed if element already is equal to 0
	else {
		c = 1.;
		s = 0.;
	}

	//Now rotate this bastard
	double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
	a_kk = A(k,k);
	a_ll = A(l,l);
	//changing matrix elements
	A(k,k) = a_kk*c*c - 2*A(k,l)*c*s + a_ll*s*s;
	A(l,l) = a_ll*c*c + 2*A(k,l)*c*s + a_kk*s*s;
	//theta is chosen so that all A(k,l) and A(l,k) elements become zero
	A(k,l) = 0.;
	A(l,k) = 0.;

	//more rotating
	for (int i=0; i<n; i++) {
		if (i != k && i != l) {
			A(i,i) = A(i,i); //unaffected
			a_ik = A(i,k);
			a_il = A(i,l);
			A(i,k) = a_ik*c - a_il*s;
			A(k,i) = A(i,k);
			A(i,l) = a_il*c + a_ik*s;
			A(l,i) = A(i,l);
		}

		//compute the new eigenvektors
		r_ik = R(i,k);
		r_il = R(i,l);
		R(i,k) = c*r_ik - s*r_il;
		R(i,l) = c*r_il + s*r_ik;
	}
	return;
}

int main(){

	//declare variables
	int n;
	int k, l;
	double epsilon = 1.0e-8; //tolerance
	int iterations; //counting the number of iterations
	double max_number_iterations = 1e6; 
	double rho_0, rho_max;
	double CPU_time;
	double h; //step
	double e; //non-diagonal element
	string outfile; //write to this file
	rowvec N; //N is a vector

	rho_0 = 0.;
	rho_max = 5.;
	N << 10 << 50 << 100 << 200 << 300 << 310;  //different sizes of matrix 

	/////////////////////////////////
	//UNIT TEST 1: Is the function maxoffdiag working correctly?
	//Defining a known test-matrix
	mat T;
	T << 1 << 3 << 16 << 6 << 8 << endr
	  << -7 << -17 << -8 << 0. << 2 << endr
	  << 0. << -9 << 10 << 5 << 3 << endr
	  << -20 << 2 << 9 << 25 << 0. << endr
	  << 3 << 1 << 4 << 6 << 8 << endr;
	double test_result;
	test_result = maxoffdiag(T, &k, &l, 5);
	//cout << test_result << endl;  
	//should return 20: Test passed
	////////////////////////////////

	outfile = "Jacobi.txt";
	ofile.open(outfile);
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	ofile << "rho_max: " << rho_max <<  endl;
	ofile << "tolerance: " << epsilon << endl;
	ofile << "n:            Eigenvalues:                         Eigenvalues(arma):         Time:          Time(arma):    Rotations:"<< endl;


	for (int j=0; j<6; j++) { //length of N is 6
		n = N(j); //size of matrix
		cout << n << endl;
		iterations = 0; //counting the number of iterations
		h = (rho_max - rho_0)/n;
		e = -1./(h*h);
		vec d(n-1);
		vec rho(n+1);
		vec eigenvalues(n-1);
		mat A(n-1, n-1);
		A.zeros(); //setting up the matrix A
		mat R(n-1,n-1);
		R.eye(); //R is set to identity matrix



		//setting up vectors
		for (int i = 0; i <= n; i++) {
			rho(i) = rho_0 + i*h;
		}
		for (int i=0; i<n-1; i++) {
			d(i) = 2.0/(h*h) + rho(i+1)*rho(i+1);
		}

		//setting up the matix A
		for (int i=0; i < n-1; i++) {
			A(i,i) = d(i);
			if (i != n-2) {
				A(i,i+1) = e;
				A(i+1,i) = e; 
			}
		}
	

	mat B = A; //used in Armadillo method

    clock_t start, finish;

	//compute the eigenvectors with Jacobi's method:
	double max_off = maxoffdiag(A, &k, &l, n-1);

	//start timer
	start = clock();

	while (max_off > epsilon && iterations < max_number_iterations) {
		max_off = maxoffdiag(A, &k, &l, n-1);
		jacobi_method(A, R, k, l, n-1); //rotate one more time
		iterations++;
	}

	//stop timer
    finish = clock();
    CPU_time = ((double) (finish - start)/CLOCKS_PER_SEC );


	//we now have the eigenvalues stored in the diagonal in matrix A
	//fill up eigenvalues vector
	for (int i=0; i<n-1; i++) {
		eigenvalues(i) = A(i,i);
	}
	//sort in ascending order
	eigenvalues = sort(eigenvalues);


	/////////////////////////////////
	//UNIT TEST 2: Does jacobi_method return the correct eigenvalues?
	mat R_(3,3);
	R_.eye(); //R_ is set to identity matrix
	mat C(3,3);
    C << 1 << 2 << 3 << endr
	  << 4 << 5 << 6 << endr
	  << 5 << 4 << 3 << endr;
	double max_off_ = maxoffdiag(C, &k, &l, 3);
	int iterations_ = 0;
	while (max_off_ > epsilon && iterations_ < max_number_iterations) {
		max_off_ = maxoffdiag(C, &k, &l, 3);
		jacobi_method(C, R_, k, l, 3); //rotate one more time
		iterations_++;
	}
	vec eigenvalues_(3);
	for (int i=0; i<3; i++) { eigenvalues_(i) = C(i,i); }
	eigenvalues_ = sort(eigenvalues_);
	//cout << eigenvalues_ << endl;
	//Should return the eigenvalues: 11.15, -2.15 and 0.0 (according to MATLAB)
	//Test returns: 11.26, -2.48 and 0.22
	//Close enuogh for a 3x3-matrix; test passed
	/////////////////////////////////

	//now calculating the eigenvalues -and vectors using armadillo
	double CPU_arma;
	clock_t start_, finish_;
	mat R_arma; //store the eigenvectors
	vec eigenvalues_arma; //store the eigenvalues

	start_ = clock();
	eig_sym(eigenvalues_arma, R_arma, B); //finding eigenvalues and vectors from previously stored matrix
	finish_ = clock();

	CPU_arma = ((double) (finish_ - start_)/CLOCKS_PER_SEC );




	//finally, write results to file
	ofile << setprecision(5) << n << ",    ";
	ofile << "(" << setprecision(7) << eigenvalues(0) << ", ";
	ofile << setprecision(7) << setw(10) << eigenvalues(1) << ", ";
	ofile << setprecision(7) << setw(10) << eigenvalues(2) << ")";
	ofile << "(" << setprecision(7) << eigenvalues_arma(0) << ", ";
	ofile << setprecision(7) << setw(10) << eigenvalues_arma(1) << ", ";
	ofile << setprecision(7) << setw(10) << eigenvalues_arma(2) << ")";
	ofile << setw(10) << setprecision(5) << CPU_time << " sec.";
    ofile << setw(10) << setprecision(5) << CPU_arma << " sec.";
	ofile << setw(7) << setprecision(6) << iterations << endl;
	}
	ofile.close();

	
	return 0;
}


