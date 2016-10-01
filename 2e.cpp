/*
Using Armadillo for finding eigenvalues eigenvectors of the symetric matrix A.
The eigenvalues of A will be on the diagonal of A, with eigenvalue i being A[i][i].
The j-th component of the i-th eigenvector is stored in R[i][j].
A: input matrix (n x n)
R: empty matrix for eigenvectors (n x n) 
n: dimention of matrices

This program  plots the eigenfunctions, with and without Coulomb-interaction, for different values of omega_R.

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


int main(){

	//declare variables
	int n = 400;
	int k, l;
	double epsilon = 1.0e-8; //tolerance
	int iterations; //counting the number of iterations
	double max_number_iterations = 1e6; 
	double rho_0, rho_max;
	double CPU_time;
	double h; //step
	double e; //non-diagonal element

	//rho_max varies as function of omega_R
	rho_0 = 0.;
	rho_max = 2.;
	double omega = 5.; //0.01, 0.5, 1.0, 5.0
	iterations = 0; //counting the number of iterations
	h = (rho_max - rho_0)/n;
	e = -1./(h*h);
	vec d(n-1);
	vec rho(n+1);
	vec eigenvalues(n-1);
	vec V(n-1); //POTENTIAL
	mat A(n-1, n-1);
	A.zeros(); //setting up the matrix A
	mat R(n-1,n-1);
	R.eye(); //R is set to identity matrix



	//setting up vectors
	for (int i = 0; i <= n; i++) {
		rho(i) = rho_0 + i*h;
	}
	for (int i = 0; i<n-1; i++) {
		V(i) = omega*omega*rho(i+1)*rho(i+1) + 1./rho(i+1); //Potetial with interaction
	}
	for (int i=0; i<n-1; i++) {
		d(i) = 2.0/(h*h) + V(i);
	}


	//setting up the matix A
	for (int i=0; i < n-1; i++) {
		A(i,i) = d(i);
		if (i != n-2) {
			A(i,i+1) = e;
			A(i+1,i) = e; 
		}
	}

	//now calculating the eigenvalues -and vectors using armadillo
	eig_sym(eigenvalues, R, A);


	//Extracting the eigenvectors:
	vec eigenvector1 = R.col(0);
	vec eigenvector2 = R.col(1);
	vec eigenvector3 = R.col(2);


	//Squaring the eigenvectors:
	for (int i=0; i<n-1; i++) {
		eigenvector1(i) *= eigenvector1(i);
		eigenvector2(i) *= eigenvector2(i);
		eigenvector3(i) *= eigenvector3(i);
	}

	double norm1 = norm(eigenvector1);
	double norm2 = norm(eigenvector2);
	double norm3 = norm(eigenvector3);

	//norming the eigenvectors:
	for (int i=0; i<n-1; i++) {
		eigenvector1(i) /= norm1; 
		eigenvector2(i) /= norm2;
		eigenvector3(i) /= norm3;  
	}
 

	//write reults to file for a specified omega_R:
	string outfile; //write to this file
	outfile = "wavefunc_omega_5_ON_not.txt";
	ofile.open(outfile);
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	ofile << "omega:    " << omega << endl;
	ofile << "   rho:            u1:                u2:                 u3:" << endl;

	for (int i=0; i<n-1; i++) {
		ofile << setw(10) << setprecision(8) << rho(i+1);
		ofile << setw(18) << setprecision(10) << eigenvector1(i);
		ofile << setw(18) << setprecision(10) << eigenvector2(i);
		ofile << setw(18) << setprecision(10) << eigenvector3(i) << endl;
	}

	ofile.close();

	return 0;

}


