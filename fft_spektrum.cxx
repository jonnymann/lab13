
Skip to content
This repository

    Pull requests
    Issues
    Gist

    @jonnymann

1
0

    15

AbdulHamada/lab13 forked from TP1-HHU/lab13
Code
Pull requests 0
Wiki
Pulse
Graphs
lab13/fft_spektrum.cxx
1d72384 6 hours ago
Abdul Hamada FFT
1 contributor
90 lines (71 sloc) 2.16 KB
#include "fftw3.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
//-------------------------------
using namespace std;
//-------------------------------
void writeData(const fftw_complex* const f, const int N, const double L,const char* const fname);
void reading (double* const S, const int N, double& L, const string fname );

//-------------------------------

int main(int argc, char** argv){

	if(argc != 3){
		cout << "Usage: " << argv[0] << " input_file \t output_file" << endl;
		exit(1);
	}

	char *in_file  = argv[1];
	char *out_file = argv[2];
	
	string fname = in_file;

	const int N = 16384;
	double L;

	// Allocate memory
	fftw_complex* f = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1));
	double* inR  = (double*) malloc(sizeof(double)*N);

  // Create plan
	fftw_plan FW  = fftw_plan_dft_r2c_1d(N, inR, f, FFTW_ESTIMATE);

	// Read input data
	reading(inR,N,L,fname);
	cout << inR[10] << '\t' << L << endl;
	// Call function which reads the data from
	// the input file into the array inR


  // Calculate FFT
  fftw_execute(FW);

  // write output file
  writeData(f, N,  L, out_file);

  // Clean up
	fftw_destroy_plan(FW);
  fftw_free(f);
  free(inR);

	return 0;
}
//------------------------------
void reading (double* const S, const int N, double& L, const string fname ){
  ifstream in(fname.c_str());
  double temp;
  for (int i=0; i<N; i++){
    in >> temp;
    in >> S[i];
  }
  L=temp; // um die ltzte Wert von x in L zuspeichern.
}
//-------------------------------
void writeData(const fftw_complex* const f, const int N, const double L,const char* const fname){
	ofstream out(fname);
	const double dk = 2*M_PI/L;
	double pk;


	for(int i=0; i<=N/2; i++){
		pk = sqrt(f[i][0]*f[i][0] + f[i][1]*f[i][1])/N;
		out << i*dk << "\t" << pk << "\t" << f[i][0] << "\t" << f[i][1] << endl;
	}

	out.close();
}
//-------------------------------
/* lab13> g++ fft_spektrum.cxx -o fft -I/.netmount/opt/fftw/include -L/.netmount/opt/fftw/lib -lfftw3
 ./fft i.txt if.txt
 gnuplot> p "i.txt" w l, "r.txt" w l
gnuplot> p "if.txt" w l, "rf.txt" w l
L druckenauf dem graphick um von normal zu logaritm zu wichseln
 */

    Status API Training Shop Blog About Pricing 

    © 2016 GitHub, Inc. Terms Privacy Security Contact Help 

