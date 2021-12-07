#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <complex>

using namespace std;


typedef vector<vector<complex<double>>> matrix;
typedef vector<complex<double>> ComplVector;

double T = 2;
double firstFrequency = 100;
double secondFrequency = 10;
double samplingFrequency = 0.001;
double phase = 0;
const double pi = 3.141592653589793238;
const double eps = 1e-6;

complex<double> imagePart(0, 1);


double sinFunction(double time, double signalFrequency) {
	return sin((2 * pi * signalFrequency * time) + phase);
}

vector<complex<double>> sampling(double time1, double time2, double samplingFreq, double (*sinSignal)(double, double), double signalFrequency) {
	vector<complex<double>> samplingSignal;
	for (double t = time1; t <= time2; t += samplingFreq) {
		samplingSignal.emplace_back(sinSignal(t, signalFrequency), 0);
	}
	return samplingSignal;
}

// vector<complex<double>> addCoef(vector<complex<double>> signal, complex<double> coef) {
// 	vector<complex<double>> tmp;
// 	for(complex<double> i : signal){
// 		tmp.emplace_back(i * coef);
// 	}
// 	return tmp;
// }

matrix createF(int N) {
	matrix F(N, vector<complex<double>>(N));
	for (int k = 0; k < N; k++) {
		for (int n = 0; n < N; n++) {
			F[k][n] = exp(imagePart * complex<double>(((2 * pi * n)/N) * k, 0));
		}
	}
	return F;
}

matrix createHermitConjugateF(int N) {
	matrix FH(N, vector<complex<double>>(N));
	for (int k = 0; k < N; k++) {
		for (int n = 0; n < N; n++) {
			FH[k][n] = conj(exp(imagePart * complex<double>(((2 * pi * n)/N) * k, 0)));
		}
	}
	return FH;
}

vector<complex<double>> product(vector<complex<double>> samplingSignal, matrix F){
	int N = samplingSignal.size();
	vector<complex<double>> restoredSamplingSignal(N);
	for(int k = 0; k < N; k++) {
		complex<double> sum = 0;
		for(int n = 0; n < N; n++) {
			sum += (F[n][k] * samplingSignal[n]);
		}
		restoredSamplingSignal[k] = sum;
	}
	return restoredSamplingSignal;
}

bool isZero(double number){
	return number < eps;
}

ComplVector spectrOffset(ComplVector spectr) {
	ComplVector offsetSpectr;
	for(complex<double> i : spectr) {
		offsetSpectr.emplace_back(i * exp((imagePart * complex<double>(2*pi * phase / spectr.size(), 0))));
	}
	return offsetSpectr;
}

int main() {

	int a = 2;
	int b = 3;

	cout << "input frequencies & T | (f1 f2 T): ";
	cin >> firstFrequency >> secondFrequency >> T;
	cout << "f1: " << firstFrequency << " f2: " << secondFrequency << " T: " << T << '\n';
	cin.get();
	// cout << "input frequencies & T & sampFreq | (f1 f2 T sF): ";
	// cin >> firstFrequency >> secondFrequency >> T >> samplingFrequency;
	// cout << "f1: " << firstFrequency << " f2: " << secondFrequency << " T: " << T << " sF: " << samplingFrequency << '\n';
	ofstream fout;
	fout << fixed << setprecision(22);


	//--signals--//
	ComplVector firstSignal = sampling(0, T, samplingFrequency, sinFunction, firstFrequency);
	ComplVector secondSignal = sampling(0, T, samplingFrequency, sinFunction, secondFrequency);

	
	//--spectr signals--//
	ComplVector spectrFirstSignal = product(firstSignal, createHermitConjugateF(firstSignal.size()));
	ComplVector spectrSecondSignal = product(secondSignal, createHermitConjugateF(secondSignal.size()));


	//--linearity signal--//
	ComplVector linearity = firstSignal;
	for (int i = 0; i < firstSignal.size(); i++) {
		linearity[i] += secondSignal[i];
	}
	
	// ComplVector spectrLinearity = product(linearity, createHermitConjugateF(linearity.size()));
	ComplVector spectrLinearity = spectrFirstSignal;
	for (int i = 0; i < spectrFirstSignal.size(); i++) {
		spectrLinearity[i] += spectrSecondSignal[i];
	}

	//--restored signals--//
	ComplVector restoredFirstSignal = product(spectrFirstSignal, createF(spectrFirstSignal.size()));
	for(int i = 0; i < restoredFirstSignal.size(); i++) {
		restoredFirstSignal[i] = restoredFirstSignal[i] / complex<double>(restoredFirstSignal.size(), 0);
	}
	ComplVector restoredLinearity = product(spectrLinearity, createF(spectrLinearity.size()));
	for(int i = 0; i < restoredLinearity.size(); i++) {
		restoredLinearity[i] = restoredLinearity[i] / complex<double>(restoredLinearity.size(), 0);
	}


	//--offset signal & spectr--//
	::phase = 0.3;
	ComplVector offsetFirstSignal = sampling(0, T, samplingFrequency, sinFunction, firstFrequency);
	ComplVector spectrFirstoffset = spectrOffset(spectrFirstSignal);
	ComplVector restoredFirstOffsetSignal = product(spectrFirstoffset, createF(spectrFirstoffset.size()));
	for(int i = 0; i < restoredFirstOffsetSignal.size(); i++) {
		restoredFirstOffsetSignal[i] = restoredFirstOffsetSignal[i] / complex<double>(restoredFirstOffsetSignal.size(), 0);
	}



	// ComplVector linearity = addCoef(firstSignal, complex<double>(a, 0));
	// ComplVector secondSignalWithCoef = addCoef(secondSignal, complex<double>(b, 0));
	// for (int i = 0; i < firstSignal.size(); i++) {
	// 	linearity[i] += secondSignalWithCoef[i];
	// }
	// ComplVector spectrLinearity = addCoef(spectrFirstSignal, complex<double>(a, 0));
	// ComplVector secondSpectrWithCoef = addCoef(spectrSecondSignal, complex<double>(b, 0));
	// for (int i = 0; i < spectrFirstSignal.size(); i++) {
	// 	spectrLinearity[i] += secondSpectrWithCoef[i];
	// }
	// ComplVector spectrLinearity = product(linearity, createHermitConjugateF(linearity.size()));




	fout.open("restoredSignal.txt");
	for(complex<double> i : restoredFirstSignal)
		fout << i.real() << '\n';
	fout.close();
	fout.open("signal1.txt");
	for(complex<double> i : firstSignal)
		fout << i.real() << '\n';
	fout.close();
	system("python3 graphDraw.py signal1.txt restoredSignal.txt &");

	fout.open("spec1.txt");
	for(complex<double> i : spectrFirstSignal)
		fout << abs(i.imag()) << '\n';
	fout.close();
	fout.open("spec2.txt");
	for(complex<double> i : spectrSecondSignal)
		fout << abs(i.imag()) << '\n';
	fout.close();
	system("python3 graphDraw.py spec1.txt spec2.txt &");

	cin.get();

	// fout.open("offsetSpec1.txt");
	// for(complex<double> i : spectrFirstoffset)
	// 	fout << abs(i.imag()) << '\n';
	// fout.close();
	// system("python3 graphDraw.py spec1.txt offsetSpec1.txt &");

	fout.open("signal1.txt");
	for(complex<double> i : firstSignal)
		fout << i.real() << '\n';
	fout.close();
	fout.open("restoredOffsetSign1.txt");
	for(complex<double> i : restoredFirstOffsetSignal)
		fout << i.real() << '\n';
	fout.close();
	system("python3 graphDraw.py signal1.txt restoredOffsetSign1.txt &");

	cin.get();

	fout.open("linearity.txt");
	for(complex<double> i : linearity)
		fout << i.real() << '\n';
	fout.close();
	fout.open("restoredLinearity.txt");
	for(complex<double> i : restoredLinearity)
		fout << i.real() << '\n';
	fout.close();
	system("python3 graphDraw.py linearity.txt  restoredLinearity.txt&");

	fout.open("spectrLinearity.txt");
	for(complex<double> i : spectrLinearity)
		fout << abs(i.imag()) << '\n';
	fout.close();
	system("python3 graphDraw.py spectrLinearity.txt &");


	
	// system("python3 graphDraw.py offsetSign1.txt &");



	for(int i = 0; i < firstSignal.size(); i++){
		if(!isZero(firstSignal[i].real() - restoredFirstSignal[i].real())) {
			cout << "\t[X] The first signal is NOT similar to the restored signal.\n";
			break;
		}
		if(i == firstSignal.size()-1){
			cout << "\t[+] The first signal is similar to the restored signal\n--------------------------------------------------\n";
		}
	}
	

	double sumSpectr = 0;
	double sumRestored = 0;
	for(int i = 0; i < spectrFirstSignal.size(); i++){
		sumSpectr += spectrFirstSignal[i].imag() * spectrFirstSignal[i].imag();
	} sumSpectr /= spectrFirstSignal.size();

	for(int i = 0; i < restoredFirstSignal.size(); i++) {
		sumRestored += restoredFirstSignal[i].real() * restoredFirstSignal[i].real();
	}

	cout << "Parseval Equality: " << sumSpectr << " -- " << sumRestored << '\n';

	return 0;
}