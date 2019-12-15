#include <iostream>
#include <vector>
#include<string>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <time.h>
using namespace std;

double f(double x) {
	return sin(x) + 0.5;
}

std::vector<double> fshum(std::vector<double> X) {
	std::vector<double>F(101);
	for (int i = 0; i < 101; i++) {
		double omega = -1 + rand() % 2;
		omega += (double)rand() / RAND_MAX;
		
		F[i] = f(X[i]) + omega/4.;
		
	}

	return F;
}
//--------------------------------методы для r=3---------------------------------------------//
std::vector<double> Alpha() {
	std::vector<double>alpha(3);
	alpha[1]= (double)rand() / RAND_MAX;
	alpha[0] = (1 - alpha[1]) / 2;
	alpha[2] = alpha[0];
	return alpha;
}

std::vector<double> fsredn(std::vector<double> fshu, std::vector<double> alpha) {
	std::vector<double> fsr(101);
	double buffer = 1;
	
	for (int i = 0; i < 101; i++) {
		if (i == 0) fsr[i] = pow(fshu[i], alpha[1]) * pow(fshu[i + 1], alpha[2]);
			
		
		if (i == 100) fsr[i]= pow(fshu[i - 1], alpha[0]) * pow(fshu[i], alpha[1]);
		if (i!=100&&i!=0) fsr[i]=pow(fshu[i-1], alpha[0]) * pow(fshu[i], alpha[1]) * pow(fshu[i+1], alpha[2]);
		
	}
	return fsr;
}
//-------------------------------------------------------------------------------------------//
double w(std::vector<double> fs) {
	double W = 0;
	for (int i = 1; i < 101; i++) W += pow((fs[i] - fs[i - 1]), 2);
	return sqrt(W);
}

double delta(std::vector<double> fshu, std::vector<double> fs) {
	double Delta = 0;
	for (int i = 0; i < 101; i++) {
		Delta += pow((fs[i] - fshu[i]), 2);
	}
	Delta = Delta / 100;
	
	return sqrt(Delta);
}

double J(double lamda, double W, double Delta) {
	return lamda * W + (1 - lamda) * Delta;
}

double distance(double W, double Delta) {
	return sqrt(W * W + Delta * Delta);
}

int veroyatnos() {
	int N = 0;
	N= log(0.05) / log(1 - (0.01 / 3.1415));
	return N;
}

void method_for_3() {	
	std::vector<double>X(101);
	std::vector<double> lamda(11);
	for (int i = 0; i < 11; i++){
		lamda[i] = i / 10.;
	}
	double Wmax = 0;
	double Deltamax = 0;
	double mint = 9999999;
	std::vector<double>alphamax(3);
	for (int i = 0; i < 101; i++) {
		X[i] = (3.1415 * i) / 100;
	}

	
	double minf = 0;
	int flag = 0;
	int N = veroyatnos();
	double Ji = 0;
	double Wi = 0;
	double Di = 0;
	double h = 0;
	for (int i = 0; i < 11; i++) {
		for (int j = 0; j < N; j++) {
			std::vector<double> fshu = fshum(X);
			std::vector<double> alpha = Alpha();
			std::vector<double> fs = fsredn(fshu, alpha);
			
			double W = w(fs);
			double Delta = delta(fshu, fs);
			if (flag == 0) {
				flag = 1;
				minf = distance(W, Delta);
				
				Wmax = W;
				Deltamax = Delta;
				alphamax[0] = alpha[0];
				alphamax[1] = alpha[1];
				alphamax[2] = alpha[2];
			}
			else {
				if (minf > distance(W, Delta)) {
					minf = distance(W, Delta);
					
					Wmax = W;
					Deltamax = Delta;
					alphamax[0] = alpha[0];
					alphamax[1] = alpha[1];
					alphamax[2] = alpha[2];
					
				}
				else continue;

			}
		}
		std::cout << fixed << setprecision(1) << lamda[i] << " | " << fixed<<setprecision(6)<<distance(Wmax, Deltamax) << "   |[" << fixed << setprecision(6) << alphamax[0] << "," <<
			fixed << setprecision(6) << alphamax[1] << "," << fixed << setprecision(6) << alphamax[2] << "]|   " << Wmax << " | " << Deltamax <<" | "<< std::endl;
		if (mint > minf) {
			mint = minf;
			h = lamda[i];
			Ji = J(h, Wmax, Deltamax);
			Wi = Wmax;
			Di = Deltamax;
		}

		flag = 0;
		minf = 0;
		Wmax = 0;
		Deltamax = 0;
		alphamax[0] = 0;
		alphamax[1] = 0;
		alphamax[2] = 0;
	}
	std::cout << std::endl;
	std::cout << "----------------------------------------------" << std::endl;
	std::cout <<  "| h*   |      J     |      w     |     d     |" << std::endl;
	std::cout << "| " << fixed << setprecision(1)<< h << "  |  " << fixed << setprecision(6) << Ji << "  |  " << Wi << "  |  " << Di << " |" << std::endl;
	std::cout << "----------------------------------------------" << std::endl;
	
}

//--------------------------------методы для r=5---------------------------------------------//
std::vector<double> Alpha2() {
	std::vector<double>alpha(5);
	alpha[2] = (double)rand() / RAND_MAX;
	alpha[1] = 9999;
	while (alpha[1] > (1 - alpha[2])) {
	alpha[1]= (double)rand() / RAND_MAX;
	}
	alpha[1] = alpha[1] / 2.;
	alpha[3] = alpha[1];


	alpha[0] = 0.5 * (1 - (alpha[1] + alpha[2] + alpha[3]));
	alpha[4] = alpha[0];

	return alpha;
}

std::vector<double> fsredn2(std::vector<double> fshu, std::vector<double> alpha) {
	std::vector<double> fsr(101);
	double buffer = 1;

	for (int i = 0; i < 101; i++) {
		if (i == 0) fsr[i] = pow(fshu[i], alpha[2]) * pow(fshu[i + 1], alpha[3]) * pow(fshu[i + 2], alpha[4]);
		if (i == 1) fsr[i] = pow(fshu[i-1], alpha[1]) * pow(fshu[i], alpha[2]) * pow(fshu[i + 1], alpha[3]) * pow(fshu[i + 2], alpha[4]);
		if (i==99) fsr[i] = pow(fshu[i - 2], alpha[0]) * pow(fshu[i-1], alpha[1]) * pow(fshu[i], alpha[2]) * pow(fshu[i + 1], alpha[3]);
		if (i == 100) fsr[i] = pow(fshu[i - 2], alpha[0]) * pow(fshu[i - 1], alpha[1]) * pow(fshu[i], alpha[2]);
		if (i != 100 && i != 0 && i!=1 && i!=99) fsr[i] = pow(fshu[i - 2], alpha[0]) * pow(fshu[i-1], alpha[1]) * pow(fshu[i], alpha[2])
			* pow(fshu[i + 1], alpha[3]) * pow(fshu[i + 2], alpha[4]);

	}
	return fsr;
}
//-------------------------------------------------------------------------------------------//

void method_for_5() {
	std::vector<double>X(101);
	std::vector<double> lamda(11);
	for (int i = 0; i < 11; i++) {
		lamda[i] = i / 10.;
	}
	double Wmax = 0;
	double Deltamax = 0;
	double mint = 9999999;
	std::vector<double>alphamax(5);
	for (int i = 0; i < 101; i++) {
		X[i] = (3.1415 * i) / 100;
	}


	double minf = 0;
	int flag = 0;
	int N = veroyatnos();
	double Ji = 0;
	double Wi = 0;
	double Di = 0;
	double h = 0;
	for (int i = 0; i < 11; i++) {
		for (int j = 0; j < N; j++) {
			std::vector<double> fshu = fshum(X);
			std::vector<double> alpha = Alpha2();
			std::vector<double> fs = fsredn2(fshu, alpha);

			double W = w(fs);
			double Delta = delta(fshu, fs);
			if (flag == 0) {
				flag = 1;
				minf = distance(W, Delta);

				Wmax = W;
				Deltamax = Delta;
				alphamax[0] = alpha[0];
				alphamax[1] = alpha[1];
				alphamax[2] = alpha[2];
				alphamax[3] = alpha[3];
				alphamax[4] = alpha[4];
			}
			else {
				if (minf > distance(W, Delta)) {
					minf = distance(W, Delta);

					Wmax = W;
					Deltamax = Delta;
					alphamax[0] = alpha[0];
					alphamax[1] = alpha[1];
					alphamax[2] = alpha[2];
					alphamax[3] = alpha[3];
					alphamax[4] = alpha[4];

				}
				else continue;

			}
		}
		std::cout << fixed << setprecision(1) << lamda[i] << " | " << fixed << setprecision(6) << distance(Wmax, Deltamax) << "   |[" << fixed << setprecision(6) << alphamax[0] << "," <<
			fixed << setprecision(6) << alphamax[1] << "," << fixed << setprecision(6) << alphamax[2] << "," << fixed << setprecision(6) << alphamax[3] <<
			"," << fixed << setprecision(6) << alphamax[4] << "]|   " << Wmax << " | " << Deltamax << " | " << std::endl;
		if (mint > minf) {
			mint = minf;
			h = lamda[i];
			Ji = J(h, Wmax, Deltamax);
			Wi = Wmax;
			Di = Deltamax;
		}

		flag = 0;
		minf = 0;
		Wmax = 0;
		Deltamax = 0;
		alphamax[0] = 0;
		alphamax[1] = 0;
		alphamax[2] = 0;
		alphamax[3] = 0;
		alphamax[4] = 0;
	}
	std::cout << std::endl ;
	std::cout << "----------------------------------------------" << std::endl;
	std::cout  << "| h*   |      J     |      w     |     d     |" << std::endl;
	std::cout << "| " << fixed << setprecision(1) << h << "  |  " << fixed << setprecision(6) << Ji << "  |  " << Wi << "  |  " << Di << " |" << std::endl;
	std::cout << "----------------------------------------------" << std::endl;

}



int main() {
	srand(time(NULL));
	std::cout << "                               r=3                                     " << std::endl;
	std::cout << " h  |  distance  |            alpha           |     W      |   Delta  |" << std::endl;
	method_for_3();
	std::cout << std::endl << std::endl;
	std::cout << "                                         r=5                                       " << std::endl;
	std::cout << " h  |  distance  |                      alpha                   |     W      |   Delta  |" << std::endl;
	method_for_5();

	return 0;
}
