#include <iostream>
#include <fstream> 
#include <cmath>
#include <cstring>
#include <string>

double *CubicInterpolation(double *src_x, double *src_y, double *new_x);
double *Jordan(double **a, double **b, int N);
using namespace std;

int k = 9;
int N = k;
int K = 10;

int main()
{
	double *y, *x, *xi, *yi;
	y = new double[9];
	yi = new double[9];
	x = new double[9];
	xi = new double[9];
	y[0] = 1;
	y[1] = 1.007568;
	y[2] = 1.031121;
	y[3] = 1.073456;
	y[4] = 1.140228;
	y[5] = 1.242129;
	y[6] = 1.400176;
	y[7] = 1.6603;
	y[8] = 2.143460;
	x[0] = 0;
	cout << x[0] <<" "<<y[0]<< endl;
	for (int i = 1; i < 9; i++) {
		x[i] = x[i - 1] + 0.15;
		cout << x[i] <<" "<<y[i]<< endl;
	} 
	cout << "______________________" << endl;
	cout << endl;
	

	//Кубический сплайн

	cout << "Cebicheski spline:" << endl;
	int a = 0;
	int b = k - 1;
	double h;
	h = 0.15;
	int iter = 0;
	double *x1 = new double[(N - 1) * K + N];
	for (int i = 0; i < ((N - 1) * K + N); i++) {
		x1[i] = a + (i * h) / (K + 1);
		iter = i;
	}
	cout << "Iteration:" << iter  << endl;
	
	double *y1 = CubicInterpolation(x, y, x1);
	for (int i = 1; i < (N - 1) * K + N; i++) {
		if (i == k * (K + N)) {
			k++;
		}
		else {
			cout << x1[i] << " " << y1[i] << endl;
		}
	}
	cout << "______________________" << endl;
	cout << endl;
	

	//Метод трапеции 
	std::ofstream fout;
	fout.open("ans1.dat");
	cout << "Metod trapeci:" << endl;
	//fout << "Metod trapeci:" << endl;
		long double shag, z=0;
		shag = x1[2] - x1[1];
		cout << "Shag:" << shag<<endl;
		z = (y[0] + y[8]) / 2;
	    //cout << z << endl;
		//cout << "y0=" << y1[0];
		for (int i = 1; i < iter; i++) {
			z += y1[i];
		}
		z = z * shag;
		cout << "Znachenie integrala methodom trapeci:" << z << endl;
		//fout << "Znachenie integrala methodom trapeci:" << z << endl<<endl;
		cout << "______________________" << endl<<endl;

		//Метод Симпсона
		std::ofstream fout1;
		fout1.open("ans2.dat");
		cout << "Metod Simpsona:" << endl;
		//fout1 << "Metod Simpsona:" << endl;
		cout << "Shag:" << shag << endl;
		long double zS = 0;
		zS = y[0]+y[8];
		for (int i = 1; i < iter; i++) {
			if(i%2==0){
				zS += 2 * y1[i];
			}
			else {
				zS += 4 * y1[i];
			}
		}
		zS = zS * (shag / 3);
		cout << "Znachenie integrala methodom Simpsona:" << zS << endl;
		//fout1 << "Znachenie integrala methodom Simpsona:" << zS << endl<<endl<<endl;
		cout << "______________________" << endl << endl;

		//Метод прямоуглоьника
		std::ofstream fout2;
		fout2.open("ans3.dat");
		cout << "Metod priamoygolnikov:" << endl;
		//fout2 << "Metod priamoygolnikov:" << endl;
		cout << "Shag:" << shag << endl;
		double zT = 0;
		for (int i = 2; i <= iter; i=i+2) {
			zT += (x1[i] - x1[i - 2])*y1[i - 1];

		}
		cout << "Znachenie integrala methodom priamoygolnikov:" << zT << endl;
		//fout2 << "Znachenie integrala methodom priamoygolnikov:" << zS << endl<<endl<<endl;
		cout << "______________________" << endl << endl;

		//Точность вычисления метода трапеции

		cout << "Tochnost vichislenia metoda trapeci:" << endl;
		int iter1 = 0;
		long double pogresh;
		for(int j =0; j<1 ; j){
			K = 2 * K;
			cout << "K=" << K<<endl;
			iter1++;
			int iter = 0;
			double *x1 = new double[(N - 1) * K + N];
			for (int i = 0; i < ((N - 1) * K + N); i++) {
				x1[i] = a + (i * h) / (K + 1);
				iter = i;
			}
			cout << "Iteration:" << iter << endl;

			double *y1 = CubicInterpolation(x, y, x1);
			for (int i = 1; i < (N - 1) * K + N; i++) {
				if (i == k * (K + N)) {
					k++;
				}
				else {
				//	cout << x1[i] << " " << y1[i] << endl;
				}
			}
				long double z1 = 0;
				//cout << "Metod trapeci:" << endl;
				double shag;
				shag = x1[2] - x1[1];
				cout << "Shag:" << shag << endl;
				z1 = (y[0] + y[8]) / 2;
				//cout << z << endl;
				//cout << "y0=" << y1[0];
				for (int i = 1; i < iter; i++) {
					z1 += y1[i];
				}
				z1 = z1 * shag;
				cout << "Posledyushie znachenia integrala methodom trapeci:" << z1 << endl;
				//fout << "Posledyushie znachenia integrala methodom trapeci:" << z1 << endl;
				pogresh = 0;
				pogresh = z - z1;
				cout << "Tochnost visheslenia:" << pogresh<< endl;
				//fout << "Tochnost visheslenia:" << pogresh << endl<<endl;
				j = 1;
				if (pogresh > 0.000001) {
					j = 0;
				}
				z = z1;
				//cout << "Iteration tochnosti:" << iter1 << endl;
				cout <<endl;
		}

		cout << "Konechnoe znachenie integrala methodom trapeci:" << z << endl;
		cout << "Konechnaia tochnost visheslenia (method trapeci):" << pogresh << endl;
		fout << z <<" " << pogresh << endl;
		//fout << "Konechnaia tochnost visheslenia (method trapeci):" << pogresh << endl;
		cout << "Iteration tochnosti (method trapeci):" << iter1 << endl;
		cout << "______________________" << endl << endl;

		//Точность вычисления метода Симпсона

		cout << "Tochnost vichislenia metoda Simpsona:" << endl;
		int iter2 = 0;
		K = 10;
		long double pogresh1;
		for (int j = 0; j < 1; j) {
			K = 2 * K;
			cout << "K=" << K << endl;
			iter2++;
			int iter = 0;
			double *x1 = new double[(N - 1) * K + N];
			for (int i = 0; i < ((N - 1) * K + N); i++) {
				x1[i] = a + (i * h) / (K + 1);
				iter = i;
			}
			cout << "Iteration:" << iter << endl;

			double *y1 = CubicInterpolation(x, y, x1);
			for (int i = 1; i < (N - 1) * K + N; i++) {
				if (i == k * (K + N)) {
					k++;
				}
				else {
					//	cout << x1[i] << " " << y1[i] << endl;
				}
			}
			//cout << "Metod Simpsona:" << endl;
			long double shag;
			shag = x1[2] - x1[1];
			cout << "Shag:" << shag << endl;
			long double zS1 = 0;
			zS1 = y[0] + y[8];
			for (int i = 1; i < iter; i++) {
				if (i % 2 == 0) {
					zS1 += 2 * y1[i];
				}
				else {
					zS1 += 4 * y1[i];
				}
			}
			zS1 = zS1 * (shag / 3);
			//cout << "Second znachenie integrala methodom Simpsona:" << zS1 << endl << endl;
			pogresh1 = 0;
			pogresh1 = abs(zS- zS1);
			//cout << "Tochnost visheslenia:" << pogresh<< endl;
			j = 1;
			if (pogresh1 > 0.000001) {
				j = 0;
			}
			zS = zS1;
			//cout << "Iteration tochnosti:" << iter1 << endl;
			cout << endl;
			
		}

		cout << "Second znachenie integrala methodom Simpsona:" << zS << endl;
		cout << "Tochnost vicheslenia (method Simpsona):" << pogresh1 << endl;
		fout1 << zS <<" " << pogresh1 << endl;
		//fout1 << "Tochnost vicheslenia (method Simpsona):" << pogresh1 << endl;
		cout << "Iteration tochnosti (method Simpsona):" << iter2 << endl;
		cout << "______________________" << endl << endl;

		//Точность вычисления метода priamoygolnika

		cout << "Tochnost vichislenia metoda priamoygolnika:" << endl;
		int iter3 = 0;
		K = 10;
		long double pogresh2;
		for (int j = 0; j < 1; j) {
			K = 2 * K;
			cout << "K=" << K << endl;
			iter3++;
			int iter = 0;
			double *x1 = new double[(N - 1) * K + N];
			for (int i = 0; i < ((N - 1) * K + N); i++) {
				x1[i] = a + (i * h) / (K + 1);
				iter = i;
			}
			cout << "Iteration:" << iter << endl;

			double *y1 = CubicInterpolation(x, y, x1);
			for (int i = 1; i < (N - 1) * K + N; i++) {
				if (i == k * (K + N)) {
					k++;
				}
				else {
					//	cout << x1[i] << " " << y1[i] << endl;
				}
			}
			//cout << "Metod priamoygolnika:" << endl;
			long double shag;
			shag = x1[2] - x1[1];
			cout << "Shag:" << shag << endl;
			long double zT1 = 0;
			for (int i = 2; i <= iter; i = i + 2) {
				zT1 += (x1[i] - x1[i - 2])*y1[i - 1];

			}
			//cout << "Second znachenie integrala methodom priamoygolnika:" << zS1 << endl << endl;
			pogresh2 = 0;
			pogresh2 = abs(zT - zT1);
			//cout << "Tochnost visheslenia:" << pogresh<< endl;
			j = 1;
			if (pogresh1 > 0.000001) {
				j = 0;
			}
			zT = zT1;
			//cout << "Iteration tochnosti:" << iter1 << endl;
			cout << endl;

		}

		cout << "Second znachenie integrala methodom priamoygolnika:" << zT << endl;
		cout << "Tochnost vicheslenia (metodom priamoygolnika):" << pogresh2 << endl;
		fout2 << zT << " " << pogresh2 << endl;
		//fout2 << "Tochnost vicheslenia (method priamoygolnika):" << pogresh1 << endl;
		cout << "Iteration tochnosti (method priamoygolnika):" << iter3 << endl;
}

double *CubicInterpolation(double *src_x, double *src_y, double *new_x) { // src_x - начальные х, src_y - начальные у, new_x - новые х для нахождения у

	double *H = new double[N - 1];
	double **gamma = new double*[1];
	*gamma = new double[N - 2];
	double **matrix = new double*[N - 2];
	double *c = new double[N - 2]; 
	double *d = new double[N - 1];
	double *b = new double[N - 1]; 
	double *new_y = new double[(N - 1) * K + N]; 
	double *ci = new double[N - 1]; 
	double *a = new double[N - 1]; 

	for (int i = 0; i < N - 2; i++) {
		matrix[i] = new double[N - 2];
	}

	for (int i = 0; i < N - 1; i++) 
	{
		H[i] = src_x[i + 1] - src_x[i];
	}

	for (int i = 0; i < N - 2; i++)  
	{
		gamma[0][i] = 3 * ((src_y[i + 2] - src_y[i + 1]) / (H[i + 1]) - (src_y[i + 1] - src_y[i]) / (H[i]));
	}

	for (int i = 0; i < N - 2; i++) {
		for (int j = 0; j < N - 2; j++) {
			matrix[i][j] = 0;
		}
	}

	for (int i = 0; i < N - 2; i++) {   
		matrix[i][i] = 2 * (H[i] + H[i + 1]);
		if (i == 0) {
			matrix[i][i + 1] = H[i + 1];
		}
		else {
			if (i + 1 == N - 2) {
				matrix[i][i - 1] = H[i];
			}
			else {
				matrix[i][i + 1] = H[i + 1];
				matrix[i][i - 1] = H[i];
			}
		}
	}

	c = Jordan(matrix, gamma, N - 2); 

	for (int i = 0; i < N - 1; i++) { 
		if (i == 0) {
			ci[0] = 0;
		}
		else {
			ci[i] = c[i - 1];
		}
	}

	for (int i = 0; i < N - 2; i++) {  
		d[i] = ((ci[i + 1] - ci[i]) / (3.0 * H[i]));
	}
	d[N - 2] = (-1 * ci[N - 2]) / (3.0 * H[N - 2]);

	for (int i = 0; i < N - 2; i++) { 
		b[i] = ((src_y[i + 1] - src_y[i]) / H[i]) - (((ci[i + 1] + 2 * ci[i]) / 3.0))* H[i];
	}
	b[N - 2] = ((src_y[N - 1] - src_y[N - 2]) / H[N - 2]) - (2.0 / 3.0)*(ci[N - 2] * H[N - 2]);

	for (int i = 0; i < N - 1; i++) {
		a[i] = src_y[i];
	}

	
	for (int k = 0; k < N - 1; k++) {
		for (int i = 1 + k * (K + 1); i < K + 2 + k * (K + 1); i++) {
			new_y[i] = a[k] + (b[k] * (new_x[i] - src_x[k])) + (ci[k] * pow((new_x[i] - src_x[k]), 2)) + (d[k] * pow((new_x[i] - src_x[k]), 3));
		}
	}
	int k = 0;
	for (int i = 0; i < ((K + 1) * N - 1); i += K + 1) {
		new_y[i] = src_y[k];
		k++;
	}

	return new_y;
}

double *Jordan(double **a, double **b, int N) {

	double *x = new double[N];

	float R;
	for (int k = 0; k < N; k++)
	{
		for (int i = 0; i < N; i++)
		{
			if (i != k)
			{
				R = a[i][k] / a[k][k];
				for (int j = k; j < N; j++)
				{
					a[i][j] -= a[k][j] * R;
				}
				b[0][i] -= b[0][k] * R;
			}
		}
	}
	for (int i = 0; i < N; i++)
	{
		x[i] = b[0][i] / a[i][i];
	}

	return x;
}
