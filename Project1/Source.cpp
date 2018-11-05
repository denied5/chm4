//Cubic spline interpolation program
//when we have two columns of data x and y in input file:
//
//x0 y0
//x1 y1
//...
//xn yn
//
//and we want to find such function f(x)  
//where f(xi) = yi
//and f(x) is cubic function on every [x_k-1, x_k] segment
//and f(x), f'(x), f''(x) are continual
//the result is four columns of cubic polinom coefficients
#define _CRT_SECURE_NO_WARNINGS
#include <math.h>
#include <stdio.h>
#include <process.h>
#include <iostream>
#include <iomanip>  
using namespace std;
float *x, *y, *h, *l, *delta, *lambda, *c, *d, *b;
int N = 24;
char filename[256];
FILE* InFile = NULL;
float dopExArrey[4]= {-0.62 , 1.12 , 2.34, 4 };

void readmatrix(FILE* InFile) {
	int i = 0;
	//read matrixes a and b from input file
	for (i = 0; i < N + 1; i++) {
		fscanf_s(InFile, "%f", &x[i]);
		fscanf_s(InFile, "%f", &y[i]);
	}
}

void allocmatrix() {
	//allocate memory for matrixes
	x = new float[N + 1];
	y = new float[N + 1];
	h = new float[N + 1];
	l = new float[N + 1];
	delta = new float[N + 1];
	lambda = new float[N + 1];
	c = new float[N + 1];
	d = new float[N + 1];
	b = new float[N + 1];
}
void freematrix() {
	delete[] x;
	delete[] y;
	delete[] h;
	delete[] l;
	delete[] delta;
	delete[] lambda;
	delete[] c;
	delete[] d;
	delete[] b;
}

void printresult() {
	int k = 0;
	printf("\nA[k]\tB[k]\tC[k]\tD[k]\n");
	for (k = 1; k <= N; k++) {
		printf("%f\t%f\t%f\t%f\n", y[k], b[k], c[k], d[k]);
	}
}

void setFile(FILE* InFile, float start = -1, float end = 3)
{

	float fy;
	float fx =start;
	float step = (end - start)/N;// h = (b - a) / n; b = 4, a = -1, n = 6
	while (fx <= end + 0.01)
	{
		fy = (1 + fx)* exp(-2 * fx);
		fprintf(InFile, "%f\t%f\n", fx, fy);
		fx += step;
		
	}
	
}

void testresult(float start = -1, float end = 3) {
	
	float step = (end - start) / (4 * N); //з кроком h = (b - a) / (4n) для відрізку[a - H, b + H]
	FILE* OutFile = fopen("out.txt", "wt");
	cout << "    Fi(x)       F(x)        []        %        x"<<endl;
	for (float s = start; s <= end + step; s += step) {
		//find k, where s in [x_k-1; x_k]
		int k;
		for (k = 1; k <= N; k++) {
			if (s >= x[k - 1] && s <= x[k]) {
				break;
			}
		}
		float fy = (1 + s)* exp(-2 * s);
		float F = y[k] + b[k] * (s - x[k]) + c[k] * pow(s - x[k], 2) + d[k] * pow(s - x[k], 3);
		cout.precision(2);
		cout.width(7);
		cout << setw(8) << F <<  "  |" << setw(8) << fy  <<"  |" << setw(8) << F - fy << "  |" << setw(8) << ((F-fy)*100)/F << "  |" << setw(8) << s <<endl;
		fprintf(OutFile, "%f\t%f\t%f\n",s, fy, F);
	}



	fclose(OutFile);
}

void accomulate()
{
	int k = 0;
	for (k = 1; k <= N; k++) {
		h[k] = x[k] - x[k - 1];
		l[k] = (y[k] - y[k - 1]) / h[k];
	}
	delta[1] = -h[2] / (2 * (h[1] + h[2]));
	lambda[1] = 1.5*(l[2] - l[1]) / (h[1] + h[2]);
	for (k = 3; k <= N; k++) {
		delta[k - 1] = -h[k] / (2 * h[k - 1] + 2 * h[k] + h[k - 1] * delta[k - 2]);
		lambda[k - 1] = (3 * l[k] - 3 * l[k - 1] - h[k - 1] * lambda[k - 2]) /
			(2 * h[k - 1] + 2 * h[k] + h[k - 1] * delta[k - 2]);
	}
	c[0] = 0;
	c[N] = 0;
	for (k = N; k >= 2; k--) {
		c[k - 1] = delta[k - 1] * c[k] + lambda[k - 1];
	}
	for (k = 1; k <= N; k++) {
		d[k] = (c[k] - c[k - 1]) / (3 * h[k]);
		b[k] = l[k] + (2 * c[k] * h[k] + h[k] * c[k - 1]) / 3;
	}
}
void main() {
	
	InFile = fopen("in.txt", "w+t");
	rewind(InFile);
	setFile(InFile);
	rewind(InFile);
	allocmatrix();
	readmatrix(InFile);
	accomulate();
	printresult();
	testresult();
	
	//Додаткове завдання 
	//Точки 1)(-1;4);(1,4;3),(2;3)(4;4)
	//n = 6 

	cout << endl << endl << endl;
	cout << "N = 6"<<endl;
	cout << "    Fi(zi)     F(zi)        []        %        zi" << endl;
	for (size_t i = 0; i < 4; i++)
	{
		int k;
		for (k = 1; k < N; k++) {
			if (dopExArrey[i] >= x[k - 1] && dopExArrey[i] <= x[k]) {
				break;
			}
		}
		float fy = (1 + dopExArrey[i])* exp(-2 * dopExArrey[i]);
		float F = y[k] + b[k] * (dopExArrey[i] - x[k]) + c[k] * pow(dopExArrey[i] - x[k], 2) + d[k] * pow(dopExArrey[i] - x[k], 3);
		cout << setw(8) << F << "  |" << setw(8) << fy << "  |" << setw(8) << F - fy << "  |" << setw(8) << ((F - fy) * 100) / F << "  |" << setw(8) << dopExArrey[i] << endl;
	}


	cout << endl << endl << endl;
	cout << "N = 24"<<endl;
	cout << "    Fi(zi)     F(zi)        []        %        zi" << endl;
	N = 24;
	FILE* InFile2 = fopen("in2.txt", "w+t");
	rewind(InFile2);
	setFile(InFile2);
	rewind(InFile2);
	allocmatrix();
	readmatrix(InFile2);
	accomulate();
	for (size_t i = 0; i < 4; i++)
	{
		int k;
		for (k = 1; k < N; k++) {
			if (dopExArrey[i] >= x[k - 1] && dopExArrey[i] <= x[k]) {
				break;
			}
		}
		float fy = (1 + dopExArrey[i])* exp(-2 * dopExArrey[i]);
		float F = y[k] + b[k] * (dopExArrey[i] - x[k]) + c[k] * pow(dopExArrey[i] - x[k], 2) + d[k] * pow(dopExArrey[i] - x[k], 3);
		cout << setw(8) << F << "  |" << setw(8) << fy << "  |" << setw(8) << F - fy << "  |" << setw(8) << ((F - fy) * 100) / F << "  |" << setw(8) << dopExArrey[i] << endl;
	}
	freematrix();
	


	system("pause");
}