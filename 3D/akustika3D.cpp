#pragma warning (disable : 4996)
#pragma warning (disable : 4101)
#include <cmath>
#include <fstream>
#include <stdio.h>
#include "iofile.h"
using namespace std;

double exp(double index)
{
	const double e = 2.7182818284590;
	return (pow(e, index));
}

double sqr (double n)
{
	return n*n;
}

double ***PNew, ***POld, ***P, ***P0, ***F;

const int N = 51;
const int T = 20;

const double PI = 3.1415926535;
const double c = (double)347;
const double L = (double)30;
const double dt = 0.01;
const double dx = L / (double)(N - 1);
const double dy = L / (double)(N - 1);
const double dz = L / (double)(N - 1);

const double dt2 = sqr(dt);
const double dx2 = sqr(dx);
const double dy2 = sqr(dy);
const double dz2 = sqr(dz);

void FillZero(double ***M)
{
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			for (int k = 0; k < N; k++){
				M[i][j][k] = (double) 0;
			}
		}
	}
}

bool check(double ***un, double ***u)
{
	double max = 0.0;
	for(int i=1; i<N-1; i++){
		for(int j=1; j<N-1; j++){
			for(int k=1; k<N-1; k++)
			{
				if(fabs(1.0 - u[i][j][k]/un[i][j][k]) > max)
					max = fabs(1.0 - u[i][j][k]/un[i][j][k]);
				else
					max = max;
			}
		}
	}
	if(max >= 1e-02)
		return false;
	else
		return true;
}

void Swap (double ***Old, double ***New)
{
	for (int i = 1; i < N-1; i++){
		for (int j = 1; j < N-1; j++){
			for (int k = 1; k < N-1; k++){
				Old[i][j][k] = New[i][j][k];
			}
		}
	}
}

void InitAllMassives()
{
	***PNew = new double **[N];
	***POld = new double **[N];
	***P = new double **[N];
	***P0 = new double **[N];
	***F = new double **[N];
	for (int i = 0; i < N; i++)
	{
		PNew[i] = new double *[N];
		POld[i] = new double *[N];
		P[i] = new double *[N];
		P0[i] = new double *[N];
		F[i] = new double *[N];

		for (int j = 0; j < N; j++)
		{
			PNew[i][j] = new double [N];
			POld[i][j] = new double [N];
			P[i][j] = new double [N];
			P0[i][j] = new double [N];
			F[i][j] = new double [N];
			for (int k = 0; k < N; k++)
			{
				F[i][j][k] = (double) 0;
				P[i][j][k] = (double) 0;
				P0[i][j][k] = (double) 0;
				PNew[i][j][k] = (double) 0;
				POld[i][j][k] = (double) 0;
			}
		}
	}
}

void delAllMassives()
{
	for(int i=0; i<N; i++)
	{
		delete F[i];
	}
}

void main()
{
	double diffnorm = (double)0;
	
	InitAllMassives();

	InitTecPlot3D ("Akustika3D.dat");
	WriteTecPlot3D ("Akustika3D.dat",PNew, 0,N,N,N,dx,dy,dz);

	for(int t = 1; t <= T; t++)
	{
		F[1][N/2][N/2] = 2.0*sin( 2.0*PI*100.0 + 90.0/360.0*2.0*PI);
		//2.0*sin( 2.*PI*0.001*(t - abs( sqrt(sqr(70) + sqr(30)) )/330. ));

		/*swap(P, POld);
		swap(PNew, P);*/
		/*Swap(P, POld);
		Swap(POld, PNew);*/

		for (int i = 1; i < N-1; i++){
			for (int j = 1; j < N-1; j++){
				for (int k = 1; k < N-1; k++){
					P[i][j][k] = POld[i][j][k];
					POld[i][j][k] = PNew[i][j][k];
				}
			}
		}

		for(int i=1; i<N-1; i++){
			for(int j=1; j<N-1; j++){
				for(int k=1; k<N-1; k++)
				{
					PNew[i][j][k] = 2.0*POld[i][j][k] - P[i][j][k] + dt2*sqr(c)*( 
						(POld[i+1][j][k] - 2*POld[i][j][k] + POld[i-1][j][k]) / dx2 + 
						(POld[i][j+1][k] - 2*POld[i][j][k] + POld[i][j-1][k]) / dy2 + 
						(POld[i][j][k+1] - 2*POld[i][j][k] + POld[i][j][k-1]) / dz2 + 
						F[i][j][k]);

				}
			}
		}

		WriteTecPlot3D("Akustika3D.dat",PNew, t,N,N,N,dx,dy,dz);
	}
	/*for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++)
	for (int k = 0; k < N; k++)
	printf("%f ",PNew[i][j][k]);
	*/
	delAllMassives();
}