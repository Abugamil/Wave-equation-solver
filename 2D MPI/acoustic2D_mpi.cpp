#pragma warning (disable : 4996)
#pragma warning (disable : 4101)
#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include "mpi.h"
using namespace std;
//Имя файла результата и путь к нему (путь\название файла.расширение)
char *File = "C:\\acoustic2D.dat";
//Размер сетки
const int n = 63;

double exp(double index)
{
	const double e = 2.7182818284590;
	return (pow(e, index));
}

double sqr (double n)
{
	return n*n;
}

void Read_TecPlot(const char *file_name,double a[n+1][n+1])
{
	FILE	*outfile = fopen(file_name,"r");
	int ii,jj;
	for(int i=0;i<=n;i++)
		for(int j=0;j<=n;j++)
			fscanf(outfile,"%li \t %li \t %+.6E \n",&ii,&jj,&a[i][j]);
}

void Read_File(const char *file_name, double a[n+1][n+1])
{
	ifstream inn(file_name,ios :: in);
	for(int i=0;i<=n;i++)
		for(int j=0;j<=n;j++)
			inn>>a[i][j];
}

void Init_TecPlot(const char *file_name)
{
	FILE	*outfile = fopen(file_name,"w");
	fprintf(outfile,"TITLE = \"Example: 2D Plot Akustika\"\n");
	fprintf(outfile,"VARIABLES = \"X\", \"Y\", \"P\"\n");
	fclose(outfile);
}

void Write_TecPlot(const char *file_name,int T, double P[n+1][n+1])
{
	FILE	*outfile = fopen(file_name,"a");
	fprintf(outfile,"ZONE T=\"%li\", I=%li, J=%li, F=POINT\n",T,n+1,n+1);
	for(int i=0;i<=n;i++)
		for(int j=0;j<=n;j++)
			fprintf(outfile,"%li \t %li \t %+.6E \n",i,j,P[i][j]);
	fclose(outfile);
}

void Write_TecPlot_B(const char *file_name,int T, double P[n+1][n+1])
{
	FILE	*outfile = fopen(file_name,"a");
	fprintf(outfile,"ZONE T=\"%li\", I=%li, J=%li, F=POINT\n",T,n+1,n+1);
	for(int i=0;i<=n;i++){
		fprintf(outfile,"%li \t %li \t %+.6E \n",i,0,P[i][0]);
	}
	for(int i=0;i<=n;i++){
		fprintf(outfile,"%li \t %li \t %+.6E \n",i,n,P[i][n]);
	}
	for(int i=0;i<=n;i++){
		fprintf(outfile,"%li \t %li \t %+.6E \n",0,i,P[0][i]);
	}
	for(int i=0;i<=n;i++){
		fprintf(outfile,"%li \t %li \t %+.6E \n",n,i,P[n][i]);
	}
	fclose(outfile);
}

int main(int argc,char **argv)
{
	//Колличество процессов
	int size;
	//индекс процесса
	int rank;
	//число ПИ
	const double Pi = 3.1415926535;
	//Колличество итераций по времени
	const int T = 400;
	//Длина сетки
	const double l=100.;
	//Скорость звука
	const double c = 330.;
	//Амплитуда волны
	const double A = 30.564319;
	//Частота волны
	const double w = 0.03;
	//Значение погрешности
	const double eps = 0.1;
	//???
	const double e = 0.1;
	//Старт расчета
	double StartTime = 0;
	//Конец расчета
	double EndTime;
	
	int k,r1,rc;
	double eps1,eps2,x0,y0,tau,h1;
	double diffnorm, dx2,dy2;
	double pn[n+1][n+1],pn1[n+1][n+1],p[n+1][n+1],pn0[n+1][n+1],pn2[n+1][n+1],f[n+1][n+1];
	double a[n+1][n+1];

	rc = MPI_Init(&argc, &argv);
	MPI_Status stat;

	if (rc != MPI_SUCCESS)
	{
		printf("\nERROR initializing MPI. Exit.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
		return 0;
	}

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	if(rank == 0)
	{
		StartTime = MPI_Wtime();
		printf("Start program\n");
	}

	if (size == 1)
	{
		printf("Process number equals to 1 \n");
	}

	k = (int)sqrt(size);
	r1 = (int)(rank - rank%k) / k;

	dx2 = dy2= sqr(1. / (double) n);
	tau = 0.0065;
	for(int i=0;i<=n;i++)
	{
		for(int j=0;j<=n;j++)
		{
			pn[i][j] = 0.0; p[i][j] = 0.0;
			pn1[i][j] = 0.0; pn0[i][j] = 0.0;
			pn2[i][j] = 0.0; f[i][j] = 0.;
			a[i][j] = 1.;
		}
	}
	x0 = 0.5;y0 = 0.5;
	//cout<<(2.*Pi + 2.*Pi*w*abs( sqrt(sqr(x0) + sqr(y0)) )/330. ) / (2.*Pi*w)<<endl;
	//**************************************************************
	Init_TecPlot(File);
	if(rank == 0)
		printf("chislo Kuranta = %f && %f\n",sqr(tau)/dx2,sqr(tau)/dy2);

	for(int t = 0; t <= T; t++)
	{
		if(t <= 100)
		{
			f[n/2][n/2] = A*sin( 2.*Pi*w*(t - abs( sqrt(sqr(x0) + sqr(y0)) )/330. ));
			//f[70][70] = A*sin( 2.*Pi*w*(t - abs( sqrt(sqr(x0) + sqr(y0)) )/330. )); 
		}
		do
		{
			diffnorm = 0.;
			for (int i = (int)(rank%k)*(n+1)/k + 1; i < (int)(rank%k + 1)*(n+1)/k - 1; i++)
			{
				for(int j = (int)(r1%k)*(n+1)/k + 1;j < (int)(r1%k + 1)*(n+1)/k - 1; j++)
				{
					p[i][j] = sqr(tau)*(
						a[i][j]*(pn[i-1][j] - 2.0*pn[i][j] + pn[i+1][j]) / dx2 + 
						a[i][j]*(pn[i][j-1] - 2.0*pn[i][j] + pn[i][j+1]) / dy2 + 
						sqr(l)*f[i][j]/c 
						) + 2.0*pn[i][j] - pn2[i][j];

					diffnorm += sqr(p[i][j] - pn[i][j]);
				}
			}

			//non-reflecting boundary
			for (int j=0; j<=n; j++)
			{
				p[n][j] = pn[n][j] - tau*(pn[n][j]-pn[n-1][j]) / sqrt(dx2);
				p[0][j] = pn[0][j] + tau*(pn[1][j]-pn[0][j]) / sqrt(dx2);
			}
			for (int i=0; i<=n; i++)
			{
				p[i][n] = pn[i][n] - tau*(pn[i][n]-pn[i][n-1]) / sqrt(dy2);
				p[i][0] = pn[i][0] + tau*(pn[i][1]-pn[i][0]) / sqrt(dy2);
			}

			for (int i = (int)(rank%k)*(n+1)/k; i < (int)(rank%k + 1)*(n+1)/k -1; i++)
			{
				for (int j = (int)(r1%k)*(n+1)/k; j < (int)(r1%k + 1)*(n+1)/k -1; j++)
				{
					pn2[i][j] = pn[i][j];
					pn[i][j] = p[i][j];
				}
			}

			diffnorm = sqrt(diffnorm);
		}while(diffnorm > eps);
		
		if(rank == 0 && t%10 == 0 )
		{
			Write_TecPlot(File,t,p);
		}
	}

	if(rank == 0)
	{
		EndTime = MPI_Wtime();
		printf("Compute time: %f\n",EndTime - StartTime);
	}
	MPI_Finalize();

	return 0;
}