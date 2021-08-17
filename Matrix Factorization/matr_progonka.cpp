#include <iostream>
#include <cmath>
#include "class_matrix.h"
#define n 10
#define exp 2.7182818284590
using namespace std;
int main()
{
	int i,j;
	double h = 1.0/ (n-1); 
	matrix A,B,C,betta,U;
	matrix *alfa=new matrix[n], *D = new matrix[n];
	
	A.SetSize(n,n);B.SetSize(n,n);C.SetSize(n,n);
	betta.SetSize(n,n); U.SetSize(n,n);
	for(i=0;i<n;i++){
		alfa[i].SetSize(n,n);
		D[i].SetSize(1,n);
	}
	cout<<"\a";
	for(i=0;i<n;i++){
			A(i,i) = 1.0;
			B(i,i) = 1.0;
			C(i,i) = 1.0;
			C(i,i+1) = C(i,i-1) = -1.0;
			alfa[1](i,i) = 1.0;
			betta[1](i,0) = h*(2.0 * cos(i*h) + 3.0 * pow(exp, 2.0*i*h));

			U[n-1](0,i) = sin(i*h +2.0) + pow(exp, 2.0*i*h + 3.0);
		}
	cout<<"\a";
	for(i=0;i<n;i++)
		for(j=0;j<n;j++){
			D[i](0,j) = (2.0* sin(i*h + 2.0*j*h) - 16.0 * pow(exp, 2.0*i*h + 3.0*j*h)) * h*h;
		}
	for(i=0;i<n;i++){
		D[1](0,i) = D[1](0,i) + sin(2.0 * i*h)+ pow(exp, 3.0 * i*h);
		D[n-2](0,i) = D[n-2](0,i) + sin(1.0 + 2.0*i*h) + pow(exp, 2.0 + 3.0*i*h);
	}
	cout<<D[1];
	for(i=2;i<n;i++){
		alfa[i] = (!(C - A*alfa[i-1]))*B;
		cout<<"\a";
		betta[i] = (!(C - A*alfa[i-1]))*(betta[i-1]*A + D[i-1]);
	}
	
	for(i=n-2;i>=0;i--)
		U[i]= U[i+1]*alfa[i+1] + betta[i+1];
	for(i=0;i<n;i++)
		cout<<U[i];
	return 0;
}