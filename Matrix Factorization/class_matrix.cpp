#include <iostream>
#include<math.h>
#include<memory.h>
#include"class_matrix.h"
using namespace std;
#define max(a,b)    (((a) > (b)) ? (a) : (b))
#define min(a,b)    (((a) < (b)) ? (a) : (b))

size_t matrix::RowNo ()
{
	return Row;
}

size_t matrix::ColNo ()
{
	return Col;
}

// constructor
matrix::matrix (size_t row, size_t col)
{
	RowSiz = Row = row;
	ColSiz = Col = col;

	Val = new double *[row];

	for(size_t i = 0; i < row; i++)
		Val[i] = new double[col];
	Null ();
}

// copy constructor
matrix::matrix (const matrix& m)
{
	RowSiz = Row = m.Row;
	ColSiz = Col = m.Col;

	Val = new double *[Row];

	size_t colsize = Col * sizeof (double);

	for(size_t i = 0; i < Row; i++)
	{
		Val[i] = new double[Col];
		memcpy (Val[i], m.Val[i], colsize);
	}
}

// destructor
matrix::~matrix ()
{
	for(size_t i = 0; i < RowSiz; i++)
		delete[]Val[i];
	delete[]Val;
}

//  reallocation method
void matrix::realloc (size_t row, size_t col)
{
	if((row == RowSiz) && (col == ColSiz))
	{
		Row = RowSiz;
		Col = ColSiz;

		return;
	}

	size_t i;

	double **Val1 = new double *[row];

	for(i = 0; i < row; i++)
		Val1[i] = new double[col];

	for(i = 0; i < RowSiz; i++)
		delete[]Val[i];
	delete[]Val;

	RowSiz = Row = row;
	ColSiz = Col = col;
	Val = Val1;

	Null ();

	return;
}

// public method for resizing matrix
void matrix::SetSize (size_t row, size_t col)
{
	size_t oldRow = Row;
	size_t oldCol = Col;

	if(row != RowSiz || col != ColSiz)
		realloc (row, col);

  return;
}

// subscript operator to get/set individual elements
double &matrix::operator() (size_t row, size_t col)
{
	return Val[row][col];
}

// input stream function
istream & operator >> (istream & istrm, matrix& m)
{
	for(size_t i = 0; i < m.Row; i++)
		for(size_t j = 0; j < m.Col; j++)
			istrm >> m.Val[i][j];
	return istrm;
}

// output stream function
ostream & operator << (ostream & ostrm, matrix& m)
{
	for(size_t i = 0; i < m.Row; i++)
	{
		for(size_t j = 0; j < m.Col; j++)
			cout << m.Val[i][j] << '\t';
		cout << endl;
	}
	return ostrm;
}

// assignment operator
matrix& matrix::operator = (const matrix& m)
{
	if(Row != m.Row || Col != m.Col)
		realloc (m.Row, m.Col);

	size_t colbyte = m.Col * sizeof (double);

	for(size_t i = 0; i < m.Row; i++)
		memcpy (Val[i], m.Val[i], colbyte);

	return *this;
}

// logical equal-to operator
bool operator == (const matrix& m1, const matrix& m2)
{
	if(m1.Row != m2.Row || m1.Col != m2.Col)
		return false;

	for(size_t i = 0; i < m1.Row; i++)
		for(size_t j = 0; j < m1.Col; i++)
			if(m1.Val[i][j] != m2.Val[i][j])
				return false;
	return true;
}

// logical no-equal-to operator
bool operator != (const matrix& m1, const matrix& m2)
{
	return ((m1 == m2) ? false : true);
}

// combined addition and assignment operator
matrix& matrix::operator += (const matrix& m)
{
	for(size_t i = 0; i < m.Row; i++)
		for(size_t j = 0; j < m.Col; j++)
			Val[i][j] += m.Val[i][j];
	return *this;
}

// combined subtraction and assignment operator
matrix& matrix::operator -= (const matrix& m)
{
	for(size_t i = 0; i < m.Row; i++)
		for(size_t j = 0; j < m.Col; j++)
			Val[i][j] -= m.Val[i][j];
	return *this;
}

// combined scalar multiplication and assignment operator
matrix& matrix::operator *= (const double &c)
{
	for(size_t i = 0; i < Row; i++)
		for(size_t j = 0; j < Col; j++)
			Val[i][j] *= c;
	return *this;
}

// combined matrix multiplication and assignment operator
matrix& matrix::operator *= (const matrix& m)
{
	*this = *this * m;
	return *this;
}

// combined scalar division and assignment operator
matrix& matrix::operator /= (const double& c)
{
	for(size_t i = 0; i < Row; i++)
		for(size_t j = 0; j < Col; j++)
			Val[i][j] /= c;
	return *this;
}

// combined power and assignment operator
matrix& matrix::operator ^= (const size_t& pow)
{
	for(size_t i = 2; i <= pow; i++)
		*this = *this * *this;
	return *this;
}

matrix matrix::operator + ()
{
	return *this;
}

// unary negation operator
matrix matrix::operator - ()
{
	matrix temp (Row, Col);

	for(size_t i = 0; i < Row; i++)
		for(size_t j = 0; j < Col; j++)
			temp.Val[i][j] = -Val[i][j];
	return temp;
}

// binary addition operator
matrix operator + (const matrix& m1, const matrix& m2)
{
	matrix temp (m1.Row, m1.Col);

	for(size_t i = 0; i < m1.Row; i++)
		for(size_t j = 0; j < m1.Col; j++)
			temp.Val[i][j] = m1.Val[i][j] + m2.Val[i][j];
	return temp;
}

// binary subtraction operator
matrix operator - (const matrix& m1, const matrix& m2)
{
	matrix temp (m1.Row, m1.Col);

	for(size_t i = 0; i < m1.Row; i++)
		for(size_t j = 0; j < m1.Col; j++)
			temp.Val[i][j] = m1.Val[i][j] - m2.Val[i][j];
	return temp;
}

// binary scalar multiplication operator
matrix operator * (const matrix& m, const double &no)
{
	matrix temp (m.Row, m.Col);

	for(size_t i = 0; i < m.Row; i++)
		for(size_t j = 0; j < m.Col; j++)
			temp.Val[i][j] = no * m.Val[i][j];
	return temp;
}

// binary matrix multiplication operator
matrix operator * (const matrix& m1, const matrix& m2)
{
	matrix temp (m1.Row, m2.Col);
	for(size_t i = 0; i < m1.Row; i++)
		for(size_t j = 0; j < m2.Col; j++)
		{
			temp.Val[i][j] = double (0);
			for(size_t k = 0; k < m1.Col; k++)
				temp.Val[i][j] += m1.Val[i][k] * m2.Val[k][j];
		}
	return temp;
}

//
matrix operator * (const double& no, const matrix& m)
{
	return (m * no);
}

matrix operator / (const matrix& m1, const matrix& m2)
{
	return (m1 * !m2);
}

matrix operator / (const matrix& m, const double& no)
{
	return (m * (1 / no));
}

matrix operator / (const double& no, const matrix& m)
{
	return (!m * no);
}

// binary power operator
matrix operator ^ (const matrix& m, const size_t & pow)
{
	matrix temp (m);

	for(size_t i = 2; i <= pow; i++)
		temp = temp * temp;

	return temp;
}

// unary transpose operator
matrix operator  ~ (const matrix& m)
{
	matrix temp (m.Col, m.Row);
	for(size_t i = 0; i < m.Row; i++)
		for(size_t j = 0; j < m.Col; j++)
			temp.Val[j][i] = m.Val[i][j];
	return temp;
}

// unary inversion operator
matrix operator ! (matrix m)
{
	size_t i, j, k;
	double a1, a2, *rowptr;
	matrix temp (m.Row, m.Col);
	temp.Unit ();

	for(k = 0; k < m.Row; k++)
	{
		int index = m.pivot (k);
		if(index != 0)
		{
			rowptr = temp.Val[k];
			temp.Val[k] = temp.Val[index];
			temp.Val[index] = rowptr;
		}
		a1 = m.Val[k][k];
		for(j = 0; j < m.Row; j++)
		{
			m.Val[k][j] /= a1;
			temp.Val[k][j] /= a1;
		}
		for(i = 0; i < m.Row; i++)
			if(i != k)
			{
				a2 = m.Val[i][k];
				for(j = 0; j < m.Row; j++)
				{
					m.Val[i][j] -= a2 * m.Val[k][j];
					temp.Val[i][j] -= a2 * temp.Val[k][j];
				}
			}
	}

	return temp;
}

// solve simultaneous equation
matrix matrix::Solve (const matrix& v) const
{
	size_t i, j, k;
	double a1;
	matrix temp (Row, Col + v.Col);

	for(i = 0; i < Row; i++)
	{
		for(j = 0; j < Col; j++)
			temp.Val[i][j] = Val[i][j];

		for(k = 0; k < v.Col; k++)
			temp.Val[i][Col + k] = v.Val[i][k];
	}

	for(k = 0; k < Row; k++)
	{
		int indx = temp.pivot (k);
		a1 = temp.Val[k][k];

		for(j = k; j < temp.Col; j++)
			temp.Val[k][j] /= a1;

		for(i = k + 1; i < Row; i++)
		{
			a1 = temp.Val[i][k];

			for(j = k; j < temp.Col; j++)
				temp.Val[i][j] -= a1 * temp.Val[k][j];
		}
	}

	matrix s (v.Row, v.Col);

	for(k = 0; k < v.Col; k++)
		for(int m = int (Row) - 1; m >= 0; m--)
		{
			s.Val[m][k] = temp.Val[m][Col + k];

			for(j = m + 1; j < Col; j++)
				s.Val[m][k] -= temp.Val[m][j] * s.Val[j][k];
		}

	return s;
}

// set zero to all elements of this matrix
void matrix::Null (const size_t& row, const size_t& col)
{
	if(row != Row || col != Col)
		realloc (row, col);

	for(size_t i = 0; i < Row; i++)
		for(size_t j = 0; j < Col; j++)
			Val[i][j] = double (0);
	return;
}

// set zero to all elements of this matrix
void matrix::Null ()
{
	for(size_t i = 0; i < Row; i++)
		for(size_t j = 0; j < Col; j++)
			Val[i][j] = double (0);
	return;
}

// set this matrix to unity
void matrix::Unit (const size_t & row)
{
	if(row != Row || row != Col)
		realloc (row, row);
	for(size_t i = 0; i < Row; i++)
		for(size_t j = 0; j < Col; j++)
			Val[i][j] = i == j ? double (1) : double (0);
	return;
}

// set this matrix to unity
void matrix::Unit ()
{
	size_t row = min (Row, Col);
	Row = Col = row;

	for(size_t i = 0; i < Row; i++)
		for(size_t j = 0; j < Col; j++)
			Val[i][j] = i == j ? double (1) : double (0);
	return;
}

// private partial pivoting method
int matrix::pivot (size_t row)
{
	int k = int (row);
	double amax, temp;
	amax = -1;

	for(size_t i = row; i < Row; i++)
		if((temp = fabs (Val[i][row])) > amax && temp != 0.0)
		{
			amax = temp;
			k = i;
		}

	if(Val[k][row] == double (0))
		return -1;

	if(k != int (row))
	{
		double *rowptr = Val[k];
		Val[k] = Val[row];
		Val[row] = rowptr;
		return k;
	}

	return 0;
}

// calculate the determinant of a matrix
double matrix::Det ()
{
	size_t i, j, k;
	double piv, detVal = double (1);
	matrix temp (*this);
	for(k = 0; k < Row; k++)
	{
		int index = temp.pivot (k);
		if(index == -1)
			return double (0);
		if(index != 0)
			detVal = -detVal;
		detVal = detVal * temp.Val[k][k];
		for(i = k + 1; i < Row; i++)
		{
			piv = temp.Val[i][k] / temp.Val[k][k];
			for(j = k + 1; j < Row; j++)
				temp.Val[i][j] -= piv * temp.Val[k][j];
		}
	}
	return detVal;
}

// calculate the norm of a matrix
double matrix::Norm ()
{
	double retVal = double (0);
	for(size_t i = 0; i < Row; i++)
		for(size_t j = 0; j < Col; j++)
			retVal += Val[i][j] * Val[i][j];
	retVal = sqrt (retVal);
	return retVal;
}

// calculate the condition number of a matrix
double matrix::Cond ()
{
	matrix inv (Row, Col);
	inv = !(*this);
	double retVal = Norm () * inv.Norm ();
	return retVal;
}

// calculate the cofactor of a matrix for a given element
double matrix::Cofact (size_t row, size_t col)
{
	size_t i, i1, j, j1;
	matrix temp (Row - 1, Col - 1);
	for(i = i1 = 0; i < Row; i++)
	{
		if(i == row)
			continue;
		for(j = j1 = 0; j < Col; j++)
		{
			if(j == col)
				continue;
			temp.Val[i1][j1] = Val[i][j];
			j1++;
		}
		i1++;
	}
	double cof = temp.Det ();
	if((row + col) % 2 == 1)
		cof = -cof;
	return cof;
}

// calculate adjoin of a matrix
matrix matrix::Adj ()
{
	matrix temp (Row, Col);
	for(size_t i = 0; i < Row; i++)
		for(size_t j = 0; j < Col; j++)
			temp.Val[i][j] = Cofact (i, j);
	temp = ~temp;
	return temp;
}

bool matrix::IsSquare ()
{
	return (Row == Col);
}

// Determine if the matrix is singular
bool matrix::IsSingular ()
{
	if(Row != Col)
		return false;
	return (Det () == double (0));
}

// Determine if the matrix is diagonal
bool matrix::IsDiagonal ()
{
	if(Row != Col)
		return false;
	for(size_t i = 0; i < Row; i++)
		for(size_t j = 0; j < Col; j++)
			if(i != j && Val[i][j] != double (0))
				return false;
	return true;
}

// Determine if the matrix is scalar
bool matrix::IsScalar ()
{
	if(!IsDiagonal ())
		return false;
	double v = Val[0][0];
	for(size_t i = 1; i < Row; i++)
		if(Val[i][i] != v)
			return false;
	return true;
}

// Determine if the matrix is a unit matrix
bool matrix::IsUnit ()
{
	if(IsScalar () && Val[0][0] == double (1))
		return true;
	return false;
}

// Determine if this is a null matrix
bool matrix::IsNull ()
{
	for(size_t i = 0; i < Row; i++)
		for(size_t j = 0; j < Col; j++)
			if(Val[i][j] != double (0))
				return false;
	return true;
}

// Determine if the matrix is symmetric
bool matrix::IsSymmetric ()
{
	if(Row != Col)
		return false;
	for(size_t i = 0; i < Row; i++)
		for(size_t j = 0; j < Col; j++)
			if(Val[i][j] != Val[j][i])
				return false;
	return true;
}

// Determine if the matrix is skew-symmetric
bool matrix::IsSkewSymmetric ()
{
	if(Row != Col)
		return false;
	for(size_t i = 0; i < Row; i++)
		for(size_t j = 0; j < Col; j++)
			if(Val[i][j] != -Val[j][i])
				return false;
	return true;
}

// Determine if the matrix is upper triangular
bool matrix::IsUpperTriangular ()
{
	if(Row != Col)
		return false;
	for(size_t i = 1; i < Row; i++)
		for(size_t j = 0; j < i - 1; j++)
			if(Val[i][j] != double (0))
				return false;
	return true;
}

// Determine if the matrix is lower triangular
bool matrix::IsLowerTriangular ()
{
	if(Row != Col)
		return false;
	for(size_t j = 1; j < Col; j++)
		for(size_t i = 0; i < j - 1; i++)
			if(Val[i][j] != double (0))
				return false;
	return true;
} 