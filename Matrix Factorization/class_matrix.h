using namespace std;
class matrix 
{
	double **Val;
	size_t Row, Col, RowSiz, ColSiz;
	void realloc (size_t row, size_t col);
	int pivot (size_t row);

public: 

//	Constructors
	matrix (const matrix& m);
	matrix (size_t row = 5, size_t col = 5);

//	Destructor
	~matrix ();

//	Value extraction method
	size_t RowNo ();
	size_t ColNo ();

//	Subscript operator
	double& operator () (size_t row, size_t col);

//	Unary operators
	matrix operator + ();
	matrix operator - ();

//	Assignment operators
	matrix& operator = (const matrix& m);

//	Combined assignment - calculation operators
	matrix& operator += (const matrix& m);
	matrix& operator -= (const matrix& m);
	matrix& operator *= (const matrix& m);
	matrix& operator *= (const double& c);
	matrix& operator /= (const double& c);
	matrix& operator ^= (const size_t& pow);

//	Logical operators
	friend bool operator == (const matrix& m1, const matrix& m2);
	friend bool operator != (const matrix& m1, const matrix& m2);

//	Calculation operators
	friend matrix operator + (const matrix& m1, const matrix& m2);
	friend matrix operator - (const matrix& m1, const matrix& m2);
	friend matrix operator * (const matrix& m1, const matrix& m2);
	friend matrix operator * (const matrix& m, const double& no);
	friend matrix operator * (const double& no, const matrix& m);
	friend matrix operator / (const matrix& m1, const matrix& m2);
	friend matrix operator / (const matrix& m, const double& no);
	friend matrix operator / (const double& no, const matrix& m);
	friend matrix operator ~ (const matrix& m);
	friend matrix operator ! (matrix m);
	friend matrix operator ^ (const matrix& m, const size_t& pow);

//	Miscellaneous -methods
	void Null (const size_t& row, const size_t& col);
	void Null ();
	void Unit (const size_t& row);
	void Unit ();
	void SetSize (size_t row, size_t col);

//	Utility methods
	matrix Solve (const matrix& v) const;
	matrix Adj ();
	double Det ();
	double Norm ();
	double Cofact (size_t row, size_t col);
	double Cond ();

//	Type of matrices
	bool IsSquare ();
	bool IsSingular ();
	bool IsDiagonal ();
	bool IsScalar ();
	bool IsUnit ();
	bool IsNull ();
	bool IsSymmetric ();
	bool IsSkewSymmetric ();
	bool IsUpperTriangular ();
	bool IsLowerTriangular ();

//	Io-stream operators
	friend istream & operator >> (istream& i, matrix& m);
	friend ostream & operator << (ostream& o, matrix& m);
};