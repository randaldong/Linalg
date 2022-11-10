#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

const double ZERO = 10e-9;

class Polynomial {
    int degree;
    double* coefficients;

    public:
    Polynomial () {};
    Polynomial (int degree);
    Polynomial (int degree, double* coefficients);
    Polynomial (const Polynomial& original);
    ~Polynomial ();
    int GetDegree () const;
    double* GetCoefficients () const;
    void SetCoefficient (int index, double value);
    double Evaluate (double point);
    double EvaluateDerivative (double point);
    double ComputeZero ();
    double* ComputeZeros ();

    friend ostream& operator << (ostream&, Polynomial&);
    friend istream& operator >> (istream&, Polynomial&);
    Polynomial& operator = (const Polynomial&);
    Polynomial operator + (const Polynomial&);
    Polynomial operator - (const Polynomial&);
    Polynomial operator * (const Polynomial&);
    Polynomial operator / (const Polynomial&);
};

class Vector {
    int size;
    double* array;

    public:
    Vector () {};
    Vector (int size);
    Vector (int size, double* array);
    Vector (const Vector& original);
    ~Vector ();
    int GetSize () const;
    double* GetArray () const;
    double GetElement (int index) const;
    void SetElement (int index, double value);
    double Norm (int p=2);
    Vector Normalized (int p=2);
    Vector ProjectedOnto (Vector axis);

    friend ostream& operator << (ostream&, Vector&);
    friend istream& operator >> (istream&, Vector&);
    Vector& operator = (const Vector&);
    Vector operator + (const Vector&);
    Vector operator - (const Vector&);
    double operator * (const Vector&); //dot product
    Vector operator * (double scalar);
    Vector operator / (double scalar);
    bool operator == (const Vector&);
    bool operator != (const Vector&);
};

class Matrix {
    int rows, cols;
    double** matrix;

    bool IsContained (int, int*, int);
    Matrix GetMatrix (int);
    Matrix GetMatrix (int*, int);
    Matrix GetCombs (int, int, int);
    Matrix Merge (const Matrix&);
	Matrix SwapRows ();
    Matrix GetPermutation (Matrix&);

    public:
    Matrix () {};
    Matrix (int rows, int cols);
    Matrix (int rows, int cols, double** matrix);
    Matrix (const Matrix& original);
    ~Matrix ();
    int GetRows () const;
    int GetCols () const;
    Vector GetRowVector (int index) const;
    Vector GetColVector (int index) const;
    double GetElement (int row_index, int col_index) const;
    void SetElement (int row_index, int col_index, double value);
    double Norm (char p='f');
    Matrix T ();
    Matrix I ();
    Matrix Triu ();
    Matrix Gauss ();
    int Rank ();
    double Trace ();
    double Determinant (char method='g');
    Matrix* Eigendecomposition ();
    Matrix* SVdecomposition ();
    Matrix* QRdecomposition ();
    Matrix* Choleskydecomposition (bool return_diag=false);
    Matrix* LUdecomposition (bool return_diag=false);
    Matrix HouseholderReflection (int index);
    Polynomial CharacteristicPol ();
    bool IsSquare ();
    bool IsSymmetric ();
    bool IsPositiveDefinite ();

    friend ostream& operator << (ostream&, Matrix&);
    friend istream& operator >> (istream&, Matrix&);
    Matrix& operator = (const Matrix&);
    Matrix operator + (const Matrix&);
    Matrix operator - (const Matrix&);
    Matrix operator * (const Matrix&); //matrix multiplication
    Matrix operator * (double scalar);
    Matrix operator / (double scalar);
    bool operator == (const Matrix&);
    bool operator != (const Matrix&);
};

double random_double (double min=0, double max=1);
int random_int (int min=0, int max=1);

Matrix Zeros (int* shape);
Matrix Ones (int* shape);
Matrix Rand (int* shape, double low=0, double high=1);
Matrix RandInt (int* shape, int low=0, int high=1);
Matrix Eye (int dim);
Vector Diag (Matrix matrix);
Matrix Diag (Vector diagonal);
double Dot (Vector vec1, Vector vec2); //dot product (see Vector::*)
Vector Dot (Matrix mat, Vector vec); //dot product
Vector Dot (Vector vec, Matrix mat); //dot product
Matrix Dot (Matrix mat1, Matrix mat2); //matrix multiplication (see Matrix::*)
Matrix Outer (Vector vec1, Vector vec2); //outer product
Vector Hadamard (Vector vec1, Vector vec2); //element-wise product
Vector* GramSchmidt (Vector* vectors, int vectors_number, int p=2); //orthonormalization
double WilkinsonShift (Matrix mat);
Vector Convolution (Vector vec1, Vector vec2); //1-dimensional
Matrix Convolution (Matrix mat1, Matrix mat2); //2-dimensional