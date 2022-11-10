#include "../hpp/linalg.hpp"

double random_double (double min, double max) {
    return min + ((double)rand()/RAND_MAX) * (max - min);
}

int random_int (int min, int max) {
    return (int)(random_double(min, max) + 0.5);
}

Matrix Zeros (int* shape) {
    int r = shape[0], c = shape[1];
    Matrix zeros (r, c);
    return zeros;
}

Matrix Ones (int* shape) {
    int r = shape[0], c = shape[1];
    double** mat = new double*[r];
    for (int i = 0; i < r; i++) {
        mat[i] = new double[c];
        for (int j = 0; j < c; j++) {
            mat[i][j] = 1;
        }
    }
    Matrix ones (r, c, mat);
    return ones;
}

Matrix Rand (int* shape, double low, double high) {
    int r = shape[0], c = shape[1];
    double** mat = new double*[r];
    for (int i = 0; i < r; i++) {
        mat[i] = new double[c];
        for (int j = 0; j < c; j++) {
            mat[i][j] = random_double(low, high);
        }
    }
    Matrix rand (r, c, mat);
    return rand;
}

Matrix RandInt (int* shape, int low, int high) {
    int r = shape[0], c = shape[1];
    double** mat = new double*[r];
    for (int i = 0; i < r; i++) {
        mat[i] = new double[c];
        for (int j = 0; j < c; j++) {
            mat[i][j] = (double)random_int(low, high);
        }
    }
    Matrix randint (r, c, mat);
    return randint;
}

Matrix Eye (int n) {
    Matrix id(n,n);
    for (int i = 0; i < n; i++) {
        id.SetElement(i,i,1);
    }
    return id;
}

Vector Diag (Matrix mat) {
    int size = (mat.GetRows() < mat.GetCols()) ? mat.GetRows() : mat.GetCols();
    Vector diag(size);
    for (int i = 0; i < size; i++) {
        diag.SetElement(i, mat.GetElement(i,i));
    }
    return diag;
}

Matrix Diag (Vector vec) {
    int size = vec.GetSize();
    Matrix diag(size, size);
    for (int i = 0; i < size; i++) {
        diag.SetElement(i,i, vec.GetElement(i));
    }
    return diag;
}

double Dot (Vector v1, Vector v2) {
    double result = v1 * v2;
    return result;
}

Vector Dot (Matrix mat, Vector vec) {
    int rows = mat.GetRows();
    int cols = mat.GetCols();
    int size = vec.GetSize();
    try {
        if (cols != size)
            throw "\033[1;31mVector Dot (Matrix, Vector) -->\n\tThe number of matrix columns and the size of the vector are different\033[0m\n";
    }
    catch (const char* err) {
        cout << "\n\033[1;31mEXCEPTION: \033[0m" << err;
        throw;
    }
    Vector result(rows);
    for (int i = 0; i < rows; i++) {
        double el = 0;
        for (int j = 0; j < cols; j++) {
            el += mat.GetElement(i,j) * vec.GetElement(j);
        }
        result.SetElement(i, el);
    }
    return result;
}

Vector Dot (Vector vec, Matrix mat) {
    int size = vec.GetSize();
    int rows = mat.GetRows();
    int cols = mat.GetCols();
    try {
        if (size != rows)
            throw "\033[1;31mVector Dot (Vector, Matrix) -->\n\tThe size of the vector and the number of matrix rows are different\033[0m\n";
    }
    catch (const char* err) {
        cout << "\n\033[1;31mEXCEPTION: \033[0m" << err;
        throw;
    }
    Vector result(cols);
    for (int i = 0; i < cols; i++) {
        double el = 0;
        for (int j = 0; j < rows; j++) {
            el += vec.GetElement(j) * mat.GetElement(i,j);
        }
        result.SetElement(i, el);
    }
    return result;
}

Matrix Dot (Matrix m1, Matrix m2) {
    Matrix result = m1 * m2;
    return result;
}

Matrix Outer (Vector v1, Vector v2) {
    int s1 = v1.GetSize();
    int s2 = v2.GetSize();
    Matrix result(s1, s2);
    for (int i = 0; i < s1; i++) {
        for (int j = 0; j < s2; j++) {
            result.SetElement(i,j, v1.GetElement(i)*v2.GetElement(j));
        }
    }
    return result;
}

Vector Hadamard (Vector v1, Vector v2) {
    int s1 = v1.GetSize();
    int s2 = v2.GetSize();
    try {
        if (s1 != s2)
            throw "\033[1;31mVector Hadamard (Vector, Vector) -->\n\tThe sizes of the two vectors are different\033[0m\n";
    }
    catch (const char* err) {
        cout << "\n\033[1;31mEXCEPTION: \033[0m" << err;
        throw;
    }
    Vector result(s1);
    for (int i = 0; i < s1; i++) {
        result.SetElement(i, v1.GetElement(i)*v2.GetElement(i));
    }
    return result;
}

Vector* GramSchmidt (Vector* vecs, int num, int p) {
    Vector* orthogonal = new Vector[num];
    for (int i = 0; i < num; i++) {
        orthogonal[i] = vecs[i];
    }
    for (int i = 1; i < num; i++) {
        for (int k = 0; k < i; k++) {
            Vector proj = vecs[i].ProjectedOnto(orthogonal[k]);
            orthogonal[i] = orthogonal[i] - proj;
        }
    }
    Vector* orthonormal = new Vector[num];
    for (int i = 0; i < num; i++) {
        orthonormal[i] = orthogonal[i].Normalized(p);
    }
    return orthonormal;
}

double WilkinsonShift (Matrix mat) {
    double a = mat.GetElement(mat.GetRows()-2, mat.GetCols()-2);
    double b = mat.GetElement(mat.GetRows()-2, mat.GetCols()-1);
    double c = mat.GetElement(mat.GetRows()-1, mat.GetCols()-1);
    double delta = (a - c) * 0.5;
    double sign = 1 ? delta > 0 : -1;
    double denominator = abs(delta) + sqrt( pow(delta,2)+pow(b,2) );
    double shift = c - ( (sign*pow(b,2)) / denominator );
    return shift;
}

Vector Convolution (Vector v1, Vector v2) {
    int s1 = v1.GetSize();
    int s2 = v2.GetSize();
    int size = s1+s2-1;
    double* arr = new double[size];
    for (int i = 0; i < s1; i++) {
        for(int j = 0; j < s2; j++) {
            arr[i+j] += v1.GetElement(i) * v2.GetElement(j);
        }
    }
    Vector result(size, arr);
    return result;
}

Matrix Convolution (Matrix m1, Matrix m2) {
    int r1 = m1.GetRows(), c1 = m1.GetCols();
    int r2 = m2.GetRows(), c2 = m2.GetCols();
    int rows = r1+r2-1, cols = c1+c2-1;
    double** arr = new double*[rows];
    for (int i = 0; i < rows; i++) {
        arr[i] = new double[cols];
    }
    for (int i1 = 0; i1 < r1; i1++) {
        for (int j1 = 0; j1 < c1; j1++) {
            for (int i2 = 0; i2 < r2; i2++) {
                for (int j2 = 0; j2 < c2; j2++) { 
                    arr[i1+i2][j1+j2] += m1.GetElement(i1,j1) * m2.GetElement(i2,j2);
                }
            }
        }
    }
    Matrix result(rows, cols, arr);
    return result;
}