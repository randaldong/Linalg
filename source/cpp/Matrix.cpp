#include "../hpp/linalg.hpp"

Matrix::Matrix (int m, int n) {
    rows = m;
    cols = n;
    matrix = new double*[rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = 0;
        }
    }
}

Matrix::Matrix (int m, int n, double** arr) {
    rows = m;
    cols = n;
    matrix = new double*[rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = arr[i][j];
        }
    }
}

Matrix::Matrix (const Matrix& mat) {
    rows = mat.rows;
    cols = mat.cols;
    matrix = new double*[rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = mat.matrix[i][j];
        }
    }
}

Matrix::~Matrix () {delete matrix;}

int Matrix::GetRows () const {return rows;}

int Matrix::GetCols () const {return cols;}

Vector Matrix::GetRowVector (int k) const {
    Vector row(cols);
    for (int i = 0; i < cols; i++) {
        row.SetElement(i, matrix[k][i]);
    }
    return row;
}

Vector Matrix::GetColVector (int k) const {
    Vector col(rows);
    for (int i = 0; i < rows; i++) {
        col.SetElement(i, matrix[k][i]);
    }
    return col;
}

double Matrix::GetElement (int i, int j) const {return matrix[i][j];}

void Matrix::SetElement (int i, int j, double val) {matrix[i][j] = val;}

bool Matrix::IsContained (int i, int* indici, int dim) {
    bool result = false;
    for (int k = 0; k < dim; k++) {
        if (i == indici[k]) {
            result = true;
        }
    }
    return result;
}

Matrix Matrix::GetMatrix (int k) {
    int riga = 0, colonna;
    Matrix mat(rows-1, cols-1);
    for (int i = 0; i < rows; i++) {
        colonna = 0;
        if (i != k) {
            for (int j = 1; j < cols; j++) {
                    mat.matrix[riga][colonna] = matrix[i][j];
                    colonna++;
            }
            riga++;
        }
    }
    return mat;
}

Matrix Matrix::GetMatrix (int* indici, int dim) {
    Matrix mat(rows-dim, cols-dim);
    int riga = 0, colonna;
    for (int i = 0; i < rows; i++) {
        colonna = 0;
        if (!IsContained(i, indici, dim)) {
            for (int j = 0; j < cols; j++) {
                if (!IsContained(j, indici, dim)) {
                    mat.matrix[riga][colonna] = matrix[i][j];
                    colonna++;
                }
            }
            riga++;
        }
    }
    return mat;
}

Matrix Matrix::GetCombs (int dim, int m, int n) {
    Matrix result(m, n);
    int* vett = new int[n];
    for (int i = 0; i < n; i++) {
        vett[i] = i + 1;
    }
    int i = 0;
    while (i < m) {
        int j = 0;
        while (j < n) {
            result.matrix[i][j] = vett[j]; 
            j++;      
        }
        int z = n - 1;
        while (z >= 0) {
            int max = dim - n + z + 1;  
            if (vett[z] < max) {     
                vett[z]++;
                int k = z;
                while (k < n - 1) {
                    vett[k + 1] = vett[k] + 1;
                    k++;
                }
                break;
            }
            z--;
        }
        i++;
    }  
    return result;
}

Matrix Matrix::Merge (const Matrix& right) {
    Matrix composta(rows, cols*2);
    for (int i = 0; i < composta.rows; i++) {
        for (int j = 0; j < composta.cols; j++) {
            if (j < cols) {
                composta.matrix[i][j] = matrix[i][j];
            } else {
                composta.matrix[i][j] = right.matrix[i][j-composta.rows];
            }
        }
    }
    return composta;
}

Matrix Matrix::SwapRows () {
    Matrix swapped = *this;
    double *temp = new double[rows];
    for (int i = 0; i < rows-1; i++) {
		if (swapped.matrix[i][i] == 0) {
			for (int p = i+1; p < rows; p++) {
				if (swapped.matrix[p][i] != 0) {
					for (int k = 0; k < cols; k++) {
						temp[k] = swapped.matrix[i][k];
						swapped.matrix[i][k] = swapped.matrix[p][k];
						swapped.matrix[p][k] = temp[k];
					}
                    break;
				}
			}
            for (int p = 0; p < i; p++) {
				if (swapped.matrix[p][i] != 0) {
					for (int k = 0; k < cols; k++) {
						temp[k] = swapped.matrix[i][k];
						swapped.matrix[i][k] = swapped.matrix[p][k];
						swapped.matrix[p][k] = temp[k];
					}
                    break;
				}
			}
		}
    }
    return swapped;
}

Matrix Matrix::GetPermutation (Matrix& original) {
    Matrix P(rows,rows);
    for (int i = 0; i < rows; i++) {
        Vector row1 = GetRowVector(i);
        for (int p = 0; p < rows; p++) {
            Vector row2 = original.GetRowVector(p);
            if (row1 == row2) {
                P.matrix[i][p] = 1;
            }
        }
    }
    return P;
}

double Matrix::Norm (char p) {
    double result;
    // maximum singular value
    if (p == '2') {
        Matrix* USV = SVdecomposition();
        Matrix sigma = USV[1];
        result = sigma.matrix[0][0];
        for (int i = 1; i < sigma.rows; i++) {
            if (sigma.matrix[i][i] > result) {
                result = sigma.matrix[i][i];
            }
        }
    }
    // frobenius norm
    if (p == 'f') {
        Matrix X = (rows < cols) ? (*this)*(this->T()) : (this->T())*(*this);
        result = sqrt(X.Trace());
    }
    // nuclear norm
    if (p == 'n') {
        Matrix* USV = SVdecomposition();
        Matrix sigma = USV[1];
        result = 0;
        for (int i = 0; i < sigma.rows; i++) {
            result += sigma.matrix[i][i];
        }
    }
    return result;
}

Matrix Matrix::T () {
    Matrix trasp(cols, rows);
    for (int i = 0; i < trasp.rows; i++) {
        for (int j = 0; j < trasp.cols; j++) {
            trasp.matrix[i][j] = matrix[j][i];
        }
    }
    return trasp;
}

Matrix Matrix::I () {
    try {
        if (Determinant() == 0)
            throw "\033[1;31mMatrix Matrix::I () -->\n\tThe matrix cannot be inverted (determinant = 0)\033[0m\n";
    }
    catch (const char* err) {
        cout << "\n\033[1;31mEXCEPTION: \033[0m" << err;
        throw;
    }
    Matrix id = Eye(rows);
    Matrix composta = this->Merge(id);
    Matrix gauss = composta.Gauss();
    for (int i = gauss.rows-1; i > 0; i--) {
		if (gauss.matrix[i][i] != 0) {
			for (int q = i-1; q >= 0; q--) {
                double k = gauss.matrix[q][i] / gauss.matrix[i][i];
				for (int j = 0; j < gauss.cols; j++) {
                    gauss.matrix[q][j] -= gauss.matrix[i][j] * k;
				}
			}
		}
	}
    for (int i = 0; i < gauss.rows; i++) {
        double q = gauss.matrix[i][i];
		if (q != 0) {
			for (int j = 0; j < gauss.cols; j++) {
			    gauss.matrix[i][j] /= q;
			}
		}
	}
    Matrix inv(rows, cols);
    for (int i = 0; i < inv.rows; i++) {
        for (int j = 0; j < inv.cols; j++) {
            inv.matrix[i][j] = gauss.matrix[i][j+cols];
        }
    }
    return inv;
}

Matrix Matrix::Triu () {
    Matrix triu = *this;
    for (int i = 1; i < triu.rows; i++) {
        for (int j = 0; j < i; j++) {
            triu.matrix[i][j] = 0;
        }
    }
    return triu;
}

Matrix Matrix::Gauss () {
    Matrix gauss = this->SwapRows();
    for (int i = 0; i < gauss.rows-1; i++) {
        if (gauss.matrix[i][i] == 0) continue;
		for (int q = i+1; q < gauss.rows; q++) {
            double k = gauss.matrix[q][i] / gauss.matrix[i][i];
			for (int j = 0; j < gauss.cols; j++) {
				gauss.matrix[q][j] -= gauss.matrix[i][j] * k;
			}
		}
	}
    return gauss;
}

int Matrix::Rank () {
    Matrix gauss = Gauss();
    int rank = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (abs(gauss.matrix[i][j]) > ZERO) {
                rank++;
                break;
            }
        }
    }
    return rank;
}

double Matrix::Trace () {
    try {
        if (!IsSquare())
            throw "\033[1;31mdouble Matrix::Trace () -->\n\tThe trace cannot be computed because the matrix is not square\033[0m\n";
    }
    catch (const char* err) {
        cout << "\n\033[1;31mEXCEPTION: \033[0m" << err;
        throw;
    }
    double result = 0;
    for (int i = 0; i < rows; i++) {
        result += matrix[i][i];
    }
    return result;
}

double Matrix::Determinant (char method) {
    try {
        if (!IsSquare())
            throw "\033[1;31mdouble Matrix::Determinant () -->\n\tThe determinant cannot be computed because the matrix is not square\033[0m\n";
    }
    catch (const char* err) {
        cout << "\n\033[1;31mEXCEPTION: \033[0m" << err;
        throw;
    }
    double result;
    if (method == 'g') {
        result = 1;
        Matrix triang = Gauss();
        for (int i = 0; i < triang.rows; i++) {
            result *= triang.matrix[i][i];
        }
    }
    if (method == 'l') {
        if (rows == 1) {
            result = matrix[0][0];
        } else {
            result = 0;
            for (int k = 0; k < rows; k++) {
                Matrix temp = GetMatrix(k);
                result += pow(-1,k)*matrix[k][0]*temp.Determinant();
            }
        }
    }
    return result;
}

Matrix* Matrix::Eigendecomposition () {
    try {
        if (!IsSquare())
            throw "\033[1;31mMatrix* Matrix::Eigendecomposition () -->\n\tThe eigen-decomposition cannot be performed for non-square matrices\033[0m\n";
        if (!IsSymmetric())
            throw "\033[1;31mMatrix* Matrix::Eigendecomposition () -->\n\tThe eigen-decomposition is not implemented for non-symmetric matrices\033[0m\n";
    }
    catch (const char* err) {
        cout << "\n\033[1;31mEXCEPTION: \033[0m" << err;
        throw;
    }
    Matrix Q = Eye(this->cols);
    Matrix A = *this;
    Matrix I = Eye(A.rows);
    for (int iter = 0; iter < 1000; iter++) {
        double shift = WilkinsonShift(A);
        A = A - (I * shift);
        Matrix* QR = A.QRdecomposition();
        Q = Q * QR[0];
        A = (QR[1] * QR[0]) + (I * shift);
        if (A == A.Triu()) break;
    }
    Matrix* QL = new Matrix[2];
    QL[0] = Q; QL[1] = Diag(Diag(A));
    return QL;
}

Matrix* Matrix::SVdecomposition () {
    Matrix* QL;
    Matrix* USV = new Matrix[3];
    if ((*this) == this->T()) {
        QL = Eigendecomposition();
        USV[0] = QL[0];
        USV[1] = QL[1];
        USV[2] = QL[0];
    } else {
        Matrix symm_pos_semidef = (rows < cols) ? (*this)*(this->T()) : (this->T())*(*this);
        QL = symm_pos_semidef.Eigendecomposition();
        Matrix Q1 = QL[0], sigma = QL[1];
        for (int i = 0; i < sigma.rows; i++) {
            sigma.matrix[i][i] = sqrt(sigma.matrix[i][i]);
        }
        Matrix Q2 = (rows < cols) ? (this->T()) * Q1 * sigma.I() : (*this) * Q1 * sigma.I();
        USV[0] = (rows < cols) ? Q1 : Q2;
        USV[1] = sigma;
        USV[2] = (rows < cols) ? Q2 : Q1;
    }
    return USV;
}

Matrix* Matrix::QRdecomposition () {
    try {
        if (!IsSquare())
            throw "\033[1;31mMatrix* Matrix::QRdecomposition () -->\n\tThe QR decomposition cannot be performed for non-square matrices\033[0m\n";
    }
    catch (const char* err) {
        cout << "\n\033[1;31mEXCEPTION: \033[0m" << err;
        throw;
    }
    Matrix Qt = Eye(rows), R = *this, H;
    for (int i = 0; i < rows; i++) {
        H = R.HouseholderReflection(i);
        Qt = H * Qt;
        R = H * R;
    }
    Matrix* QR = new Matrix[2];
    QR[0] = Qt.T();
    QR[1] = R;
    return QR;
}

Matrix* Matrix::Choleskydecomposition (bool diag) {
    try {
        if (!IsSquare())
            throw "\033[1;31mMatrix* Matrix::Choleskydecomposition (bool) -->\n\tThe Cholesky decomposition cannot be performed for non-square matrices\033[0m\n";
        if (!IsSymmetric())
            throw "\033[1;31mMatrix* Matrix::Choleskydecomposition (bool) -->\n\tThe Cholesky decomposition cannot be performed for non-symmetric matrices\033[0m\n";
        if (!IsPositiveDefinite())
            throw "\033[1;31mMatrix* Matrix::Choleskydecomposition (bool) -->\n\tThe Cholesky decomposition cannot be performed for non-positive definite matrices\033[0m\n";
    }
    catch (const char* err) {
        cout << "\n\033[1;31mEXCEPTION: \033[0m" << err;
        throw;
    }
    Matrix C(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j <= i; j++) {
            int sum = 0;
            if (j == i) {
                for (int k = 0; k < j; k++)
                    sum += pow(C.matrix[j][k], 2);
                C.matrix[j][j] = sqrt(matrix[j][j] - sum);
            } else {
                for (int k = 0; k < j; k++)
                    sum += C.matrix[i][k] * C.matrix[j][k];
                C.matrix[i][j] = (matrix[i][j] - sum) / C.matrix[j][j];
            }
        }
    }
    Matrix* factors;
    if (!diag) {
        factors = new Matrix[2];
        factors[0] = C;
        factors[1] = C.T();
    } else {
        factors = new Matrix[3];
        Matrix S = Diag(Diag(C));
        factors[0] = C * S.I();
        factors[1] = S * S;
        factors[2] = factors[0].T();
    }
    return factors;
}

Matrix* Matrix::LUdecomposition (bool diag) {
    try {
        if (!IsSquare())
            throw "\033[1;31mMatrix* Matrix::LUdecomposition (bool) -->\n\tThe LU decomposition cannot be performed for non-square matrices\033[0m\n";
    }
    catch (const char* err) {
        cout << "\n\033[1;31mEXCEPTION: \033[0m" << err;
        throw;
    }
    if (IsSymmetric()) return Choleskydecomposition(diag);
    Matrix L = Eye(rows);
    Matrix U = this->SwapRows();
    Matrix P = U.GetPermutation(*this);
    for (int i = 0; i < U.rows-1; i++) {
        if (U.matrix[i][i] == 0) continue;
        for (int q = i+1; q < U.rows; q++) {
            double k = U.matrix[q][i] / U.matrix[i][i];
            L.matrix[q][i] = k;
            for (int j = 0; j < U.cols; j++) {
                U.matrix[q][j] -= U.matrix[i][j] * k;
            }
        }
    }
    Matrix* factors;
    if (!diag) {
        factors = new Matrix[3];
        factors[0] = L;
        factors[1] = U;
        factors[2] = P;
    } else {
        factors = new Matrix[4];
        Matrix D = Diag(Diag(U));
        factors[0] = L;
        factors[1] = D;
        factors[2] = D.I() * L.I() * (*this);
        factors[3] = P;
    }
    return factors;
}

Matrix Matrix::HouseholderReflection(int i) {
    Vector _x = GetColVector(i);
    Vector x(rows-i);
    for (int j = 0; j < rows-i; j++) {
        x.SetElement(j, _x.GetElement(j+i));
    }
    Vector xp = Eye(rows-i).GetColVector(0) * x.Norm();
    Vector u = x - xp;
    double dot = Dot(u,u);
    Matrix H = Eye(rows);
    if (abs(dot) > ZERO) {
        Matrix _H = Eye(rows-i) - ( Outer(u,u) * (2./dot) );
        for (int r = 0; r < _H.rows; r++) {
            for (int c = 0; c < _H.cols; c++) {
                H.matrix[r+i][c+i] = _H.matrix[r][c];
            }
        }
    }
    return H;
}

Polynomial Matrix::CharacteristicPol () {
    auto Fattoriale = [] (int n) {
        int res = 1;
        for (int i = 2; i <= n; i++) {
            res = res * i;
        }
        return res;
    };
    Polynomial car (rows);
    int iCoeff = car.GetDegree();
    car.SetCoefficient(iCoeff, pow(-1, iCoeff));
    car.SetCoefficient(iCoeff-1, pow(-1, iCoeff-1)*Trace());
    for (int k = iCoeff-2; k > 0; k--) {
        int segno = pow(-1 ,k);   
        int comb = Fattoriale(iCoeff)/(Fattoriale(iCoeff-k)*Fattoriale(k));    
        Matrix temp1 = GetCombs(rows, comb, k);
        double minore = 0;
        for (int t = 0; t < comb; t++) {
            int* indici = new int[k];
            for (int z = 0; z < k; z++) {
                indici[z] = (int)temp1.matrix[t][z]-1;
            }
            Matrix temp2 = GetMatrix(indici, k);
            minore += temp2.Determinant();
        }
        car.SetCoefficient(k, segno*minore);
    }
    car.SetCoefficient(0, Determinant());
    return car;
}

bool Matrix::IsSquare () {
    return (rows == cols) ? true : false;
}

bool Matrix::IsSymmetric () {
    if (!IsSquare()) return false;
    return (*this == this->T()) ? true : false;
}

bool Matrix::IsPositiveDefinite () {
    if (!IsSymmetric()) return false;
    Matrix gauss = Gauss();
    for (int i = 0; i < gauss.rows; i++) {
        if (gauss.matrix[i][i] < 0 || abs(gauss.matrix[i][i]) < ZERO) return false;
    }
    return true;
}

ostream& operator << (ostream& os, Matrix& mat) {
    for (int i = 0; i < mat.rows; i++) {
        os << "[";
        for (int j = 0; j < mat.cols; j++) {
            double entry = abs(mat.matrix[i][j]) >= 0.005 ? mat.matrix[i][j] : 0;
            os << setprecision(2) << setw(8) << entry << ' ';
        }
        os << "]\n";
    }
    return os;
}

istream& operator >> (istream& is, Matrix& mat) {
    if (cin) {
        clog << "\nDigita il numero di righe: ";
        is >> mat.rows;
    }
    if (cin) {
        clog << "\nDigita il numero di colonne: ";
        is >> mat.cols;
    }
    clog << "\n";
    mat.matrix = new double*[mat.rows];
    for (int i = 0; i < mat.rows; i++) {
        mat.matrix[i] = new double[mat.cols];
    }
    if (cin) {
        for (int i = 0; i < mat.rows; i++) {
            for (int j = 0; j < mat.cols; j++) {
                clog << "Mat[" << i+1 << "][" << j+1 << "] = ";
                is >> mat.matrix[i][j];
            }
            clog << "\n";
        }
    }
    return is;
}

Matrix& Matrix::operator = (const Matrix& mat) {
    rows = mat.rows;
    cols = mat.cols;
    matrix = new double*[rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = mat.matrix[i][j];
        }
    }
    return *this;
}

Matrix Matrix::operator + (const Matrix& mat) {
    try {
        if (rows != mat.rows || cols != mat.cols)
            throw "\033[1;31mMatrix Matrix::operator + (Matrix) -->\n\tThe matrices have different dimensions\033[0m\n";
    }
    catch (const char* err) {
        cout << "\n\033[1;31mEXCEPTION: \033[0m" << err;
        throw;
    }
    Matrix somma(rows, cols);
    for (int i = 0; i < somma.rows; i++) {
        for (int j = 0; j < somma.cols; j++) {
            somma.matrix[i][j] = matrix[i][j] + mat.matrix[i][j];
        }
    }
    return somma;
}

Matrix Matrix::operator - (const Matrix& mat) {
    try {
        if (rows != mat.rows || cols != mat.cols)
            throw "\033[1;31mMatrix Matrix::operator - (Matrix) -->\n\tThe matrices have different dimensions\033[0m\n";
    }
    catch (const char* err) {
        cout << "\n\033[1;31mEXCEPTION: \033[0m" << err;
        throw;
    }
    Matrix diff(rows, cols);
    for (int i = 0; i < diff.rows; i++) {
        for (int j = 0; j < diff.cols; j++) {
            diff.matrix[i][j] = matrix[i][j] - mat.matrix[i][j];
        }
    }
    return diff;
}

Matrix Matrix::operator * (const Matrix& mat) {
    try {
        if (cols != mat.rows)
            throw "\033[1;31mMatrix Matrix::operator * (Matrix) -->\n\tThe matrices have incompatible dimensions for multiplication\033[0m\n";
    }
    catch (const char* err) {
        cout << "\n\033[1;31mEXCEPTION: \033[0m" << err;
        throw;
    }
    Matrix prod(rows, mat.cols);
    for (int i = 0; i < rows; i++) {
        for (int k = 0; k < cols; k++) {
            for (int j = 0; j < mat.cols; j++) {
            prod.matrix[i][j] += matrix[i][k] * mat.matrix[k][j];
            }
        }
    }
    return prod;
}

Matrix Matrix::operator * (double k) {
    Matrix res(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            res.matrix[i][j] = matrix[i][j] * k;
        }
    }
    return res;
}

Matrix Matrix::operator / (double k) {
    try {
        if (k == 0)
            throw "\033[1;31mMatrix Matrix::operator / (double) -->\n\tDivision by zero\033[0m\n";
    }
    catch (const char* err) {
        cout << "\n\033[1;31mEXCEPTION: \033[0m" << err;
        throw;
    }
    Matrix res(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            res.matrix[i][j] = matrix[i][j] / k;
        }
    }
    return res;
}

bool Matrix::operator == (const Matrix& mat) {
    if (rows != mat.rows) return false;
    if (cols != mat.cols) return false;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (abs(matrix[i][j] - mat.matrix[i][j]) > ZERO) return false;
        }
    }
    return true;
}

bool Matrix::operator != (const Matrix& mat) {
    return !(*this == mat);
}