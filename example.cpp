#include "source/hpp/linalg.hpp"

int main () {

    srand((unsigned int)time(NULL));


    cout << "\n############################################################\n";
    cout << "SINGULAR VALUES DECOMPOSITION OF A RANDOM RECTANGULAR MATRIX\n";
    cout << "############################################################\n";
    int shape[2] = {4,6};
    Matrix A = Rand(shape, -2, 3);
    cout << "Matrix A =\n" << A;

    Matrix* USV = A.SVdecomposition();
    cout << "\nU =\n" << USV[0];
    cout << "\nSigma =\n" << USV[1];
    cout << "\nV =\n" << USV[2];
    Matrix A_again = USV[0] * USV[1] * USV[2].T();
    cout << "\nA = U * Sigma * V.T =\n" << A_again;


    cout << "\n\n\n###############################################\n";
    cout << "EIGENDECOMPOSITION OF A RANDOM SYMMETRIC MATRIX\n";
    cout << "###############################################\n";
    shape[0] = 5; shape[1] = 5;
    Matrix B = Rand(shape, -5, 5);
    Matrix S = B.Triu() + (B.Triu()).T() - Diag(Diag(B));
    cout << "Matrix S =\n" << S;

    double norm = S.Norm('f');
    cout << "Frobenius Norm = " << norm;

    double tr = S.Trace();
    cout << "\nTrace = " << tr;

    double det = S.Determinant();
    cout << "\nDeterminant = " << det << "\n";

    Polynomial polchar = S.CharacteristicPol();
    cout << "Characteristic polynomial = " << polchar << "\n";

    Matrix* QL = S.Eigendecomposition();
    cout << "\nQ =\n" << QL[0];
    cout << "\nLambda =\n" << QL[1];
    Matrix S_again = QL[0] * QL[1] * QL[0].T();
    cout << "\nS = Q * Lambda * Q.T =\n" << S_again;


    cout << "\n\n\n#############################################################\n";
    cout << "LU DECOMPOSITION (WITH PERMUTATION) OF A RANDOM SQUARE MATRIX\n";
    cout << "#############################################################\n";
    shape[0] = 4; shape[1] = 4;
    Matrix C = Rand(shape, -3, 5);
    C.SetElement(0,0, 0);
    C.SetElement(3,3, 0);
    cout << "Matrix C =\n" << C;

    Matrix* dec = C.LUdecomposition();
    cout << "\nL =\n" << dec[0];
    cout << "\nU =\n" << dec[1];
    cout << "\nP =\n" << dec[2];
    Matrix C_again = dec[2].T() * dec[0] * dec[1];
    cout << "\nC = P.T * L * U =\n" << C_again;

    return 0;
}