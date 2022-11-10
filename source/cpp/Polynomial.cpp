#include "../hpp/linalg.hpp"

Polynomial::Polynomial (int k) {
    degree = k;
    coefficients = new double[degree+1];
    for (int i = 0; i <= degree; i++) {
        coefficients[i] = 0;
    }
}

Polynomial::Polynomial (int k, double *coeff) {
    degree = k;
    coefficients = new double[degree+1];
    for (int i = 0; i <= degree; i++) {
        coefficients[i] = coeff[i];
    }
}

Polynomial::Polynomial (const Polynomial& pol) {
    degree = pol.degree;
    coefficients = new double[degree+1];
    for (int i = 0; i <= degree; i++) {
        coefficients[i] = pol.coefficients[i];
    }
}

Polynomial::~Polynomial () {delete coefficients;}

int Polynomial::GetDegree () const {return degree;}

double* Polynomial::GetCoefficients () const {return coefficients;}

void Polynomial::SetCoefficient (int indice, double valore) {
    coefficients[indice] = valore;
}

double Polynomial::Evaluate (double x) {
    double valore = 0;  
    for (int i = degree; i >= 0; i--) {
        valore += pow(x, i)*coefficients[i];
    }
    return valore;
}

double Polynomial::EvaluateDerivative (double x) {
    double valore = 0;
    for (int i = degree; i > 0; i--) {
        valore += i*pow(x, i-1)*coefficients[i];
    }
    return valore;
}

double Polynomial::ComputeZero () {
    double rad = static_cast <double> (rand());
    int i = 0;
    int MAXITER = 1000;
    try {
        while (abs(Evaluate(rad)) > ZERO) {
            if (i > MAXITER)
                throw "\033[1;31mdouble Polynomial::ComputeZero () -->\n\tThe polynomial has at least one complex zero solution\033[0m\n";
            double val = EvaluateDerivative(rad);
            if (val == 0)
                throw "\033[1;31mdouble Polynomial::ComputeZero () -->\n\tDivision by zero\033[0m\n";
            rad -= Evaluate(rad) / val;
            i++;
        }
    }
    catch (const char* err) {
        cout << "\n\033[1;31mEXCEPTION: \033[0m" << err;
        throw;
    }
    return rad;
}

double* Polynomial::ComputeZeros () {
    double* radici = new double[degree];
    Polynomial div, temp;
    double rad;
    if (coefficients[0] == 0) {
        rad = 0;
    } else {      
        rad = ComputeZero();
    }
    radici[0] = rad;
    int curIndex = 1;
    if (degree > 1) {
        Polynomial div(1);
        div.SetCoefficient(1, 1);
        div.SetCoefficient(0, -rad);
        temp = *this / div;
        double* tempRadici = temp.ComputeZeros();
        for (int i = 0; i < temp.GetDegree(); i++) {
            radici[curIndex++] = tempRadici[i];
        }
    }
    return radici;
}

ostream& operator << (ostream& os, Polynomial& pol) {
    for (int k = pol.degree; k > 0; k--) {   
        double abs_coeff = abs(pol.coefficients[k]) > 0.005 ? abs(pol.coefficients[k]) : 0;
        if (abs_coeff != 0) {
            os << (pol.coefficients[k] > 0 ? '+' : '-');
            if (abs_coeff != 1) {
                os << abs_coeff;
            }
            if (k != 1) {
                os << "x^" << k;
            } else {
                os << "x";
            }
            os << ' ';
        }
    }
    if (pol.coefficients[0] != 0) {
        os << (pol.coefficients[0] > 0 ? '+' : '-') << abs(pol.coefficients[0]);
    }
    return os;
}

istream& operator >> (istream& is, Polynomial& pol) {
    pol.coefficients = new double[pol.degree+1];
    if (cin) {
        clog << "Enter the " << pol.degree+1 << " coefficients of the polynomial:\n";
        for (int k = pol.degree; k >= 0; k--) {
            is >> pol.coefficients[k];
        }
    }
    return is;
}

Polynomial& Polynomial::operator = (const Polynomial& pol) {
    degree = pol.degree;
    coefficients = new double[degree+1];
    for (int i = 0; i <= degree; i++) {
        coefficients[i] = pol.coefficients[i];
    }
    return *this;
}

Polynomial Polynomial::operator + (const Polynomial& pol) {
    Polynomial max, min;
    if (degree > pol.degree) {
        max = *this;
        min = pol;
    } else {
        max = pol;
        min = *this;  
    }
    Polynomial somma(max.degree);
    for (int i = 0; i <= somma.degree; i++) { 
        somma.coefficients[i] = max.coefficients[i] + (i <= min.degree ? min.coefficients[i] : 0);
    }       
    return somma;
}

Polynomial Polynomial::operator - (const Polynomial& pol) {
    Polynomial diff(pol.degree <= degree ? degree : pol.degree);
    for (int i = 0; i <= diff.degree; i++) { 
        diff.coefficients[i] = coefficients[i] - pol.coefficients[i];
    }       
    return diff;
}

Polynomial Polynomial::operator * (const Polynomial& pol) {
    Vector v1(degree+1, coefficients);
    Vector v2(pol.degree+1, pol.coefficients);
    Vector conv = Convolution(v1, v2);
    Polynomial prod(degree+pol.degree, conv.GetArray());
    return prod;
}

Polynomial Polynomial::operator / (const Polynomial& pol) {
    Polynomial ret(degree-pol.degree);
    ret.coefficients[ret.degree] = coefficients[degree] / pol.coefficients[pol.degree];
    for (int i = 1; i <= ret.degree; i++) {
        double* u = new double[degree-i+2];
        for (int k = 0; k <= degree-i+1; k++) {
            u[k] = coefficients[k];
        }
        Polynomial temp1(degree-i);
        Polynomial temp2(degree-i+1, u);
        temp1.coefficients[degree-i] = ret.coefficients[ret.degree-i+1];
        temp1 = temp1 * pol,
        temp2 = temp2 - temp1;
        ret.coefficients[ret.degree-i] = temp2.coefficients[degree-i] / pol.coefficients[pol.degree];
    }
    return ret;
}