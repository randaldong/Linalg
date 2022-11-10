#include "../hpp/linalg.hpp"

Vector::Vector (int s) {
    size = s;
    array = new double[size];
    for (int i = 0; i < size; i++) {
        array[i] = 0;
    }
}

Vector::Vector (int s, double* arr) {
    size = s;
    array = new double[size];
    for (int i = 0; i < size; i++) {
        array[i] = arr[i];
    }
}

Vector::Vector (const Vector& vec) {
    size = vec.size;
    array = new double[size];
    for (int i = 0; i < size; i++) {
        array[i] = vec.array[i];
    }
}

Vector::~Vector () {delete array;}

int Vector::GetSize () const {return size;}

double* Vector::GetArray () const {return array;}

double Vector::GetElement (int i) const {return array[i];}

void Vector::SetElement (int i, double val) {array[i] = val;}

double Vector::Norm (int p) {
    double result = 0;
    for (int i = 0; i < size; i++) {
        result += pow(abs(array[i]), p);
    }
    result = pow(result, 1/(double)p);
    return result;
}

Vector Vector::Normalized (int p) {
    Vector vec = *this / Norm(p);
    return vec;
}

Vector Vector::ProjectedOnto (Vector vec) {
    double k = Dot(*this,vec) / Dot(vec,vec);
    Vector proj = vec*k;
    return proj;
}

ostream& operator << (ostream& os, Vector& vec) {
    os << "[";
    for (int i = 0; i < vec.size; i++) {
        double entry = abs(vec.array[i]) >= 0.005 ? vec.array[i] : 0;
        os << setprecision(2) << entry;
        if (i != vec.size-1) {os << "  ";}
    }
    os << "]";
    return os;
}

istream& operator >> (istream& is, Vector& vec) {
    vec.array = new double[vec.size];
    if (cin) {
        clog << "Enter the " << vec.size << " components of the vector:\n";
        for (int k = 0; k < vec.size; k++) { 
            is >> vec.array[k];
        }
    }
    return is;
}

Vector& Vector::operator = (const Vector& vec) {
    size = vec.size;
    array = new double[size];
    for (int i = 0; i < size; i++) {
        array[i] = vec.array[i];
    }
    return *this;
}

Vector Vector::operator + (const Vector& vec) {
    try {
        if (size != vec.size)
            throw "\033[1;31mVector Vector::operator + (Vector) -->\n\tThe vectors have different sizes\033[0m\n";
    }
    catch (const char* err) {
        cout << "\n\033[1;31mEXCEPTION: \033[0m" << err;
        throw;
    }
    Vector sum(size);
    for (int i = 0; i < size; i++) {
        sum.array[i] = array[i] + vec.array[i];
    }
    return sum;
}

Vector Vector::operator - (const Vector& vec) {
    try {
        if (size != vec.size)
            throw "\033[1;31mVector Vector::operator - (Vector) -->\n\tThe vectors have different sizes\033[0m\n";
    }
    catch (const char* err) {
        cout << "\n\033[1;31mEXCEPTION: \033[0m" << err;
        throw;
    }
    Vector diff(size);
    for (int i = 0; i < size; i++) {
        diff.array[i] = array[i] - vec.array[i];
    }
    return diff;
}

double Vector::operator * (const Vector& vec) {
    try {
        if (size != vec.size)
            throw "\033[1;31mVector Vector::operator * (Vector) -->\n\tThe vectors have different sizes\033[0m\n";
    }
    catch (const char* err) {
        cout << "\n\033[1;31mEXCEPTION: \033[0m" << err;
        throw;
    }
    double result = 0;
    for (int i = 0; i < size; i++) {
        result += array[i] * vec.array[i];
    }
    return result;
}

Vector Vector::operator * (double k) {
    Vector res(size);
    for (int i = 0; i < size; i++) {
        res.array[i] = array[i] * k;
    }
    return res;
}

Vector Vector::operator / (double k) {
    try {
        if (k == 0)
            throw "\033[1;31mVector Vector::operator / (double) -->\n\tDivision by zero\033[0m\n";
    }
    catch (const char* err) {
        cout << "\n\033[1;31mEXCEPTION: \033[0m" << err;
        throw;
    }
    Vector res(size);
    for (int i = 0; i < size; i++) {
        res.array[i] = array[i] / k;
    }
    return res;
}

bool Vector::operator == (const Vector& vec) {
    if (size != vec.size) {return false;}
    for (int i = 0; i < size; i++) {
        if (abs(array[i] - vec.array[i]) > ZERO) return false;
    }
    return true;
}

bool Vector::operator != (const Vector& vec) {
    return !(*this == vec);
}