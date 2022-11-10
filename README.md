# Linalg
> C++ Implementation of A Linear Algebra Library

## Table of contents

- [LINear ALGebra library](#linear-algebra-library)
  - [Overview](#overview)
  - [Usage](#usage)
  - [References](#references)


## Overview

This library contains a basic implementation of three different classes on the real numbers domain: Polynomial, Vector and Matrix. Here, a brief overview of the functionalities for each kind of object is provided.
To have a more precise idea of all the available methods just have a look to the header file *linalg.hpp*.

* Polynomial
  * basic algebraic operations (addition, subtraction, multiplication and division)
  * evaluation in a specific point
  * evaluation of the derivative in a specific point
  * computation of the zeros (if all real)

* Vector
  * basic algebraic operations (addition, subtraction)
  * multiplication by a scalar
  * Hadamard product (element-wise)
  * dot product
  * outer product
  * computation of the Lp norm (default: L2)
  * projection onto another Vector

* Matrix
  * basic algebraic operations (addition, subtraction)
  * multiplication by a scalar
  * matrix multiplication
  * computation of the norm (L2, nuclear, default: Frobenius)
  * transposition
  * gaussian elimination
  * inversion
  * determination of the rank
  * computation of the trace
  * computation of the determinant by different methods (Laplace, default: Gauss)
  * eigendecomposition (only for symmetric matrices)
  * singular values decomposition (SVD) for rectangular matrices
  * QR decomposition
  * Cholesky decomposition for positive definite matrices
  * LU decomposition
  * computation of the characteristic Polynomial

* Others
  * dot product between Matrix and Vector
  * Gram-Schmidt orthonormalization
  * 1D convolution between vectors
  * 2D convolution between matrices


## Usage

First of all, you have to clone the repository and change directory:

```bash
git clone https://github.com/ryandoren/Linalg.git
cd linalg
```

Then, you can easily build the `linalg` library (on Linux system only) using *CMake*:

```bash
cmake .
make
```

Finally, you can run as an example the executable *example.out* simply typing:
```
./example.out
```
