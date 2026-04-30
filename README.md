# Linear ODE System Solver

A self-contained Python script for solving linear systems of ordinary differential equations of the form **X' = AX + F(t)**. No external libraries required — all arithmetic uses exact rational fractions internally.

```bash
python3 ode_solver.py
```

> Supports 2×2 and 3×3 systems. Handles real, repeated, and complex eigenvalues. Works on any Python 3.7+ installation with no pip installs needed.

## Table of Contents

- [Requirements](#requirements)
- [Main Menu](#main-menu)
- [Functions](#functions)
- [How to Enter Numbers](#how-to-enter-numbers)
- [Mathematical Background](#mathematical-background)
- [Example Problems](#example-problems)
- [Limitations](#limitations)

## Requirements

Python 3.7 or higher. No third-party packages.

## Main Menu

When you run the script you will see:

```
1. Scalar equations  →  Matrix form
2. Matrix form       →  Scalar equations
3. Eigenvalues & Eigenvectors
4. Solve Homogeneous System   X' = AX
5. Solve Initial Value Problem
6. Solve Nonhomogeneous System
7. Verify a proposed solution
8. Fundamental Set check
9. Example Library
0. Quit
```

## Functions

### 1 — Scalar Equations → Matrix Form

**What it does:** Takes equations written out as `dx/dt = ...` and rewrites the system as the matrix equation X' = AX + F(t), identifying the matrix A and forcing vector F(t).

**When to use it:** You have a system written in scalar form and want to put it in matrix form before solving.

**How to input:**

1. Choose option `1`
2. Enter dimension (2 or 3)
3. Type the right-hand side of each equation — variables are `x`, `y`, `z`

Supported terms:

| Term | Example |
|------|---------|
| Variable with coefficient | `3x`, `-2y`, `5z` |
| Constant | `7`, `-4` |
| Linear in t | `2t`, `-t` |
| Quadratic in t | `t^2`, `3t^2` |
| Exponential | `e^t`, `2e^t`, `e^(3t)` |

**Example:**

```
dx/dt = 8x - 4y
dy/dt = 2x + 3y

→  A = [ 8  -4 ]    F(t) = [0]
       [ 2   3 ]            [0]
```

### 2 — Matrix Form → Scalar Equations

**What it does:** The reverse of option 1. Given matrix A and an optional forcing vector, writes out the individual scalar equations.

**How to input:**

1. Choose option `2`
2. Enter dimension and matrix A
3. Choose forcing type: none, constant vector, or b·e^(kt)

**Example:**

```
A = [ 6   4 ]    b·e^t = [-5]·e^t
    [-1   2 ]             [ 6]

→  dx/dt = 6x + 4y - 5e^t
   dy/dt = -x + 2y + 6e^t
```

### 3 — Eigenvalues & Eigenvectors

**What it does:** Computes the characteristic polynomial, eigenvalues λ, and eigenvectors K for a 2×2 or 3×3 matrix.

**How to input:**

1. Choose option `3`
2. Enter dimension and matrix A

**Example:**

```
A = [ 3  -2 ]
    [ 4  -1 ]

Characteristic polynomial:  λ² - 2λ + 5 = 0
λ₁ = 1 + 2i,   K₁ = [1/2 + 1/2i,  1]
λ₂ = 1 - 2i,   K₂ = [1/2 - 1/2i,  1]
```

### 4 — Homogeneous System (X' = AX)

**What it does:** Fully solves X' = AX, producing the general solution with arbitrary constants C₁, C₂ (and C₃ for 3×3). Automatically handles distinct real, repeated, and complex eigenvalues.

**How to input:**

1. Choose option `4`
2. Enter dimension and matrix A

**Example — real eigenvalues:**

```
x' = 3x + 2y - 8z
y' = y + 4z
z' = 5z

A[1][1]=3   A[1][2]=2   A[1][3]=-8
A[2][1]=0   A[2][2]=1   A[2][3]=4
A[3][1]=0   A[3][2]=0   A[3][3]=5

→  X(t) = C1·[v₁]e^(3t) + C2·[v₂]e^t + C3·[v₃]e^(5t)
```

**Example — complex eigenvalues:**

```
A = [ 3  -2 ]
    [ 4  -1 ]

→  X(t) = C1·e^t([a]cos(2t) - [b]sin(2t))
         + C2·e^t([a]sin(2t) + [b]cos(2t))
```

### 5 — Initial Value Problem (X' = AX, X(0) = X₀)

**What it does:** Solves X' = AX with given initial conditions, finding the specific constants C₁, C₂ (C₃) that satisfy X(0) = X₀.

**How to input:**

1. Choose option `5`
2. Enter dimension and matrix A
3. Enter each component of X(0)

**Example:**

```
A = [ 1/2   0  ]    X(0) = [4]
    [  1   -1/2]            [6]

→  C₁ = 2,  C₂ = 4
   X(t) = 2[0,1]e^(-1/2·t) + 4[1,1]e^(1/2·t)
```

### 6 — Nonhomogeneous System (X' = AX + F(t))

**What it does:** Solves X' = AX + F(t) by the method of undetermined coefficients. Finds a particular solution Xₚ and adds it to the complementary solution Xc.

**Supported forcing types:**

| F(t) | Trial solution |
|------|----------------|
| Constant vector **b** | constant vector **a** |
| **b**·e^(kt), k not an eigenvalue | **a**·e^(kt) |
| **b**·e^(kt), k equals an eigenvalue (resonance) | **a**·t·e^(kt) |
| Polynomial degree ≤ 2 | matching polynomial vector |

**How to input:**

1. Choose option `6`
2. Enter dimension and matrix A
3. Choose forcing type and enter the vector components

**Example:**

```
A = [ 5   9 ]    F(t) = [ 7 ]
    [-1  11 ]            [21]

→  X(t) = Xc + Xp
```

### 7 — Verify a Solution

**What it does:** Checks whether a proposed solution actually satisfies X' = AX by computing both sides symbolically and comparing.

**Supported types:** Exponential X = v·e^(λt), or trig components with sin(t) and cos(t).

**How to input:**

1. Choose option `7`
2. Select solution type
3. Enter matrix A, eigenvector v, and eigenvalue λ

**Example:**

```
A = [ 3  -3 ]    Proposed: X = [1, 3]·e^(-6t)
    [ 3  -7 ]

→  AX = λX  ✓
```

### 8 — Fundamental Set Check

**What it does:** Computes the Wronskian of n proposed solutions at t = 0 to verify they are linearly independent and form a complete fundamental set.

**How to input:**

1. Choose option `8`
2. Enter dimension n
3. For each solution, enter eigenvector v and eigenvalue λ

**Example:**

```
X₁ = [1,  1]·e^(-9t)
X₂ = [1, -1]·e^(-3t)

W(0) = det([ 1   1 ]) = -2 ≠ 0  →  Fundamental set ✓
           [ 1  -1 ]
```

### 9 — Example Library

13 built-in worked examples. Select option `9` and choose a number.

| # | Problem |
|---|---------|
| 1 | Scalar → matrix: `dx/dt = 8x - 4y, dy/dt = 2x + 3y` |
| 2 | 3×3 polynomial forcing → matrix form |
| 3 | Matrix → scalar: `X' = [[6,4],[-1,2]]X + [-5,6]e^t` |
| 4 | Verify `X = [1,3]e^(-6t)` satisfies a 2×2 system |
| 5 | Fundamental set: `X₁=[1,1]e^(-9t)`, `X₂=[1,-1]e^(-3t)` |
| 6 | 3×3 homogeneous: `dx/dt = x+y-z, dy/dt = 5y, dz/dt = y-z` |
| 7 | IVP: `A=[[1/2,0],[1,-1/2]]`, `X(0)=[4,6]` |
| 8 | Repeated eigenvalue: `dx/dt = 6x-y, dy/dt = 36x-6y` |
| 9 | 3×3 homogeneous with three distinct eigenvalues |
| 10 | Nonhomogeneous with constant forcing |
| 11 | Nonhomogeneous with exponential forcing |
| 12 | Complex eigenvalues: `A=[[3,-2],[4,-1]]` giving λ = 1 ± 2i |
| 13 | Complex IVP: `A=[[0,1],[-1,0]]`, `X(0)=[1,0]` |

## How to Enter Numbers

The script accepts integers, negatives, and fractions anywhere a number is expected.

| Value | Type this |
|-------|-----------|
| Integer | `3` or `-7` |
| Fraction | `1/2` or `-3/4` |
| Zero | `0` |

Fractions are kept **exact** throughout — no floating-point rounding on rational inputs.

**Entering a matrix** — you are prompted entry by entry, row by row:

```
[ 3  -2 ]   →   A[1][1] = 3    A[1][2] = -2
[ 4  -1 ]       A[2][1] = 4    A[2][2] = -1
```

## Mathematical Background

### Eigenvalues and Eigenvectors

An eigenvalue λ and eigenvector K of matrix A satisfy **AK = λK**. To find eigenvalues, solve the characteristic equation `det(A − λI) = 0`, which gives a quadratic (2×2) or cubic (3×3) polynomial. Each eigenvalue is substituted back into `(A − λI)K = 0` and solved by row reduction.

### Homogeneous Systems

For X' = AX with n distinct real eigenvalues λ₁ … λₙ and eigenvectors K₁ … Kₙ, the general solution is a linear combination:

```
X(t) = C₁K₁e^(λ₁t) + C₂K₂e^(λ₂t) + ... + CₙKₙe^(λₙt)
```

### Repeated Eigenvalues

When eigenvalue λ is repeated but the matrix is defective (only one independent eigenvector K), a **generalized eigenvector** P is found by solving `(A − λI)P = K`. This yields two independent solutions:

```
X₁ = K·e^(λt)
X₂ = (Kt + P)·e^(λt)
```

The script detects defectiveness automatically via rank analysis.

### Complex Eigenvalues

When eigenvalues come as conjugate pairs λ = α ± βi, the complex eigenvector K = **a** + i**b** yields two real solutions:

```
X₁(t) = e^(αt) [ a·cos(βt) − b·sin(βt) ]
X₂(t) = e^(αt) [ a·sin(βt) + b·cos(βt) ]
```

The script computes exact eigenvectors when β is rational. When β is irrational, it falls back to a decimal approximation.

### Initial Value Problems

Evaluating the general solution at t = 0 (where e^0 = 1, cos 0 = 1, sin 0 = 0) reduces to a linear system `C₁K₁ + C₂K₂ + ... = X₀`, solved by Gaussian elimination.

### Nonhomogeneous Systems (Undetermined Coefficients)

For X' = AX + F(t) the general solution is **X = Xc + Xp**. A trial form for Xp is guessed based on F(t), substituted into the ODE, and solved by matching coefficients. If the forcing rate matches an eigenvalue (**resonance**), the trial is multiplied by t.

### Fundamental Sets and the Wronskian

A set of n solutions is a fundamental set if and only if the Wronskian is nonzero:

```
W(t) = det[ X₁ | X₂ | ... | Xₙ ] ≠ 0
```

The script evaluates W at t = 0 for simplicity.

## Example Problems

### x' = 22x − 30y,  y' = 35x − 33y

The eigenvalues are −5.5 ± 17.14i (β = 5√47/2, irrational, so decimal output):

```
Option 4  →  dimension 2

A[1][1] = 22    A[1][2] = -30
A[2][1] = 35    A[2][2] = -33

X(t) = C1·e^(-5.5t)([0.7857, 1]cos(17.14t) − [0.4897, 0]sin(17.14t))
     + C2·e^(-5.5t)([0.7857, 1]sin(17.14t) + [0.4897, 0]cos(17.14t))
```

### x' = 3x + 2y − 8z,  y' = y + 4z,  z' = 5z

Upper-triangular matrix — eigenvalues are just the diagonal entries (3, 1, 5):

```
Option 4  →  dimension 3

A[1][1]=3   A[1][2]=2   A[1][3]=-8
A[2][1]=0   A[2][2]=1   A[2][3]=4
A[3][1]=0   A[3][2]=0   A[3][3]=5
```

## Limitations

- Supports **2×2 and 3×3** systems only
- Polynomial forcing up to **degree 2** (t²)
- Nonhomogeneous solver handles **one forcing type at a time** — mixed forcing requires combining results manually
- Complex eigenvalues with **irrational β** produce decimal approximations only
- No support for **non-constant coefficient** matrices (e.g. x' = t·x + y)
