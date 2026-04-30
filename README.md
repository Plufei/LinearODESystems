\# Linear ODE System Solver



A self-contained Python script for solving linear systems of ordinary differential equations of the form \*\*X' = AX + F(t)\*\*. No external libraries required — all arithmetic is done with exact rational fractions internally.



```bash

python3 ode\_solver.py

```



\---



\## Table of Contents



\- \[Features](#features)

\- \[Requirements](#requirements)

\- \[Main Menu Overview](#main-menu-overview)

\- \[Functions \& How to Use Each](#functions--how-to-use-each)

&#x20; - \[1. Convert Scalar Equations → Matrix Form](#1-convert-scalar-equations--matrix-form)

&#x20; - \[2. Convert Matrix Form → Scalar Equations](#2-convert-matrix-form--scalar-equations)

&#x20; - \[3. Eigenvalues \& Eigenvectors](#3-eigenvalues--eigenvectors)

&#x20; - \[4. Homogeneous System](#4-homogeneous-system)

&#x20; - \[5. Initial Value Problem (IVP)](#5-initial-value-problem-ivp)

&#x20; - \[6. Nonhomogeneous System](#6-nonhomogeneous-system)

&#x20; - \[7. Verify a Solution](#7-verify-a-solution)

&#x20; - \[8. Fundamental Set Check](#8-fundamental-set-check)

&#x20; - \[9. Example Library](#9-example-library)

\- \[How to Enter Numbers](#how-to-enter-numbers)

\- \[Mathematical Background](#mathematical-background)

&#x20; - \[Eigenvalues \& Eigenvectors](#eigenvalues--eigenvectors-math)

&#x20; - \[Homogeneous Systems](#homogeneous-systems-math)

&#x20; - \[Repeated Eigenvalues](#repeated-eigenvalues-math)

&#x20; - \[Complex Eigenvalues](#complex-eigenvalues-math)

&#x20; - \[Initial Value Problems](#initial-value-problems-math)

&#x20; - \[Nonhomogeneous Systems](#nonhomogeneous-systems-math)

&#x20; - \[Fundamental Sets \& the Wronskian](#fundamental-sets--the-wronskian-math)

\- \[Example Problems Walkthrough](#example-problems-walkthrough)

\- \[Limitations](#limitations)



\---



\## Features



\- Solves 2×2 and 3×3 linear ODE systems

\- Exact rational arithmetic (no floating-point errors on rational inputs)

\- Handles all eigenvalue cases:

&#x20; - Distinct real eigenvalues

&#x20; - Repeated eigenvalues (defective matrices, generalized eigenvectors)

&#x20; - Complex conjugate pairs — exact when β is rational, decimal approximation otherwise

\- Nonhomogeneous forcing by undetermined coefficients: constant, exponential, and polynomial F(t)

\- Resonance detection for nonhomogeneous systems

\- IVP solver — applies initial conditions to find specific constants C₁, C₂, C₃

\- Solution verification — checks if a proposed vector satisfies X' = AX

\- Wronskian check to confirm a fundamental set of solutions

\- Bidirectional conversion between scalar equations and matrix form

\- Step-by-step output with copy-paste answer block



\---



\## Requirements



\- Python 3.7 or higher

\- No third-party packages (no numpy, scipy, or sympy)



```bash

python3 ode\_solver.py

```



\---



\## Main Menu Overview



```

1\. Scalar equations  →  Matrix form

2\. Matrix form       →  Scalar equations

3\. Eigenvalues \& Eigenvectors

4\. Solve Homogeneous System   X' = AX

5\. Solve Initial Value Problem

6\. Solve Nonhomogeneous System

7\. Verify a proposed solution

8\. Fundamental Set check

9\. Example Library

0\. Quit

```



\---



\## Functions \& How to Use Each



\---



\### 1. Convert Scalar Equations → Matrix Form



\*\*What it does:\*\* Takes a system written as individual scalar equations and rewrites it as the matrix equation X' = AX + F(t).



\*\*When to use it:\*\* You have equations like `dx/dt = 3x - 2y + t` and want to identify A and F(t) before solving.



\*\*How to input:\*\*



1\. Choose option `1` from the main menu.

2\. Enter the system dimension (2 or 3).

3\. For each equation, type the right-hand side. Variables are `x`, `y`, `z`.



Supported terms in each equation:

| Term | Example input |

|------|---------------|

| Coefficient times variable | `3x`, `-2y`, `5z` |

| Constant | `7`, `-4` |

| Linear in t | `2t`, `-t` |

| Quadratic in t | `t^2`, `3t^2` |

| Exponential | `e^t`, `2e^t`, `e^(3t)` |



\*\*Example:\*\*

```

System:  dx/dt = 8x - 4y

&#x20;        dy/dt = 2x + 3y



Input:

&#x20; dx/dt = 8x - 4y

&#x20; dy/dt = 2x + 3y



Output:

&#x20; A = \[ 8  -4 ]     F(t) = \[0]

&#x20;     \[ 2   3 ]             \[0]

```



\---



\### 2. Convert Matrix Form → Scalar Equations



\*\*What it does:\*\* The reverse — given a matrix A and optional forcing vector, writes out the individual scalar equations.



\*\*When to use it:\*\* You have a matrix A from a textbook and want to see what the equations look like written out.



\*\*How to input:\*\*



1\. Choose option `2`.

2\. Enter dimension, then matrix A row by row.

3\. Choose forcing type: none, constant vector b, or b·e^(kt).



\*\*Example:\*\*

```

A = \[ 6   4 ]    b = \[-5]·e^t

&#x20;   \[-1   2 ]        \[ 6]



Output:

&#x20; dx/dt = 6x + 4y - 5e^t

&#x20; dy/dt = -x + 2y + 6e^t

```



\---



\### 3. Eigenvalues \& Eigenvectors



\*\*What it does:\*\* Computes the characteristic polynomial, eigenvalues, and corresponding eigenvectors for a 2×2 or 3×3 matrix.



\*\*When to use it:\*\* As a standalone computation, or to check your work before solving the full system.



\*\*How to input:\*\*



1\. Choose option `3`.

2\. Enter dimension and matrix A.



\*\*Example:\*\*

```

A = \[ 3  -2 ]

&#x20;   \[ 4  -1 ]



Output:

&#x20; Characteristic polynomial: λ² - 2λ + 5 = 0

&#x20; λ₁ = 1 + 2i,   K₁ = \[1/2 + 1/2i, 1]

&#x20; λ₂ = 1 - 2i,   K₂ = \[1/2 - 1/2i, 1]

```



\---



\### 4. Homogeneous System



\*\*What it does:\*\* Fully solves X' = AX, producing the general solution with arbitrary constants C₁, C₂, (C₃). Handles distinct real, repeated, and complex eigenvalues automatically.



\*\*When to use it:\*\* Any time you need the general solution to a homogeneous linear system.



\*\*How to input:\*\*



1\. Choose option `4`.

2\. Enter dimension and matrix A.



\*\*Example — distinct real eigenvalues:\*\*

```

System:  x' = 3x + 2y - 8z

&#x20;        y' = y + 4z

&#x20;        z' = 5z



A = \[ 3  2  -8 ]

&#x20;   \[ 0  1   4 ]

&#x20;   \[ 0  0   5 ]



Eigenvalues: 3, 1, 5  (diagonal of upper-triangular matrix)

Output: X(t) = C1\[v₁]e^(3t) + C2\[v₂]e^t + C3\[v₃]e^(5t)

```



\*\*Example — complex eigenvalues:\*\*

```

A = \[ 3  -2 ]

&#x20;   \[ 4  -1 ]



Eigenvalues: 1 ± 2i

Output:

&#x20; X(t) = C1·e^t(\[a]cos(2t) - \[b]sin(2t)) + C2·e^t(\[a]sin(2t) + \[b]cos(2t))

```



\---



\### 5. Initial Value Problem (IVP)



\*\*What it does:\*\* Solves X' = AX with X(0) = X₀ — applies initial conditions to find the specific constants C₁, C₂, (C₃), giving a unique particular solution.



\*\*When to use it:\*\* You have both a system and starting values at t = 0.



\*\*How to input:\*\*



1\. Choose option `5`.

2\. Enter dimension and matrix A.

3\. Enter the initial condition vector X(0), one component at a time.



\*\*Example:\*\*

```

A = \[ 1/2   0  ]     X(0) = \[4]

&#x20;   \[  1   -1/2]             \[6]



Output:

&#x20; C₁ = 2,  C₂ = 4

&#x20; X(t) = 2\[0,1]e^(-1/2·t) + 4\[1,1]e^(1/2·t)

```



\---



\### 6. Nonhomogeneous System



\*\*What it does:\*\* Solves X' = AX + F(t) using the \*\*method of undetermined coefficients\*\*. Finds a particular solution Xₚ(t) and adds it to the homogeneous general solution Xc.



\*\*Supported forcing functions F(t):\*\*

| Type | Example |

|------|---------|

| Constant vector | `b = \[7, 21]` |

| Exponential vector | `b·e^(kt)` e.g. `\[-5, 6]e^t` |

| Polynomial vector | `a + bt + ct²` (per component) |



\*\*How to input:\*\*



1\. Choose option `6`.

2\. Enter dimension and matrix A.

3\. Choose forcing type and enter the components.



\*\*Resonance:\*\* If the exponential rate k matches an eigenvalue of A, the script automatically applies the resonance correction (multiplies the trial solution by t).



\*\*Example:\*\*

```

A = \[ 5   9 ]    F(t) = \[ 7 ]

&#x20;   \[-1  11 ]            \[21]



Output: X(t) = Xc + Xp

&#x20; Xc = C1\[v₁]e^(λ₁t) + C2\[v₂]e^(λ₂t)

&#x20; Xp = \[particular constant vector]

```



\---



\### 7. Verify a Solution



\*\*What it does:\*\* Checks whether a proposed vector solution actually satisfies X' = AX by computing both sides and comparing.



\*\*When to use it:\*\* To double-check a textbook answer or your own result.



\*\*Supported solution types:\*\*

\- Exponential: X = v·e^(λt)

\- Trigonometric: X = \[a·sin(t) + b·cos(t), ...] per component



\*\*How to input:\*\*



1\. Choose option `7`.

2\. Select solution type (exponential or trig).

3\. Enter matrix A, then the proposed vector v and eigenvalue λ.



\*\*Example:\*\*

```

A = \[ 3  -3 ]    Proposed: X = \[1, 3]·e^(-6t)

&#x20;   \[ 3  -7 ]



Output: AX = λX  ✓  (solution verified)

```



\---



\### 8. Fundamental Set Check



\*\*What it does:\*\* Given n proposed solution vectors X₁, X₂, (X₃), computes the \*\*Wronskian\*\* at t = 0 to verify they form a fundamental set (i.e., are linearly independent).



\*\*When to use it:\*\* To confirm that your set of solutions is complete and can represent the general solution.



\*\*How to input:\*\*



1\. Choose option `8`.

2\. Enter dimension n.

3\. For each solution, enter the eigenvector v and eigenvalue λ (solution is v·e^(λt)).



\*\*Example:\*\*

```

X₁ = \[1, 1]·e^(-9t)    X₂ = \[1, -1]·e^(-3t)



Wronskian at t=0:

&#x20; W = det(\[1, 1], \[1, -1]) = -2 ≠ 0  →  Fundamental set ✓

```



\---



\### 9. Example Library



13 built-in worked examples covering every feature. Run the script and select option `9` to browse them.



| # | Problem |

|---|---------|

| 1 | Scalar → matrix: `dx/dt = 8x - 4y, dy/dt = 2x + 3y` |

| 2 | 3×3 polynomial forcing → matrix form |

| 3 | Matrix → scalar: `X' = \[\[6,4],\[-1,2]]X + \[-5,6]e^t` |

| 4 | Verify: `X = \[1,3]e^(-6t)` for a 2×2 system |

| 5 | Fundamental set: `X₁=\[1,1]e^(-9t)`, `X₂=\[1,-1]e^(-3t)` |

| 6 | 3×3 homogeneous: `dx/dt = x+y-z, dy/dt = 5y, dz/dt = y-z` |

| 7 | IVP: `A=\[\[1/2,0],\[1,-1/2]]`, `X(0)=\[4,6]` |

| 8 | Repeated eigenvalue: `dx/dt = 6x-y, dy/dt = 36x-6y` |

| 9 | 3×3 homogeneous with three distinct eigenvalues |

| 10 | Nonhomogeneous constant forcing |

| 11 | Nonhomogeneous exponential forcing |

| 12 | Complex eigenvalues: `A=\[\[3,-2],\[4,-1]]` → λ = 1 ± 2i |

| 13 | Complex IVP: `A=\[\[0,1],\[-1,0]]`, `X(0)=\[1,0]` → `\[cos(t), -sin(t)]` |



\---



\## How to Enter Numbers



The script accepts fractions and negatives everywhere a number is expected:



| Value | Input |

|-------|-------|

| Integer | `3`, `-7`, `0` |

| Fraction | `1/2`, `-3/4`, `5/3` |

| Zero | `0` |

| Negative fraction | `-1/2` |



\*\*Matrix entry:\*\* Entered row by row, element by element. For the matrix

```

\[ 3  -2 ]

\[ 4  -1 ]

```

you type: `3`, `-2`, `4`, `-1` at the prompts.



\*\*Fractions are kept exact throughout\*\* — the solver uses rational arithmetic internally, so `1/3 + 2/3` gives exactly `1`, not `0.9999...`.



\---



\## Mathematical Background



\### Eigenvalues \& Eigenvectors (math)



For a square matrix A, an \*\*eigenvalue\*\* λ and \*\*eigenvector\*\* K satisfy:



```

AK = λK

```



To find eigenvalues, solve the \*\*characteristic equation\*\*:



```

det(A - λI) = 0

```



For a 2×2 matrix this is a quadratic in λ. For 3×3 it's a cubic. Each eigenvalue λᵢ is then substituted back into `(A - λᵢI)K = 0` to find the corresponding eigenvector Kᵢ by row reduction.



\---



\### Homogeneous Systems (math)



A homogeneous linear system has the form:



```

X' = AX

```



If A has \*\*n distinct real eigenvalues\*\* λ₁, λ₂, ..., λₙ with eigenvectors K₁, K₂, ..., Kₙ, the general solution is:



```

X(t) = C₁K₁e^(λ₁t) + C₂K₂e^(λ₂t) + ... + CₙKₙe^(λₙt)

```



Each term is an independent solution; the Cᵢ are arbitrary constants determined by initial conditions.



\---



\### Repeated Eigenvalues (math)



When an eigenvalue λ is repeated (multiplicity 2) but the matrix is \*\*defective\*\* (only one linearly independent eigenvector K), a second independent solution is constructed using a \*\*generalized eigenvector\*\* P:



```

(A - λI)P = K

```



The two independent solutions are then:



```

X₁ = K·e^(λt)

X₂ = (Kt + P)·e^(λt)

```



The script detects defectiveness automatically via rank analysis and switches solvers accordingly.



\---



\### Complex Eigenvalues (math)



When the characteristic polynomial has a negative discriminant, eigenvalues come in \*\*conjugate pairs\*\*: λ = α ± βi. Rather than using complex exponentials, the two real solutions are extracted from the complex eigenvector K = \*\*a\*\* + i\*\*b\*\*:



```

X₁(t) = e^(αt) \[ a·cos(βt) − b·sin(βt) ]

X₂(t) = e^(αt) \[ a·sin(βt) + b·cos(βt) ]

```



where \*\*a\*\* is the real part and \*\*b\*\* is the imaginary part of the eigenvector. The general contribution of this pair to the solution is `C₁X₁ + C₂X₂`.



When β is rational, the script computes exact eigenvectors. When β is irrational (e.g. β = 5√47/2), it falls back to a decimal approximation.



\---



\### Initial Value Problems (math)



Applying X(0) = X₀ to the general solution at t = 0:



```

X(0) = C₁K₁ + C₂K₂ + ... = X₀

```



Since e^(λᵢ·0) = 1 and cos(0) = 1, sin(0) = 0, this reduces to a linear algebraic system in C₁, C₂, .... The script solves it by Gaussian elimination to find the unique constants.



\---



\### Nonhomogeneous Systems (math)



For X' = AX + F(t), the \*\*general solution\*\* is:



```

X(t) = Xc(t) + Xp(t)

```



where \*\*Xc\*\* is the complementary (homogeneous) solution and \*\*Xp\*\* is a \*\*particular solution\*\* found by undetermined coefficients:



| Forcing F(t) | Trial form for Xp |

|---|---|

| Constant vector \*\*b\*\* | \*\*a\*\* (constant vector) |

| \*\*b\*\*·e^(kt), k not an eigenvalue | \*\*a\*\*·e^(kt) |

| \*\*b\*\*·e^(kt), k equals eigenvalue (resonance) | \*\*a\*\*t·e^(kt) |

| Polynomial \*\*b\*\*t² + \*\*c\*\*t + \*\*d\*\* | \*\*a\*\*t² + \*\*b\*\*t + \*\*c\*\* |



The undetermined coefficient vectors are found by substituting the trial form into X' = AX + F(t) and matching coefficients.



\---



\### Fundamental Sets \& the Wronskian (math)



A set of n solutions {X₁, X₂, ..., Xₙ} is a \*\*fundamental set\*\* if and only if they are linearly independent, verified by the \*\*Wronskian\*\*:



```

W(t) = det\[ X₁ | X₂ | ... | Xₙ ]

```



If W(t₀) ≠ 0 for any single point t₀, the set is fundamental. The script evaluates W at t = 0 (where e^(λt) = 1 simplifies the computation).



\---



\## Example Problems Walkthrough



\### Solving x' = 22x − 30y, y' = 35x − 33y



The eigenvalues turn out to be −5.5 ± 17.14i (irrational β):



```

Main menu → 4 (Homogeneous System) → dimension 2



A\[1]\[1] = 22    A\[1]\[2] = -30

A\[2]\[1] = 35    A\[2]\[2] = -33



Output:

&#x20; λ ≈ −5.5 ± 17.14i

&#x20; X(t) = C1·e^(-5.5t)(\[0.7857,1]cos(17.14t) − \[0.4897,0]sin(17.14t))

&#x20;      + C2·e^(-5.5t)(\[0.7857,1]sin(17.14t) + \[0.4897,0]cos(17.14t))

```



\---



\### Solving x' = 3x + 2y − 8z, y' = y + 4z, z' = 5z



Upper-triangular matrix → eigenvalues are the diagonal entries (3, 1, 5):



```

Main menu → 4 → dimension 3



A\[1]\[1]=3  A\[1]\[2]=2  A\[1]\[3]=-8

A\[2]\[1]=0  A\[2]\[2]=1  A\[2]\[3]=4

A\[3]\[1]=0  A\[3]\[2]=0  A\[3]\[3]=5

```



\---



\### Applying an Initial Condition



```

Main menu → 5 (IVP) → dimension 2

Enter A, then X(0) component by component.

The script solves for C₁ and C₂ and prints the specific solution.

```



\---



\## Limitations



\- Supports \*\*2×2 and 3×3\*\* systems only.

\- Polynomial forcing up to \*\*degree 2\*\* (t²).

\- Nonhomogeneous solver handles one forcing term at a time; mixed forcing (e.g. constant + exponential simultaneously) requires combining results manually.

\- For complex eigenvalues with irrational β, eigenvectors and the solution are given as \*\*decimal approximations\*\* (4 significant figures).

\- No support for non-constant coefficient matrices (e.g. x' = t·x + y).

