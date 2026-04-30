#!/usr/bin/env python3
"""
Stepwise Wronskian Linear Independence Checker
No external dependencies.

Works with simple functions built from:
  - sin(t), cos(t), tan(t)
  - e^(kt), exp(kt)
  - t^n, constants, and linear combinations
  - products like e^(2t)*sin(t), t^2*e^t, 3*t*cos(t)

Main purpose:
  Given n functions f1(t), ..., fn(t), compute the Wronskian matrix,
  determinant, and decide whether the functions are linearly independent
  when the Wronskian is not identically zero.

Important math note:
  If W(t) is not identically zero, the functions are linearly independent.
  If W(t) is identically zero, this script reports that the Wronskian test
  is inconclusive or suggests dependence for common simple cases.
"""

import math

EPS = 1e-9

# ============================================================
# Expression tree
# ============================================================

class Expr:
    def eval(self, t):
        raise NotImplementedError

    def deriv(self):
        raise NotImplementedError

    def simplify(self):
        return self

    def __add__(self, other):
        return Add([self, to_expr(other)]).simplify()

    def __sub__(self, other):
        return Add([self, Mul([Const(-1), to_expr(other)])]).simplify()

    def __mul__(self, other):
        return Mul([self, to_expr(other)]).simplify()

    def __neg__(self):
        return Mul([Const(-1), self]).simplify()

class Const(Expr):
    def __init__(self, value):
        self.value = float(value)

    def eval(self, t):
        return self.value

    def deriv(self):
        return Const(0)

    def simplify(self):
        if abs(self.value) < EPS:
            return Const(0)
        return self

    def __str__(self):
        if abs(self.value - round(self.value)) < EPS:
            return str(int(round(self.value)))
        return str(self.value)

class Var(Expr):
    def eval(self, t):
        return t

    def deriv(self):
        return Const(1)

    def __str__(self):
        return "t"

class Add(Expr):
    def __init__(self, terms):
        self.terms = terms

    def eval(self, t):
        return sum(term.eval(t) for term in self.terms)

    def deriv(self):
        return Add([term.deriv() for term in self.terms]).simplify()

    def simplify(self):
        new_terms = []
        const_sum = 0.0

        for term in self.terms:
            term = term.simplify()
            if isinstance(term, Add):
                for sub in term.terms:
                    new_terms.append(sub)
            elif isinstance(term, Const):
                const_sum += term.value
            else:
                new_terms.append(term)

        if abs(const_sum) > EPS:
            new_terms.append(Const(const_sum))

        if not new_terms:
            return Const(0)

        if len(new_terms) == 1:
            return new_terms[0]

        return Add(new_terms)

    def __str__(self):
        out = ""
        for i, term in enumerate(self.terms):
            s = str(term)
            if i == 0:
                out += s
            else:
                if s.startswith("-"):
                    out += " - " + s[1:]
                else:
                    out += " + " + s
        return out

class Mul(Expr):
    def __init__(self, factors):
        self.factors = factors

    def eval(self, t):
        result = 1.0
        for factor in self.factors:
            result *= factor.eval(t)
        return result

    def deriv(self):
        # Product rule:
        # (f1*f2*...*fn)' = sum over i of f1*...*fi'*...*fn
        terms = []
        for i in range(len(self.factors)):
            factors = []
            for j, factor in enumerate(self.factors):
                factors.append(factor.deriv() if i == j else factor)
            terms.append(Mul(factors))
        return Add(terms).simplify()

    def simplify(self):
        new_factors = []
        const_product = 1.0

        for factor in self.factors:
            factor = factor.simplify()
            if isinstance(factor, Mul):
                for sub in factor.factors:
                    new_factors.append(sub)
            elif isinstance(factor, Const):
                const_product *= factor.value
            else:
                new_factors.append(factor)

        if abs(const_product) < EPS:
            return Const(0)

        if abs(const_product - 1.0) > EPS:
            new_factors.insert(0, Const(const_product))

        if not new_factors:
            return Const(1)

        if len(new_factors) == 1:
            return new_factors[0]

        return Mul(new_factors)

    def __str__(self):
        parts = []
        for factor in self.factors:
            if isinstance(factor, Add):
                parts.append("(" + str(factor) + ")")
            else:
                parts.append(str(factor))
        return "*".join(parts)

class Pow(Expr):
    def __init__(self, base, power):
        self.base = base
        self.power = int(power)

    def eval(self, t):
        return self.base.eval(t) ** self.power

    def deriv(self):
        # Only supports integer powers.
        if self.power == 0:
            return Const(0)
        return Mul([Const(self.power), Pow(self.base, self.power - 1), self.base.deriv()]).simplify()

    def simplify(self):
        base = self.base.simplify()
        if self.power == 0:
            return Const(1)
        if self.power == 1:
            return base
        if isinstance(base, Const):
            return Const(base.value ** self.power)
        return Pow(base, self.power)

    def __str__(self):
        if isinstance(self.base, Add):
            return "(" + str(self.base) + ")^" + str(self.power)
        return str(self.base) + "^" + str(self.power)

class Sin(Expr):
    def __init__(self, arg):
        self.arg = arg

    def eval(self, t):
        return math.sin(self.arg.eval(t))

    def deriv(self):
        return Mul([Cos(self.arg), self.arg.deriv()]).simplify()

    def simplify(self):
        self.arg = self.arg.simplify()
        return self

    def __str__(self):
        return "sin(" + str(self.arg) + ")"

class Cos(Expr):
    def __init__(self, arg):
        self.arg = arg

    def eval(self, t):
        return math.cos(self.arg.eval(t))

    def deriv(self):
        return Mul([Const(-1), Sin(self.arg), self.arg.deriv()]).simplify()

    def simplify(self):
        self.arg = self.arg.simplify()
        return self

    def __str__(self):
        return "cos(" + str(self.arg) + ")"

class Tan(Expr):
    def __init__(self, arg):
        self.arg = arg

    def eval(self, t):
        return math.tan(self.arg.eval(t))

    def deriv(self):
        # d/dt tan(u) = sec^2(u)u' = u' / cos^2(u)
        return Mul([Pow(Cos(self.arg), -2), self.arg.deriv()]).simplify()

    def simplify(self):
        self.arg = self.arg.simplify()
        return self

    def __str__(self):
        return "tan(" + str(self.arg) + ")"

class Exp(Expr):
    def __init__(self, arg):
        self.arg = arg

    def eval(self, t):
        return math.exp(self.arg.eval(t))

    def deriv(self):
        return Mul([Exp(self.arg), self.arg.deriv()]).simplify()

    def simplify(self):
        self.arg = self.arg.simplify()
        return self

    def __str__(self):
        return "e^(" + str(self.arg) + ")"

def to_expr(x):
    if isinstance(x, Expr):
        return x
    return Const(x)

# ============================================================
# Parser
# ============================================================

class Parser:
    def __init__(self, text):
        self.text = text.replace(" ", "")
        self.i = 0

    def peek(self):
        if self.i >= len(self.text):
            return ""
        return self.text[self.i]

    def eat(self, ch=None):
        if ch is None:
            c = self.peek()
            self.i += 1
            return c
        if self.peek() == ch:
            self.i += 1
            return True
        return False

    def parse(self):
        if self.text == "":
            raise ValueError("Empty expression")
        expr = self.parse_sum()
        if self.i != len(self.text):
            raise ValueError("Unexpected text near: " + self.text[self.i:])
        return expr.simplify()

    def parse_sum(self):
        expr = self.parse_product()
        while True:
            if self.eat("+"):
                expr = Add([expr, self.parse_product()]).simplify()
            elif self.eat("-"):
                expr = Add([expr, Mul([Const(-1), self.parse_product()])]).simplify()
            else:
                break
        return expr

    def starts_factor(self):
        c = self.peek()
        return c.isdigit() or c == "." or c == "t" or c == "e" or c == "(" or c.isalpha() or c == "-"

    def parse_product(self):
        factors = [self.parse_power()]

        while True:
            if self.eat("*"):
                factors.append(self.parse_power())
            elif self.eat("/"):
                denom = self.parse_power()
                factors.append(Pow(denom, -1))
            else:
                # Allow implicit multiplication:
                # 2t, 3sin(t), e^t cos(t), etc.
                c = self.peek()
                if c and (c == "t" or c == "(" or c.isalpha()):
                    factors.append(self.parse_power())
                else:
                    break

        return Mul(factors).simplify()

    def parse_power(self):
        expr = self.parse_unary()
        while self.eat("^"):
            power_expr = self.parse_unary()

            if isinstance(power_expr, Const):
                power = int(round(power_expr.value))
                expr = Pow(expr, power).simplify()
            elif isinstance(expr, Const) and abs(expr.value - math.e) < EPS:
                # e^(anything)
                expr = Exp(power_expr).simplify()
            else:
                raise ValueError("Only integer powers are supported, except e^(...)")
        return expr

    def parse_unary(self):
        if self.eat("+"):
            return self.parse_unary()
        if self.eat("-"):
            return Mul([Const(-1), self.parse_unary()]).simplify()
        return self.parse_atom()

    def parse_atom(self):
        c = self.peek()

        if c == "(":
            self.eat("(")
            expr = self.parse_sum()
            if not self.eat(")"):
                raise ValueError("Missing closing parenthesis")
            return expr

        if c.isdigit() or c == ".":
            return self.parse_number()

        if self.text.startswith("sin", self.i):
            self.i += 3
            return Sin(self.parse_function_argument()).simplify()

        if self.text.startswith("cos", self.i):
            self.i += 3
            return Cos(self.parse_function_argument()).simplify()

        if self.text.startswith("tan", self.i):
            self.i += 3
            return Tan(self.parse_function_argument()).simplify()

        if self.text.startswith("exp", self.i):
            self.i += 3
            return Exp(self.parse_function_argument()).simplify()

        if self.text.startswith("e", self.i):
            self.i += 1
            return Const(math.e)

        if self.text.startswith("t", self.i):
            self.i += 1
            return Var()

        raise ValueError("Unexpected character near: " + self.text[self.i:])

    def parse_function_argument(self):
        if self.eat("("):
            expr = self.parse_sum()
            if not self.eat(")"):
                raise ValueError("Missing closing parenthesis")
            return expr
        return self.parse_atom()

    def parse_number(self):
        start = self.i
        while self.peek().isdigit() or self.peek() == ".":
            self.i += 1
        return Const(float(self.text[start:self.i]))

def parse_expr(text):
    # Normalize common inputs.
    text = text.strip()
    text = text.replace(" ", "")
    text = text.replace("π", str(math.pi))
    text = text.replace("pi", str(math.pi))
    return Parser(text).parse()

# ============================================================
# Determinants and Wronskian
# ============================================================

def determinant_numeric(M):
    n = len(M)
    A = [[float(M[i][j]) for j in range(n)] for i in range(n)]
    det = 1.0

    for col in range(n):
        pivot = col
        for row in range(col + 1, n):
            if abs(A[row][col]) > abs(A[pivot][col]):
                pivot = row

        if abs(A[pivot][col]) < EPS:
            return 0.0

        if pivot != col:
            A[col], A[pivot] = A[pivot], A[col]
            det *= -1

        pivot_value = A[col][col]
        det *= pivot_value

        for row in range(col + 1, n):
            factor = A[row][col] / pivot_value
            for j in range(col, n):
                A[row][j] -= factor * A[col][j]

    return det

def determinant_symbolic(M):
    n = len(M)

    if n == 1:
        return M[0][0].simplify()

    if n == 2:
        return (M[0][0] * M[1][1] - M[0][1] * M[1][0]).simplify()

    total = Const(0)
    for col in range(n):
        minor = []
        for i in range(1, n):
            row = []
            for j in range(n):
                if j != col:
                    row.append(M[i][j])
            minor.append(row)

        cofactor = determinant_symbolic(minor)
        if col % 2 == 1:
            cofactor = -cofactor

        total = total + M[0][col] * cofactor

    return total.simplify()

def wronskian_matrix(functions):
    n = len(functions)
    rows = []
    current = functions[:]

    for row_index in range(n):
        rows.append(current[:])
        current = [f.deriv().simplify() for f in current]

    return rows

def print_matrix(M):
    for row in M:
        print("  [" + " , ".join(str(x) for x in row) + "]")

def numeric_zero_test(expr):
    # Test many points. If any point is nonzero, W is not identically zero.
    test_points = [-3, -2, -1, -0.5, 0, 0.25, 0.5, 1, 2, 3]
    values = []
    nonzero_points = []

    for t in test_points:
        try:
            val = expr.eval(t)
            values.append((t, val))
            if abs(val) > 1e-7 and math.isfinite(val):
                nonzero_points.append((t, val))
        except Exception:
            pass

    return values, nonzero_points

def check_wronskian(function_texts, verbose=True):
    funcs = [parse_expr(s) for s in function_texts]
    n = len(funcs)

    if verbose:
        print("\n" + "=" * 64)
        print("STEP 1: Functions")
        print("=" * 64)
        for i, f in enumerate(funcs):
            print("f" + str(i + 1) + "(t) = " + str(f))

    Wmat = wronskian_matrix(funcs)

    if verbose:
        print("\n" + "=" * 64)
        print("STEP 2: Wronskian matrix")
        print("=" * 64)
        print("Rows are derivatives from order 0 to order " + str(n - 1) + ".")
        print_matrix(Wmat)

    W = determinant_symbolic(Wmat).simplify()

    if verbose:
        print("\n" + "=" * 64)
        print("STEP 3: Determinant")
        print("=" * 64)
        print("W(t) = " + str(W))

    values, nonzero_points = numeric_zero_test(W)

    if verbose:
        print("\n" + "=" * 64)
        print("STEP 4: Test whether W(t) is identically zero")
        print("=" * 64)
        for t, val in values:
            if abs(val - round(val)) < 1e-7:
                shown = str(int(round(val)))
            else:
                shown = str(round(val, 8))
            print("W(" + str(t) + ") = " + shown)

    if nonzero_points:
        t0, val0 = nonzero_points[0]

        if verbose:
            print("\n" + "=" * 64)
            print("CONCLUSION")
            print("=" * 64)
            print("Since W(" + str(t0) + ") is not 0, W(t) is not identically zero.")
            print("Therefore the functions are linearly independent.")

        return True, W

    else:
        if verbose:
            print("\n" + "=" * 64)
            print("CONCLUSION")
            print("=" * 64)
            print("The sampled values of W(t) are all 0.")
            print("For many homework problems, this means the functions are dependent.")
            print("More carefully: W(t) = 0 makes the Wronskian test inconclusive unless")
            print("the functions are known solutions of the same linear homogeneous ODE.")

        return False, W

# ============================================================
# Simple interface
# ============================================================

def example_cos_sin():
    print("\nExample: cos(t), sin(t)")
    check_wronskian(["cos(t)", "sin(t)"], verbose=True)

def menu():
    while True:
        print("\n" + "=" * 64)
        print("Stepwise Wronskian Linear Independence Checker")
        print("=" * 64)
        print("1. Enter functions manually")
        print("2. Example: cos(t), sin(t)")
        print("3. Quit")

        choice = input("Choose an option: ").strip()

        if choice == "1":
            try:
                n = int(input("How many functions? ").strip())
                if n < 1:
                    print("Enter at least 1 function.")
                    continue

                funcs = []
                print("\nEnter functions using t as the variable.")
                print("Examples: cos(t), sin(t), e^t, e^(2t), t^2, t*e^t")
                for i in range(n):
                    funcs.append(input("f" + str(i + 1) + "(t) = ").strip())

                check_wronskian(funcs, verbose=True)

            except Exception as e:
                print("Error:", e)

        elif choice == "2":
            example_cos_sin()

        elif choice == "3" or choice.lower() in ("q", "quit", "exit"):
            print("Goodbye.")
            break

        else:
            print("Invalid option.")

if __name__ == "__main__":
    menu()
