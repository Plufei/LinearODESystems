#!/usr/bin/env python3
"""
Linear Systems of ODE Solver
Solves X' = AX + F(t) with exact fraction arithmetic.
No external dependencies (no numpy, scipy, sympy).
"""
import sys
import math

# ============================================================
# SECTION 1: FRACTION ARITHMETIC
# Fractions represented as tuples (numerator, denominator)
# Invariant: denominator > 0 always, reduced form
# ============================================================

def gcd(a, b):
    a, b = abs(int(a)), abs(int(b))
    while b:
        a, b = b, a % b
    return a if a else 1

def frac_simplify(num, den):
    num, den = int(num), int(den)
    if den == 0:
        raise ZeroDivisionError("Denominator cannot be zero")
    if num == 0:
        return (0, 1)
    sign = -1 if (num < 0) ^ (den < 0) else 1
    num, den = abs(num), abs(den)
    g = gcd(num, den)
    return (sign * (num // g), den // g)

def add_frac(a, b):
    an, ad = a; bn, bd = b
    return frac_simplify(an * bd + bn * ad, ad * bd)

def sub_frac(a, b):
    bn, bd = b
    return add_frac(a, (-bn, bd))

def mul_frac(a, b):
    an, ad = a; bn, bd = b
    return frac_simplify(an * bn, ad * bd)

def div_frac(a, b):
    bn, bd = b
    if bn == 0:
        raise ZeroDivisionError("Division by zero")
    return mul_frac(a, (bd, bn))

def neg_frac(a):
    an, ad = a
    return (-an, ad)

def pow_frac(a, n):
    if n == 0:
        return (1, 1)
    result = (1, 1)
    base = a if n > 0 else div_frac((1, 1), a)
    for _ in range(abs(n)):
        result = mul_frac(result, base)
    return result

def frac_from_int(n):
    return (int(n), 1)

def frac_zero():
    return (0, 1)

def frac_one():
    return (1, 1)

def frac_eq(a, b):
    an, ad = a; bn, bd = b
    return an * bd == bn * ad

def frac_is_zero(a):
    return a[0] == 0

def frac_to_float(a):
    return a[0] / a[1]

def parse_number(s):
    """Parse string into exact fraction tuple. Handles int, fraction, decimal."""
    s = s.strip()
    if not s:
        raise ValueError("Empty string")
    if '/' in s:
        parts = s.split('/')
        if len(parts) != 2:
            raise ValueError(f"Invalid fraction: {s}")
        return frac_simplify(int(parts[0].strip()), int(parts[1].strip()))
    elif '.' in s:
        neg = s.startswith('-')
        s2 = s.lstrip('-').lstrip('+')
        int_part, dec_part = s2.split('.', 1)
        dec_part = dec_part or '0'
        den = 10 ** len(dec_part)
        num = int(int_part or '0') * den + int(dec_part)
        return frac_simplify(-num if neg else num, den)
    else:
        return frac_simplify(int(s), 1)

def format_frac(x):
    """Format fraction as clean string."""
    n, d = x
    if n == 0:
        return "0"
    if d == 1:
        return str(n)
    return f"{n}/{d}"

def format_expr_term(coef, variable):
    """Format coefficient * variable cleanly."""
    n, d = coef
    if n == 0:
        return ""
    if n == d:
        return variable
    if n == -d:
        return f"-{variable}"
    return f"{format_frac(coef)}{variable}"

# ============================================================
# SECTION 1B: COMPLEX FRACTION ARITHMETIC
# A complex fraction is (real_frac, imag_frac), each a frac tuple.
# Used for finding eigenvectors of complex eigenvalues exactly.
# ============================================================

def cf_make(real_part=None, imag_part=None):
    if real_part is None: real_part = frac_zero()
    if imag_part is None: imag_part = frac_zero()
    if isinstance(real_part, int): real_part = frac_from_int(real_part)
    if isinstance(imag_part, int): imag_part = frac_from_int(imag_part)
    return (real_part, imag_part)

def cf_zero():
    return (frac_zero(), frac_zero())

def cf_one():
    return (frac_one(), frac_zero())

def cf_from_frac(f):
    return (f, frac_zero())

def cf_add(a, b):
    return (add_frac(a[0], b[0]), add_frac(a[1], b[1]))

def cf_sub(a, b):
    return (sub_frac(a[0], b[0]), sub_frac(a[1], b[1]))

def cf_mul(a, b):
    # (a + bi)(c + di) = (ac - bd) + (ad + bc)i
    real = sub_frac(mul_frac(a[0], b[0]), mul_frac(a[1], b[1]))
    imag = add_frac(mul_frac(a[0], b[1]), mul_frac(a[1], b[0]))
    return (real, imag)

def cf_div(a, b):
    # (a+bi)/(c+di) = ((ac+bd) + (bc-ad)i) / (c²+d²)
    denom = add_frac(mul_frac(b[0], b[0]), mul_frac(b[1], b[1]))
    if frac_is_zero(denom):
        raise ZeroDivisionError("Division by complex zero")
    num_r = add_frac(mul_frac(a[0], b[0]), mul_frac(a[1], b[1]))
    num_i = sub_frac(mul_frac(a[1], b[0]), mul_frac(a[0], b[1]))
    return (div_frac(num_r, denom), div_frac(num_i, denom))

def cf_neg(a):
    return (neg_frac(a[0]), neg_frac(a[1]))

def cf_is_zero(a):
    return frac_is_zero(a[0]) and frac_is_zero(a[1])

def cf_eq(a, b):
    return frac_eq(a[0], b[0]) and frac_eq(a[1], b[1])

def cf_format(a):
    """Format complex fraction as 'a + bi' or 'a - bi'."""
    re_s = format_frac(a[0])
    im_n = a[1][0]
    if im_n == 0:
        return re_s
    abs_im = (abs(im_n), a[1][1])
    abs_im_s = format_frac(abs_im)
    if a[0][0] == 0:
        prefix = "" if im_n > 0 else "-"
        return f"{prefix}i" if abs_im_s == "1" else f"{prefix}{abs_im_s}i"
    sign = " + " if im_n > 0 else " - "
    return f"{re_s}{sign}i" if abs_im_s == "1" else f"{re_s}{sign}{abs_im_s}i"


# ============================================================
# SECTION 2: SYMBOLIC EXPRESSION ARITHMETIC
# Expressions are lists of term dicts:
#   {'coef': frac, 'type': 'const'}
#   {'coef': frac, 'type': 'exp',     'k': frac}         -> coef*e^(k*t)
#   {'coef': frac, 'type': 'sin',     'k': frac}         -> coef*sin(k*t)
#   {'coef': frac, 'type': 'cos',     'k': frac}         -> coef*cos(k*t)
#   {'coef': frac, 'type': 'poly',    'n': int}           -> coef*t^n
#   {'coef': frac, 'type': 'exp_sin', 'k': frac, 'm': frac}
#   {'coef': frac, 'type': 'exp_cos', 'k': frac, 'm': frac}
# ============================================================

def _term_key(term):
    t = term['type']
    if t == 'const':  return ('const',)
    if t == 'exp':    return ('exp',     term['k'])
    if t == 'sin':    return ('sin',     term['k'])
    if t == 'cos':    return ('cos',     term['k'])
    if t == 'poly':   return ('poly',    term['n'])
    if t == 'exp_sin':return ('exp_sin', term['k'], term['m'])
    if t == 'exp_cos':return ('exp_cos', term['k'], term['m'])
    return (t,)

def expr_simplify(expr):
    combined = {}
    for term in expr:
        key = _term_key(term)
        if key in combined:
            new_coef = add_frac(combined[key]['coef'], term['coef'])
            combined[key] = {**term, 'coef': new_coef}
        else:
            combined[key] = dict(term)
    return [v for v in combined.values() if not frac_is_zero(v['coef'])]

def expr_add(a, b):
    return expr_simplify(a + b)

def expr_sub(a, b):
    neg_b = [{**t, 'coef': neg_frac(t['coef'])} for t in b]
    return expr_simplify(a + neg_b)

def expr_mul_const(c, expr):
    if frac_is_zero(c):
        return []
    return [{**t, 'coef': mul_frac(c, t['coef'])} for t in expr]

def expr_equal(a, b):
    diff = expr_sub(a, b)
    return len(expr_simplify(diff)) == 0

def make_const_expr(c):
    if frac_is_zero(c): return []
    return [{'coef': c, 'type': 'const'}]

def make_exp_expr(c, k):
    if frac_is_zero(c): return []
    if frac_is_zero(k): return make_const_expr(c)
    return [{'coef': c, 'type': 'exp', 'k': k}]

def make_sin_expr(c, k):
    if frac_is_zero(c): return []
    return [{'coef': c, 'type': 'sin', 'k': k}]

def make_cos_expr(c, k):
    if frac_is_zero(c): return []
    return [{'coef': c, 'type': 'cos', 'k': k}]

def derivative_expr(expr):
    """Differentiate expression list with respect to t."""
    result = []
    for term in expr:
        t, c = term['type'], term['coef']
        if t == 'const':
            pass
        elif t == 'exp':
            k = term['k']
            nc = mul_frac(c, k)
            if not frac_is_zero(nc):
                result.append({'coef': nc, 'type': 'exp', 'k': k})
        elif t == 'sin':
            k = term['k']
            nc = mul_frac(c, k)
            if not frac_is_zero(nc):
                result.append({'coef': nc, 'type': 'cos', 'k': k})
        elif t == 'cos':
            k = term['k']
            nc = neg_frac(mul_frac(c, k))
            if not frac_is_zero(nc):
                result.append({'coef': nc, 'type': 'sin', 'k': k})
        elif t == 'poly':
            n = term['n']
            if n > 0:
                nc = mul_frac(c, frac_from_int(n))
                result.append({'coef': nc, 'type': 'poly' if n > 1 else 'const',
                               **({'n': n - 1} if n > 1 else {})})
        elif t == 'exp_sin':
            k, m = term['k'], term['m']
            ck = mul_frac(c, k); cm = mul_frac(c, m)
            if not frac_is_zero(ck):
                result.append({'coef': ck, 'type': 'exp_sin', 'k': k, 'm': m})
            if not frac_is_zero(cm):
                result.append({'coef': cm, 'type': 'exp_cos', 'k': k, 'm': m})
        elif t == 'exp_cos':
            k, m = term['k'], term['m']
            ck = mul_frac(c, k); cm_neg = neg_frac(mul_frac(c, m))
            if not frac_is_zero(ck):
                result.append({'coef': ck, 'type': 'exp_cos', 'k': k, 'm': m})
            if not frac_is_zero(cm_neg):
                result.append({'coef': cm_neg, 'type': 'exp_sin', 'k': k, 'm': m})
    return expr_simplify(result)

def _format_coef_prefix(c):
    """Return coefficient string for inline use (suppress 1/-1)."""
    n, d = c
    if n == d:   return ""
    if n == -d:  return "-"
    return format_frac(c)

def format_term(term):
    t, c = term['type'], term['coef']
    cs = format_frac(c)
    cp = _format_coef_prefix(c)

    if t == 'const':
        return cs
    elif t == 'exp':
        k = term['k']; ks = format_frac(k)
        exp_str = f"e^({ks}t)" if ks not in ('1', '0') else ("e^t" if ks == '1' else '1')
        if frac_is_zero(k): return cs
        return f"{cp}{exp_str}"
    elif t == 'sin':
        k = term['k']; ks = format_frac(k)
        inner = "t" if ks == "1" else f"{ks}t"
        return f"{cp}sin({inner})"
    elif t == 'cos':
        k = term['k']; ks = format_frac(k)
        inner = "t" if ks == "1" else f"{ks}t"
        return f"{cp}cos({inner})"
    elif t == 'poly':
        n = term['n']
        if n == 0: return cs
        if n == 1: return f"{cp}t"
        return f"{cp}t^{n}"
    elif t == 'exp_sin':
        k, m = term['k'], term['m']
        return f"{cp}e^({format_frac(k)}t)sin({format_frac(m)}t)"
    elif t == 'exp_cos':
        k, m = term['k'], term['m']
        return f"{cp}e^({format_frac(k)}t)cos({format_frac(m)}t)"
    return cs

def format_expr(expr):
    if not expr: return "0"
    parts = []
    for i, term in enumerate(expr):
        s = format_term(term)
        if i == 0:
            parts.append(s)
        else:
            if s.startswith('-'):
                parts.append(f" - {s[1:]}")
            else:
                parts.append(f" + {s}")
    return "".join(parts)

# ============================================================
# SECTION 3: MATRIX AND VECTOR ARITHMETIC
# Matrices: list of lists of fractions
# Vectors:  list of fractions
# ============================================================

def zero_matrix(rows, cols):
    return [[frac_zero() for _ in range(cols)] for _ in range(rows)]

def identity_matrix(n):
    M = zero_matrix(n, n)
    for i in range(n): M[i][i] = frac_one()
    return M

def matrix_add(A, B):
    r, c = len(A), len(A[0])
    return [[add_frac(A[i][j], B[i][j]) for j in range(c)] for i in range(r)]

def matrix_sub(A, B):
    r, c = len(A), len(A[0])
    return [[sub_frac(A[i][j], B[i][j]) for j in range(c)] for i in range(r)]

def matrix_mul(A, B):
    ra, ca, cb = len(A), len(A[0]), len(B[0])
    C = zero_matrix(ra, cb)
    for i in range(ra):
        for j in range(cb):
            s = frac_zero()
            for k in range(ca):
                s = add_frac(s, mul_frac(A[i][k], B[k][j]))
            C[i][j] = s
    return C

def matrix_vector_mul(A, v):
    rows = len(A)
    result = []
    for i in range(rows):
        s = frac_zero()
        for j in range(len(v)):
            s = add_frac(s, mul_frac(A[i][j], v[j]))
        result.append(s)
    return result

def scalar_matrix_mul(c, A):
    return [[mul_frac(c, A[i][j]) for j in range(len(A[0]))] for i in range(len(A))]

def scalar_vector_mul(c, v):
    return [mul_frac(c, x) for x in v]

def vector_add(u, v):
    return [add_frac(u[i], v[i]) for i in range(len(u))]

def vector_sub(u, v):
    return [sub_frac(u[i], v[i]) for i in range(len(u))]

def zero_vector(n):
    return [frac_zero() for _ in range(n)]

# ============================================================
# SECTION 4: DETERMINANTS AND CHARACTERISTIC POLYNOMIALS
# ============================================================

def determinant_2x2(A):
    return sub_frac(mul_frac(A[0][0], A[1][1]), mul_frac(A[0][1], A[1][0]))

def determinant_3x3(A):
    a,b,c = A[0]; d,e,f = A[1]; g,h,k = A[2]
    t1 = mul_frac(a, sub_frac(mul_frac(e, k), mul_frac(f, h)))
    t2 = mul_frac(b, sub_frac(mul_frac(d, k), mul_frac(f, g)))
    t3 = mul_frac(c, sub_frac(mul_frac(d, h), mul_frac(e, g)))
    return sub_frac(add_frac(t1, t3), t2)

def determinant(A):
    n = len(A)
    if n == 1: return A[0][0]
    if n == 2: return determinant_2x2(A)
    if n == 3: return determinant_3x3(A)
    det = frac_zero()
    for j in range(n):
        minor = [[A[i][kk] for kk in range(n) if kk != j] for i in range(1, n)]
        sign = frac_from_int(1 if j % 2 == 0 else -1)
        det = add_frac(det, mul_frac(mul_frac(sign, A[0][j]), determinant(minor)))
    return det

def trace_2x2(A):
    return add_frac(A[0][0], A[1][1])

def characteristic_polynomial_2x2(A):
    """Returns [a,b,c] for a*λ² + b*λ + c = 0"""
    tr = add_frac(A[0][0], A[1][1])
    det = determinant_2x2(A)
    return [frac_one(), neg_frac(tr), det]

def characteristic_polynomial_3x3(A):
    """Returns [a,b,c,d] for a*λ³ + b*λ² + c*λ + d = 0"""
    tr = add_frac(add_frac(A[0][0], A[1][1]), A[2][2])
    m00 = sub_frac(mul_frac(A[1][1], A[2][2]), mul_frac(A[1][2], A[2][1]))
    m11 = sub_frac(mul_frac(A[0][0], A[2][2]), mul_frac(A[0][2], A[2][0]))
    m22 = sub_frac(mul_frac(A[0][0], A[1][1]), mul_frac(A[0][1], A[1][0]))
    sum_minors = add_frac(add_frac(m00, m11), m22)
    det = determinant_3x3(A)
    return [frac_one(), neg_frac(tr), sum_minors, neg_frac(det)]

# ============================================================
# SECTION 5: POLYNOMIAL ROOT FINDING
# ============================================================

def _int_sqrt(n):
    if n < 0: return None
    s = int(math.isqrt(n))
    for c in [s - 1, s, s + 1]:
        if c >= 0 and c * c == n:
            return c
    return None

def quadratic_roots_exact(a, b, c):
    """Solve a*x² + b*x + c = 0 exactly. Returns list of roots (fracs or complex tuples)."""
    b2 = mul_frac(b, b)
    four_ac = mul_frac(mul_frac(frac_from_int(4), a), c)
    disc = sub_frac(b2, four_ac)
    disc_n, disc_d = disc

    if disc_n == 0:
        root = div_frac(neg_frac(b), mul_frac(frac_from_int(2), a))
        return [root, root]

    two_a = mul_frac(frac_from_int(2), a)
    neg_b = neg_frac(b)

    if disc_n > 0:
        # Positive discriminant: check for perfect rational square
        # sqrt(disc_n/disc_d) = sqrt(disc_n * disc_d) / disc_d
        product = disc_n * disc_d
        sqrt_prod = _int_sqrt(product)
        if sqrt_prod is not None:
            sqrt_disc = frac_simplify(sqrt_prod, disc_d)
            r1 = div_frac(add_frac(neg_b, sqrt_disc), two_a)
            r2 = div_frac(sub_frac(neg_b, sqrt_disc), two_a)
            return [r1, r2]
        else:
            # Irrational: approximate
            disc_f = frac_to_float(disc)
            sqrt_disc_f = math.sqrt(disc_f)
            af = frac_to_float(a)
            bf = frac_to_float(b)
            r1_f = (-bf + sqrt_disc_f) / (2 * af)
            r2_f = (-bf - sqrt_disc_f) / (2 * af)
            def approx(f):
                for den in range(1, 200):
                    num = round(f * den)
                    if abs(num / den - f) < 1e-9:
                        return frac_simplify(num, den)
                return frac_simplify(round(f * 100), 100)
            return [approx(r1_f), approx(r2_f)]
    else:
        # Negative discriminant: complex roots α ± βi
        # α = -b/(2a)
        alpha = div_frac(neg_b, two_a)
        # β = sqrt(|disc|) / (2|a|)
        # |disc| = -disc_n / disc_d (disc_n < 0, disc_d > 0)
        product = -disc_n * disc_d  # positive integer
        sqrt_prod = _int_sqrt(product)
        if sqrt_prod is not None:
            sqrt_abs_disc = frac_simplify(sqrt_prod, disc_d)
            beta = div_frac(sqrt_abs_disc, two_a)
            # Force beta positive (we'll add ±)
            if beta[0] < 0:
                beta = neg_frac(beta)
            return [('complex', alpha, beta),
                    ('complex', alpha, neg_frac(beta))]
        else:
            # Irrational imaginary part — fall back to floats
            af = frac_to_float(a)
            bf = frac_to_float(b)
            disc_f = frac_to_float(disc)
            real_f = frac_to_float(alpha)
            imag_f = math.sqrt(abs(disc_f)) / (2 * abs(af))
            return [('complex_f', real_f, imag_f),
                    ('complex_f', real_f, -imag_f)]

def eval_polynomial(coeffs, x):
    """Evaluate polynomial at x (frac). coeffs = [leading, ..., constant]."""
    val = frac_zero()
    for c in coeffs:
        val = add_frac(mul_frac(val, x), c)
    return val

def divide_polynomial(coeffs, root):
    """Synthetic division of polynomial by (x - root)."""
    result = []
    carry = frac_zero()
    for c in coeffs:
        carry = add_frac(mul_frac(carry, root), c)
        result.append(carry)
    return result[:-1]  # drop remainder

def integer_roots_of_polynomial(coeffs):
    """Find rational roots using rational root theorem."""
    # Multiply by LCM of denominators to get integer coefficients
    dens = [c[1] for c in coeffs]
    lcm_den = dens[0]
    for d in dens[1:]:
        lcm_den = lcm_den * d // gcd(lcm_den, d)
    int_coeffs = []
    for c in coeffs:
        val = mul_frac(c, frac_from_int(lcm_den))
        n, d = val
        if d != 1:
            return []  # Can't make integer
        int_coeffs.append(n)
    const = abs(int_coeffs[-1]) if int_coeffs[-1] != 0 else 1
    lead  = abs(int_coeffs[0])
    fc = [i for i in range(1, const + 1) if const % i == 0]
    fl = [i for i in range(1, lead  + 1) if lead  % i == 0]
    candidates = set()
    for p in fc:
        for q in fl:
            candidates.add((p, q))
            candidates.add((-p, q))
    roots = []
    for p, q in candidates:
        x = frac_simplify(p, q)
        if frac_is_zero(eval_polynomial(coeffs, x)):
            if x not in roots:
                roots.append(x)
    return roots

def factor_polynomial_roots(coeffs):
    """Find all rational eigenvalues of a polynomial."""
    roots = []
    remaining = list(coeffs)
    while len(remaining) > 1:
        deg = len(remaining) - 1
        if deg == 0:
            break
        if deg == 1:
            root = div_frac(neg_frac(remaining[1]), remaining[0])
            roots.append(root)
            break
        if deg == 2:
            quad = quadratic_roots_exact(remaining[0], remaining[1], remaining[2])
            roots.extend(quad)
            break
        int_roots = integer_roots_of_polynomial(remaining)
        if not int_roots:
            # Try quadratic if we have 2 remaining
            break
        root = int_roots[0]
        roots.append(root)
        remaining = divide_polynomial(remaining, root)
    return roots

# ============================================================
# SECTION 6: ROW REDUCTION AND LINEAR SYSTEM SOLVING
# ============================================================

def row_reduce(M):
    """Reduced row echelon form using exact fraction arithmetic."""
    M = [list(row) for row in M]
    rows, cols = len(M), len(M[0])
    pivot_row = 0
    for col in range(cols - 1):
        found = next((r for r in range(pivot_row, rows) if not frac_is_zero(M[r][col])), -1)
        if found == -1:
            continue
        M[pivot_row], M[found] = M[found], M[pivot_row]
        pivot = M[pivot_row][col]
        M[pivot_row] = [div_frac(x, pivot) for x in M[pivot_row]]
        for r in range(rows):
            if r != pivot_row and not frac_is_zero(M[r][col]):
                factor = M[r][col]
                M[r] = [sub_frac(M[r][j], mul_frac(factor, M[pivot_row][j])) for j in range(cols)]
        pivot_row += 1
    return M

def solve_linear_system(A, b):
    """
    Solve Ax = b.
    Returns ('unique', vector) | ('free', (particular, null_basis)) | ('no_solution', None)
    """
    rows, cols = len(A), len(A[0])
    aug = [A[i] + [b[i]] for i in range(rows)]
    rref = row_reduce(aug)
    for row in rref:
        if all(frac_is_zero(row[j]) for j in range(cols)) and not frac_is_zero(row[cols]):
            return ('no_solution', None)
    pivot_cols = {}
    for r in range(rows):
        for c in range(cols):
            if not frac_is_zero(rref[r][c]):
                pivot_cols[c] = r
                break
    free_cols = [c for c in range(cols) if c not in pivot_cols]
    if not free_cols:
        solution = [frac_zero()] * cols
        for c, r in pivot_cols.items():
            solution[c] = rref[r][cols]
        return ('unique', solution)
    else:
        particular = [frac_zero()] * cols
        for c, r in pivot_cols.items():
            particular[c] = rref[r][cols]
        null_basis = []
        for fc in free_cols:
            vec = [frac_zero()] * cols
            vec[fc] = frac_one()
            for pc, r in pivot_cols.items():
                vec[pc] = neg_frac(rref[r][fc])
            null_basis.append(vec)
        return ('free', (particular, null_basis))

def nullspace(A):
    """Basis for the nullspace of A."""
    rows, cols = len(A), len(A[0])
    b = [frac_zero()] * rows
    stype, result = solve_linear_system(A, b)
    if stype == 'unique':
        return []
    elif stype == 'free':
        _, basis = result
        return basis
    return []

def row_reduce_complex(M):
    """Reduced row echelon form over the complex rationals."""
    M = [list(row) for row in M]
    rows, cols = len(M), len(M[0])
    pivot_row = 0
    for col in range(cols - 1):
        found = next((r for r in range(pivot_row, rows) if not cf_is_zero(M[r][col])), -1)
        if found == -1:
            continue
        M[pivot_row], M[found] = M[found], M[pivot_row]
        pivot = M[pivot_row][col]
        M[pivot_row] = [cf_div(x, pivot) for x in M[pivot_row]]
        for r in range(rows):
            if r != pivot_row and not cf_is_zero(M[r][col]):
                factor = M[r][col]
                M[r] = [cf_sub(M[r][j], cf_mul(factor, M[pivot_row][j])) for j in range(cols)]
        pivot_row += 1
    return M

def nullspace_complex(A):
    """Basis for the nullspace of a complex-fraction matrix."""
    rows, cols = len(A), len(A[0])
    aug = [row + [cf_zero()] for row in A]
    rref = row_reduce_complex(aug)
    pivot_cols = {}
    for r in range(rows):
        for c in range(cols):
            if not cf_is_zero(rref[r][c]):
                pivot_cols[c] = r
                break
    free_cols = [c for c in range(cols) if c not in pivot_cols]
    basis = []
    for fc in free_cols:
        vec = [cf_zero()] * cols
        vec[fc] = cf_one()
        for pc, r in pivot_cols.items():
            vec[pc] = cf_neg(rref[r][fc])
        basis.append(vec)
    return basis

def complex_eigenvector(A, alpha, beta):
    """
    Find a complex eigenvector for the eigenvalue λ = α + βi (α, β fracs).
    Returns the complex vector V = a + bi as a list of (real, imag) frac pairs,
    or None if no eigenvector found.
    """
    n = len(A)
    lam = (alpha, beta)
    M = []
    for i in range(n):
        row = []
        for j in range(n):
            elem = cf_from_frac(A[i][j])
            if i == j:
                elem = cf_sub(elem, lam)
            row.append(elem)
        M.append(row)
    ns = nullspace_complex(M)
    return ns[0] if ns else None

def split_complex_vector(v):
    """Split complex vector V = a + bi into (a, b) with frac components."""
    a = [comp[0] for comp in v]
    b = [comp[1] for comp in v]
    return a, b

def complex_eigenvector_float(A, alpha, beta):
    """
    Compute a float complex eigenvector for eigenvalue alpha + beta*i.
    Returns list of (real_float, imag_float) pairs, or None.
    Uses Python's built-in complex arithmetic.
    """
    n = len(A)
    lam = complex(alpha, beta)
    # Build (A - λI)
    M = [[complex(frac_to_float(A[i][j])) - (lam if i == j else 0)
          for j in range(n)] for i in range(n)]
    # Row reduction over complex floats
    pivot_row = 0
    for col in range(n):
        found = next((r for r in range(pivot_row, n) if abs(M[r][col]) > 1e-10), -1)
        if found == -1:
            continue
        M[pivot_row], M[found] = M[found], M[pivot_row]
        piv = M[pivot_row][col]
        M[pivot_row] = [x / piv for x in M[pivot_row]]
        for r in range(n):
            if r != pivot_row and abs(M[r][col]) > 1e-10:
                f = M[r][col]
                M[r] = [M[r][j] - f * M[pivot_row][j] for j in range(n)]
        pivot_row += 1
    # Find free column
    pivot_cols = {}
    for r in range(n):
        for c in range(n):
            if abs(M[r][c]) > 0.5:
                pivot_cols[c] = r
                break
    free_cols = [c for c in range(n) if c not in pivot_cols]
    if not free_cols:
        return None
    fc = free_cols[0]
    vec = [complex(0)] * n
    vec[fc] = complex(1)
    for pc, r in pivot_cols.items():
        vec[pc] = -M[r][fc]
    return [(v.real, v.imag) for v in vec]

def _fmt_float_vec(v):
    """Format a list of (real, imag) float pairs as a vector string."""
    parts = []
    for r, im in v:
        r = round(r, 6); im = round(im, 6)
        if abs(im) < 1e-9:
            parts.append(f"{r:.4g}")
        elif abs(r) < 1e-9:
            parts.append(f"{im:.4g}i")
        else:
            sign = "+" if im >= 0 else "-"
            parts.append(f"{r:.4g} {sign} {abs(im):.4g}i")
    return "[" + ", ".join(parts) + "]"

def _split_float_vec(v):
    """Split float complex vector into (a_real, b_imag) as float lists."""
    a = [round(comp[0], 10) for comp in v]
    b = [round(comp[1], 10) for comp in v]
    return a, b

def _fmt_float_rv(v):
    return "[" + ", ".join(f"{x:.4g}" for x in v) + "]"


# ============================================================
# SECTION 7: EIGENVALUES AND EIGENVECTORS
# ============================================================

def eigenvalues(A):
    """Eigenvalues of 2x2 or 3x3 matrix."""
    n = len(A)
    if n == 2:
        return quadratic_roots_exact(*characteristic_polynomial_2x2(A))
    elif n == 3:
        return factor_polynomial_roots(characteristic_polynomial_3x3(A))
    raise ValueError("Only 2x2 and 3x3 supported")

def eigenvectors(A, lambdas):
    """List of (lambda, eigenvector) for each eigenvalue.
    For real eigenvalues: vector is a list of fracs.
    For 'complex' eigenvalues (rational α,β): vector is a list of (real,imag) frac pairs.
    For 'complex_f' eigenvalues (irrational): vector is None.
    """
    n = len(A)
    I = identity_matrix(n)
    pairs = []
    for lam in lambdas:
        if isinstance(lam, tuple) and lam[0] == 'complex_f':
            _, alpha, beta = lam
            cvec = complex_eigenvector_float(A, alpha, beta)
            pairs.append((lam, cvec))
            continue
        if isinstance(lam, tuple) and lam[0] == 'complex':
            _, alpha, beta = lam
            cvec = complex_eigenvector(A, alpha, beta)
            pairs.append((lam, cvec))
            continue
        # Real eigenvalue
        lam_I = scalar_matrix_mul(lam, I)
        B = matrix_sub(A, lam_I)
        ns = nullspace(B)
        pairs.append((lam, ns[0] if ns else zero_vector(n)))
    return pairs

def sort_eigenpairs(epairs):
    """Sort eigenvalues smallest to largest. Real first, then complex pairs (+β before -β)."""
    def key(pair):
        lam = pair[0]
        if isinstance(lam, tuple):
            if lam[0] == 'complex':
                return (1, frac_to_float(lam[1]), -frac_to_float(lam[2]))
            elif lam[0] == 'complex_f':
                return (1, lam[1], -lam[2])
        return (0, frac_to_float(lam))
    return sorted(epairs, key=key)

def generalized_eigenvector(A, lam, v):
    """Solve (A - λI)w = v for generalized eigenvector."""
    n = len(A)
    I = identity_matrix(n)
    B = matrix_sub(A, scalar_matrix_mul(lam, I))
    stype, result = solve_linear_system(B, v)
    if stype == 'unique':
        return result
    elif stype == 'free':
        particular, _ = result
        return particular
    return zero_vector(n)

# ============================================================
# SECTION 8: FORMATTING
# ============================================================

def format_vector(v):
    """Format as [a, b, c]."""
    if v is None: return "[?]"
    return "[" + ", ".join(format_frac(x) for x in v) + "]"

def format_matrix(A):
    """Aligned matrix display."""
    cols = len(A[0])
    widths = [max(len(format_frac(A[i][j])) for i in range(len(A))) for j in range(cols)]
    lines = []
    for row in A:
        inner = "  ".join(format_frac(row[j]).rjust(widths[j]) for j in range(cols))
        lines.append(f"  [ {inner} ]")
    return "\n".join(lines)

def _exp_str(lam):
    """Return e^(λt) string; empty string if λ=0."""
    if frac_is_zero(lam): return ""
    ls = format_frac(lam)
    if ls == "1":  return "e^t"
    if ls == "-1": return "e^(-t)"
    return f"e^({ls}t)"

def format_complex_eigenvalue(lam):
    """Format a complex eigenvalue tuple as 'α + βi' or 'α - βi'."""
    if not _is_complex_lam(lam):
        return format_frac(lam)
    if lam[0] == 'complex':
        _, alpha, beta = lam
        return cf_format((alpha, beta))
    elif lam[0] == 'complex_f':
        _, a, b = lam
        sign = "+" if b >= 0 else "-"
        return f"{a:.4f} {sign} {abs(b):.4f}i"
    return str(lam)

def format_complex_vector(v):
    """Format a complex vector (list of (real,imag) frac pairs)."""
    if v is None: return "[?]"
    return "[" + ", ".join(cf_format(c) for c in v) + "]"

def format_complex_pair_solution(c1_label, c2_label, alpha, beta, a, b):
    """
    Format the real solution pair from a complex eigenvalue:
      C1·e^(αt)·([a]cos(βt) - [b]sin(βt)) + C2·e^(αt)·([a]sin(βt) + [b]cos(βt))
    """
    a_str = format_vector(a)
    b_str = format_vector(b)
    bs = format_frac(beta)
    arg = "t" if bs == "1" else f"{bs}t"
    inner1 = f"({a_str}cos({arg}) - {b_str}sin({arg}))"
    inner2 = f"({a_str}sin({arg}) + {b_str}cos({arg}))"
    exp_str = _exp_str(alpha)
    return f"{c1_label}{exp_str}{inner1} + {c2_label}{exp_str}{inner2}"

def format_complex_pair_specific(alpha, beta, cos_vec, sin_vec):
    """
    Format the specific (no-constants) version of a complex-pair contribution
    after IC fitting:
      e^(αt)·([cos_vec]cos(βt) + [sin_vec]sin(βt))
    """
    bs = format_frac(beta)
    arg = "t" if bs == "1" else f"{bs}t"
    cos_zero = all(frac_is_zero(x) for x in cos_vec)
    sin_zero = all(frac_is_zero(x) for x in sin_vec)
    if cos_zero and sin_zero:
        return ""
    parts = []
    if not cos_zero:
        parts.append(f"{format_vector(cos_vec)}cos({arg})")
    if not sin_zero:
        parts.append(f"{format_vector(sin_vec)}sin({arg})")
    inner = " + ".join(parts)
    if len(parts) > 1:
        inner = f"({inner})"
    exp_str = _exp_str(alpha)
    return f"{exp_str}{inner}" if exp_str else inner

def _is_complex_lam(lam):
    return isinstance(lam, tuple) and lam[0] in ('complex', 'complex_f')

def _find_conjugate_index(epairs, i, alpha, beta):
    """Locate index of conjugate eigenvalue α - βi in epairs (returns -1 if absent)."""
    for j in range(len(epairs)):
        if j == i: continue
        lam_j = epairs[j][0]
        if isinstance(lam_j, tuple) and lam_j[0] == 'complex':
            if frac_eq(lam_j[1], alpha) and frac_eq(lam_j[2], neg_frac(beta)):
                return j
    return -1


def format_general_solution(terms):
    """terms = list of strings; returns X(t) = term1 + term2 + ..."""
    return "X(t) = " + " + ".join(terms)

def print_steps(title, steps):
    """Print step-by-step solution."""
    print(f"\n{'='*62}")
    print(f"  {title}")
    print('='*62)
    for step in steps:
        print(step)
    print('='*62)

def copy_paste_answer_only(result):
    """Print clean final answer."""
    print(f"\n{'─'*55}")
    print("  ANSWER:")
    print(f"  {result}")
    print('─'*55)

# ============================================================
# SECTION 9: HOMOGENEOUS SYSTEM SOLVER
# ============================================================

def _get_eigenpairs(A):
    """Compute and sort eigenpairs (no printing)."""
    evals = eigenvalues(A)
    epairs = eigenvectors(A, evals)
    return sort_eigenpairs(epairs)

def _detect_defect(A, epairs):
    """Check for defective (repeated) eigenvalues."""
    # Group by eigenvalue float
    groups = {}
    for lam, vec in epairs:
        if _is_complex_lam(lam): continue
        key = round(frac_to_float(lam), 8)
        groups.setdefault(key, []).append((lam, vec))
    for key, pairs in groups.items():
        if len(pairs) > 1:
            lam = pairs[0][0]
            n = len(A)
            I = identity_matrix(n)
            B = matrix_sub(A, scalar_matrix_mul(lam, I))
            ns = nullspace(B)
            if len(ns) < len(pairs):
                return True, lam, ns[0] if ns else zero_vector(n)
    return False, None, None

def solve_homogeneous_system(A, verbose=True):
    """
    Solve X' = AX.
    Returns (solution_string, sorted_epairs).
    Handles real eigenvalues, complex conjugate pairs (rational α,β), and defective repeated eigenvalues.
    """
    n = len(A)
    epairs = _get_eigenpairs(A)
    steps = [f"Matrix A:\n{format_matrix(A)}"]

    defect, def_lam, def_vec = _detect_defect(A, epairs)
    if defect:
        return solve_repeated_eigenvalue_system(A, verbose=verbose)

    # Char poly display
    if n == 2:
        coeffs = characteristic_polynomial_2x2(A)
        tr = format_frac(add_frac(A[0][0], A[1][1]))
        det = format_frac(determinant_2x2(A))
        steps.append(f"Characteristic polynomial: λ² - ({tr})λ + ({det}) = 0")
    else:
        coeffs = characteristic_polynomial_3x3(A)
        cstrs = [format_frac(c) for c in coeffs]
        steps.append(f"Characteristic polynomial:\n  ({cstrs[0]})λ³ + ({cstrs[1]})λ² + ({cstrs[2]})λ + ({cstrs[3]}) = 0")

    steps.append("Eigenvalues and eigenvectors:")
    for i, (lam, vec) in enumerate(epairs):
        if isinstance(lam, tuple) and lam[0] == 'complex':
            steps.append(f"  λ{i+1} = {format_complex_eigenvalue(lam)}")
            if vec is not None:
                steps.append(f"     K{i+1} = {format_complex_vector(vec)}")
                a_part, b_part = split_complex_vector(vec)
                steps.append(f"           = {format_vector(a_part)} + i·{format_vector(b_part)}")
        elif isinstance(lam, tuple) and lam[0] == 'complex_f':
            steps.append(f"  λ{i+1} ≈ {lam[1]:.4f} + {lam[2]:.4f}i  (irrational imaginary part)")
            if vec is not None:
                steps.append(f"     K{i+1} ≈ {_fmt_float_vec(vec)}")
        else:
            steps.append(f"  λ{i+1} = {format_frac(lam)},  K{i+1} = {format_vector(vec)}")

    # Build solution: real eigenvalues first, then complex pairs.
    term_parts = []
    constant_idx = 1

    # Real eigenvalues
    for lam, vec in epairs:
        if _is_complex_lam(lam): continue
        if vec is None: continue
        exp = _exp_str(lam)
        term_parts.append(
            f"C{constant_idx}{format_vector(vec)}{exp}" if exp
            else f"C{constant_idx}{format_vector(vec)}"
        )
        constant_idx += 1

    # Complex eigenvalue pairs
    processed = set()
    for i, (lam, vec) in enumerate(epairs):
        if i in processed: continue
        if not _is_complex_lam(lam): continue
        if vec is None:
            processed.add(i); continue

        if lam[0] == 'complex':
            _, alpha, beta = lam
            a_part, b_part = split_complex_vector(vec)
            conj = _find_conjugate_index(epairs, i, alpha, beta)
            if conj >= 0: processed.add(conj)
            processed.add(i)
            sol_str = format_complex_pair_solution(
                f"C{constant_idx}", f"C{constant_idx+1}",
                alpha, beta, a_part, b_part
            )
            steps.append(f"  Complex pair λ = {format_complex_eigenvalue(lam)}:")
            steps.append(f"    a = {format_vector(a_part)},  b = {format_vector(b_part)}")
        else:  # complex_f — float approximation
            _, alpha, beta = lam
            a_part, b_part = _split_float_vec(vec)
            # Find conjugate
            for j in range(len(epairs)):
                if j == i or j in processed: continue
                lj = epairs[j][0]
                if isinstance(lj, tuple) and lj[0] == 'complex_f':
                    if abs(lj[1] - alpha) < 1e-9 and abs(lj[2] + beta) < 1e-9:
                        processed.add(j); break
            processed.add(i)
            af = f"{alpha:.4g}"
            bf = f"{abs(beta):.4g}"
            arg = "t" if bf == "1" else f"{bf}t"
            exp_str = "" if abs(alpha) < 1e-9 else (f"e^({af}t)" if abs(alpha - 1) > 1e-9 else "e^t")
            a_str = _fmt_float_rv(a_part)
            b_str = _fmt_float_rv(b_part)
            sol_str = (f"C{constant_idx}{exp_str}({a_str}cos({arg}) - {b_str}sin({arg})) + "
                       f"C{constant_idx+1}{exp_str}({a_str}sin({arg}) + {b_str}cos({arg}))")
            steps.append(f"  Complex pair λ ≈ {alpha:.4g} ± {abs(beta):.4g}i  (approx):")
            steps.append(f"    a ≈ {a_str},  b ≈ {b_str}")

        term_parts.append(sol_str)
        constant_idx += 2

    sol = "X(t) = " + " + ".join(term_parts) if term_parts else "X(t) = 0"
    steps.append(f"General solution:\n  {sol}")

    if verbose:
        print_steps("Homogeneous System  X' = AX", steps)
    return sol, epairs

def solve_repeated_eigenvalue_system(A, verbose=True):
    """Handle repeated (possibly defective) eigenvalues."""
    n = len(A)
    steps = [f"Matrix A:\n{format_matrix(A)}"]

    epairs = _get_eigenpairs(A)

    # Group eigenvalues
    groups = {}
    for lam, vec in epairs:
        if _is_complex_lam(lam): continue
        key = round(frac_to_float(lam), 8)
        groups.setdefault(key, {'lam': lam, 'alg': 0, 'geom_vecs': []})
        groups[key]['alg'] += 1

    # Find geometric multiplicities
    I = identity_matrix(n)
    for key in groups:
        lam = groups[key]['lam']
        B = matrix_sub(A, scalar_matrix_mul(lam, I))
        ns = nullspace(B)
        groups[key]['geom_vecs'] = ns
        groups[key]['geom'] = len(ns)

    term_parts = []
    ci = 1
    for key, g in sorted(groups.items()):
        lam = g['lam']
        alg  = g['alg']
        ns   = g['geom_vecs']
        geom = g['geom']
        ls = format_frac(lam)
        exp = _exp_str(lam)

        steps.append(f"Eigenvalue λ = {ls}  (algebraic mult={alg}, geometric mult={geom})")

        if geom == alg:
            # Enough eigenvectors
            for v in ns:
                steps.append(f"  Eigenvector: {format_vector(v)}")
                part = f"C{ci}{format_vector(v)}{exp}" if exp else f"C{ci}{format_vector(v)}"
                term_parts.append(part)
                ci += 1
        else:
            # Defective: need generalized eigenvectors
            v = ns[0] if ns else zero_vector(n)
            steps.append(f"  Eigenvector v = {format_vector(v)}")
            w = generalized_eigenvector(A, lam, v)
            steps.append(f"  Generalized eigenvector w (solve (A-λI)w = v): {format_vector(w)}")
            if exp:
                part1 = f"C{ci}{format_vector(v)}{exp}"
                part2 = f"C{ci+1}(t{format_vector(v)} + {format_vector(w)}){exp}"
            else:
                part1 = f"C{ci}{format_vector(v)}"
                part2 = f"C{ci+1}(t{format_vector(v)} + {format_vector(w)})"
            term_parts.append(part1)
            term_parts.append(part2)
            ci += 2

    sol = "X(t) = " + " + ".join(term_parts)
    steps.append(f"General solution:\n  {sol}")

    if verbose:
        print_steps("Repeated Eigenvalue System", steps)
    return sol, epairs

def solve_constants_from_initial_conditions(vectors, X0):
    """Solve C1*v1 + C2*v2 + ... = X0 for constants."""
    n = len(X0)
    m = len(vectors)
    A = [[vectors[j][i] for j in range(m)] for i in range(n)]
    stype, result = solve_linear_system(A, X0)
    if stype == 'unique':
        return result
    elif stype == 'free':
        particular, _ = result
        return particular
    return [frac_zero()] * m

def solve_ivp_homogeneous(A, X0, verbose=True):
    """Solve X' = AX, X(0) = X0. Handles real and complex eigenvalues."""
    n = len(A)
    steps = [
        f"Matrix A:\n{format_matrix(A)}",
        f"Initial condition: X(0) = {format_vector(X0)}",
    ]

    defect, _, _ = _detect_defect(A, _get_eigenpairs(A))
    if defect:
        sol_str, _ = solve_repeated_eigenvalue_system(A, verbose=False)
        steps.append(f"General homogeneous solution:\n  {sol_str}")
        steps.append("⚠  Defective repeated eigenvalue case — returning general form.")
        if verbose:
            print_steps("Initial Value Problem  X' = AX,  X(0) = X0", steps)
        return sol_str

    epairs = sort_eigenpairs(_get_eigenpairs(A))

    # Build basis for fitting initial condition.
    # Each entry: ('real', lam, vec, ic_vector)  — contributes C·v at t=0
    # or         ('complex', alpha, beta, a, b, which) where which='a' or 'b'
    #            'a' contributes C·a at t=0, 'b' contributes C·b at t=0
    basis = []  # in order matching constants C1, C2, ...

    # Real eigenvalues first
    for lam, vec in epairs:
        if _is_complex_lam(lam): continue
        if vec is None: continue
        basis.append(('real', lam, vec, vec))

    # Complex pairs
    processed = set()
    for i, (lam, vec) in enumerate(epairs):
        if i in processed: continue
        if not (isinstance(lam, tuple) and lam[0] == 'complex'): continue
        if vec is None:
            processed.add(i); continue
        _, alpha, beta = lam
        a_part, b_part = split_complex_vector(vec)
        conj = _find_conjugate_index(epairs, i, alpha, beta)
        if conj >= 0: processed.add(conj)
        processed.add(i)
        # Two real solutions: at t=0 the first contributes a (cos(0)=1, sin(0)=0)
        #                     and the second contributes b
        basis.append(('complex', alpha, beta, a_part, b_part, 'a'))
        basis.append(('complex', alpha, beta, a_part, b_part, 'b'))

    if not basis:
        steps.append("⚠  No usable eigenvectors found.")
        if verbose:
            print_steps("Initial Value Problem", steps)
        return "X(t) = 0"

    steps.append("Eigenvalues and eigenvectors:")
    seen_pair = set()
    for entry in basis:
        if entry[0] == 'real':
            _, lam, vec, _ = entry
            steps.append(f"  λ = {format_frac(lam)},  K = {format_vector(vec)}")
        else:
            _, alpha, beta, a, b, which = entry
            key = (alpha, beta)
            if key in seen_pair: continue
            seen_pair.add(key)
            steps.append(f"  λ = {format_complex_eigenvalue(('complex', alpha, beta))}")
            steps.append(f"    a = {format_vector(a)},  b = {format_vector(b)}")

    # Solve for constants: stack ic_vectors as columns of M, solve M·C = X0
    column_vecs = []
    for entry in basis:
        if entry[0] == 'real':
            column_vecs.append(entry[3])
        else:
            _, alpha, beta, a, b, which = entry
            column_vecs.append(a if which == 'a' else b)

    M = [[column_vecs[j][i] for j in range(len(column_vecs))] for i in range(n)]
    stype, result = solve_linear_system(M, X0)
    if stype == 'unique':
        constants = result
    elif stype == 'free':
        constants, _ = result
    else:
        constants = [frac_zero()] * len(column_vecs)

    steps.append("Solve for constants from X(0) = X0:")
    for i, c in enumerate(constants):
        steps.append(f"  C{i+1} = {format_frac(c)}")

    # Build specific solution
    term_parts = []
    constants_idx = 0
    pair_acc = {}  # key=(alpha,beta) -> [c1, c2, a, b]
    pair_order = []

    for entry in basis:
        c = constants[constants_idx]
        constants_idx += 1
        if entry[0] == 'real':
            _, lam, vec, _ = entry
            if frac_is_zero(c): continue
            exp = _exp_str(lam)
            cn, cd = c
            prefix = "" if cn == cd else ("-" if cn == -cd else format_frac(c))
            part = f"{prefix}{format_vector(vec)}{exp}" if exp else f"{prefix}{format_vector(vec)}"
            term_parts.append(part)
        else:
            _, alpha, beta, a, b, which = entry
            key = (alpha, beta)
            if key not in pair_acc:
                pair_acc[key] = {'c1': frac_zero(), 'c2': frac_zero(), 'a': a, 'b': b}
                pair_order.append(key)
            if which == 'a':
                pair_acc[key]['c1'] = c
            else:
                pair_acc[key]['c2'] = c

    # Now combine each complex pair contribution
    # X_pair(t) = C1·e^(αt)[a·cos(βt) - b·sin(βt)] + C2·e^(αt)[a·sin(βt) + b·cos(βt)]
    #          = e^(αt)·[(C1·a + C2·b)·cos(βt) + (C2·a - C1·b)·sin(βt)]
    for key in pair_order:
        info = pair_acc[key]
        alpha, beta = key
        c1, c2, a, b = info['c1'], info['c2'], info['a'], info['b']
        cos_vec = vector_add(scalar_vector_mul(c1, a), scalar_vector_mul(c2, b))
        sin_vec = vector_sub(scalar_vector_mul(c2, a), scalar_vector_mul(c1, b))
        s = format_complex_pair_specific(alpha, beta, cos_vec, sin_vec)
        if s:
            term_parts.append(s)

    sol = "X(t) = " + " + ".join(term_parts) if term_parts else "X(t) = 0"
    steps.append(f"Specific solution:\n  {sol}")

    if verbose:
        print_steps("Initial Value Problem  X' = AX,  X(0) = X0", steps)
    return sol

# ============================================================
# SECTION 10: NONHOMOGENEOUS SYSTEMS (UNDETERMINED COEFFICIENTS)
# ============================================================

def check_resonance(A, k):
    """Is k an eigenvalue of A?"""
    evals = eigenvalues(A)
    kf = frac_to_float(k)
    for lam in evals:
        if not _is_complex_lam(lam) and abs(frac_to_float(lam) - kf) < 1e-9:
            return True
    return False



def _homogeneous_parts_from_epairs(epairs):
    """Return printable homogeneous real-eigenvalue terms and component pieces."""
    hom_parts = []
    component_pieces = None
    ci = 1
    for lam, vec in epairs:
        if _is_complex_lam(lam) or vec is None:
            continue
        exp = _exp_str(lam)
        term = f"C{ci}{format_vector(vec)}{exp}" if exp else f"C{ci}{format_vector(vec)}"
        hom_parts.append(term)
        if component_pieces is None:
            component_pieces = [[] for _ in range(len(vec))]
        for r, coef in enumerate(vec):
            if frac_is_zero(coef):
                continue
            cs = format_frac(coef)
            if cs == "1":
                piece = f"C{ci}{exp}"
            elif cs == "-1":
                piece = f"-C{ci}{exp}"
            else:
                piece = f"{cs}C{ci}{exp}"
            component_pieces[r].append(piece)
        ci += 1
    return hom_parts, component_pieces

def _join_component_pieces(pieces):
    if not pieces:
        return "0"
    out = pieces[0]
    for p in pieces[1:]:
        if p.startswith('-'):
            out += " - " + p[1:]
        else:
            out += " + " + p
    return out

def _add_vector_exp_to_components(component_pieces, vec, exp_str):
    if component_pieces is None:
        component_pieces = [[] for _ in range(len(vec))]
    for i, coef in enumerate(vec):
        if frac_is_zero(coef):
            continue
        cs = format_frac(coef)
        if cs == "1":
            piece = exp_str
        elif cs == "-1":
            piece = "-" + exp_str
        else:
            piece = cs + exp_str
        component_pieces[i].append(piece)
    return component_pieces

def _add_vector_constant_to_components(component_pieces, vec):
    if component_pieces is None:
        component_pieces = [[] for _ in range(len(vec))]
    for i, coef in enumerate(vec):
        if not frac_is_zero(coef):
            component_pieces[i].append(format_frac(coef))
    return component_pieces

def _component_form_lines(component_pieces):
    if component_pieces is None:
        return []
    names = ['x(t)', 'y(t)', 'z(t)', 'w(t)']
    lines = []
    for i, pieces in enumerate(component_pieces):
        nm = names[i] if i < len(names) else 'x' + str(i+1) + '(t)'
        lines.append(f"  {nm} = {_join_component_pieces(pieces)}")
    return lines

def _describe_homogeneous_part(A, epairs):
    """Detailed homogeneous steps for the nonhomogeneous solvers."""
    n = len(A)
    steps = []
    steps.append("STEP 1: Solve the associated homogeneous system X' = AX")
    if n == 2:
        tr = add_frac(A[0][0], A[1][1])
        det = determinant_2x2(A)
        steps.append(f"  trace(A) = {format_frac(tr)}")
        steps.append(f"  det(A) = {format_frac(det)}")
        steps.append(f"  Characteristic equation: λ² - ({format_frac(tr)})λ + ({format_frac(det)}) = 0")
    elif n == 3:
        coeffs = characteristic_polynomial_3x3(A)
        cs = [format_frac(c) for c in coeffs]
        steps.append(f"  Characteristic equation: λ³ + ({cs[1]})λ² + ({cs[2]})λ + ({cs[3]}) = 0")
    steps.append("  Eigenvalues and eigenvectors:")
    for i, (lam, vec) in enumerate(epairs):
        if isinstance(lam, tuple) and lam[0] == 'complex':
            steps.append(f"    λ{i+1} = {format_complex_eigenvalue(lam)},  K{i+1} = {format_complex_vector(vec)}")
        elif isinstance(lam, tuple) and lam[0] == 'complex_f':
            steps.append(f"    λ{i+1} ≈ {lam[1]:.4f} + {lam[2]:.4f}i")
        else:
            steps.append(f"    λ{i+1} = {format_frac(lam)},  K{i+1} = {format_vector(vec)}")
    hom_parts, _ = _homogeneous_parts_from_epairs(epairs)
    if hom_parts:
        steps.append("  Homogeneous solution:")
        steps.append("    Xh(t) = " + " + ".join(hom_parts))
    else:
        steps.append("  Homogeneous solution contains complex terms or could not be formatted here.")
    return steps

def solve_nonhomogeneous_constant(A, b, verbose=True):
    """
    Solve X' = AX + b (constant b).
    Particular solution: Ap = -b.
    Returns (full_solution_string, Xp_vector).
    """
    n = len(A)
    steps = [
        "Forcing: constant vector",
        f"b = {format_vector(b)}",
        "Trial: Xp = p  (constant vector)",
        "Substituting: 0 = Ap + b  ⟹  Ap = -b",
    ]
    neg_b = [neg_frac(x) for x in b]
    stype, result = solve_linear_system(A, neg_b)
    if stype == 'unique':
        Xp = result
    elif stype == 'free':
        Xp, _ = result
    else:
        steps.append("⚠  No particular solution (system inconsistent).")
        Xp = zero_vector(n)

    steps.append(f"Particular solution: Xp = {format_vector(Xp)}")

    hom_sol, epairs = solve_homogeneous_system(A, verbose=False)
    epairs = sort_eigenpairs(epairs)

    hom_parts = []
    for i, (lam, vec) in enumerate(epairs):
        if _is_complex_lam(lam) or vec is None: continue
        exp = _exp_str(lam)
        hom_parts.append(f"C{i+1}{format_vector(vec)}{exp}" if exp else f"C{i+1}{format_vector(vec)}")

    hom_part = " + ".join(hom_parts)
    Xp_any = any(not frac_is_zero(x) for x in Xp)
    full = f"X(t) = {hom_part} + {format_vector(Xp)}" if Xp_any else f"X(t) = {hom_part}"
    steps.append(f"Full solution:\n  {full}")

    if verbose:
        print_steps("Nonhomogeneous  X' = AX + b  (constant)", steps)
    return full, Xp

def solve_nonhomogeneous_exp(A, b, k, verbose=True):
    """
    Solve X' = AX + b·e^(kt).
    Returns (full_solution_string, Xp_vector).
    """
    n = len(A)
    ks = format_frac(k)
    steps = [
        f"Forcing: b·e^({ks}t),  b = {format_vector(b)}",
    ]

    I = identity_matrix(n)
    hom_sol, epairs = solve_homogeneous_system(A, verbose=False)
    epairs = sort_eigenpairs(epairs)

    if check_resonance(A, k):
        steps.append("⚠  Resonance: k is an eigenvalue. Trial: Xp = (tp + q)e^(kt)")
        # Eigenvector with eigenvalue k
        kI = scalar_matrix_mul(k, I)
        B = matrix_sub(A, kI)
        ns = nullspace(B)
        p = ns[0] if ns else b
        # Solve (A - kI)q = p - b
        rhs = vector_sub(p, b)
        stype, result = solve_linear_system(B, rhs)
        q = result if stype == 'unique' else (result[0] if stype == 'free' else zero_vector(n))
        steps.append(f"  p = {format_vector(p)}")
        steps.append(f"  q = {format_vector(q)}")
        Xp_str = f"t{format_vector(p)}e^({ks}t) + {format_vector(q)}e^({ks}t)"
        steps.append(f"Particular solution: Xp = {Xp_str}")
    else:
        steps.append(f"No resonance. Trial: Xp = p·e^({ks}t)")
        steps.append(f"Substituting: k·p·e^(kt) = A·p·e^(kt) + b·e^(kt)")
        steps.append(f"⟹  (kI - A)p = b")
        kI = scalar_matrix_mul(k, I)
        B = matrix_sub(kI, A)
        steps.append(f"(kI - A) =\n{format_matrix(B)}")
        stype, result = solve_linear_system(B, b)
        if stype == 'unique':
            Xp = result
        elif stype == 'free':
            Xp, _ = result
        else:
            steps.append("⚠  No particular solution.")
            Xp = zero_vector(n)
        steps.append(f"p = {format_vector(Xp)}")
        Xp_str = f"{format_vector(Xp)}e^({ks}t)"
        steps.append(f"Particular solution: Xp = {Xp_str}")

    hom_parts = []
    for i, (lam, vec) in enumerate(epairs):
        if _is_complex_lam(lam) or vec is None: continue
        exp = _exp_str(lam)
        hom_parts.append(f"C{i+1}{format_vector(vec)}{exp}" if exp else f"C{i+1}{format_vector(vec)}")
    hom_part = " + ".join(hom_parts)
    full = f"X(t) = {hom_part} + {Xp_str}"
    steps.append(f"Full solution:\n  {full}")

    if verbose:
        print_steps(f"Nonhomogeneous  X' = AX + b·e^({ks}t)", steps)
    return full, (Xp_str if check_resonance(A, k) else Xp)

def solve_nonhomogeneous_polynomial(A, polynomial_vectors, verbose=True):
    """
    Solve X' = AX + b_m·t^m + ... + b_0.
    polynomial_vectors: list of (power, vector).
    """
    n = len(A)
    max_power = max(pw for pw, _ in polynomial_vectors)
    bvecs = {pw: v for pw, v in polynomial_vectors}

    steps = [f"Forcing: polynomial up to t^{max_power}"]
    for pw, v in sorted(polynomial_vectors, reverse=True):
        steps.append(f"  b_{pw} = {format_vector(v)}")

    steps.append(f"Trial: Xp = p_{max_power}·t^{max_power} + ... + p_0")
    steps.append(f"Matching t^{max_power}: A·p_{max_power} = -b_{max_power}")

    pvecs = {}
    b_max = bvecs.get(max_power, zero_vector(n))
    neg_bm = [neg_frac(x) for x in b_max]
    stype, result = solve_linear_system(A, neg_bm)
    pvecs[max_power] = result if stype == 'unique' else (result[0] if stype == 'free' else zero_vector(n))
    steps.append(f"  p_{max_power} = {format_vector(pvecs[max_power])}")

    for k in range(max_power - 1, -1, -1):
        factor = frac_from_int(k + 1)
        pk1 = pvecs[k + 1]
        bk = bvecs.get(k, zero_vector(n))
        rhs = [sub_frac(mul_frac(factor, pk1[i]), bk[i]) for i in range(n)]
        steps.append(f"Matching t^{k}: A·p_{k} = {k+1}·p_{k+1} - b_{k}")
        stype, result = solve_linear_system(A, rhs)
        pvecs[k] = result if stype == 'unique' else (result[0] if stype == 'free' else zero_vector(n))
        steps.append(f"  p_{k} = {format_vector(pvecs[k])}")

    Xp_parts = []
    for pw in range(max_power, -1, -1):
        if pw in pvecs and any(not frac_is_zero(x) for x in pvecs[pw]):
            v_str = format_vector(pvecs[pw])
            if pw == 0:   Xp_parts.append(v_str)
            elif pw == 1: Xp_parts.append(f"{v_str}t")
            else:         Xp_parts.append(f"{v_str}t^{pw}")
    Xp_str = " + ".join(Xp_parts) if Xp_parts else "0"

    hom_sol, epairs = solve_homogeneous_system(A, verbose=False)
    epairs = sort_eigenpairs(epairs)
    hom_parts = []
    for i, (lam, vec) in enumerate(epairs):
        if _is_complex_lam(lam) or vec is None: continue
        exp = _exp_str(lam)
        hom_parts.append(f"C{i+1}{format_vector(vec)}{exp}" if exp else f"C{i+1}{format_vector(vec)}")
    hom_part = " + ".join(hom_parts)
    full = f"X(t) = {hom_part} + {Xp_str}"
    steps.append(f"Particular solution: Xp = {Xp_str}")
    steps.append(f"Full solution:\n  {full}")

    if verbose:
        print_steps("Nonhomogeneous Polynomial System", steps)
    return full, pvecs

def solve_undetermined_coefficients(A, forcing, verbose=True):
    """
    Main undetermined coefficients dispatcher.
    forcing: list of ('const', b) | ('exp', b, k) | ('poly', [(pw, v), ...])

    This version is intentionally verbose. It prints the same work a student
    normally writes: homogeneous solution, trial form, substitution, matrix
    equation for the unknown coefficients, particular solution, full vector
    solution, and component form when possible.
    """
    n = len(A)
    steps = []
    steps.append("Original system:")
    steps.append("  X' = AX + forcing")
    steps.append("Matrix A:")
    steps.append(format_matrix(A))

    hom_sol, epairs = solve_homogeneous_system(A, verbose=False)
    epairs = sort_eigenpairs(epairs)
    steps.extend(_describe_homogeneous_part(A, epairs))

    hom_parts, component_pieces = _homogeneous_parts_from_epairs(epairs)
    hom_part = " + ".join(hom_parts) if hom_parts else hom_sol.replace("X(t) = ", "")

    Xp_parts = []

    for force_index, force in enumerate(forcing, start=1):
        ftype = force[0]

        if ftype == 'const':
            _, b = force
            steps.append("")
            steps.append(f"STEP 2.{force_index}: Find a particular solution for constant forcing")
            steps.append(f"  b = {format_vector(b)}")
            steps.append("  Trial form: Xp = p, where p is a constant vector")
            steps.append("  Since Xp' = 0, substitute into X' = AX + b:")
            steps.append("    0 = Ap + b")
            steps.append("    Ap = -b")
            neg_b = [neg_frac(x) for x in b]
            steps.append(f"  Right side -b = {format_vector(neg_b)}")
            steps.append("  Solve the linear system Ap = -b.")
            stype, result = solve_linear_system(A, neg_b)
            if stype == 'unique':
                Xp = result
            elif stype == 'free':
                Xp, _ = result
                steps.append("  The system has free variables; using one valid particular vector.")
            else:
                Xp = zero_vector(n)
                steps.append("  No particular solution was found by this method.")
            steps.append(f"  p = {format_vector(Xp)}")
            steps.append(f"  Particular solution: Xp(t) = {format_vector(Xp)}")
            if any(not frac_is_zero(x) for x in Xp):
                Xp_parts.append(format_vector(Xp))
                component_pieces = _add_vector_constant_to_components(component_pieces, Xp)

        elif ftype == 'exp':
            _, b, k = force
            ks = format_frac(k)
            steps.append("")
            steps.append(f"STEP 2.{force_index}: Find a particular solution for exponential forcing")
            steps.append(f"  Forcing term: b·e^({ks}t), where b = {format_vector(b)}")
            steps.append(f"  Check resonance: compare k = {ks} with the eigenvalues of A.")
            eval_strings = []
            for lam, _vec in epairs:
                eval_strings.append(format_complex_eigenvalue(lam) if _is_complex_lam(lam) else format_frac(lam))
            steps.append("  Eigenvalues of A: " + ", ".join(eval_strings))

            I = identity_matrix(n)
            if check_resonance(A, k):
                steps.append("  Resonance occurs because k is an eigenvalue of A.")
                steps.append("  Use the resonant trial form:")
                steps.append("    Xp(t) = (t p + q)e^(kt)")
                steps.append("  The script solves for a valid p and q using generalized-eigenvector style equations.")
                kI = scalar_matrix_mul(k, I)
                B = matrix_sub(A, kI)
                steps.append("  A - kI =")
                steps.append(format_matrix(B))
                ns = nullspace(B)
                pvec = ns[0] if ns else b
                rhs = vector_sub(pvec, b)
                steps.append(f"  Choose p = {format_vector(pvec)} from the nullspace of A - kI.")
                steps.append(f"  Solve (A - kI)q = p - b = {format_vector(rhs)}")
                stype, result = solve_linear_system(B, rhs)
                qvec = result if stype == 'unique' else (result[0] if stype == 'free' else zero_vector(n))
                steps.append(f"  q = {format_vector(qvec)}")
                Xp_str = f"t{format_vector(pvec)}e^({ks}t) + {format_vector(qvec)}e^({ks}t)"
                steps.append(f"  Particular solution: Xp(t) = {Xp_str}")
                Xp_parts.append(Xp_str)
            else:
                steps.append("  No resonance, so use the trial form:")
                steps.append(f"    Xp(t) = p e^({ks}t)")
                steps.append("  Differentiate the trial solution:")
                steps.append(f"    Xp'(t) = {ks}p e^({ks}t)")
                steps.append("  Substitute into X' = AX + b e^(kt):")
                steps.append(f"    {ks}p e^({ks}t) = A p e^({ks}t) + b e^({ks}t)")
                steps.append("  Divide out the common exponential factor:")
                steps.append(f"    {ks}p = Ap + b")
                steps.append("  Move Ap to the left:")
                steps.append("    (kI - A)p = b")
                kI = scalar_matrix_mul(k, I)
                B = matrix_sub(kI, A)
                steps.append("  kI - A =")
                steps.append(format_matrix(B))
                steps.append("  Solve:")
                steps.append(f"    (kI - A)p = {format_vector(b)}")
                stype, result = solve_linear_system(B, b)
                if stype == 'unique':
                    pvec = result
                elif stype == 'free':
                    pvec, _ = result
                    steps.append("  The system has free variables; using one valid particular vector.")
                else:
                    pvec = zero_vector(n)
                    steps.append("  No particular solution was found by this method.")
                steps.append(f"  p = {format_vector(pvec)}")
                Xp_str = f"{format_vector(pvec)}e^({ks}t)"
                steps.append(f"  Particular solution: Xp(t) = {Xp_str}")
                Xp_parts.append(Xp_str)
                exp_piece = "e^t" if ks == "1" else f"e^({ks}t)"
                component_pieces = _add_vector_exp_to_components(component_pieces, pvec, exp_piece)

        elif ftype == 'poly':
            _, poly_vecs = force
            steps.append("")
            steps.append(f"STEP 2.{force_index}: Find a particular solution for polynomial forcing")
            max_power = max(pw for pw, _ in poly_vecs)
            bvecs = {pw: v for pw, v in poly_vecs}
            steps.append(f"  Highest power: t^{max_power}")
            for pw, v in sorted(poly_vecs, reverse=True):
                steps.append(f"  b_{pw} = {format_vector(v)}")
            steps.append(f"  Trial form: Xp = p_{max_power}t^{max_power} + ... + p_1t + p_0")
            steps.append("  Match coefficients of equal powers of t.")

            pvecs = {}
            b_max = bvecs.get(max_power, zero_vector(n))
            neg_bm = [neg_frac(x) for x in b_max]
            steps.append(f"  Highest power equation: A p_{max_power} = -b_{max_power} = {format_vector(neg_bm)}")
            stype, result = solve_linear_system(A, neg_bm)
            pvecs[max_power] = result if stype == 'unique' else (result[0] if stype == 'free' else zero_vector(n))
            steps.append(f"  p_{max_power} = {format_vector(pvecs[max_power])}")

            for pw in range(max_power - 1, -1, -1):
                factor = frac_from_int(pw + 1)
                pk1 = pvecs[pw + 1]
                bk = bvecs.get(pw, zero_vector(n))
                rhs = [sub_frac(mul_frac(factor, pk1[i]), bk[i]) for i in range(n)]
                steps.append(f"  Power t^{pw}: A p_{pw} = {pw+1}p_{pw+1} - b_{pw} = {format_vector(rhs)}")
                stype, result = solve_linear_system(A, rhs)
                pvecs[pw] = result if stype == 'unique' else (result[0] if stype == 'free' else zero_vector(n))
                steps.append(f"  p_{pw} = {format_vector(pvecs[pw])}")

            pp = []
            for pw in range(max_power, -1, -1):
                if pw in pvecs and any(not frac_is_zero(x) for x in pvecs[pw]):
                    vs = format_vector(pvecs[pw])
                    if pw == 0:
                        pp.append(vs)
                        component_pieces = _add_vector_constant_to_components(component_pieces, pvecs[pw])
                    elif pw == 1:
                        pp.append(f"{vs}t")
                    else:
                        pp.append(f"{vs}t^{pw}")
            Xp_poly = " + ".join(pp) if pp else "0"
            steps.append(f"  Particular solution: Xp(t) = {Xp_poly}")
            if pp:
                Xp_parts.append(Xp_poly)

    full = "X(t) = " + hom_part
    if Xp_parts:
        full += " + " + " + ".join(Xp_parts)

    steps.append("")
    steps.append("STEP 3: Combine homogeneous and particular parts")
    steps.append("  X(t) = Xh(t) + Xp(t)")
    steps.append(f"  {full}")

    comp_lines = _component_form_lines(component_pieces)
    if comp_lines:
        steps.append("")
        steps.append("STEP 4: Component form")
        steps.extend(comp_lines)

    if verbose:
        print_steps("Undetermined Coefficients — Detailed Solution", steps)
    return full

def parse_forcing_terms(n):
    """Interactive: collect forcing term from user."""
    print("\n  Forcing term types:")
    print("    1. Constant vector b")
    print("    2. Exponential  b·e^(kt)")
    print("    3. Polynomial   b₂t² + b₁t + b₀")
    choice = input("  Forcing type (1/2/3): ").strip()
    if choice == '1':
        print("  Enter constant vector b:")
        b = input_vector(n)
        return [('const', b)]
    elif choice == '2':
        k_str = input("  Enter k (exponent, e.g. 1 for e^t, -2 for e^(-2t)): ").strip()
        k = parse_number(k_str)
        print(f"  Enter vector b for b·e^({format_frac(k)}t):")
        b = input_vector(n)
        return [('exp', b, k)]
    elif choice == '3':
        deg_str = input("  Highest power (1 or 2): ").strip()
        max_deg = int(deg_str)
        poly_vecs = []
        for d in range(max_deg, -1, -1):
            print(f"  Enter vector for t^{d} term (enter zeros to skip):")
            v = input_vector(n)
            if any(not frac_is_zero(x) for x in v):
                poly_vecs.append((d, v))
        return [('poly', poly_vecs)]
    return []

# ============================================================
# SECTION 11: VERIFICATION
# ============================================================

def verify_homogeneous_vector_solution(A, v, lam, verbose=True):
    """Verify X = v·e^(λt) solves X' = AX by checking Av = λv."""
    Av = matrix_vector_mul(A, v)
    lv = scalar_vector_mul(lam, v)
    steps = [
        f"X = v·e^(λt),  v = {format_vector(v)},  λ = {format_frac(lam)}",
        f"X' = λ·v·e^(λt) = {format_frac(lam)}{format_vector(v)}·e^(λt)",
        f"AX = A·v·e^(λt)",
        f"Compute Av: {format_vector(Av)}",
        f"Compute λv: {format_vector(lv)}",
    ]
    ok = all(frac_eq(Av[i], lv[i]) for i in range(len(v)))
    if ok:
        steps.append("✓  Av = λv  ⟹  X = v·e^(λt) is a solution of X' = AX.")
    else:
        steps.append("✗  Av ≠ λv  ⟹  X does NOT satisfy X' = AX.")
    if verbose:
        print_steps("Verify Exponential Solution", steps)
    return ok

def verify_solution(A, X_components, verbose=True):
    """
    Verify a vector function X(t) satisfies X' = AX.
    X_components: list of symbolic expression lists.
    """
    n = len(A)
    X_prime = [derivative_expr(comp) for comp in X_components]
    AX = []
    for i in range(n):
        row_sum = []
        for j in range(n):
            coef = A[i][j]
            if not frac_is_zero(coef):
                row_sum = expr_add(row_sum, expr_mul_const(coef, X_components[j]))
        AX.append(row_sum)
    steps = ["Components of X(t):"]
    for i, comp in enumerate(X_components):
        steps.append(f"  x{i+1} = {format_expr(comp)}")
    steps.append("Components of X'(t):")
    for i, comp in enumerate(X_prime):
        steps.append(f"  x{i+1}' = {format_expr(comp)}")
    steps.append("Components of AX:")
    for i, comp in enumerate(AX):
        steps.append(f"  (AX){i+1} = {format_expr(comp)}")
    steps.append("Comparison X' vs AX:")
    all_match = True
    for i in range(n):
        if expr_equal(X_prime[i], AX[i]):
            steps.append(f"  ✓ Component {i+1} matches.")
        else:
            steps.append(f"  ✗ Component {i+1} differs.")
            all_match = False
    steps.append("✓  X(t) is a valid solution." if all_match else "✗  X(t) does NOT satisfy X' = AX.")
    if verbose:
        print_steps("Verify Solution", steps)
    return all_match

def wronskian_matrix(solution_vectors):
    """Build Wronskian matrix (each vector is a column)."""
    n = len(solution_vectors[0])
    return [[solution_vectors[j][i] for j in range(len(solution_vectors))] for i in range(n)]

def fundamental_set_check(solution_pairs, verbose=True):
    """
    Check if solutions form a fundamental set.
    solution_pairs: [(vector, lambda), ...] where Xi = vi·e^(λi·t)
    """
    n = len(solution_pairs)
    vecs = [v for v, _ in solution_pairs]
    lams = [l for _, l in solution_pairs]
    W = wronskian_matrix(vecs)
    det_V = determinant(W)

    steps = [
        "Build Wronskian matrix W from solution vectors (as columns):",
        format_matrix(W),
        f"det(V) = {format_frac(det_V)}",
    ]
    if not frac_is_zero(det_V):
        sum_lam = frac_zero()
        for l in lams:
            if not isinstance(l, tuple): sum_lam = add_frac(sum_lam, l)
        exp_factor = f"e^({format_frac(sum_lam)}t)"
        steps.append(f"W(t) = {format_frac(det_V)}·{exp_factor}")
        steps.append("Since det(V) ≠ 0, W(t) ≠ 0 for all t.")
        steps.append("✓  Yes — the solution set is linearly independent on (-∞, ∞).")
        steps.append("    The vectors form a fundamental set of solutions.")
        result = True
    else:
        steps.append("det(V) = 0  ⟹  solutions are linearly DEPENDENT.")
        steps.append("✗  No — the vectors do NOT form a fundamental set.")
        result = False
    if verbose:
        print_steps("Fundamental Set Check", steps)
    return result

# ============================================================
# SECTION 12: MATRIX ↔ SCALAR EQUATION CONVERSION
# ============================================================

def write_matrix_form_from_equations(n, equations):
    """
    Convert list of parsed equations into matrix form display.
    equations: list of (A_row, forcing_terms) from parse_scalar_equation.
    """
    A = []
    const_vec = zero_vector(n)
    t_vec = zero_vector(n)
    t2_vec = zero_vector(n)
    exp_vecs = {}  # k -> vector

    for i, (coefs, forcing) in enumerate(equations):
        A.append(coefs)
        for ftype, val in forcing:
            if ftype == 'const':
                const_vec[i] = add_frac(const_vec[i], val)
            elif ftype == 't':
                t_vec[i] = add_frac(t_vec[i], val)
            elif ftype == 't2':
                t2_vec[i] = add_frac(t2_vec[i], val)
            elif isinstance(ftype, tuple) and ftype[0] == 'exp':
                k = ftype[1]
                kf = frac_to_float(k)
                key = round(kf, 8)
                if key not in exp_vecs:
                    exp_vecs[key] = (k, zero_vector(n))
                k2, v = exp_vecs[key]
                v[i] = add_frac(v[i], val)

    print(f"\nMatrix form:  X' = AX + F(t)")
    print(f"\nA =\n{format_matrix(A)}")

    forcing_parts = []
    if any(not frac_is_zero(x) for x in t2_vec):
        forcing_parts.append(f"t²{format_vector(t2_vec)}")
    if any(not frac_is_zero(x) for x in t_vec):
        forcing_parts.append(f"t{format_vector(t_vec)}")
    if any(not frac_is_zero(x) for x in const_vec):
        forcing_parts.append(format_vector(const_vec))
    for key, (k, v) in exp_vecs.items():
        if any(not frac_is_zero(x) for x in v):
            forcing_parts.append(f"{format_vector(v)}e^({format_frac(k)}t)")

    if forcing_parts:
        print(f"F(t) = {' + '.join(forcing_parts)}")
    else:
        print("F(t) = 0  (homogeneous system)")

    return A, const_vec, t_vec, t2_vec, exp_vecs

def parse_scalar_equation(s, var_names, n):
    """
    Parse a scalar ODE RHS. Returns (coefs_list, forcing_terms).
    Very simplified tokenizer — handles: ±coef*var, ±coef*t^n, ±coef*e^(kt), constants.
    """
    import re
    coefs = [frac_zero()] * n
    forcing = []

    s = s.strip()
    # Normalize spacing around +/-
    s = re.sub(r'\s*\+\s*', '+', s)
    s = re.sub(r'\s*-\s*', '-', s)
    s = s.replace(' ', '')

    # Split into signed tokens
    tokens = []
    current = ''
    for i, ch in enumerate(s):
        if ch in '+-' and i > 0 and s[i-1] not in '(^*':
            if current:
                tokens.append(current)
            current = ch
        else:
            current += ch
    if current:
        tokens.append(current)

    for tok in tokens:
        if not tok:
            continue
        # Parse sign
        sign = frac_one()
        if tok.startswith('-'):
            sign = neg_frac(frac_one()); tok = tok[1:]
        elif tok.startswith('+'):
            tok = tok[1:]

        matched_var = False
        for j, vname in enumerate(var_names):
            # Check for patterns like: vname, coef*vname, coefvname
            patterns = [vname, f'*{vname}']
            for pat in patterns:
                if tok.endswith(pat):
                    coef_str = tok[:-len(pat)].rstrip('*')
                    if not coef_str:
                        coef = frac_one()
                    else:
                        try: coef = parse_number(coef_str)
                        except: coef = frac_one()
                    coefs[j] = add_frac(coefs[j], mul_frac(sign, coef))
                    matched_var = True
                    break
            if matched_var: break

        if not matched_var:
            tl = tok.lower()
            # Exponential
            if 'e^' in tl:
                e_pos = tl.index('e^')
                coef_str = tok[:e_pos]
                exp_str = tok[e_pos+2:]
                try: coef = parse_number(coef_str) if coef_str else frac_one()
                except: coef = frac_one()
                coef = mul_frac(sign, coef)
                # Parse exponent like (kt) or t or (t) or -t or (-2t)
                exp_str = exp_str.strip('()')
                if exp_str.endswith('t'):
                    k_str = exp_str[:-1]
                    if not k_str or k_str == '+': k = frac_one()
                    elif k_str == '-': k = neg_frac(frac_one())
                    else:
                        try: k = parse_number(k_str)
                        except: k = frac_one()
                else:
                    k = frac_one()
                forcing.append((('exp', k), coef))
            elif tok == 't^2' or tok.endswith('t^2'):
                coef_str = tok[:-3]
                try: coef = parse_number(coef_str) if coef_str else frac_one()
                except: coef = frac_one()
                forcing.append(('t2', mul_frac(sign, coef)))
            elif tok == 't' or tok.endswith('t'):
                coef_str = tok[:-1]
                try: coef = parse_number(coef_str) if coef_str else frac_one()
                except: coef = frac_one()
                forcing.append(('t', mul_frac(sign, coef)))
            else:
                try:
                    coef = parse_number(tok)
                    forcing.append(('const', mul_frac(sign, coef)))
                except:
                    pass  # unknown token

    return coefs, forcing

def write_equations_from_matrix_form(A, b_forcing, k):
    """Display X' = AX + b·e^(kt) as scalar equations."""
    n = len(A)
    var_names   = ['x', 'y', 'z'][:n]
    deriv_names = ['dx/dt', 'dy/dt', 'dz/dt'][:n]
    print("\nScalar equations:")
    for i in range(n):
        parts = []
        for j in range(n):
            coef = A[i][j]
            if not frac_is_zero(coef):
                parts.append(format_expr_term(coef, var_names[j]))
        if b_forcing:
            bi = b_forcing[i]
            if not frac_is_zero(bi):
                bi_s = format_frac(bi)
                if k is None or frac_is_zero(k):
                    parts.append(bi_s)
                else:
                    ks = format_frac(k)
                    exp = "e^t" if ks == "1" else f"e^({ks}t)"
                    if bi_s == "1":   parts.append(exp)
                    elif bi_s == "-1":parts.append(f"-{exp}")
                    else:             parts.append(f"{bi_s}{exp}")
        if not parts:
            rhs = "0"
        else:
            rhs = parts[0]
            for p in parts[1:]:
                rhs += (" - " + p[1:]) if p.startswith('-') else (" + " + p)
        print(f"  {deriv_names[i]} = {rhs}")

# ============================================================
# SECTION 13: INPUT HELPERS
# ============================================================

def ask_dimension():
    while True:
        s = input("  System size — 2 or 3 variables? ").strip()
        if s in ('2', '3'):
            return int(s)
        print("  Please enter 2 or 3.")

def ask_solution_type():
    opts = {
        '1': 'Homogeneous', '2': 'Nonhomogeneous constant',
        '3': 'Nonhomogeneous exponential', '4': 'Nonhomogeneous polynomial',
        '5': 'Initial value problem', '6': 'Verification', '7': 'Matrix form conversion',
    }
    for k, v in opts.items():
        print(f"    {k}. {v}")
    while True:
        s = input("  Choice: ").strip()
        if s in opts:
            return int(s)
        print("  Please enter 1–7.")

def input_matrix(rows, cols):
    print(f"\n  Enter {rows}×{cols} matrix — one row per line, entries space-separated.")
    print("  Supports fractions (1/2), decimals (0.5), negatives (-3).")
    A = []
    for i in range(rows):
        while True:
            row_str = input(f"  Row {i+1}: ").strip()
            parts = row_str.split()
            if len(parts) == cols:
                try:
                    A.append([parse_number(p) for p in parts])
                    break
                except ValueError as e:
                    print(f"  ✗ {e}  — please re-enter.")
            else:
                print(f"  Need {cols} entries, got {len(parts)}.")
    return A

def input_vector(size):
    while True:
        row_str = input(f"  Enter {size}-component vector (space-separated): ").strip()
        parts = row_str.split()
        if len(parts) == size:
            try:
                return [parse_number(p) for p in parts]
            except ValueError as e:
                print(f"  ✗ {e}  — please re-enter.")
        else:
            print(f"  Need {size} entries, got {len(parts)}.")

def safe_input_loop(prompt, parser):
    while True:
        try:
            s = input(prompt).strip()
            if not s:
                print("  Empty input — please try again.")
                continue
            return parser(s)
        except (ValueError, ZeroDivisionError) as e:
            print(f"  ✗ {e}  — please try again.")

# ============================================================
# SECTION 14: MENU HANDLERS
# ============================================================

def menu_matrix_form():
    print("\n--- Convert Scalar Equations → Matrix Form ---")
    n = ask_dimension()
    var_names   = ['x', 'y', 'z'][:n]
    deriv_names = ['dx/dt', 'dy/dt', 'dz/dt'][:n]
    print(f"\n  Enter the RHS of each equation.")
    print(f"  Variables: {', '.join(var_names)}")
    print(f"  Supported: constants, t, t^2, e^t, e^(kt)  e.g.  3x - 2y + t - 1")
    equations = []
    for i in range(n):
        while True:
            s = input(f"  {deriv_names[i]} = ").strip()
            if s:
                row = parse_scalar_equation(s, var_names, n)
                equations.append(row)
                break
    write_matrix_form_from_equations(n, equations)

def menu_equations_from_matrix():
    print("\n--- Matrix System → Scalar Equations ---")
    n = ask_dimension()
    print("  Enter matrix A:")
    A = input_matrix(n, n)
    print("\n  Forcing type?")
    print("    1. None (homogeneous)")
    print("    2. Constant vector b")
    print("    3. b·e^(kt)")
    ftype = input("  Choice: ").strip()
    if ftype == '1':
        write_equations_from_matrix_form(A, None, None)
    elif ftype == '2':
        print("  Enter forcing vector b:")
        b = input_vector(n)
        write_equations_from_matrix_form(A, b, frac_zero())
    elif ftype == '3':
        k = safe_input_loop("  Enter k (exponent): ", parse_number)
        print("  Enter forcing vector b:")
        b = input_vector(n)
        write_equations_from_matrix_form(A, b, k)

def menu_verify_solution():
    print("\n--- Verify Proposed Solution ---")
    print("  Solution types:")
    print("    1. Exponential  X = v·e^(λt)")
    print("    2. Trig components  (sin/cos)")
    stype = input("  Type (1/2): ").strip()
    n = ask_dimension()
    print("  Enter matrix A:")
    A = input_matrix(n, n)

    if stype == '1':
        print("  Enter eigenvector v:")
        v = input_vector(n)
        lam = safe_input_loop("  Enter eigenvalue λ: ", parse_number)
        verify_homogeneous_vector_solution(A, v, lam)

    elif stype == '2':
        print("\n  Enter each component of X as a·sin(t) + b·cos(t)")
        X_components = []
        for i in range(n):
            print(f"\n  Component {i+1}:")
            a = safe_input_loop("    Coef of sin(t): ", parse_number)
            b = safe_input_loop("    Coef of cos(t): ", parse_number)
            comp = []
            if not frac_is_zero(a): comp.append({'coef': a, 'type': 'sin', 'k': frac_one()})
            if not frac_is_zero(b): comp.append({'coef': b, 'type': 'cos', 'k': frac_one()})
            X_components.append(comp)
        verify_solution(A, X_components)

def menu_fundamental_set():
    print("\n--- Fundamental Set Check ---")
    n = ask_dimension()
    print(f"  Enter {n} solution vectors of the form Xi = vi·e^(λi·t):")
    solution_pairs = []
    for i in range(n):
        print(f"\n  Solution X{i+1}:")
        v = input_vector(n)
        lam = safe_input_loop(f"  Exponent λ{i+1}: ", parse_number)
        solution_pairs.append((v, lam))
    fundamental_set_check(solution_pairs)

def menu_eigenvalues():
    print("\n--- Eigenvalues and Eigenvectors ---")
    n = ask_dimension()
    print("  Enter matrix A:")
    A = input_matrix(n, n)

    epairs = sort_eigenpairs(_get_eigenpairs(A))
    steps = [f"Matrix A:\n{format_matrix(A)}"]
    if n == 2:
        tr = add_frac(A[0][0], A[1][1])
        det = determinant_2x2(A)
        steps.append(f"trace = {format_frac(tr)},  det = {format_frac(det)}")
        steps.append(f"Characteristic polynomial: λ² - ({format_frac(tr)})λ + ({format_frac(det)}) = 0")
    else:
        coeffs = characteristic_polynomial_3x3(A)
        cs = [format_frac(c) for c in coeffs]
        steps.append(f"Characteristic polynomial:\n  λ³ + ({cs[1]})λ² + ({cs[2]})λ + ({cs[3]}) = 0")
    steps.append("Eigenvalues and eigenvectors (sorted smallest to largest):")
    for i, (lam, vec) in enumerate(epairs):
        if isinstance(lam, tuple) and lam[0] == 'complex':
            steps.append(f"  λ{i+1} = {format_complex_eigenvalue(lam)}")
            if vec is not None:
                steps.append(f"     K{i+1} = {format_complex_vector(vec)}")
        elif isinstance(lam, tuple) and lam[0] == 'complex_f':
            steps.append(f"  λ{i+1} ≈ {lam[1]:.4f} + {lam[2]:.4f}i  (irrational complex)")
        else:
            steps.append(f"  λ{i+1} = {format_frac(lam)},  K{i+1} = {format_vector(vec)}")
    print_steps("Eigenvalues and Eigenvectors", steps)

    ans_lines = []
    for i, (l, v) in enumerate(epairs):
        if isinstance(l, tuple) and l[0] == 'complex':
            ans_lines.append(f"λ{i+1} = {format_complex_eigenvalue(l)},  K{i+1} = {format_complex_vector(v)}")
        elif isinstance(l, tuple) and l[0] == 'complex_f':
            ans_lines.append(f"λ{i+1} ≈ {l[1]:.4f} + {l[2]:.4f}i  (irrational)")
        else:
            ans_lines.append(f"λ{i+1} = {format_frac(l)},  K{i+1} = {format_vector(v)}")
    copy_paste_answer_only("\n  ".join(ans_lines))

def menu_homogeneous():
    print("\n--- Solve Homogeneous System  X' = AX ---")
    n = ask_dimension()
    print("  Enter matrix A:")
    A = input_matrix(n, n)
    sol, _ = solve_homogeneous_system(A)
    copy_paste_answer_only(sol)

def menu_nonhomogeneous():
    print("\n--- Solve Nonhomogeneous System (Undetermined Coefficients) ---")
    n = ask_dimension()
    print("  Enter matrix A:")
    A = input_matrix(n, n)
    forcing = parse_forcing_terms(n)
    if not forcing:
        print("  No forcing term entered.")
        return
    sol = solve_undetermined_coefficients(A, forcing)
    copy_paste_answer_only(sol)

def menu_ivp():
    print("\n--- Solve Initial Value Problem  X' = AX,  X(0) = X0 ---")
    n = ask_dimension()
    print("  Enter matrix A:")
    A = input_matrix(n, n)
    print("  Enter initial condition X(0):")
    X0 = input_vector(n)
    sol = solve_ivp_homogeneous(A, X0)
    copy_paste_answer_only(sol)

# ============================================================
# SECTION 15: EXAMPLE LIBRARY
# ============================================================

def run_example_library():
    print("\n" + "="*62)
    print("  EXAMPLE LIBRARY")
    print("="*62)
    examples = [
        " 1. Scalar → matrix form:  dx/dt=8x-4y, dy/dt=2x+3y",
        " 2. 3×3 polynomial forcing matrix form",
        " 3. Matrix → scalar:  X'=[[6,4],[-1,2]]X + [1,-1]e^t",
        " 4. Verify exponential solution:  X=[1,3]e^(-6t)",
        " 5. Fundamental set:  X₁=[1,1]e^(-9t), X₂=[1,-1]e^(-3t)",
        " 6. 3×3 eigenvalue system: dx/dt=x+y-z, dy/dt=5y, dz/dt=y-z",
        " 7. 2×2 IVP:  A=[[1/2,0],[1,-1/2]], X(0)=[4,6]",
        " 8. Repeated eigenvalue:  dx/dt=6x-y, dy/dt=36x-6y",
        " 9. 3×3 homogeneous:  dx/dt=3x-y-z, ...",
        "10. Nonhomogeneous constant:  5x+9y+7, -x+11y+21",
        "11. Nonhomogeneous exponential:  X'=AX+[-5,6]e^t",
        "12. Complex eigenvalues (2×2):  A=[[3,-2],[4,-1]]  (λ=1±2i)",
        "13. Complex IVP (pure imaginary):  A=[[0,1],[-1,0]], X(0)=[1,0]",
        " 0. Back to main menu",
    ]
    for ex in examples:
        print(f"  {ex}")
    choice = input("\n  Choice: ").strip()

    F = frac_from_int
    NEG = neg_frac
    Z = frac_zero; ONE = frac_one

    if choice == '0':
        return

    elif choice == '1':
        A = [[F(8), F(-4)], [F(2), F(3)]]
        equations = [
            ([F(8), F(-4)], []),
            ([F(2), F(3)],  []),
        ]
        write_matrix_form_from_equations(2, equations)
        copy_paste_answer_only(f"A =\n{format_matrix(A)}\nF(t) = [0, 0]")

    elif choice == '2':
        A = [[ONE(), NEG(ONE()), ONE()],
             [F(3),  ONE(),       NEG(ONE())],
             [ONE(), ONE(),       ONE()]]
        equations = [
            ([ONE(), NEG(ONE()), ONE()],  [('const', F(-1)), ('t', ONE())]),
            ([F(3),  ONE(),      NEG(ONE())], [('t2', F(-2))]),
            ([ONE(), ONE(),      ONE()],   [('const', F(3)), ('t2', ONE()), ('t', NEG(ONE()))]),
        ]
        write_matrix_form_from_equations(3, equations)

    elif choice == '3':
        A = [[F(6), F(4)], [NEG(ONE()), F(2)]]
        b = [ONE(), NEG(ONE())]
        k = ONE()
        write_equations_from_matrix_form(A, b, k)

    elif choice == '4':
        # dx/dt = 3x - 3y,  dy/dt = 3x - 7y
        A = [[F(3), F(-3)], [F(3), F(-7)]]
        v = [ONE(), F(3)]
        lam = F(-6)
        verify_homogeneous_vector_solution(A, v, lam)

    elif choice == '5':
        solution_pairs = [
            ([ONE(), ONE()],       F(-9)),
            ([ONE(), NEG(ONE())],  F(-3)),
        ]
        fundamental_set_check(solution_pairs)

    elif choice == '6':
        A = [[ONE(), ONE(), NEG(ONE())],
             [Z(), F(5), Z()],
             [Z(), ONE(), NEG(ONE())]]
        sol, _ = solve_homogeneous_system(A)
        copy_paste_answer_only(sol)

    elif choice == '7':
        A = [[parse_number('1/2'), Z()],
             [ONE(), parse_number('-1/2')]]
        X0 = [F(4), F(6)]
        sol = solve_ivp_homogeneous(A, X0)
        copy_paste_answer_only(sol)

    elif choice == '8':
        A = [[F(6), NEG(ONE())], [F(36), F(-6)]]
        sol, _ = solve_repeated_eigenvalue_system(A)
        copy_paste_answer_only(sol)

    elif choice == '9':
        A = [[F(3), NEG(ONE()), NEG(ONE())],
             [ONE(), ONE(),      NEG(ONE())],
             [ONE(), NEG(ONE()), ONE()]]
        sol, _ = solve_homogeneous_system(A)
        copy_paste_answer_only(sol)

    elif choice == '10':
        A = [[F(5), F(9)], [NEG(ONE()), F(11)]]
        b = [F(7), F(21)]
        sol, Xp = solve_nonhomogeneous_constant(A, b)
        copy_paste_answer_only(sol)

    elif choice == '11':
        A = [[F(6), F(4)], [NEG(ONE()), F(2)]]
        b = [F(-5), F(6)]
        k = ONE()
        sol, Xp = solve_nonhomogeneous_exp(A, b, k)
        copy_paste_answer_only(sol)

    elif choice == '12':
        # A = [[3, -2], [4, -1]]  →  trace=2, det=-3+8=5  →  λ²-2λ+5=0  →  λ=1±2i
        A = [[F(3), F(-2)], [F(4), F(-1)]]
        sol, _ = solve_homogeneous_system(A)
        copy_paste_answer_only(sol)

    elif choice == '13':
        # A = [[0, 1], [-1, 0]]  →  λ=±i,  pure rotation
        A = [[Z(), ONE()], [NEG(ONE()), Z()]]
        X0 = [ONE(), Z()]
        sol = solve_ivp_homogeneous(A, X0)
        copy_paste_answer_only(sol)

    else:
        print("  Invalid choice.")
        return

    input("\n  [Press Enter to continue]")

# ============================================================
# SECTION 16: MAIN MENU
# ============================================================

def main_menu():
    while True:
        print(f"\n{'='*62}")
        print("       LINEAR SYSTEMS OF ODE SOLVER")
        print("       Exact arithmetic · No dependencies")
        print('='*62)
        print("  1. Write system in matrix form")
        print("  2. Write matrix system as scalar equations")
        print("  3. Verify a proposed vector solution")
        print("  4. Check if solution vectors form a fundamental set")
        print("  5. Find eigenvalues and eigenvectors")
        print("  6. Solve homogeneous linear system  X' = AX")
        print("  7. Solve nonhomogeneous system (undetermined coefficients)")
        print("  8. Solve initial value problem")
        print("  9. Example library (built-in problem templates)")
        print("  0. Exit")
        print()
        choice = input("  Choice: ").strip()

        try:
            if   choice == '1': menu_matrix_form()
            elif choice == '2': menu_equations_from_matrix()
            elif choice == '3': menu_verify_solution()
            elif choice == '4': menu_fundamental_set()
            elif choice == '5': menu_eigenvalues()
            elif choice == '6': menu_homogeneous()
            elif choice == '7': menu_nonhomogeneous()
            elif choice == '8': menu_ivp()
            elif choice == '9': run_example_library()
            elif choice == '0':
                print("\n  Goodbye!\n")
                sys.exit(0)
            else:
                print("  Please enter 0–9.")
        except KeyboardInterrupt:
            print("\n  (Interrupted — returning to menu)")
        except Exception as e:
            print(f"\n  ⚠  Unexpected error: {e}")
            import traceback; traceback.print_exc()

# ============================================================
# ENTRY POINT
# ============================================================

if __name__ == '__main__':
    print("\nLinear Systems of ODE Solver")
    print("Pure Python · Exact fraction arithmetic · No external packages\n")
    main_menu()
