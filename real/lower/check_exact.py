"""Small exact checks for the analytic identities and one specified matrix."""
from fractions import Fraction as Q
from pathlib import Path
import json


def add(a, b):
    result = a.copy()
    for degree, c in b.items():
        result[degree] = result.get(degree, 0)+c
    return {k:v for k,v in result.items() if v}


def scale(a, c):
    return {k:v*c for k,v in a.items() if v*c}


def mul(a, b):
    result = {}
    for (i, j), c in a.items():
        for (k, l), d in b.items():
            degree = (i+k, j+l)
            result[degree] = result.get(degree, 0)+c*d
    return {k:v for k,v in result.items() if v}


def power(a, n):
    result = {(0, 0): 1}
    for _ in range(n):
        result = mul(result, a)
    return result


def arctan_bounds(q, pairs=16):
    # The series alternates with strictly decreasing term magnitudes.
    x = Q(1, q)
    lower = sum(((-1)**j*x**(2*j+1)/Q(2*j+1) for j in range(2*pairs)), Q(0))
    upper = lower+x**(4*pairs+1)/Q(4*pairs+1)
    return lower, upper


def main():
    one = {(0, 0): 1}
    u, v = {(1, 0): 1}, {(0, 1): 1}
    u2, v2 = power(u, 2), power(v, 2)
    U = add(scale(power(u, 3), 4), scale(u, -3))
    V = add(scale(power(v, 3), 4), scale(v, -3))
    z1 = add(one, scale(add(u2, scale(mul(u, v), -1)), 4))
    z3 = add(one, scale(add(power(U, 2), scale(mul(U, V), -1)), 4))
    lhs = add(scale(one, 9), scale(add(scale(z1, 8), z3), -1))
    bracket = add(scale(one, 8), mul(add(scale(add(add(u2, mul(u, v)), v2), 4), scale(one, -3)),
                                    add(scale(u2, 4), scale(one, -3))))
    rhs = scale(mul(mul(u, add(v, scale(u, -1))), bracket), 4)
    assert lhs == rhs
    sos = add(mul(power(add(scale(u, 2), scale(one, -1)), 2),
                  add(add(scale(u2, 4), scale(u, 8)), scale(one, 5))),
              scale(mul(mul(add(one, scale(v, -1)), add(add(one, v), u)),
                        add(scale(one, 3), scale(u2, -4))), 4))
    assert bracket == sos
    l5, u5 = arctan_bounds(5)
    l239, u239 = arctan_bounds(239, 4)
    pi_lo, pi_hi = 16*l5-4*u239, 16*u5-4*l239
    simple_pi_lo = Q(3141592653589793, 10**15)
    simple_pi_hi = Q(3141592653589794, 10**15)
    assert pi_lo > simple_pi_lo and pi_hi < simple_pi_hi
    # Machin identity: tan(4 atan(1/5))=120/119, and subtraction
    # of atan(1/239) gives tangent one in the first quadrant.
    assert (Q(120, 119)-Q(1, 239))/(1+Q(120, 119*239)) == 1
    a, C = Q(13, 8), Q(2917, 2000)
    circle_gap = C-(Q(35, 8)/simple_pi_lo+Q(47, 72)*Q(11, 50)**2)
    assert circle_gap > 0
    A_gaps = [C-(a*x*x+Q(3, 2)*(1-x)+Q(26, 125)) for x in (Q(11, 50), Q(2, 3))]
    assert min(A_gaps) > 0
    # Target B: C - expression = y/20-y^2/8 >= y/120 on [0,1/3].
    assert Q(1, 20)-Q(1, 3)/8 == Q(1, 120)
    assert a-Q(1, 6)+Q(1, 6000) == C
    n, q, R = 100000, 40, 20
    # 3^(-n/2) < 2^(-100), exp(-R^2/2)<2^(-200).
    eps2 = Q(n, 2**(2*q))+Q(12, 2**100)+Q(4*n, 2**200)
    assert eps2 < Q(1, 10**18)
    vector_lower = (a+1)*(1-Q(1, n))**2-1-2*a*Q(1, 10**9)
    ratio_lower = vector_lower*simple_pi_lo*1000/Q(2917)
    assert ratio_lower > Q(7, 4)
    result = dict(circle_polynomial_identities=True,
                   pi_lower=str(simple_pi_lo), pi_upper=str(simple_pi_hi),
                   circle_endpoint_gap=str(circle_gap),
                   A_endpoint_gaps=list(map(str, A_gaps)),
                   B_gap_lower="y/120 for 0 <= y <= 1/3",
                   finite_dimension=n, cutoff=R, mesh_exponent=q,
                   squared_discretization_error_upper=str(eps2),
                   vector_value_lower=str(vector_lower),
                   finite_ratio_lower=str(ratio_lower),
                   finite_ratio_gap_above_7_over_4=str(ratio_lower-Q(7, 4)),
                   conditional_on_universal_fiber_targets=True)
    Path(__file__).with_name("exact_checks.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
