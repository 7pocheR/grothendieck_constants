"""Exact algebra and final finite-matrix margins; no numerical integration."""
from fractions import Fraction as Q
from pathlib import Path
import json
import sys

if not __debug__:
    raise RuntimeError("Assertions must be enabled")
sys.dont_write_bytecode = True
from check_exact import add, scale, mul, power, arctan_bounds
from check_fibers import matrix_checks


def q30_mul(a, b):
    return (a[0]*b[0] + 30*a[1]*b[1], a[0]*b[1] + a[1]*b[0])


def main():
    one, u, v = {(0, 0): Q(1)}, {(1, 0): Q(1)}, {(0, 1): Q(1)}
    u2, v2 = power(u, 2), power(v, 2)
    U = add(scale(power(u, 3), 4), scale(u, -3))
    V = add(scale(power(v, 3), 4), scale(v, -3))
    z1 = add(one, scale(add(u2, scale(mul(u, v), -1)), 4))
    z3 = add(one, scale(add(power(U, 2), scale(mul(U, V), -1)), 4))
    bracket = add(scale(one, 8), mul(add(scale(u2, 4), scale(one, -3)),
                  add(scale(add(add(u2, mul(u, v)), v2), 4), scale(one, -3))))
    assert add(scale(one, 9), scale(add(scale(z1, 8), z3), -1)) == scale(mul(mul(u, add(v, scale(u, -1))), bracket), 4)
    assert bracket == add(mul(power(add(scale(u, 2), scale(one, -1)), 2),
                              add(add(scale(u2, 4), scale(u, 8)), scale(one, 5))),
                          scale(mul(mul(add(one, scale(v, -1)), add(add(one, v), u)),
                                    add(scale(one, 3), scale(u2, -4))), 4))
    l5, h5 = arctan_bounds(5)
    l239, h239 = arctan_bounds(239, 4)
    pi_lo, pi_hi = Q(3141592653589793, 10**15), Q(3141592653589794, 10**15)
    assert 16*l5 - 4*h239 > pi_lo and 16*h5 - 4*l239 < pi_hi
    assert (Q(120,119)-Q(1,239))/(1+Q(120,119*239)) == 1
    a, C = Q(13,8), Q(2917,2000)
    circle_gap = C - Q(35,8)/pi_lo - Q(47,72)*Q(11,50)**2
    assert circle_gap > 0
    A_gaps = [C - (a*x*x + Q(3,2)*(1-x) + Q(26,125)) for x in (Q(11,50), Q(2,3))]
    assert A_gaps == [Q(37,20000), Q(509,18000)]
    assert a - Q(1,6) + Q(1,6000) == C and Q(1,20) - Q(1,3)/8 == Q(1,120)
    quadrature_moments = []
    for k in range(4):
        total = (Q(0), Q(0))
        for sign in (-1, 1):
            square, weight = (Q(3,7), sign*Q(2,35)), (Q(1,2), -sign*Q(1,36))
            power_k = (Q(1), Q(0))
            for _ in range(k):
                power_k = q30_mul(power_k, square)
            term = q30_mul(weight, power_k)
            total = tuple(total[j] + 2*term[j] for j in range(2))
        assert total == (Q(2,2*k+1), Q(0))
        quadrature_moments.append(str(total[0]))
    polynomial = [Q(3,35), Q(0), -Q(6,7), Q(0), Q(1)]
    norm = sum((co*cp*Q(2,i+j+1) for i,co in enumerate(polynomial)
                for j,cp in enumerate(polynomial) if (i+j)%2 == 0), Q(0))
    assert norm == Q(128,11025) and norm/2**9 == Q(1,44100)
    n, R, q = 100000, 20, 40
    eps2 = Q(n,2**80) + Q(12,2**100) + Q(4*n,2**200)
    assert eps2 < Q(1,10**18)
    W = (a+1)*(1-Q(1,n))**2 - 1 - 2*a*Q(1,10**9)
    D_upper = Q(2917,1000)/pi_lo
    ratio = W/D_upper
    assert ratio > Q(7,4)
    entry_error = Q(1,10**6)
    rational_ratio = (W-entry_error)/(D_upper+entry_error)
    assert rational_ratio > Q(7,4)
    assert Q(1625)*pi_lo/2917 > Q(7,4)
    result = dict(pi_lower=str(pi_lo), pi_upper=str(pi_hi), circle_endpoint_gap=str(circle_gap),
                  A_endpoint_gaps=list(map(str,A_gaps)), matrix_determinants=matrix_checks(),
                  quadrature_even_moments=quadrature_moments, quadrature_node_polynomial_norm=str(norm),
                  quadrature_scaled_error_coefficient="1/44100", finite_dimension=n, cutoff=R,
                  mesh_exponent=q, squared_discretization_error_upper=str(eps2),
                  vector_value_lower=str(W), scalar_norm_upper=str(D_upper),
                  finite_ratio_lower=str(ratio), finite_ratio_gap_above_7_over_4=str(ratio-Q(7,4)),
                  total_matrix_entry_error=str(entry_error), rational_matrix_ratio_lower=str(rational_ratio),
                  rational_matrix_ratio_gap_above_7_over_4=str(rational_ratio-Q(7,4)))
    Path(__file__).with_name("final_exact_02.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
