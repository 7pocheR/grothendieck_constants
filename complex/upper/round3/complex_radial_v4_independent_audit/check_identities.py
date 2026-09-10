"""Small exact symbolic checks; no numerical boundary computation."""
import argparse
import ast
from fractions import Fraction
from hashlib import sha256
from itertools import permutations
import json
from math import comb, factorial
from pathlib import Path
import sympy as sp

from audit_exact import load_json, phase_modes, require


def determinant(matrix):
    n = len(matrix)
    result = 0
    for p in permutations(range(n)):
        inversions = sum(p[i] > p[j] for i in range(n) for j in range(i+1, n))
        term = (-1)**inversions
        for i in range(n):
            term *= matrix[i][p[i]]
        result += term
    return sp.expand(result)


def multiply(A, B):
    return [[sp.expand(sum(A[i][k]*B[k][j] for k in range(len(B))))
             for j in range(len(B[0]))] for i in range(len(A))]


def independent_primitives():
    a, b, h, r, s, alpha, beta = symbols = sp.symbols("a b h r s alpha beta")
    S = [[1,0,a,0],[0,1,0,b],[a,0,1,0],[0,b,0,1]]
    J = [[-alpha,-r,0,0],[r,-alpha,0,0],[0,0,-beta,-s],[0,0,s,-beta]]
    U = [[1,h,0,0],[0,0,1,-h]]
    SJ = multiply(S, J)
    M = [[int(i==j)-SJ[i][j] for j in range(4)] for i in range(4)]
    D = determinant(M)
    adj = [[(-1)**(i+j)*determinant([[M[k][l] for l in range(4) if l != i]
                                    for k in range(4) if k != j])
            for j in range(4)] for i in range(4)]
    UT = list(map(list, zip(*U)))
    N = multiply(multiply(multiply(U, adj), S), UT)
    return symbols, (D, N[0][0], N[1][1], N[0][1], N[1][0])


def isolated_functions(source):
    tree = ast.parse(Path(source).read_text())
    functions = []
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in ("polys", "exact_modes"):
            functions.append(node)
    module = ast.fix_missing_locations(ast.Module(body=functions, type_ignores=[]))
    namespace = dict(F=Fraction, factorial=factorial)
    exec(compile(module, str(source), "exec"), namespace)
    series_node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "series")
    prefix = []
    for stmt in series_node.body:
        prefix.append(stmt)
        if isinstance(stmt, ast.Assign) and isinstance(stmt.targets[0], ast.Name) and stmt.targets[0].id == "Z":
            break
    prefix.append(ast.parse("return D, N1, N2, V+W, V-W, Z").body[0])
    series_node.body = prefix
    module = ast.fix_missing_locations(ast.Module(body=[series_node], type_ignores=[]))
    namespace["arb_series"] = lambda values: sp.Symbol("y")
    exec(compile(module, str(source), "exec"), namespace)
    return namespace


def run(source, output):
    source, output = Path(source), Path(output)
    output.mkdir(parents=True, exist_ok=False)
    symbols, primitives = independent_primitives()
    a,b,h,r,s,alpha,beta = symbols
    names = ("D", "N1", "N2", "P0", "Q0")
    author = isolated_functions(source / "certify_radial.py")
    compact = author["polys"](a,b,h,r,s,alpha,beta)
    checks = {}
    for name, actual, claimed in zip(names, primitives, compact):
        require(sp.Poly(actual-claimed, *symbols).is_zero, f"Failed generic determinant identity: {name}")
        checks[name] = "exact polynomial identity"
    swap = {r:s, s:r, alpha:beta, beta:alpha}
    expected_swap = (0,2,1,3,4)
    for i,name in enumerate(names):
        require(sp.Poly(primitives[i].xreplace(swap)-primitives[expected_swap[i]], *symbols).is_zero,
                f"Failed row/column symmetry: {name}")
    A,B,y,z = sp.symbols("A B y z")
    claimed = author["series"](A,B,h,r,s,alpha,beta,4)
    for i,name in enumerate(names):
        actual = primitives[i].subs({a:z*A,b:z*B})
        if i >= 3:
            actual = sp.cancel(actual/z)
        require(sp.Poly(actual-claimed[i].subs(y,z*z), z,A,B,h,r,s,alpha,beta).is_zero,
                f"Failed parity substitution: {name}")
    require(sp.cancel(claimed[5]-y*claimed[3]*claimed[4]/(claimed[1]*claimed[2])) == 0,
            "Hypergeometric argument mismatch")
    data = load_json(source / "radial_low_candidate.json")
    modes = phase_modes(data, data["phase_order"])
    require(modes == author["exact_modes"](data,data["phase_order"]), "Mode weight identity failed")
    coefficient = Fraction(1)
    for j in range(421):
        if j:
            coefficient *= Fraction((2*j-1)**2,4*j*(j+1))
        require(coefficient == Fraction(comb(2*j,j)**2,16**j*(j+1)), f"Hypergeometric coefficient {j}")
    polynomial_data = {name:[dict(powers=list(powers),coefficient=str(c))
                             for powers,c in sp.Poly(value,*symbols).terms()]
                       for name,value in zip(names,primitives)}
    formula = dict(variables=list(map(str,symbols)),
                   definition="D=det(I-SJ); [[N1,P0],[Q0,N2]]=U adj(I-SJ) S U^T",
                   polynomials=polynomial_data)
    with (output / "primitive_polynomials.json").open("x") as f:
        json.dump(formula,f,indent=2);f.write("\n")
    result = dict(generic_primitive_identities=checks, exact_row_column_symmetry=True,
                  exact_even_odd_series_substitution=True, independently_reconstructed_modes=len(modes),
                  exact_mode_weight_agreement=True, hypergeometric_coefficients_checked=421,
                  primitive_term_counts={name:len(v) for name,v in polynomial_data.items()},
                  source_sha256=sha256((source / "certify_radial.py").read_bytes()).hexdigest())
    with (output / "identities.json").open("x") as f:
        json.dump(result,f,indent=2);f.write("\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--source",required=True)
    parser.add_argument("--output",required=True)
    args=parser.parse_args()
    run(args.source,args.output)
