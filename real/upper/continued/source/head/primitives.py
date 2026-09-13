"""Independent interval primitives. Physical execution belongs on Linux only."""
from fractions import Fraction as Q
from math import comb, factorial
from flint import arb, acb, arb_poly, arb_series, fmpz_poly, ctx
from exact import require, hermite_coefficients, laurent

N, F, M = 403, 768, 31
NODES = (1, 13, 127, 255, 384, 512, 767, 768)
CHECK_DEGREES = (0, 1, 17, 127, 402, 403)


def configure():
    ctx.prec, ctx.cap, ctx.threads = 512, N+1, 1
    require(ctx.prec == 512 and ctx.cap == 404 and ctx.threads == 1, 'Wrong context')


def a(q):
    return arb(str(Q(q)))


def endpoint(x, upper):
    m, e = map(int, (x.upper() if upper else x.lower()).man_exp())
    return Q(m)*Q(2)**e


def encode(x):
    require(x.is_finite(), 'Nonfinite interval')
    return [str(endpoint(x, False)), str(endpoint(x, True))]


def decode(v):
    require(isinstance(v, list) and len(v) == 2, 'Malformed dyadic endpoints')
    lo, hi = map(Q, v)
    require(lo <= hi, 'Reversed interval')
    return a((lo+hi)/2)+arb(0, a((hi-lo)/2).upper())


def producer_ball(v):
    require(set(v) == {'m', 'r', 'e'} and int(v['r']) >= 0, 'Malformed producer ball')
    scale = Q(10)**int(v['e'])
    return a(int(v['m'])*scale)+arb(0, a(int(v['r'])*scale).upper())


def enlarged(x, error):
    require(error >= 0, 'Negative analytic error')
    return x+arb(0, a(error).upper())


def compare(left, right, label):
    require(len(left) == len(right), 'Comparison length mismatch')
    gaps, widths = [], []
    for k, (x, y) in enumerate(zip(left, right)):
        require(x.is_finite() and y.is_finite(), 'Nonfinite comparison')
        # An explicit endpoint test, used only as a consistency diagnostic.
        xl, xu, yl, yu = endpoint(x, False), endpoint(x, True), endpoint(y, False), endpoint(y, True)
        require(max(xl, yl) <= min(xu, yu), label+': disjoint coefficient '+str(k)
                +'; left='+str([str(xl), str(xu)])+'; right='+str([str(yl), str(yu)]))
        gaps.append(max(abs(xl-yu), abs(xu-yl)))
        widths.append(xu-xl)
    return {'coefficients_checked': len(left), 'maximum_difference_upper': str(max(gaps, default=Q(0))),
            'maximum_left_width': str(max(widths, default=Q(0))),
            'meaning': 'consistency only; no endpoint is replaced by an intersection'}


def series(c, degree=N):
    require(ctx.cap >= degree+1, 'Global cap too small')
    return checked(arb_series(c, prec=degree+1), degree)


def checked(s, degree=N):
    require(ctx.cap >= degree+1 and s.prec >= degree+1, 'Silent series truncation')
    return s


def compose_block(coeff, P, degree=N, block=16):
    """Evaluate sum_q P^(block*q) sum_r coeff[block*q+r] P^r."""
    require(len(coeff) == degree+1 and P[0].is_zero(), 'Invalid composition')
    powers = [series([arb(1)], degree)]
    for j in range(block):
        powers.append(checked(powers[-1]*P, degree))
    result = series([], degree)
    for q in reversed(range((degree+block)//block)):
        part = series([], degree)
        for r in range(min(block, degree+1-q*block)):
            part = checked(part+checked(powers[r]*coeff[q*block+r], degree), degree)
        result = checked(checked(result*powers[block], degree)+part, degree)
    return result


def pair(pu, pv, vu, vv, P, degree=N, alternate=False):
    require(all(len(v) == degree+1 for v in (pu, pv, vu, vv)), 'Incomplete pair inputs')
    W = [vu[k]*vv[k] for k in range(degree+1)]
    K = series([(-1)**k*pu[k]*pv[k] for k in range(degree+1)], degree)
    if alternate:
        C = checked(series(W, degree)(P), degree)
        K4 = checked(K**4, degree)
    else:
        C = compose_block(W, P, degree)
        K2 = checked(K*K, degree)
        K4 = checked(K2*K2, degree)
    out = checked(C*K4, degree)
    return [out[k] for k in range(degree+1)]


def feature_errors(data, par):
    """Independent positive spatial-error evaluation, at the maximum degree/frequency."""
    h, y, t = a(par['x_spacing']), a(par['strip_height_x']), a(par['Mehler_t_real'])
    T = F*2*arb.pi()/a(par['period'])
    DN = sum((a(Q(comb(N, r), factorial(r)))*y**(2*r) for r in range(N+1)), arb(0))
    den = (2*arb.pi()*y/h).exp()-1
    gx = 2*(y*y/2).exp()/den
    cs = [abs(a(c)) for c in data['sine_coefficients']]
    ws = [abs(a(Q(data['frequency_step'])*(r+1))) for r in range(len(cs))]
    imag = abs(a(data['linear']))*y+sum((c*(w*y).sinh() for c, w in zip(cs, ws)), arb(0))
    strip = 2*DN.sqrt()*(y*y/2+T*imag).exp()/den
    cutoff = h*par['x_cutoff']
    tail = (1-t*t)**(-a(Q(1, 4)))*t**(-a(Q(N, 2)))*(1+t).sqrt()*(cutoff/(2*(1+t)).sqrt()).erfc()
    mass = ((1+gx)*(1+gx*DN)).sqrt()
    j = par['bin_Taylor_degree']
    remainder = mass*(T*a(par['bin_radius']))**(j+1)/factorial(j+1)
    return {key: str(endpoint(v, True)) for key, v in
            [('strip', strip), ('tail', tail), ('bin', remainder)]}


def direct_features(data, par, nodes, tick=lambda: None):
    """Signed lattice sums, all degrees; explicit polynomial checks at every x node.

    Unnormalized integer Hermite coefficients and a single sqrt(k!) normalization
    differ from the producer's normalized recurrence with repeated square roots.
    Selected degrees use explicit polynomial evaluations for the actual sum.
    """
    require(1 <= len(nodes) <= 8, 'Feature batch exceeds eight nodes')
    h, delta = a(par['x_spacing']), 2*arb.pi()/a(par['period'])
    ell = a(data['linear'])
    cs = list(map(a, data['sine_coefficients']))
    ws = [a(Q(data['frequency_step'])*(r+1)) for r in range(len(cs))]
    roots = [arb(factorial(k)).sqrt() for k in range(N+1)]
    explicit = {k: arb_poly(hermite_coefficients(k)) for k in CHECK_DEGREES}
    result = {n: [arb(0) for k in range(N+1)] for n in nodes}
    maximum_hermite_gap = Q(0)
    for ix in range(par['x_cutoff']+1):
        if ix % 16 == 0:
            tick()
        x = ix*h
        B = ell*x+sum((c*(w*x).sin() for c, w in zip(cs, ws)), arb(0))
        phi = h*(-x*x/2).exp()/(2*arb.pi()).sqrt()*(1 if ix == 0 else 2)
        trig = {n: ((n*delta*B).cos(), (n*delta*B).sin()) for n in nodes}
        previous, current = arb(0), arb(1)
        for k in range(N+1):
            hk = current
            if k in explicit:
                hk = explicit[k](x)
                require((hk-current).contains(0), 'Explicit Hermite polynomial disagrees')
                maximum_hermite_gap = max(maximum_hermite_gap,
                                         endpoint(abs((hk-current)/roots[k]), True))
            weight = phi*hk/roots[k]*(-1)**(k//2)
            for n in nodes:
                result[n][k] += weight*trig[n][k % 2]
            previous, current = current, x*current-k*previous
    return result, {'positive_lattice_nodes': par['x_cutoff']+1,
                    'symmetric_lattice_nodes': 2*par['x_cutoff']+1,
                    'explicit_polynomial_degrees': list(CHECK_DEGREES),
                    'explicit_polynomial_checks': (par['x_cutoff']+1)*len(CHECK_DEGREES),
                    'maximum_normalized_hermite_difference_upper': str(maximum_hermite_gap)}


def profile(u, rows, degree=N):
    phi0 = 1/(2*arb.pi()).sqrt()
    components = []
    for s in sorted(set().union(*rows)):
        x = u-2*s
        # Regularize before evaluation; never divide a shifted polynomial by x.
        v = arb(0)
        for m in range(len(rows)-1, 0, -1):
            v = v*x+a(rows[m].get(s, 0))
        components.append((x, phi0*(-x*x/2).exp()*v))
    gaussian = phi0*(-u*u/2).exp()
    out = [arb(0)]*(degree+1)
    if not u.contains(0):
        out[0] = gaussian/u+sum((v for x, v in components), arb(0))
    else:
        require(u.is_zero(), 'Unresolved principal-value node')
        regular = sum((v for x, v in components), arb(0))
        require(regular.contains(0), 'V0 regular part parity at zero failed')
    for j in range(1, degree+1):
        # Fresh integer exponentiation and one factorial normalization per j.
        value = gaussian*u**(j-1)+sum((v*x**j for x, v in components), arb(0))
        out[j] = value/arb(factorial(j)).sqrt()
    return out


def signed_reflection_checks(data, par, fresh, tick=lambda: None):
    """Separate unsymmetrized complex sums at both extreme frequencies and zero."""
    frequencies = (0, 1, -1, 768, -768)
    delta, h = 2*arb.pi()/a(par['period']), a(par['x_spacing'])
    polynomials = {k: arb_poly(hermite_coefficients(k)) for k in CHECK_DEGREES}
    roots = {k: arb(factorial(k)).sqrt() for k in CHECK_DEGREES}
    values = {(n, k): acb(0) for n in frequencies for k in CHECK_DEGREES}
    cs = list(map(a, data['sine_coefficients']))
    ws = [a(Q(data['frequency_step'])*(r+1)) for r in range(len(cs))]
    for ix in range(-par['x_cutoff'], par['x_cutoff']+1):
        if ix % 16 == 0:
            tick()
        x = ix*h
        B = a(data['linear'])*x+sum((c*(w*x).sin() for c, w in zip(cs, ws)), arb(0))
        weight = h*(-x*x/2).exp()/(2*arb.pi()).sqrt()
        exp = {n: acb(0, n*delta*B).exp() for n in frequencies}
        for k in CHECK_DEGREES:
            hk = weight*polynomials[k](x)/roots[k]
            for n in frequencies:
                values[n, k] += hk*exp[n]
    error = feature_errors(data, par)
    allowance = Q(error['strip'])+Q(error['tail'])
    records = []
    for n in frequencies:
        for k in CHECK_DEGREES:
            rotated = values[n, k]/acb(0, 1)**k
            require(rotated.imag.contains(0), 'Unsymmetrized sum violates Hermite phase parity')
            if n:
                target = fresh[abs(n)][k]*((-1)**k if n < 0 else 1)
                diagnostic = compare([rotated.real], [target], 'Signed lattice/reflection')
            else:
                diagnostic = compare([enlarged(rotated.real, allowance)],
                                     [arb(1 if k == 0 else 0)], 'Continuous zero-frequency normalization')
            records.append({'node': n, 'degree': k, 'real': encode(rotated.real),
                            'imaginary': encode(rotated.imag), 'consistency': diagnostic})
    return {'signed_lattice_nodes': 2*par['x_cutoff']+1, 'records': records,
            'zero_frequency_analytic_allowance': str(allowance)}


def zero_vectors(data, rows, delta, degree=N):
    v0 = profile(arb(0), rows, degree)
    primary = [delta*v0[k] if k % 2 else arb(0) for k in range(degree+1)]
    for k in range(2, degree+1, 2):
        require(v0[k].contains(0), 'Primary zero parity failed')
    phi0 = 1/(2*arb.pi()).sqrt()
    cs = list(map(a, data['sine_coefficients']))
    ws = [a(Q(data['frequency_step'])*(r+1)) for r in range(len(cs))]
    auxiliary = [arb(0)]*(degree+1)
    for k in range(1, degree+1, 2):
        derivative = sum((c*w**k*(-w*w/2).exp() for c, w in zip(cs, ws)), arb(0))/arb(factorial(k)).sqrt()
        if k == 1:
            derivative += a(data['linear'])
        auxiliary[k] = delta*phi0*derivative
    return primary, auxiliary


def corrections(p, v, z, t, delta, P, degree=N):
    require(set(p) == set(v), 'Correction node sets differ')
    primary, auxiliary = [arb(0)]*(degree+1), [arb(0)]*(degree+1)
    for k in range(1, degree+1, 2):
        aa = 2*delta*sum((v[n][k]*p[n][0]**4 for n in sorted(p)), arb(0))
        bb = 2*delta*sum((v[n][0]*p[n][k]*p[n][0]**3 for n in sorted(p)), arb(0))
        primary[k] = 2*aa*z[k]+z[k]*z[k]
        auxiliary[k] = -4*(2*bb*t[k]+t[k]*t[k])
    out = checked(compose_block(primary, P, degree)+series(auxiliary, degree), degree)
    return [out[k] for k in range(degree+1)]
