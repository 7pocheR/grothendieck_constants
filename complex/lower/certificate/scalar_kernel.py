"""Exact scalar endpoint upper bounds for one rational auxiliary tuple."""
from fractions import Fraction as F
from math import isqrt
from arithmetic import IV, SCALE, exp_negative_step, tail_polynomials
MIDPOINT_ROOT_SCALE=10**8
ENDPOINT_ROOT_SCALE=10**24
MIXED_SUBDIVISIONS=16
_SUBWEIGHTS=None


def mixed_cell_upper(row,u,tau,lam):
    """Monotone upper integration on 16 exact rational subintervals.

    The original step weights remain constant. Subinterval mass divided by
    whole-cell mass is e^(-k*h/16)*(1-e^(-h/16))/(1-e^(-h)), so it is
    enclosed once with exact directed interval arithmetic.
    """
    global _SUBWEIGHTS
    fl,fr,an2,adiff,rp2,qr0,qi0,fq,vl,lm0,m0,qr,qi,denom,vr,cache=row
    if _SUBWEIGHTS is None:
        step=exp_negative_step(F(1,500*MIXED_SUBDIVISIONS))
        whole=exp_negative_step(F(1,500))
        base=(1-step)/(1-whole)
        factor=IV(1);_SUBWEIGHTS=[]
        for _ in range(MIXED_SUBDIVISIONS):
            _SUBWEIGHTS.append(factor*base);factor=factor*step
    if cache['values'] is None:
        values=[];b=cache['b'];A=1-2*b;index=cache['index']
        for k in range(1,MIXED_SUBDIVISIONS+1):
            s=F(index*MIXED_SUBDIVISIONS+k,500*MIXED_SUBDIVISIONS)
            f=s*(A+b*s)**2
            v=ceil_sqrt_fraction(f.numerator,f.denominator,ENDPOINT_ROOT_SCALE)
            values.append((f,v))
        cache['values']=values
    result=IV(0);L=tau+u/tau
    for mass,(f,v) in zip(_SUBWEIGHTS,cache['values']):
        if tau*tau*f<=denom*denom:right=max(qr,qi+u*f/denom)
        else:right=max(qr,L*v-2*lam)
        result+=mass*right
    return m0*result

def floor_sqrt_scaled(numerator, denominator, scale):
    assert numerator >= 0 and denominator > 0 and scale > 0
    root = isqrt(numerator*scale*scale//denominator)
    assert root*root*denominator <= numerator*scale*scale
    assert (root+1)*(root+1)*denominator > numerator*scale*scale
    return root


def ceil_sqrt_fraction(numerator, denominator, scale):
    root = floor_sqrt_scaled(numerator, denominator, scale)
    if root*root*denominator < numerator*scale*scale:
        root += 1
    assert root*root*denominator >= numerator*scale*scale
    return F(root, scale)


def prepare(data):
    grid, cutoff, qden = data['grid'], data['cutoff'], data['q_denominator']
    rn, inn = data['qR_numerators'], data['qI_numerators']
    lam, beta, b = (F(*data[key]) for key in ['lambda', 'beta', 'b'])
    assert (lam, beta) == (F(53326, 10**6), F(499883, 10**6))
    assert 0 < b < F(1,2)
    assert (grid, cutoff, qden) == (500, 24, 10**8)
    assert len(rn) == len(inn) == grid*cutoff+1
    assert all(type(x) is int and x > 0 for x in rn+inn)
    assert all(r >= i for r, i in zip(rn, inn))
    assert 0 < lam and 0 < beta and lam+beta < 1 and 2*b*b < beta
    A = 1-2*b
    assert A > 0 and b > 0

    # V(i/grid)^2 has one common exact integer denominator.
    bn, bd = data['b']
    An = bd-2*bn
    fden = bd*bd*grid**3
    fnums = [i*(An*grid+bn*i)**2 for i in range(grid*cutoff+1)]
    assert all(x < y for x, y in zip(fnums, fnums[1:]))
    vright = [ceil_sqrt_fraction(x, fden, ENDPOINT_ROOT_SCALE) for x in fnums]
    twice_lambda_qden = 2*lam*qden
    assert twice_lambda_qden.denominator == 1
    twice_lambda_qden = twice_lambda_qden.numerator

    exponential = exp_negative_step(F(1, grid))
    endpoint_exp = IV(1)
    endpoint_polys = tail_polynomials(F(0))
    rows = []
    tail = None
    for index, (nr, ni) in enumerate(zip(rn, inn)):
        qr, qi = F(nr, qden), F(ni, qden)
        left_tails = [endpoint_exp*p for p in endpoint_polys]
        if index == len(rn)-1:
            moments = left_tails
        else:
            endpoint_exp = endpoint_exp*exponential
            endpoint_polys = tail_polynomials(F(index+1, grid))
            moments = [left_tails[j]-endpoint_exp*endpoint_polys[j] for j in range(4)]
        m0, m1, m2, m3 = moments
        assert all(term.above(0) for term in moments)
        f_integral = A*A*m1+2*A*b*m2+b*b*m3
        assert f_integral.above(0)
        denom = 2*lam+qi
        if index == len(rn)-1:
            # On [24,infinity), max(qR,qI+Phi)<=qR+qI+u V^2/denom.
            # This bound is affine in u and is used in place of a tangent.
            tail = ((qr+qi)*m0, f_integral/denom)
            break

        root_int = floor_sqrt_scaled(2*index+1, 2*grid, MIDPOINT_ROOT_SCALE)
        assert root_int > 0
        r0 = F(root_int, MIDPOINT_ROOT_SCALE)
        # sqrt(s)<=(s+r0^2)/(2r0), multiplied by positive A+b*s.
        linear_majorant = (A*r0/2)*m0+(A/(2*r0)+b*r0/2)*m1+(b/(2*r0))*m2
        assert linear_majorant.above(0)
        an = ni+twice_lambda_qden
        rp = nr+twice_lambda_qden
        rows.append((
            fnums[index], fnums[index+1], an*an, an*(nr-ni), rp*rp,
            qr*m0, qi*m0, f_integral/denom, linear_majorant, -2*lam*m0,
            m0, qr, qi, denom, vright[index+1],
            {'index':index,'b':b,'values':None}
        ))
    assert len(rows) == grid*cutoff and tail is not None
    return lam, beta, b, fden, qden, rows, tail


def endpoint_upper(u, tau, lam, c, fden, qden, rows, tail):
    """Enclose an explicit upper bound for one convex majorant endpoint.

    The returned interval encloses the integrated polynomial/right-endpoint
    majorant, not the original scalar function from below.
    """
    assert 0 <= u <= 1 and tau > 0
    linear_coefficient = tau+u/tau
    u_iv, linear_iv = IV(u), IV(linear_coefficient)
    tau_left = tau.numerator**2*qden*qden
    tau_right = tau.denominator**2*fden
    quadratic_left = u.numerator*qden*qden
    quadratic_right = u.denominator*fden
    linear_left = linear_coefficient.numerator**2*qden*qden
    linear_right = linear_coefficient.denominator**2*fden
    total = IV(lam-c*u)+tail[0]+u_iv*tail[1]
    counts = {'inactive': 0, 'quadratic': 0, 'linear': 0, 'mixed': 0}
    for row in rows:
        fl, fr, an2, adiff, rp2, qr0, qi0, fq, vl, lm0, m0, qr, qi, denom, vr, cache = row
        if tau_left*fr <= tau_right*an2:
            # The whole cell is on the quadratic tangent branch.
            if quadratic_left*fr <= quadratic_right*adiff:
                value = qr0
                counts['inactive'] += 1
            elif quadratic_left*fl >= quadratic_right*adiff:
                value = qi0+u_iv*fq
                counts['quadratic'] += 1
            else:
                value = mixed_cell_upper(row,u,tau,lam)
                counts['mixed'] += 1
        elif tau_left*fl >= tau_right*an2:
            # The whole cell is on the linear tangent branch.
            if linear_left*fr <= linear_right*rp2:
                value = qr0
                counts['inactive'] += 1
            elif linear_left*fl >= linear_right*rp2:
                value = linear_iv*vl+lm0
                counts['linear'] += 1
            else:
                value = mixed_cell_upper(row,u,tau,lam)
                counts['mixed'] += 1
        else:
            # The tangent branch changes in this cell. Both branches, and
            # their maximum with qR, are nondecreasing functions of V.
            # The right endpoint lies on the linear branch.
            assert tau_left*fr > tau_right*an2
            value = mixed_cell_upper(row,u,tau,lam)
            counts['mixed'] += 1
        total += value
    assert sum(counts.values()) == len(rows)
    return total, counts

