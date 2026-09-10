"""Exact finite-family Gram certificate, with complete exponential tails."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,sys
sys.dont_write_bytecode=True
from arithmetic import IV,SCALE,exp_negative_step,tail_polynomials

DATA_HASH='d0d7b6795988bf9846f3dffd55007ddaa906db45edeef1b830cba22570c1f9fe'


def load():
    raw=Path(__file__).with_name('family.json').read_bytes()
    assert hashlib.sha256(raw).hexdigest()==DATA_HASH
    data=json.loads(raw)
    lam,beta=F(*data['lambda']),F(*data['beta'])
    assert (lam,beta)==(F(53326,10**6),F(499883,10**6))
    assert (data['grid'],data['cutoff'],data['q_denominator'])==(500,24,10**8)
    assert len(data['members'])==16 and len(data['intervals'])==64
    ids=[row['id'] for row in data['members']]
    assert len(set(ids))==len(ids)
    for member in data['members']:
        b=F(*member['b'])
        assert 0<b<F(1,2) and 2*b*b<beta
        rn,inn=member['qR_numerators'],member['qI_numerators']
        assert len(rn)==len(inn)==12001
        assert all(type(x) is int and x>0 for x in rn+inn)
        assert all(r>=i for r,i in zip(rn,inn))
    for j,interval in enumerate(data['intervals']):
        assert interval['interval']==j and interval['u_interval']==[j,64,j+1,64]
        assert interval['selected_candidate'] in ids
        assert F(interval['tau_numerator'],interval['tau_denominator'])>0
    return data


def moments(data):
    grid=data['grid'];rows=[]
    exponential=exp_negative_step(F(1,grid))
    endpoint_exp=IV(1);endpoint_poly=tail_polynomials(F(0))
    for index in range(12001):
        left=[endpoint_exp*p for p in endpoint_poly]
        if index==12000:m=left
        else:
            endpoint_exp=endpoint_exp*exponential
            endpoint_poly=tail_polynomials(F(index+1,grid))
            m=[left[j]-endpoint_exp*endpoint_poly[j] for j in range(4)]
        m0,m1,m2,m3=m
        p2=m2/2;p3=m2-2*m1+m0;p4=m3/6;p5=(m3-4*m2+4*m1)/2
        for term in [m0,m1,p2,p3,p4,p5]:assert term.above(0)
        rows.append((m0,m1,p2,p3,p4,p5,m2-m1,m3-2*m2))
    return rows


def verify(emit):
    data=load();emit('data_sha256',DATA_HASH)
    emit('number_of_feasible_auxiliary_tuples',len(data['members']))
    beta=F(*data['beta']);den=data['q_denominator'];rows=moments(data)
    names=['p0','r1','i1','g2a','g2b','c2','r5','i5','g3a','g3b','c3']
    largest_diagonal=IV(0);min_det2=None;min_det3=None
    for member in data['members']:
        vals={key:IV(0) for key in names}
        for nr,ni,row in zip(member['qR_numerators'],member['qI_numerators'],rows):
            qr,qi=F(nr,den),F(ni,den)
            aa=(1/qr+1/qi)/2;dd=(1/qr-1/qi)/2
            m0,m1,p2,p3,p4,p5,c2,c3=row
            vals['p0']+=beta*aa*m0
            vals['r1']+=(beta/qr)*m1;vals['i1']+=(beta/qi)*m1
            vals['g2a']+=beta*aa*p3;vals['g2b']+=beta*aa*p2;vals['c2']+=beta*dd*c2
            vals['r5']+=(beta/qr)*p5;vals['i5']+=(beta/qi)*p5
            vals['g3a']+=beta*aa*p5;vals['g3b']+=beta*aa*p4;vals['c3']+=beta*dd*c3
        for key in ['p0','r1','i1','r5','i5','g2a','g2b','g3a','g3b']:
            assert vals[key].below(1),(member['id'],key,vals[key].text())
            if vals[key].hi>largest_diagonal.hi:largest_diagonal=vals[key]
        det2=(1-vals['g2a'])*(1-vals['g2b'])-vals['c2']*vals['c2']/2
        det3=(1-vals['g3a'])*(1-vals['g3b'])-vals['c3']*vals['c3']/12
        assert det2.above(0) and det3.above(0)
        if min_det2 is None or det2.lo<min_det2.lo:min_det2=det2
        if min_det3 is None or det3.lo<min_det3.lo:min_det3=det3
        emit('member',member['id'],'b',F(*member['b']))
        for key in names:emit(key,vals[key])
        emit('det(I-beta G2)',det2,'det(I-beta G3)',det3)
    emit('largest_diagonal_enclosure',largest_diagonal)
    emit('minimum_determinant2_enclosure',min_det2,'minimum_determinant3_enclosure',min_det3)
    emit('PASS: every frozen member satisfies all seven Gram conditions.')
    emit('The scalar cover is verified separately using the same frozen family hash.')


if __name__=='__main__':
    with Path(__file__).with_suffix('.log').open('w') as log:
        def emit(*items):
            text=' '.join(item.text() if isinstance(item,IV) else str(item) for item in items)
            print(text,flush=True);print(text,file=log,flush=True)
        verify(emit)
