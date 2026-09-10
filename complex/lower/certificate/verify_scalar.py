"""Exact all-u scalar certificate for a finite family of auxiliary weights.

The operator parameters lambda,beta are common to all family members. Each
u interval uses exactly one member and one tangent at BOTH endpoints.
"""
from pathlib import Path
from fractions import Fraction as F
import sys
sys.dont_write_bytecode=True
from arithmetic import IV,SCALE
from verify_gram import load,DATA_HASH
from scalar_kernel import prepare,endpoint_upper

SCALAR_TARGET=F(3447,5000)
RATIO_TARGET=F(1373,1000)


def verify(emit):
    data=load();emit('data_sha256',DATA_HASH)
    emit('proof_domain_u','[0,1]','cover_intervals',64,'scalar_target',SCALAR_TARGET)
    lam,beta=F(*data['lambda']),F(*data['beta'])
    best=None;seen=set()
    aggregate={'inactive':0,'quadratic':0,'linear':0,'mixed':0}
    for member in data['members']:
        member_id=member['id']
        intervals=[row for row in data['intervals'] if row['selected_candidate']==member_id]
        assert intervals
        tuple_data={key:data[key] for key in ['lambda','beta','grid','cutoff','q_denominator']}
        tuple_data.update({key:member[key] for key in ['b','qR_numerators','qI_numerators']})
        emit('preparing_member',member_id,'assigned_intervals',[row['interval'] for row in intervals])
        lv,bv,b,fden,qden,rows,tail=prepare(tuple_data)
        assert (lv,bv)==(lam,beta)
        c=1-2*b*b/beta
        for interval in intervals:
            index=interval['interval']
            tau=F(interval['tau_numerator'],interval['tau_denominator'])
            for side in [0,1]:
                u=F(index+side,64)
                assert (index,side) not in seen;seen.add((index,side))
                bound,counts=endpoint_upper(u,tau,lam,c,fden,qden,rows,tail)
                assert bound.below(SCALAR_TARGET),(member_id,index,side,bound.text())
                for key in aggregate:aggregate[key]+=counts[key]
                emit('member',member_id,'interval',index,'endpoint_u',u,'tangent_tau',tau,
                     'upper_majorant_enclosure',bound,'cell_counts',counts)
                if best is None or bound.hi>best[0].hi:best=(bound,member_id,index,side,u,tau)
    assert seen=={(index,side) for index in range(64) for side in [0,1]}
    bound,member_id,index,side,u,tau=best
    emit('largest_certified_endpoint_upper',F(bound.hi,SCALE))
    emit('largest_upper_outward_decimal',bound)
    emit('largest_location member interval side u tau',member_id,index,side,u,tau)
    emit('aggregate_cell_counts',aggregate)
    emit('scalar_target_margin_lower',IV(SCALAR_TARGET)-IV.raw(bound.hi,bound.hi))
    emit('limiting_ratio_from_coarse_scalar_target',(1-lam)/SCALAR_TARGET)
    d=10**9
    Q=F(2*d,2*d+1)*(1-beta/F(4*(d+1)))-lam-F(2,10**12)
    assert max(abs(1-lam),abs(lam+beta),abs(lam))<1
    assert Q>RATIO_TARGET*SCALAR_TARGET
    emit('finite_numerator_lower_for_explicit_mesh',Q)
    emit('finite_comparison_margin_Q_minus_1373_over_1000_times_D',Q-RATIO_TARGET*SCALAR_TARGET)
    emit('PASS: the complete parameter range is covered using one feasible tuple per interval.')
    emit('Scalar norm <= 3447/5000; combined with the separate Gram certificate, finite ratio > 1373/1000.')


if __name__=='__main__':
    with Path(__file__).with_suffix('.log').open('w') as log:
        def emit(*items):
            text=' '.join(item.text() if isinstance(item,IV) else str(item) for item in items)
            print(text,flush=True);print(text,file=log,flush=True)
        verify(emit)
