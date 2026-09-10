"""Arb certificate of all forward coefficients through N for the fixed candidate.
Composite Gauss--Legendre with one uniform proved Bernstein-ellipse error for
all Hermite integrals. Matrix products and all final sums use Arb.
"""
from flint import arb,arb_mat,ctx
from pathlib import Path
from math import factorial,comb
import math,json,time,argparse
from scipy.special import roots_legendre
HERE=Path(__file__).resolve().parent
DATA=json.loads((HERE/'candidate_rational.json').read_text())
ctx.prec=384
C=[arb(v) for v in DATA['boundary_coefficients']]
N=251;NG=112;PANELS=128;X=arb(16);RHO=arb(3)/2;TAU=arb(99)/100
SQ=[arb(j).sqrt() for j in range(1001)]
NU=(2/arb.pi()).sqrt();PHI=1/(2*arb.pi()).sqrt()
def hermites(x,M):
 H=[arb(1),x]
 for j in range(1,M):H.append((x*H[-1]-SQ[j]*H[-2])/SQ[j+1])
 return H[:M+1]
def Lder(x,d):
 H=hermites(x,9)
 return sum((C[(n-1)//2]*(arb(factorial(n))/factorial(n-d)).sqrt()*H[n-d] for n in range(1,10,2) if n>=d),arb(0))
def legendre(x,n):
 pm=arb(1);p=x
 for k in range(1,n):pm,p=p,((2*k+1)*x*p-k*pm)/(k+1)
 dp=n*(x*p-pm)/(x*x-1)
 return p,dp

def gauss_nodes():
 oldprec=ctx.prec;ctx.prec=640
 nodes=[];weights=[];records=[];previous_hi=None
 approx,_=roots_legendre(NG)
 assert len(approx)==NG
 for z in approx:
  lo=arb(float(z))-arb(2)**-38;hi=arb(float(z))+arb(2)**-38
  assert -1<lo and lo<hi and hi<1
  if previous_hi is not None:assert previous_hi<lo
  previous_hi=hi;initial_lo=lo;initial_hi=hi
  flo=legendre(lo,NG)[0];fhi=legendre(hi,NG)[0]
  assert flo*fhi<0
  # The disjoint sign brackets contain all NG roots of P_NG.
  for _ in range(260):
   mid=(lo+hi)/2;fm=legendre(mid,NG)[0]
   if fm*flo>0:lo=mid;flo=fm
   else:
    assert fm*fhi>0
    hi=mid;fhi=fm
  zz=arb.union(lo,hi);mid=(lo+hi)/2;rad=(hi-lo)/2
  # Markov endpoint value bounds |P_NG second derivative| on [-1,1].
  derivative_lipschitz=arb((NG-1)*NG*(NG+1)*(NG+2))/8
  dp=legendre(mid,NG)[1]+arb(0,(derivative_lipschitz*rad).upper())
  weight=2/((1-zz*zz)*dp*dp)
  assert weight>0
  nodes.append(zz);weights.append(weight)
  records.append({"initial_lo":save_ball(initial_lo),"initial_hi":save_ball(initial_hi),"root":save_ball(zz),"weight":save_ball(weight)})
 assert len(records)==NG and sum(weights,arb(0)).contains(2)
 (HERE/f"gauss_nodes_{NG}.json").write_text(json.dumps({"degree":NG,"precision_bits":ctx.prec,"all_initial_brackets_in_minus1_1":True,"all_initial_brackets_disjoint":True,"all_sign_changes_checked":True,"all_weights_positive":True,"weights_sum":save_ball(sum(weights,arb(0))),"records":records},indent=2)+"\n")
 ctx.prec=oldprec
 return nodes,weights

def abs_lower(x):
 if x>0:return arb(x.lower())
 if x<0:return arb((-x).lower())
 return arb(0)
def ellipse_error():
 hh=X/(2*PANELS);xr=hh*(RHO+1/RHO)/2;yi=hh*(RHO-1/RHO)/2
 radius=(xr*xr+yi*yi).sqrt();sumM=arb(0);bounds=[]
 factor=TAU**(-arb(N)/2)*(1-TAU*TAU)**(-arb(1)/4)
 for panel in range(PANELS):
  center=(2*panel+1)*hh
  ders=[Lder(center,j) for j in range(10)]
  dev=sum((abs(ders[j])*radius**j/factorial(j) for j in range(1,10)),arb(0))
  lum=abs(ders[0])-dev;reU=arb(lum.lower()) if lum>0 else arb(0)
  derbound=sum((abs(ders[j+1])*radius**j/factorial(j) for j in range(9)),arb(0))
  imU=yi*derbound
  qa=NU*factor*(-reU*reU/(2*(1+TAU))+imU*imU/(2*(1-TAU))).exp()
  qzero=1+NU*imU*((imU*imU-reU*reU)/2).exp()
  qbound=arb(max(qa.upper(),qzero.upper()))
  reX=center-xr;reX=arb(reX.lower()) if reX>0 else arb(0)
  hphi=PHI*factor*(-reX*reX/(2*(1+TAU))+yi*yi/(2*(1-TAU))).exp()
  m=arb((qbound*hphi).upper());sumM+=m;bounds.append(str(m))
 err=2*8*hh*sumM*RHO**(-2*NG)/(1-1/RHO)
 tail=(X/SQ[2]).erfc().sqrt()
 return arb((err+tail).upper()),bounds

def save_ball(z):
 m,r,e=z.mid_rad_10exp();return {'m':str(m),'r':str(r),'e':int(e)}
def load_ball(z):return arb(z['m']+'e'+str(z['e']),z['r']+'e'+str(z['e']))
def run():
 t0=time.time();nodes,weights=gauss_nodes();err,bounds=ellipse_error()
 print('N',N,'nodes',NG,'panels',PANELS,'uniform_A_error',err,flush=True)
 assert err<arb('1e-14')
 evens=list(range(0,N+1,2));odds=list(range(1,N+1,2));Ae=arb_mat(len(evens),len(odds));Ao=arb_mat(len(odds),len(evens))
 chunk=4;hh=X/(2*PANELS)
 for begin in range(0,PANELS,chunk):
  mnode=NG*min(chunk,PANELS-begin)
  Qe=arb_mat(len(evens),mnode);Qo=arb_mat(len(odds),mnode);He=arb_mat(mnode,len(evens));Ho=arb_mat(mnode,len(odds))
  ni=0
  for panel in range(begin,min(begin+chunk,PANELS)):
   center=(2*panel+1)*hh
   for zz,ww in zip(nodes,weights):
    xx=center+hh*zz;H=hermites(xx,N);u=Lder(xx,0);Q=[(u/SQ[2]).erf(),NU*(-u*u/2).exp()]
    for a in range(1,N):Q.append((-u*SQ[a]*Q[a]-(a-1)*Q[a-1])/(SQ[a]*SQ[a+1]))
    weight=2*hh*ww*PHI*(-xx*xx/2).exp()
    for j,a in enumerate(evens):Qe[j,ni]=weight*Q[a];He[ni,j]=H[a]
    for j,a in enumerate(odds):Qo[j,ni]=weight*Q[a];Ho[ni,j]=H[a]
    ni+=1
  Ae+=Qe*Ho;Ao+=Qo*He
  print('panels',begin+chunk,'elapsed',round(time.time()-t0,2),flush=True)
 out={}
 for i,a in enumerate(evens):
  for j,b in enumerate(odds):
   if a+b<=N:out[f'{a},{b}']=save_ball(Ae[i,j]+arb(0,err))
 for i,a in enumerate(odds):
  for j,b in enumerate(evens):
   if a+b<=N:out[f'{a},{b}']=save_ball(Ao[i,j]+arb(0,err))
 output={'N':N,'NG':NG,'PANELS':PANELS,'precision_bits':ctx.prec,'uniform_A_error':str(err),'ellipse_M':bounds,'A':out,'seconds':time.time()-t0}
 (HERE/f'head_moments_{N}.json').write_text(json.dumps(output,indent=2)+'\n')
 print('CERTIFIED_MOMENTS_SAVED',len(out),'seconds',time.time()-t0,flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--N',type=int,default=251);ap.add_argument('--nodes',type=int,default=112);args=ap.parse_args();N=args.N;NG=args.nodes;run()
