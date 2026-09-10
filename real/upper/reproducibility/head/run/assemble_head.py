"""Assemble a rigorous forward correlation head from certified Hermite moments."""
from flint import arb,arb_series,ctx
from pathlib import Path
import json,time
HERE=Path(__file__).resolve().parent
ctx.prec=160
DATA=json.loads((HERE/'candidate_rational.json').read_text())
GRID=json.loads((HERE/'head_moments_501.json').read_text());N=GRID['N'];ctx.cap=N+1

def read(z):return arb(z['m']+'e'+str(z['e']),z['r']+'e'+str(z['e']))
def save(z):
 m,r,e=z.mid_rad_10exp();return {'m':str(m),'r':str(r),'e':int(e)}

def run():
 t0=time.time();A={tuple(map(int,k.split(','))):read(v) for k,v in GRID['A'].items()}
 p=[arb(0)]*(N+1)
 for n,c in enumerate(DATA['W_preprocessing_odd_coefficients']):p[2*n+1]=arb(c)
 P=arb_series(p);powP=arb_series([1]);H=arb_series([0])
 for a in range(N+1):
  row=[arb(0)]*(N+1)
  for b in range(N+1-a):
   if (a,b) in A:row[b]=(-1)**b*A[a,b]**2
  H+=powP*arb_series(row);powP*=P
 H*=arb.pi()/2
 head=sum((abs(H[n]) for n in range(3,N+1,2)),arb(0))
 for n in range(0,N+1,2):assert H[n].contains(0)
 out={'N':N,'b1':save(H[1]),'nonlinear_head_upper':str(head.upper()),'coefficients':{str(n):save(H[n]) for n in range(1,N+1,2)},'seconds':time.time()-t0}
 (HERE/'forward_head_501.json').write_text(json.dumps(out,indent=2)+'\n')
 print('b1',H[1]);print('nonlinear_head',head)
 for n in [3,5,7,9,11,31,61,81,83,101,151,201,251]:print(n,H[n])
 print('seconds',time.time()-t0)
if __name__=='__main__':run()
