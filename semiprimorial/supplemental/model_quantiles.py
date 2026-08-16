import math
from sympy import primerange, isprime, nextprime, primepi

GAMMA = 0.5772156649015329
EG = math.exp(GAMMA)

LIM = 400000
spf = list(range(LIM+1))
for i in range(2,int(LIM**0.5)+1):
    if spf[i]==i:
        for j in range(i*i,LIM+1,i):
            if spf[j]==j: spf[j]=i
def Omega(x):
    c=0
    while x>1:
        p=spf[x]
        while x%p==0:
            x//=p; c+=1
    return c
semis=[x for x in range(4,LIM+1) if Omega(x)==2]
logQ=[0.0]
for s in semis:
    logQ.append(logQ[-1]+math.log(s))

def sp(n): return semis[n-1]

def setup(n):
    S=sp(n); yy=S/2.0; L=logQ[n]
    p=2
    while not 2*p>S: p=nextprime(p)
    Q0=p
    rho=EG*math.log(yy)
    u=L/math.log(yy)
    return S,yy,L,Q0,rho,u

def vp_is_one(p,S):
    # v_p(Q_n)=1  iff  pi(S/p)+[p^2<=S] == 1
    return int(primepi(S//p)) + (1 if p*p<=S else 0) == 1

def hazards(n, tmax_mult=8.0):
    """return sorted list of (m, prob, channel) for offsets m>=2 up to tmax"""
    S,yy,L,Q0,rho,u = setup(n)
    base = rho/L
    sem  = base*math.log(u-1)
    tmax = int(tmax_mult*yy)
    smallprimes=set(primerange(2,int(yy)+1))
    out=[]
    for m in range(2,tmax+1):
        if isprime(m):
            if m<=yy: out.append((m,base,'I'))
            else:      out.append((m,sem,'II'))
            continue
        # composite: factor
        x=m; f={}
        while x>1:
            p=spf[x]; c=0
            while x%p==0: x//=p; c+=1
            f[p]=c
        small=[p for p in f if p<=yy]
        if len(small)!=1: continue
        p=small[0]; a=f[p]
        rest=m//(p**a)
        if a==1:
            if rest>1 and isprime(rest) and rest>=Q0 and m<Q0*Q0:
                out.append((m,base,'III'))
        elif a==2 and rest==1:
            if vp_is_one(p,S) and m<Q0*Q0:
                out.append((m,base,'III2'))
    return out,(S,yy,L,Q0,rho,u)

def analyse(n):
    ev,(S,yy,L,Q0,rho,u)=hazards(n)
    H=0.0; surv=1.0
    q25=q50=q75=None
    Hy=None; H2y=None
    pcomp=0.0
    for m,p,ch in ev:
        if Hy is None and m>yy: Hy=H
        if H2y is None and m>2*yy: H2y=H
        if ch.startswith('III'):
            pcomp += math.exp(-H)*p
        H+=p
        if q25 is None and H>=-math.log(0.75): q25=m/yy
        if q50 is None and H>=math.log(2): q50=m/yy
        if q75 is None and H>=-math.log(0.25): q75=m/yy
    return dict(n=n,S=S,y=yy,L=L,q0=Q0,u=u,logu1=math.log(u-1),
                HI_y=Hy, H_2y=H2y, En=H2y-Hy,
                q25=q25,q50=q50,q75=q75,
                surv2y=math.exp(-H2y), surv2y_noI=math.exp(-(H2y-Hy)),
                pcomp=pcomp)

print("=== per-n model quantiles of a(n)/y (calibration range) ===")
print(f"{'n':>4} {'H_I(y)':>7} {'q25':>6} {'med':>6} {'q75':>6} {'P(a>y)':>7}")
rows=[]
for n in range(66,94):
    r=analyse(n); rows.append(r)
    print(f"{n:>4} {r['HI_y']:>7.3f} {r['q25']:>6.3f} {r['q50']:>6.3f} {r['q75']:>6.3f} {math.exp(-r['HI_y']):>7.3f}")
import statistics as st
print("\nmedian across n of model q25/q50/q75 :",
      round(st.median([r['q25'] for r in rows]),3),
      round(st.median([r['q50'] for r in rows]),3),
      round(st.median([r['q75'] for r in rows]),3))
print("observed (paper)                     : 0.50  1.10  1.71")
print("mean H_I(y) over 66..93:", round(st.mean([r['HI_y'] for r in rows]),3))

print("\n=== corrected composite-term probabilities ===")
print(f"{'n':>5} {'E_n':>6} {'e^-En':>7} {'H_I(y)':>7} {'e^-(H_I+En)':>11} {'paper P':>8} {'exact P':>8}")
paper={100:0.0095,250:0.0066,1000:0.0041,3000:0.0028}
for n in [100,250,1000,3000]:
    r=analyse(n)
    print(f"{n:>5} {r['En']:>6.3f} {r['surv2y_noI']:>7.4f} {r['HI_y']:>7.3f} {r['surv2y']:>11.4f} {paper[n]:>8.4f} {r['pcomp']:>8.5f}")
