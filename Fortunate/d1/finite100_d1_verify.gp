\\ finite100_d1_verify.gp
\\ Independent re-verification of the saved d1(n) certificate package.

is_semiprime_small(k)=
{
  my(F=factor(k), s=0, r);
  r=matsize(F)[1];
  for(i=1,r, s += F[i,2]);
  return(s==2);
};

S=List(); k=4;
while(#S<100, if(is_semiprime_small(k),listput(S,k)); k++);
Q=1; QS=vector(100);
for(i=1,100, Q*=S[i]; QS[i]=Q);

W=[3,5,7,19,13,17,17,13,17,23,17,37,31,41,43,97,139,47,109,37,61,73,37,109,59,79,227,61,127,71,59,131,269,79,73,79,149,71,61,97,181,73,241,193,197,163,101,139,241,541,109,149,229,149,139,293,191,227,97,659,173,823,241,109,103,139,353,193,431,229,149,229,607,337,653,467,421,173,271,337,139,139,181,283,449,967,829,157,181,199,911,163,281,331,181,229,227,157,251,251];

verify_N(N,tag)=
{
  my(c);
  if(N < 2^64, if(!isprime(N),error(Str(tag,": small-N primality check failed"))); return());
  c=read(Str(tag,".cert.gp"));
  if(!primecertisvalid(c),error(Str(tag,": INVALID SAVED CERTIFICATE")));
  if(c[1][1] != N,error(Str(tag,": saved certificate is for the wrong integer")));
};

verifyrow(n)=
{
  my(Q=QS[n], sn=S[n], y=sn\2, q0=nextprime(y+1), p=W[n], m, N, tag);
  tag=Str("d1_n",n,"_m",p);
  if(!isprime(p),error(Str(tag,": displacement is not prime")));
  if(p < q0,error(Str(tag,": displacement lies below q0")));
  if(p >= q0^2,error(Str(tag,": displacement lies outside protected window")));
  for(m=2,p-1, if(isprime(Q+m),error(Str(tag,": earlier prime hit at m=",m))));
  N=Q+p;
  verify_N(N,tag);
  print("VALID n=",n,"  d1=",p);
};

for(n=1,100,verifyrow(n));
print("ALL 100 EXACT d1 VALUES AND SAVED CERTIFICATES VALID");
quit;
