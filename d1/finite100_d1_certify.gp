\\ finite100_d1_certify.gp
\\ Rigorous exact verification of d1(n) for every 1 <= n <= 100.
\\ d1(n) = min { m > 1 : Q_n + m is prime }.
\\
\\ The script independently reconstructs the first 100 semiprimes and Q_n,
\\ checks every smaller offset for primality, verifies the recorded d1(n),
\\ verifies d1(n) < q0(n)^2, proves d1(n) itself prime, and archives an
\\ ECPP certificate for each large prime Q_n + d1(n).

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

certify_N(N,tag)=
{
  my(c,f,s);
  if(N < 2^64, if(!isprime(N),error(Str(tag,": small-N primality check failed"))); print("PROVED SMALL ",tag); return());
  print("CERTIFY ",tag,"  (",#digits(N)," digits)");
  c=primecert(N);
  if(c==0,error(Str(tag,": primecert() returned 0")));
  if(!primecertisvalid(c),error(Str(tag,": in-memory certificate invalid")));
  if(c[1][1] != N,error(Str(tag,": certificate is for the wrong integer")));
  f=fileopen(Str(tag,".cert.gp"),"w");
  filewrite(f,c);
  fileclose(f);
  c=read(Str(tag,".cert.gp"));
  if(!primecertisvalid(c),error(Str(tag,": on-disk certificate invalid")));
  if(c[1][1] != N,error(Str(tag,": on-disk certificate is for the wrong integer")));
  s=primecertexport(c,1);
  f=fileopen(Str(tag,".primo"),"w");
  filewrite(f,s);
  fileclose(f);
  print("PROVED ",tag);
};

checkrow(n)=
{
  my(Q=QS[n], sn=S[n], y=sn\2, q0=nextprime(y+1), p=W[n], m, N, tag);
  tag=Str("d1_n",n,"_m",p);
  if(!isprime(p),error(Str(tag,": recorded displacement is not prime")));
  if(p < q0,error(Str(tag,": recorded displacement lies below q0")));
  if(p >= q0^2,error(Str(tag,": recorded displacement is outside protected window")));
  for(m=2,p-1, if(isprime(Q+m),error(Str(tag,": earlier prime hit at m=",m))));
  N=Q+p;
  certify_N(N,tag);
  print("ROW PROVED n=",n,"  d1=",p,"  q0=",q0,"  q0^2=",q0^2);
};

for(n=1,100,checkrow(n));
print("ALL 100 EXACT d1 VALUES PROVED; EVERY d1(n) IS PRIME FOR 1 <= n <= 100");
quit;
