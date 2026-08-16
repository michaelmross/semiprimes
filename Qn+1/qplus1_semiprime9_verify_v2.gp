\\ qplus1_semiprime9_verify.gp
\\ Independent verification of the saved certificates for the nine Q_n+1 hits.

is_semiprime_small(k)=
{
  my(F=factor(k), s=0, r);
  r=matsize(F)[1];
  for(i=1,r, s += F[i,2]);
  return(s==2);
};

S=List(); k=4;
while(#S<65, if(is_semiprime_small(k),listput(S,k)); k++);
Q=1; QS=vector(65);
for(i=1,65, Q*=S[i]; QS[i]=Q);

verify_prime(p,tag)=
{
  my(c);
  if(p < 2^64,
    if(!isprime(p),error(Str(tag,": small primality proof failed")));
    return()
  );
  c=read(Str(tag,".cert.gp"));
  if(!primecertisvalid(c),error(Str(tag,": INVALID SAVED CERTIFICATE")));
  if(c[1][1] != p,error(Str(tag,": saved certificate is for wrong integer")));
};

verifyrow(n,a,b)=
{
  my(N=QS[n]+1,tag=Str("qplus1_n",n));
  if(a*b != N,error(Str(tag,": factor product mismatch")));
  verify_prime(a,Str(tag,"_f1"));
  verify_prime(b,Str(tag,"_f2"));
  print("VALID n=",n," SEMIPRIME");
};

verifyrow(15,75727279,104685162319);
verifyrow(16,254329,1433835837549769);
verifyrow(19,11949359,4194494712496025039);
verifyrow(21,929,178365721756162223150467169);
verifyrow(23,2113,316033163611567942445025802177);
verifyrow(27,417414594349,51576268034339176788146858149);
verifyrow(29,269,585036066791804707525369245229345356133829);
verifyrow(43,17862581810753,257799112548968383148018972589842706489167603522885861972417);
verifyrow(65,194098513184483,143758554825735157140023967361174230521984398298852971994665174174841083313703928737290083574855805263507147);
print("ALL 9 Q_n+1 SEMIPRIME FACTORIZATIONS AND SAVED CERTIFICATES VALID");
quit;
