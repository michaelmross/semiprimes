\ qplus1_semiprime9_certify.gp
\ Rigorous certification of the nine Q_n+1 semiprime candidates
\ n = 15,16,19,21,23,27,29,43,65.
\ Reconstructs the first 65 semiprimes and Q_n independently.

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

certify_prime(p,tag)=
{
  my(c,f,s);
  if(p < 2^64,
    if(!isprime(p),error(Str(tag,": small primality proof failed")));
    print("PROVED SMALL ",tag);
    return()
  );
  print("CERTIFY ",tag,"  (",#digits(p)," digits)");
  c=primecert(p);
  if(c==0,error(Str(tag,": primecert() returned 0")));
  if(!primecertisvalid(c),error(Str(tag,": in-memory certificate invalid")));
  if(c[1][1] != p,error(Str(tag,": certificate is for wrong integer")));
  f=fileopen(Str(tag,".cert.gp"),"w");
  filewrite(f,c);
  fileclose(f);
  c=read(Str(tag,".cert.gp"));
  if(!primecertisvalid(c),error(Str(tag,": on-disk certificate invalid")));
  if(c[1][1] != p,error(Str(tag,": on-disk certificate is for wrong integer")));
  s=primecertexport(c,1);
  f=fileopen(Str(tag,".primo"),"w");
  filewrite(f,s);
  fileclose(f);
  print("PROVED ",tag);
};

checkrow(n,a,b)=
{
  my(N=QS[n]+1,tag=Str("qplus1_n",n));
  if(a*b != N,error(Str(tag,": factor product mismatch")));
  certify_prime(a,Str(tag,"_f1"));
  certify_prime(b,Str(tag,"_f2"));
  print("ROW PROVED n=",n," SEMIPRIME");
};

checkrow(15,75727279,104685162319);
checkrow(16,254329,1433835837549769);
checkrow(19,11949359,4194494712496025039);
checkrow(21,929,178365721756162223150467169);
checkrow(23,2113,316033163611567942445025802177);
checkrow(27,417414594349,51576268034339176788146858149);
checkrow(29,269,585036066791804707525369245229345356133829);
checkrow(43,17862581810753,257799112548968383148018972589842706489167603522885861972417);
checkrow(65,194098513184483,143758554825735157140023967361174230521984398298852971994665174174841083313703928737290083574855805263507147);
print("ALL 9 Q_n+1 SEMIPRIME CANDIDATES RIGOROUSLY PROVED");
quit;
