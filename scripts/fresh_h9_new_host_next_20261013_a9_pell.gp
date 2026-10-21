\p 100
default(parisizemax,1000000000);
f=y^4-4*y^2+2;t=Mod(y,f);d=10*(t-1);
nf=nfinit(f);rn=rnfinit(nf,x^2-d);ab=rnfequation(f,x^2-d);print("ABSOLUTE ",ab);b=bnfinit(ab,1);print("CERTIFY ",bnfcertify(b));print("CLASS ",b.no);print("SIGNATURE ",b.sign);
enc(v)=vector(4,j,Str(polcoef(lift(v),j-1,y)));
rel(v)={my(w=lift(rnfeltabstorel(rn,v)));[enc(polcoef(w,0,x)),enc(polcoef(w,1,x))]};
for(i=1,#b.fu,print("UNIT ",i," ",rel(b.fu[i])));
rni=rnfisnorminit(f,x^2-d,1);ss=rnfisnorm(rni,20);print("PARTICULAR ",[enc(polcoef(lift(ss[1]),0,x)),enc(polcoef(lift(ss[1]),1,x))]," ",ss[2]);
ss=bnfisintnorm(b,160000);print("INTEGRAL_COUNT ",#ss);for(j=1,#ss,print("INTEGRAL ",rel(ss[j])));
quit;
