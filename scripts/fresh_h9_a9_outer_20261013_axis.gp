\p 100
default(parisizemax,1000000000);
f=y^4-4*y^2+2;t=Mod(y,f);bk=bnfinit(f,1);print("BASE_CERTIFY ",bnfcertify(bk));
enc(v)=vector(4,j,Str(polcoef(lift(v),j-1,y)));
for(i=1,#bk.fu,print("BASE_UNIT ",i," ",enc(bk.fu[i])));
d=10*(t-1);nf=nfinit(f);rn=rnfinit(nf,x^2-d);ab=rnfequation(f,x^2-d);print("ABSOLUTE ",ab);b=bnfinit(ab,1);print("CERTIFY ",bnfcertify(b));print("CLASS ",b.no);print("SIGNATURE ",b.sign);
rel(v)={my(w=lift(rnfeltabstorel(rn,v)));[enc(polcoef(w,0,x)),enc(polcoef(w,1,x))]};
for(i=1,#b.fu,v=lift(rnfeltabstorel(rn,b.fu[i]));aa=polcoef(v,0,x);bb=polcoef(v,1,x);print("UNIT ",i," ",rel(b.fu[i])," MAP ",Vec(bnfisunit(bk,aa^2-d*bb^2))));
ns=[10000,40000,160000,640000,2560000,6250000,12500000,12960000,25000000,50000000,100000000,200000000,400000000,800000000,1600000000,3906250000,7812500000,15625000000,31250000000,62500000000,125000000000,250000000000,500000000000,1000000000000];
for(i=1,#ns,ss=bnfisintnorm(b,ns[i]);print("NORM ",ns[i]," ",#ss);for(j=1,#ss,print("SOLUTION ",ns[i]," ",rel(ss[j]))));
quit;
