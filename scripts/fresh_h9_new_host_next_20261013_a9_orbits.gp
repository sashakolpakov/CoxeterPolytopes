\p 100
default(parisizemax,1000000000);
f=y^4-4*y^2+2;t=Mod(y,f);d=10*(t-1);
nf=nfinit(f);rn=rnfinit(nf,x^2-d);ab=rnfequation(f,x^2-d);print("ABSOLUTE ",ab);b=bnfinit(ab,1);print("CERTIFY ",bnfcertify(b));print("CLASS ",b.no);print("SIGNATURE ",b.sign);
enc(v)=vector(4,j,Str(polcoef(lift(v),j-1,y)));
rel(v)={my(w=lift(rnfeltabstorel(rn,v)));[enc(polcoef(w,0,x)),enc(polcoef(w,1,x))]};
for(i=1,#b.fu,print("UNIT ",i," ",rel(b.fu[i])));
ns=[16,64,256,1024,4096,10000,20000,20736,40000,80000,160000,320000,640000,1280000,2560000,6250000,12500000,25000000,50000000,100000000,200000000,400000000,800000000,1600000000];
for(i=1,#ns,ss=bnfisintnorm(b,ns[i]);print("NORM ",ns[i]," ",#ss);for(j=1,#ss,print("SOLUTION ",ns[i]," ",rel(ss[j]))));
quit;
