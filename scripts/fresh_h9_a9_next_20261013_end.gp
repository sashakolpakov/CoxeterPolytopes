\p 100
default(parisizemax,1000000000);
f=y^4-4*y^2+2;t=Mod(y,f);bk=bnfinit(f,1);print("BASE_CERTIFY ",bnfcertify(bk));
enc(v)=vector(4,j,Str(polcoef(lift(v),j-1,y)));
for(i=1,#bk.fu,print("BASE_UNIT ",i," ",enc(bk.fu[i])));
d=1*(t-1);nf=nfinit(f);rn=rnfinit(nf,x^2-d);ab=rnfequation(f,x^2-d);print("ABSOLUTE ",ab);b=bnfinit(ab,1);print("CERTIFY ",bnfcertify(b));print("CLASS ",b.no);print("SIGNATURE ",b.sign);
rel(v)={my(w=lift(rnfeltabstorel(rn,v)));[enc(polcoef(w,0,x)),enc(polcoef(w,1,x))]};
for(i=1,#b.fu,v=lift(rnfeltabstorel(rn,b.fu[i]));aa=polcoef(v,0,x);bb=polcoef(v,1,x);print("UNIT ",i," ",rel(b.fu[i])," MAP ",Vec(bnfisunit(bk,aa^2-d*bb^2))));
ns=[1,256,625,1024,1296,2500,4096,6561,10000,16384,38416,40000,50625,65536,160000,320000,331776,390625,640000,810000,1280000,2560000,5120000,10240000,20480000,40960000,100000000,200000000,400000000,800000000,1600000000,3200000000,6400000000,12800000000,25600000000];
for(i=1,#ns,ss=bnfisintnorm(b,ns[i]);print("NORM ",ns[i]," ",#ss);for(j=1,#ss,print("SOLUTION ",ns[i]," ",rel(ss[j]))));
quit;
