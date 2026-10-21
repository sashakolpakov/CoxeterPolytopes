\p 100
default(parisizemax,1000000000);
f=y^4-4*y^2+2;t=Mod(y,f);bk=bnfinit(f,1);print("BASE_CERTIFY ",bnfcertify(bk));
enc(v)=vector(4,j,Str(polcoef(lift(v),j-1,y)));
for(i=1,#bk.fu,print("BASE_UNIT ",i," ",enc(bk.fu[i])));
d=6*(t-1);nf=nfinit(f);rn=rnfinit(nf,x^2-d);ab=rnfequation(f,x^2-d);print("ABSOLUTE ",ab);b=bnfinit(ab,1);print("CERTIFY ",bnfcertify(b));print("CLASS ",b.no);print("SIGNATURE ",b.sign);
rel(v)={my(w=lift(rnfeltabstorel(rn,v)));[enc(polcoef(w,0,x)),enc(polcoef(w,1,x))]};
for(i=1,#b.fu,v=lift(rnfeltabstorel(rn,b.fu[i]));aa=polcoef(v,0,x);bb=polcoef(v,1,x);print("UNIT ",i," ",rel(b.fu[i])," MAP ",Vec(bnfisunit(bk,aa^2-d*bb^2))));
ns=[1,256,625,1296,5184,6561,10000,20736,40000,50625,65536,82944,160000,194481,202500,331776,390625,640000,810000,1620000,1679616,2560000,3240000,4100625,6480000,12960000,25920000,51840000,103680000,207360000,506250000,1012500000,2025000000,4050000000,8100000000,16200000000,32400000000,64800000000,129600000000];
for(i=1,#ns,ss=bnfisintnorm(b,ns[i]);print("NORM ",ns[i]," ",#ss);for(j=1,#ss,print("SOLUTION ",ns[i]," ",rel(ss[j]))));
quit;
