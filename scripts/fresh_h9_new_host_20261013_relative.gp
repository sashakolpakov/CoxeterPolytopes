\p 100
default(parisizemax,1000000000);
f=y^4-4*y^2+2;t=Mod(y,f);d=6*(t-1)*(t+2);
nf=nfinit(f);rn=rnfinit(nf,x^2-d);ab=rnfequation(f,x^2-d);print("ABSOLUTE ",ab);b=bnfinit(ab,1);print("CERTIFY ",bnfcertify(b));print("CLASS ",b.no);print("SIGNATURE ",b.sign);
enc(v)=vector(4,j,Str(polcoef(lift(v),j-1,y)));
rel(v)={my(w=lift(rnfeltabstorel(rn,v)));[enc(polcoef(w,0,x)),enc(polcoef(w,1,x))]};
for(i=1,#b.fu,print("UNIT ",i," ",rel(b.fu[i])));
ns=[1,16,64,81,256,324,625,1296,2592,4096,5184,6561,10368,20736];
for(i=1,#ns,ss=bnfisintnorm(b,ns[i]);print("NORM ",ns[i]," ",#ss);for(j=1,#ss,print("SOLUTION ",ns[i]," ",rel(ss[j]))));
quit;
