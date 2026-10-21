\p 100
default(parisizemax,1000000000);
f=y^4-4*y^2+2;t=Mod(y,f);d=1*(t-1);
nf=nfinit(f);rn=rnfinit(nf,x^2-d);ab=rnfequation(f,x^2-d);
print("ABSOLUTE ",ab);b=bnfinit(ab,1);print("CERTIFY ",bnfcertify(b));print("CLASS ",b.no);print("SIGNATURE ",b.sign);
enc(v)=vector(4,j,Str(polcoef(lift(v),j-1,y)));
rel(v)={my(w=lift(rnfeltabstorel(rn,v)));[enc(polcoef(w,0,x)),enc(polcoef(w,1,x))]};
bk=bnfinit(f,1);print("BASE_CERTIFY ",bnfcertify(bk));
for(i=1,#bk.fu,print("BASE_UNIT ",i," ",enc(bk.fu[i])));
for(i=1,#b.fu,v=lift(rnfeltabstorel(rn,b.fu[i]));nu=polcoef(v,0,x)^2-d*polcoef(v,1,x)^2;print("UNIT ",i," ",rel(b.fu[i]));print("UNIT_NORM ",i," ",enc(nu)," ",bnfisunit(bk,nu)~));
ns=[1, 4, 16, 64, 81, 162, 256, 324, 648, 1296, 2592, 5184, 6561, 10368, 13122, 20736, 26244, 41472, 52488, 82944, 104976, 165888, 209952, 331776, 419904, 839808, 1679616, 3359232, 6718464, 13436928, 26873856];
for(i=1,#ns,ss=bnfisintnorm(b,ns[i]);print("NORM ",ns[i]," ",#ss);for(j=1,#ss,print("SOLUTION ",ns[i]," ",rel(ss[j]))));
quit;
