\p 80
default(parisizemax,1000000000);
f=y^4-4*y^2+2;t=Mod(y,f);d=(t-1)*(2+t);
for(k=1,2,c=if(k==1,1,3);print("FACTOR ",c);rn=rnfinit(nfinit(f),x^2-c*d);ab=rnfequation(f,x^2-c*d);print("ABSOLUTE ",ab);b=bnfinit(ab,1);for(i=1,#b.fu,u=lift(rnfeltabstorel(rn,b.fu[i]));v=subst(u,x,-x);w=lift(Mod(u,x^2-c*d)/Mod(v,x^2-c*d));print("UNIT ",i," ",w)));
quit;
