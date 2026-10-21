\\ Optional rational relative-norm probe; does not solve the integral problem.
default(parisize,100000000);
p=y^2-y-1;
a=Mod(2*y,p);b=Mod(82+107*y,p);
T=rnfisnorminit(p,x^2-a*b,1);
for(q=1,2,print(q," ",rnfisnorm(T,q*b)));
