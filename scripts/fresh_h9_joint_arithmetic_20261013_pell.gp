\\ Proposal only. The saved integral root is verified without PARI.
f=y^3-3*y-1; t=Mod(y,f); rn=rnfisnorminit(f,x^2-30*t,1);
print(rnfisnorm(rn,-5*t));
quit;
