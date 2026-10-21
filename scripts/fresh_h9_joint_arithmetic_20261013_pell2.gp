\\ PARI is a proposal mechanism only; the Python artifact checks the answer
\\ with exact field arithmetic and all prior Coxeter pairings.
f=y^3-3*y-1; t=Mod(y,f);
delta=120-50*t;
rn=rnfisnorminit(f,x^2-t*delta,1);
print(rnfisnorm(rn,-5*t));
