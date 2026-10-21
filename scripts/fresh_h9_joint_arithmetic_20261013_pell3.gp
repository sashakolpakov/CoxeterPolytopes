\\ PARI proposes a solution; the saved root is replayed independently in Python.
default(parisizemax,1000000000);
f=y^3-3*y-1;t=Mod(y,f);
delta=Mod((5820)+(-3760)*y+(430)*y^2,f);
rn=rnfisnorminit(f,x^2-t*delta,1);
sol=rnfisnorm(rn,-5*t); print(sol);
