\\ Exact Hasse-Minkowski triage only; no integral lattice identification.
f=x^3-3*x-1; nf=nfinit(f); t=Mod(x,f); b=4-t^2;
hasse(v,pr)={my(s=1);for(i=1,#v,for(j=i+1,#v,s*=nfhilbert(nf,v[i],v[j],pr)));s};
for(k=0,1,a=2*t*b^k; olda=2*t*b^(1-k); v=concat([-a],concat(vector(6,i,if(i%2,2,b/2)),[1,1,1])); w=concat([-olda],concat(vector(8,i,if(i%2,2,b/2)),[1])); ws=if(k,b*w,w); for(p=2,3,prs=idealprimedec(nf,p);for(j=1,#prs,print([k,p,j,hasse(v,prs[j]),hasse(w,prs[j]),hasse(ws,prs[j])]))));
quit;
