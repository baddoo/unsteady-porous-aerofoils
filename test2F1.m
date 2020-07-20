z = linspace(0,1); z(1) = []; z(end) = [];

a = .5+.4i; b = .3+.5i; c = .5+.3i;
z = -1+eps
tic
f1 = hypergeom([a,b],c,z);
toc
tic
f2 = hypergeom2F1(a,b,c,z);
toc

err = norm(f1-f2,'inf')