function [Fx,Fy,Lz,Sij] = get_forces(Mij,xs,ys,Ux,Uy,Omega,Eij)

np = length(xs);

u = zeros(2*np,1);
u(1:np) = Ux - Omega*ys.' + Eij(1,1)*(xs.') + Eij(1,2)*(ys.');  % the x velocity
u(np+1:end) = Uy + Omega*(xs.') + Eij(2,2)*(ys.') + Eij(2,1)*(xs.');       % y velocity
%        u(1:np) = Ux-vxext-Omega*ys.';  
%        u(np+1:end) = Uy-vyext+Omega*xs.';
f = gmres(Mij,u,[],[],500,[],[],randn(size(u)));
fx = f(1:np).';
fy = f(np+1:end).';

Fx = sum(fx);
Fy = sum(fy);
Lz = sum(xs.*fy -ys.*fx);
Sij = 0.5*[sum(xs.*fx)*2 sum(xs.*fy+ys.*fx) ; sum(ys.*fx+xs.*fy) sum(ys.*fy)*2];