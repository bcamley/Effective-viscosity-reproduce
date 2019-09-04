function [xss,yss,thetas] = integrate_faxen(cmstart,dt,steps,lsd,epsilonfactor,spacings,vext)

clf;
xss = zeros(1,steps);
yss = zeros(1,steps);
thetas = zeros(1,steps);

cmloc = cmstart;
theta = 0;

[xsgrid,ysgrid] = meshgrid(linspace(-3,3,25));
vbackground = vext(xsgrid(:),ysgrid(:));
vbackx = reshape(vbackground(:,1),size(xsgrid));
vbacky = reshape(vbackground(:,2),size(ysgrid));

thf = linspace(0,2*pi,100);

for s = 1:steps
xss(s) = cmloc(1);
yss(s) = cmloc(2);
thetas(s) = theta;
[Fx_ext,Fy_ext,Lz_ext,Sij_ext,Fxs,Fys,Lzs,Sijs,xs,ys,ux,uy,fx,fy,Mij,Vxs,Vys,Omegas,Vx_ext,Vy_ext,Omega_ext] = regstokes_circle_refiner_compute_V_Omega(lsd,epsilonfactor,spacings,vext,cmloc);
cmloc = cmloc + dt*[Vx_ext Vy_ext];
theta = theta + dt*Omega_ext;

clf
subplot(1,2,1)
quiver(xsgrid,ysgrid,vbackx,vbacky);
hold on
plot(xss(1:s),yss(1:s),'k','LineWidth',2)
plot(xss(s)+cos(thf),yss(s)+sin(thf),'r','LineWidth',3)
axis equal
subplot(1,2,2)
plot(thetas(1:s))
drawnow

end

