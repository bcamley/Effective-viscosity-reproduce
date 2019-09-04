function [xss,yss,thetas] = integrate_faxen_DO_theory(cmstart,dt,steps,vext_sym)
syms x y
vext_x = vext_sym(1)
vext_y = vext_sym(2)

clf;
xss = zeros(1,steps);
yss = zeros(1,steps);
thetas = zeros(1,steps);

cmloc = cmstart;
theta = 0;

[xsgrid,ysgrid] = meshgrid(linspace(-3,3,25));
vext = matlabFunction(vext_sym);
vbackground = vext(xsgrid(:),ysgrid(:));
vbackx = reshape(vbackground(:,1),size(xsgrid));
vbacky = reshape(vbackground(:,2),size(ysgrid));

thf = linspace(0,2*pi,100);


% test velocity divergence
div_v = simplify(diff(vext_x,x)+diff(vext_y,y))
if(~((simplify(div_v==0))))
    error('Velocity field is not divergence-free!!!')
end

vlap_x = diff(diff(vext_x,x),x)+diff(diff(vext_x,y),y);
vlap_y = diff(diff(vext_y,x),x)+diff(diff(vext_y,y),y);

vpred_x_sym = vext_x + (1/4)*vlap_x;
vpred_y_sym = vext_y + (1/4)*vlap_y;

%vpred_x = vext_x_zero + (1/4)*vlap_x_zero;
%vpred_y = vext_y_zero + (1/4)*vlap_y_zero;

curl_v = diff(vext_y,x)-diff(vext_x,y);
lap_curl_v = diff(diff(curl_v,x),x)+diff(diff(curl_v,y),y);

omegapred_sym = (1/2)*curl_v + (1/2)*(1/8)*lap_curl_v;

vpred_x = matlabFunction(vpred_x_sym,'Vars',[x y]);
vpred_y = matlabFunction(vpred_y_sym,'Vars',[x y]);
omegapred = matlabFunction(omegapred_sym,'Vars',[x y]);

for s = 1:steps
xss(s) = cmloc(1);
yss(s) = cmloc(2);
thetas(s) = theta;
%vpred_x = subs(subs(vpred_x_sym,x,cmloc(1)),y,cmloc(2));
%vpred_y = subs(subs(vpred_y_sym,x,cmloc(1)),y,cmloc(2));
%omegapred = subs(subs(omegapred_sym,x,cmloc(1)),y,cmloc(2));

%[Fx_ext,Fy_ext,Lz_ext,Sij_ext,Fxs,Fys,Lzs,Sijs,xs,ys,ux,uy,fx,fy,Mij,Vxs,Vys,Omegas,Vx_ext,Vy_ext,Omega_ext] = regstokes_circle_refiner_compute_V_Omega(lsd,epsilonfactor,spacings,vext,cmloc);
cmloc = cmloc + dt*[vpred_x(cmloc(1),cmloc(2)) vpred_y(cmloc(1),cmloc(2))];
theta = theta + dt*omegapred(cmloc(1),cmloc(2));

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

