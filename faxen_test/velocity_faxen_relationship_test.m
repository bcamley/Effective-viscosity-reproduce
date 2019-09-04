% generate velocity faxen relationship ... 
% 
wd = pwd;
libloc = [wd '/../commoncode']; % put the folder "commoncode" on the path
addpath(libloc)

lsds = logspace(-2,2,25)
epsilonfactor = 0.5
%spacings = 0.03:0.02:0.2 % hires
spacings = 0.05:0.05:0.2

syms x y
%vext_x = y.^2-y;
%vext_y = x.^2+x;

vext_x = -(1/24)*(y-1).^3+(1/8)*(y-1).^2-(y-1);
vext_y = (1/24)*x.^3+(1/8)*x.^2-x;

%vext_x = +(1/8)*(y-1).^2-(y-1);
%vext_y = +(1/8)*x.^2-x;

% test velocity divergence
div_v = simplify(diff(vext_x,x)+diff(vext_y,y))
if(~((simplify(div_v==0))))
    error('Velocity field is not divergence-free!!!')
end

vlap_x = diff(diff(vext_x,x),x)+diff(diff(vext_x,y),y);
vlap_y = diff(diff(vext_y,x),x)+diff(diff(vext_y,y),y);

vlap_x_zero = subs(subs(vlap_x,x,0),y,0);
vlap_y_zero = subs(subs(vlap_y,x,0),y,0);

vext_sym = [vext_x , vext_y];
vext = matlabFunction(vext_sym);

vext_x_zero = subs(subs(vext_x,x,0),y,0);
vext_y_zero = subs(subs(vext_y,x,0),y,0);

vpred_x = vext_x_zero + (1/4)*vlap_x_zero;
vpred_y = vext_y_zero + (1/4)*vlap_y_zero;

curl_v = diff(vext_y,x)-diff(vext_x,y);
lap_curl_v = diff(diff(curl_v,x),x)+diff(diff(curl_v,y),y);

omegapred_sym = (1/2)*curl_v + (1/2)*(1/8)*lap_curl_v;
omegapred = subs(subs(omegapred_sym,x,0),y,0);

[Fx_ext,Fy_ext,Lz_ext,Sij_ext,Fxs,Fys,Lzs,Sijs,xs,ys,ux,uy,fx,fy,Mij,Vxs,Vys,Omegas,Vx_ext,Vy_ext,Omega_ext] = regstokes_circle_refiner_compute_V_Omega(lsds,epsilonfactor,spacings,vext);

%% Plotting of results + extrapolation
clf
ms = 15; lw = 3;
fs = 24;

subplot(1,3,1)
plot(1./lsds,vpred_x*ones(size(lsds)),'k--','LineWidth',lw);
hold on
plot(1./lsds,Vx_ext,'bo','MarkerSize',ms,'LineWidth',lw);
legend({'Prediction (Oppenheimer-Diamant)','Regularized Stokeslets'})
set(gca,'xscale','log')
set(gca,'FontSize',fs)
xlabel('a/L_{sd}'); ylabel('U_x')
subplot(1,3,2)
plot(1./lsds,vpred_y*ones(size(lsds)),'k--','LineWidth',lw);
hold on
plot(1./lsds,Vy_ext,'bo','MarkerSize',ms,'LineWidth',lw);
set(gca,'xscale','log')
set(gca,'FontSize',fs)
xlabel('a/L_{sd}'); ylabel('U_y')
subplot(1,3,3)
plot(1./lsds,omegapred*ones(size(lsds)),'k--','LineWidth',lw);
hold on
plot(1./lsds,Omega_ext,'bo','MarkerSize',ms,'LineWidth',lw);
set(gca,'xscale','log')
set(gca,'FontSize',fs)
xlabel('a/L_{sd}'); ylabel('\Omega')
figure
numsub = ceil(sqrt(length(Omega_ext)));
% Do a check that extrapolation on Omega appears to be good - if these
% linear fits are not great, prediction as s->0 will also not be good!
for jj = 1:numsub*numsub
    subplot(numsub,numsub,jj);
    p = polyfit(spacings,Omegas(:,jj).',1);
    plot(spacings,Omegas(:,jj),'ko');
    hold on
    plot(spacings,polyval(p,spacings));
end


figure
[xsgrid,ysgrid] = meshgrid(linspace(-2,2,25));
vbackground = vext(xsgrid(:),ysgrid(:));
vbackx = reshape(vbackground(:,1),size(xsgrid));
vbacky = reshape(vbackground(:,2),size(ysgrid));

qq=quiver(xsgrid,ysgrid,vbackx,vbacky);
set(qq,'LineWidth',1);
hold on
thf = linspace(0,2*pi,100);
plot(cos(thf),sin(thf),'k-','LineWidth',3)
axis equal
axis tight
set(gca,'FontSize',fs);