% flow around circle plot
wd = pwd;
libloc = [wd '/../commoncode']; % put the folder "commoncode" on the path
addpath(libloc)

lsd = 100;
epsilonfactor = 0.5;
spacing = 0.1;
%vext = @(x,y) [x,-y];
%vext = @(x,y) [y.^2 x.^2]
% vext = @(x,y) [y.^2-y,x.^2+x]   % some ranges of flow patterns
vext = @(x,y) [y,x]

[Fx_ext,Fy_ext,Lz_ext,Sij_ext,Fxs,Fys,Lzs,Sijs,xs,ys,ux,uy,fx,fy,Mij,Vxs,Vys,Omegas,Vx_ext,Vy_ext,Omega_ext] = regstokes_circle_refiner_compute_V_Omega(lsd,epsilonfactor,spacing,vext);
Ux = Vx_ext;
Uy = Vy_ext;
Omega = Omega_ext;
%fprintf('Ux = %3.3g, Uy = %3.3g, Omega = %3.3g \n',Ux,Uy,Omega');
numlines = 20;
thetas = linspace(0,2*pi,numlines);
Lh = 2;



[xsgrid,ysgrid] = meshgrid(linspace(-Lh,Lh,25));
[xs,ys,fx,fy,vxgrid,vygrid] = get_flow_around_circle_specifymore(lsd,epsilonfactor,spacing,xsgrid,ysgrid,xs,ys,fx,fy);

fprintf('Ux = %3.3g, Uy = %3.3g, Omega = %3.3g \n',Ux,Uy,Omega);


clf;


vbackground = vext(xsgrid(:),ysgrid(:));
vbackx = reshape(vbackground(:,1),size(xsgrid));
vbacky = reshape(vbackground(:,2),size(ysgrid));

vscale = 0.13/max(sqrt(vbackx(:).^2+vbacky(:).^2));

titlefontsize = 28;



linestartsx = [linspace(-Lh,Lh,numlines) linspace(-Lh,Lh,numlines)];
linestartsy = [ones(1,numlines)*Lh*0.95 ones(1,numlines)*(-Lh)*0.95];

subplot(1,2,1)
q1 = quiver(xsgrid,ysgrid,vscale*vbackx,vscale*vbacky,0);
set(q1,'LineWidth',1);
hstream = streamline(xsgrid,ysgrid,vscale*(vbackx),vscale*(vbacky),linestartsx,linestartsy);
set(hstream,'color','red','linewidth',3);
axis equal
axis tight
%axis off
set(gca,'XTick',[]); set(gca,'YTick',[]);
title('Background flow','FontSize',titlefontsize)

subplot(1,2,2)

q2 = quiver(xsgrid,ysgrid,vscale*(vxgrid+vbackx),vscale*(vygrid+vbacky),0);
set(q2,'LineWidth',1);
hold on


hstream = streamline(xsgrid,ysgrid,vscale*(vxgrid+vbackx),vscale*(vygrid+vbacky),linestartsx,linestartsy);
set(hstream,'color','red','linewidth',3);
plot(xs,ys,'ko');
axis equal
axis tight
set(gca,'XTick',[]); set(gca,'YTick',[]);
title('Flow with inclusion','FontSize',titlefontsize)