% compare traces
wd = pwd;
libloc = [wd '/../commoncode'];
addpath(libloc)
syms x y
%%%% TO CHANGE THE FLOW FIELD VEXT, ALTER THESE LINES DEFINING vext_sym
%vext_sym = [-y-(1/8)*y.^2 x-(1/8)*x.^2];
vext_sym = [-(y.^3)./(1+y.^4) (x.^3)./(1+x.^4)]
%vext_sym = [-y x]
vext = matlabFunction(vext_sym);

% how long to integrate for
dt = 0.1;
steps = 150;

% initial (x,y) starting position
cmstart = [0 2];

[xss_DO,yss_DO,thetas_DO] = integrate_faxen_DO_theory(cmstart,dt,steps,vext_sym);

epsilonfactor = 0.5;
spacings = 0.1;
%spacings = [0.1 0.2];

% choose lsd = 100a, integrate faxen equations with this assumption -- this
% allows us to match to a << Lsd Oppenheimer-Diamant theory
[xss100,yss100,thetas100] = integrate_faxen(cmstart,dt,steps,100,epsilonfactor,spacings,vext);


%% plotting section

figure
%subplot(1,2,1)
maxpoint = max([max(abs(xss100)) max(abs(yss100)) max(abs(xss_DO)) max(abs(yss_DO))]);
[xsgrid,ysgrid] = meshgrid(linspace(-maxpoint,maxpoint,25));
vbackground = vext(xsgrid(:),ysgrid(:));
vbackx = reshape(vbackground(:,1),size(xsgrid));
vbacky = reshape(vbackground(:,2),size(ysgrid));

qq=quiver(xsgrid,ysgrid,vbackx,vbacky);
set(qq,'LineWidth',1,'color','r');
hold on

thf = linspace(0,2*pi,100);

plot(xss_DO,yss_DO,'-','LineWidth',8,'color',[0.8 0.8 0.8]);
% patch(xss_DO(end)+cos(thf),yss_DO(end)+sin(thf),'-','LineWidth',8,'color',[0.8 0.8 0.8])
patch(xss_DO(end)+cos(thf),yss_DO(end)+sin(thf),[0.8 0.8 0.8]);
%hold on
plot(xss100,yss100,'b-.','LineWidth',2);
%plot(xss100(end)+cos(thf),yss100(end)+sin(thf),'b','LineWidth',2)
patch(xss100(end)+cos(thf),yss100(end)+sin(thf),'b');

axis equal
axis tight
set(gca,'FontSize',32)
%plot(xss0p01,yss0p01,'r.-','LineWidth',2);
%subplot(1,2,2)
figure

plot(thetas_DO,'k--','LineWidth',3);
hold on
plot(thetas100,'b','LineWidth',4);
%plot(thetas0p01,'r.-','LineWidth',2);