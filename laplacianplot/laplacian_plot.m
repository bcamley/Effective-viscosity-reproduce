%vext_sym = [-y-(1/8)*y.^2 x-(1/8)*x.^2];
syms x y
vext_sym = [-(y.^3)./(1+y.^4) (x.^3)./(1+x.^4)]
%vext_sym = [-y x]
vext = matlabFunction(vext_sym);

mr = 2.5;

[xsgrid,ysgrid] = meshgrid(linspace(-mr,mr,25));

vext_x_sym = vext_sym(1)
vext_y_sym = vext_sym(2)

vext_x = matlabFunction(vext_x_sym,'Vars',[x y]);
vext_y = matlabFunction(vext_y_sym,'Vars',[x y]);

vext_x_grid = vext_x(xsgrid,ysgrid);
vext_y_grid = vext_y(xsgrid,ysgrid);

% test velocity divergence
div_v = simplify(diff(vext_x_sym,x)+diff(vext_y_sym,y))
if(~((simplify(div_v==0))))
    error('Velocity field is not divergence-free!!!')
end

vlap_x_sym = diff(diff(vext_x_sym,x),x)+diff(diff(vext_x_sym,y),y);
vlap_y_sym = diff(diff(vext_y_sym,x),x)+diff(diff(vext_y_sym,y),y);

vlap_x = matlabFunction(vlap_x_sym,'Vars',[x y]);
vlap_y = matlabFunction(vlap_y_sym,'Vars',[x y]);

vlap_x_grid = vlap_x(xsgrid,ysgrid);
vlap_y_grid = vlap_y(xsgrid,ysgrid);

vpred_x_sym = vext_x_sym + (1/4)*vlap_x_sym;
vpred_y_sym = vext_y_sym + (1/4)*vlap_y_sym;

vpred_x = matlabFunction(vpred_x_sym,'Vars',[x y]);
vpred_y = matlabFunction(vpred_y_sym,'Vars',[x y]);

vpred_x_grid = vpred_x(xsgrid,ysgrid);
vpred_y_grid = vpred_y(xsgrid,ysgrid);

qw = 2;
fs = 22;

subplot(1,3,1)
q=quiver(xsgrid,ysgrid,vext_x_grid,vext_y_grid);
set(q,'LineWidth',qw);
set(gca,'FontSize',fs);
axis equal
xlim([-mr mr]);
ylim([-mr mr]);
title('Background velocity')

subplot(1,3,2)
q=quiver(xsgrid,ysgrid,vlap_x_grid,vlap_y_grid);
set(q,'LineWidth',qw);
set(gca,'FontSize',fs);
axis equal
xlim([-mr mr]);
ylim([-mr mr]);
title('Laplacian of velocity')

subplot(1,3,3)
q=quiver(xsgrid,ysgrid,vpred_x_grid,vpred_y_grid);
set(q,'LineWidth',qw);
set(gca,'FontSize',fs);
axis equal
xlim([-mr mr]);
ylim([-mr mr]);
title('Predicted by Faxen')
