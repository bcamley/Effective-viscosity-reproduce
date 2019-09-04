function [xs,ys,fx,fy,vxgrid,vygrid] = get_flow_around_circle_specifymore(lsd,epsilonfactor,spacing,xsgrid,ysgrid,xs,ys,fx,fy)


eta_m = 1; % set membrane viscosity to 1 

epsilon = epsilonfactor*spacing;

np = length(xs);

xsall = [xs xsgrid(:).'];
ysall = [ys ysgrid(:).'];
npall = length(xsall);
Mijall = reg_stokeslet_matrix(xsall',ysall',eta_m,lsd,epsilon);
fall = zeros(1,2*npall);
fall(1:np) = fx;
fall((npall+1):(npall+np)) = fy;

vgrid = Mijall*(fall.');
vxgrid = vgrid((np+1):npall);
vygrid = vgrid((npall+1+np):end);
vxgrid = reshape(vxgrid,size(xsgrid));
vygrid = reshape(vygrid,size(ysgrid));
