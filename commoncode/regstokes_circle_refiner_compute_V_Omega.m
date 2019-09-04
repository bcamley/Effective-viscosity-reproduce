function [Fx_ext,Fy_ext,Lz_ext,Sij_ext,Fxs,Fys,Lzs,Sijs,xs,ys,ux,uy,fx,fy,Mij,Vxs,Vys,Omegas,Vx_ext,Vy_ext,Omega_ext] = regstokes_circle_refiner_compute_V_Omega(lsd,epsilonfactor,spacings,vext,cmloc)
% [Fx_ext,Fy_ext,Lz_ext,Sij_ext,Fxs,Fys,Lzs,Sijs,xs,ys,ux,uy,fx,fy,Mij,Vxs,Vys,Omegas,Vx_ext,Vy_ext,Omega_ext] = regstokes_circle_refiner_compute_V_Omega(lsd,epsilonfactor,spacings,vext,cmloc)

if(nargin<5)
    cmloc = [0 0];
end
Fxs = zeros(length(spacings),length(lsd));
Fys = zeros(length(spacings),length(lsd));
Lzs = zeros(length(spacings),length(lsd));

Vxs = zeros(length(spacings),length(lsd));
Vys = zeros(length(spacings),length(lsd));
Omegas = zeros(length(spacings),length(lsd));

Sijs = zeros(4,length(spacings),length(lsd));

Fx_ext = zeros(1,length(lsd));
Fy_ext = zeros(1,length(lsd));
Lz_ext = zeros(1,length(lsd));
Sij_ext = zeros(4,length(lsd));

Vx_ext = zeros(1,length(lsd));
Vy_ext = zeros(1,length(lsd));
Omega_ext = zeros(1,length(lsd));

eta_m = 1;
R = 1;
%if(nargin < 5)
%    squeeze = 1.42/1.73;
%end

for s = 1:length(spacings)
    spacing = spacings(s);
    epsilon = epsilonfactor*spacing;
    [xs,ys] = smoother_circle(spacing,R);
    xs = xs-mean(xs);  % watch out - this part is super sensitive if there is a zero rotation... 
    ys = ys-mean(ys);
    %xs = xs + cmloc(1);
    %ys = ys + cmloc(2);
    np = length(xs);
    ve = vext(xs.'+cmloc(1),ys.'+cmloc(2)); % we put in the center of mass location here because this allows us to still have CM = 0 to apply our definitions of torque that assume this
    vxext = ve(:,1);
    vyext = ve(:,2);
    
    for i = 1:length(lsd)
        %xs = cos(th(1:np)); ys = sin(th(1:np));
        Mij = reg_stokeslet_matrix(xs',ys',eta_m,lsd(i),epsilon);
        %Mijs{i} = Mij;
        % FIRST PART: compute only the force required to cancel off the
        % external velocity
        u = zeros(2*np,1);
        u(1:np) = -vxext;  
        u(np+1:end) = -vyext;        
%        u(1:np) = Ux-vxext-Omega*ys.';  
%        u(np+1:end) = Uy-vyext+Omega*xs.';
        f = gmres(Mij,u,[],[],500,[],[],randn(size(u)));
        
        fxcancel = f(1:np).';
        fycancel = f(np+1:end).';
        
        u = zeros(2*np,1);
        u(1:np) = 1;  
        u(np+1:end) = 0;        
%        u(1:np) = Ux-vxext-Omega*ys.';  
%        u(np+1:end) = Uy-vyext+Omega*xs.';
        f = gmres(Mij,u,[],[],500,[],[],randn(size(u)));
        
        fxtransx = f(1:np).';
        fytransx = f(np+1:end).';        
        
        u = zeros(2*np,1);
        u(1:np) = 0;  
        u(np+1:end) = 1;        
%        u(1:np) = Ux-vxext-Omega*ys.';  
%        u(np+1:end) = Uy-vyext+Omega*xs.';
        f = gmres(Mij,u,[],[],500,[],[],randn(size(u)));
        
        fxtransy = f(1:np).';
        fytransy = f(np+1:end).';

        % now only rotational motion
        u = zeros(2*np,1);
        u(1:np) = -ys.';  
        u(np+1:end) = xs.';        
        f = gmres(Mij,u,[],[],500,[],[],randn(size(u)));
        
        fxrot = f(1:np).';
        fyrot = f(np+1:end).';        
        
        % We need to have 3 requirements: sum of x forces = 0, sum of y
        % forces = 0, sum of torques = 0. This sets Vx, Vy, Vz
        
        Amatrix = [sum(fxtransx) sum(fxtransy) sum(fxrot) ; ...
                   sum(fytransx) sum(fytransy) sum(fyrot) ; ...
                   sum(xs.*fytransx-ys.*fxtransx) sum(xs.*fytransy-ys.*fxtransy) sum(xs.*fyrot-ys.*fxrot)];
        Bmatrix = [-sum(fxcancel) ; -sum(fycancel) ; -sum(xs.*fycancel-ys.*fxcancel)];
        
        VV = linsolve(Amatrix,Bmatrix)
        Vx = VV(1); Vy = VV(2); Omega = VV(3);
        fx = fxcancel+Vx*fxtransx+Vy*fxtransy +Omega*fxrot;
        fy = fycancel+Vx*fytransx+Vy*fytransy +Omega*fyrot;
        
        Vxs(s,i) = Vx;
        Vys(s,i) = Vy;
        Omegas(s,i) = Omega;
        
        Fxs(s,i) = sum(fx);
        Fys(s,i) = sum(fy);
        Lzs(s,i) = sum(xs.*fy -ys.*fx);
        Fxs(s,i)
        Fys(s,i)
        Lzs(s,i)
        Sij = 0.5*[sum(xs.*fx)*2 sum(xs.*fy+ys.*fx) sum(ys.*fx+xs.*fy) sum(ys.*fy)*2].';
        Sijs(:,s,i) = Sij;
        ux = u(1:np);
        uy = u(np+1:end);
    end
    
end

for i = 1:length(lsd)
    try
        p = polyfit(spacings,Fxs(:,i).',1);
        Fx_ext(i) = p(2);
        p = polyfit(spacings,Fys(:,i).',1);
        Fy_ext(i) = p(2);
        p = polyfit(spacings,Lzs(:,i).',1);
        Lz_ext(i) = p(2);
        p = polyfit(spacings,Vxs(:,i).',1);
        Vx_ext(i) = p(2);
        p = polyfit(spacings,Vys(:,i).',1);
        Vy_ext(i) = p(2);
        p = polyfit(spacings,Omegas(:,i).',1);
        Omega_ext(i) = p(2);
        for k = 1:4
                p = polyfit(spacings,squeeze(Sijs(k,:,i)),1);
                Sij_ext(k,i) = p(2);
        end
    catch err
        getReport(err)
    end
    
end
