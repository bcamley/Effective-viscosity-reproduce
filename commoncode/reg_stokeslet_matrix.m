function [Mij,Dtref,Xtref,bs] = reg_stokeslet_matrix(xs,ys,eta_m,lsd,epsilon,n_interpolate)
% function Mij = reg_stokeslet_matrix(xs,ys,eta_m,lsd,epsilon,n_interpolate)

if(nargin < 6)
   n_interpolate = 1000; 
end

    N = length(xs);
    XIM = xs*ones(1,N);
    Rx = XIM - XIM'; 
    YIM = ys*ones(1,N);
    Ry = YIM - YIM';
    RIJ = sqrt(Rx.^2+Ry.^2);
    Rhatx = Rx./RIJ;
    Rhaty = Ry./RIJ;
Mij = zeros(2*N,2*N);


% we know that the non-isotropic part vanishes at R = 0, so this is convenient:
Rhatx(isnan(Rhatx)) = 0;
Rhaty(isnan(Rhaty)) = 0;

bs = linspace(min(0.1*RIJ(RIJ > eps)/lsd),max(RIJ(:)/lsd),n_interpolate);
min(bs)

[Dtref,Xtref] = membrane_regularized_BH(bs,epsilon/lsd);

[Dt,Xt] = membrane_regularized_BH_interp(RIJ(:)/lsd,epsilon/lsd,NaN,Dtref,Xtref,bs);
Dt = reshape(Dt',N,N);
Xt = reshape(Xt',N,N);
Dt = Dt/eta_m;
Xt = Xt/eta_m;  
% the xx components

alphaxx = Dt + (Rhatx.*Rhatx).*Xt;
Mij(1:N,1:N) = alphaxx;

alphaxy = Rhatx.*Rhaty.*Xt;
Mij((N+1):(2*N),1:N) = alphaxy;
Mij(1:N,(N+1):(2*N)) = alphaxy;

alphayy = Dt +Rhaty.*Rhaty.*Xt;
Mij((N+1):(2*N),(N+1):(2*N)) = alphayy;



%for i = 1:N
%    for j = 1:N
%        Mij(i,j) = 
%    end
%end


%for i = 1:N
%    for j = 1:N
%        RIJ
%        ri = 
%        rhat_i = 
%        Mij(i,j) = 
%    end
%end
