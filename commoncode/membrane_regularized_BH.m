function [Dt,Xt] = membrane_regularized_BH(beta,c)
% function [Dt,Xt] = membrane_regularized_BH(beta,c)
% beta = R/lsd, c = epsilon/lsd
% where Sij = Dt(delta_ij) + Xt (xi xj / r^2)
%
% NOTE THESE FUNCTIONS DT,XT are NOT quite the same as the longitudinal and transverse response T_L, T_T defined in Camley, Brown 2013
% Instead, this code was built on breaking up the response function into
% T_{ij} = D (delta_{ij}) + X*(rhat_i rhat_j)
% This is entirely equivalent when you go through your bessel function identities, but a little confusing.
%
Dt = zeros(size(c));
Xt = zeros(size(c));

%c
if(length(c)==1)
    c = c * ones(size(beta));
end


if(length(beta)==1)
    beta = beta* ones(size(c));
end

length(beta)

f = @(x,beta,c) ((besselj(2,x)-besselj(0,x)).*exp(-0.5*(x.*c./beta).^2)./(x+beta))/(4*pi);
% the integral of f is B''

g = @(x,beta,c) -((besselj(1,x)./x).*exp(-0.5*(x.*c./beta).^2)./(x+beta))/(2*pi);
% the integral of g is B'/r

h = @(x,beta,c) -((besselj(0,x)).*exp(-0.5*(x.*c./beta).^2)./(x+beta))/(2*pi);
% the integral of h is H

perp_asympt = @(x) (1/(4*pi))*(pi*struve0(x)-(pi./x).*struve1(x) + 2./(x.^2) - (pi/2)*(bessely(0,x)-bessely(2,x)));
par_asympt = @(x) (1/(4*pi))*((pi./x).*struve1(x) - 2./(x.^2)-(pi/ ...
						  2).*(bessely(0,x)+bessely(2,x)));

lastwarn('');

scf = 1e5;  % small cutoff

for i = 1:length(c)
    
    if( (real(beta(i)) < 200*real(c(i))) && (real(beta(i)) > real(c(i))/scf))
      Bpr = quadgk(@(x) g(x,beta(i),c(i)),0,Inf);
      Bpp = quadgk(@(x) f(x,beta(i),c(i)),0,Inf);
      H = quadgk(@(x) h(x,beta(i),c(i)),0,Inf);
      Dt(i) = (Bpr-H);
      Xt(i) = (Bpp-Bpr);
      if(~isempty(lastwarn))
          beta(i)
          beta(i)/c(i)
          lastwarn('');
      end
    elseif(real(beta(i)) > real(c(i))/scf)   % much larger than cutoff, apply
                                  % asymptotic result
      Dt(i) = perp_asympt(beta(i));
      Xt(i) = par_asympt(beta(i))-perp_asympt(beta(i));
    else    % much smaller than cutoff, 
      disp('using small cutoff')
      beta(i)
      c(i)
      Bpr = quadgk(@(x) g(x,c(i)/scf,c(i)),0,Inf);
      H = quadgk(@(x) h(x,c(i)/scf,c(i)),0,Inf);
      Dt(i) = (Bpr-H);
      Xt(i) = 0;  % this is guaranteed
    end

   if(rem(i,100)==0)
     i
   end 
end

