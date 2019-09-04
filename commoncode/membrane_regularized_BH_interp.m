function [Dt,Xt,Dtref,Xtref,bs] = membrane_regularized_BH_interp(beta,c,interp_points,Dtref,Xtref,bs)
% function [Dt,Xt] = membrane_regularized_BH(beta,c)
% where Sij = Dt(delta_ij) + Xt (xi xj / r^2)

if(nargin < 3)
    interp_points = 400;
    Dtref = 0;
    Xtref = 0;
end

if(nargin < 4)
    Dtref = 0; Xtref = 0;
end

Dt = zeros(size(beta));
Xt = zeros(size(beta));


%if(length(c)==1)
%    c = c * ones(size(beta));
%end

if(length(beta)==1)
    beta = beta* ones(size(c));
end


bsmax = max(abs(beta));
bsmin = max(min(abs(beta(beta>eps)))/10,c/1000);

if(Dtref==0)
    disp('using logarithmically spaced betas ...')
    bs = logspace( log(bsmin)/log(10),log(bsmax)/log(10),interp_points);
    [Dtref,Xtref] = membrane_regularized_BH(bs,c);
end

%f = @(x,beta,c) ((besselj(2,x)-besselj(0,x)).*exp(-x*c./beta).*(1+ ...
%						  x*c./beta)./(x+beta))/(4*pi);
% the integral of f is B''

g = @(x,beta,c) -((besselj(1,x)./x).*exp(-x*c./beta).*(1+ x*c./beta)./(x+beta))/(2*pi);
% the integral of g is B'/r

h = @(x,beta,c) -((besselj(0,x)).*exp(-x*c./beta).*(1+ x*c./beta)./(x+beta))/(2*pi);
% the integral of h is H

%perp_asympt = @(x) (1/(4*pi))*(pi*struve0(x)-(pi./x).*struve1(x) + 2./(x.^2) - (pi/2)*(bessely(0,x)-bessely(2,x)));
%par_asympt = @(x) (1/(4*pi))*((pi./x).*struve1(x) - 2./(x.^2)-(pi/ ...%
%						  2).*(bessely(0,x)+bessely(2,x)));


if(imag(beta) < eps)
% old default was 'spline'    
%Dt(beta > bsmin) = interp1(bs,Dtref,beta(beta>bsmin),'spline');
%Xt(beta > bsmin) = interp1(bs,Xtref,beta(beta>bsmin),'spline');
Dt(beta > bsmin) = interp1(bs,Dtref,beta(beta>bsmin),'linear');
Xt(beta > bsmin) = interp1(bs,Xtref,beta(beta>bsmin),'linear');

else   % do complex beta
    betagood = abs(beta)>bsmin;
    Dt(betagood) = interp2(real(bs),imag(bs),Dtref,real(beta(betagood)),imag(beta(betagood)));
    Xt(betagood) = interp2(real(bs),imag(bs),Xtref,real(beta(betagood)),imag(beta(betagood)));
end


if(sum( (abs(beta)>bsmax)) > 0)
    disp('have beta too large, interpolation will be bad');
end

%for i = 1:length(c)

%    if( (beta(i) < 200*c(i)) && (beta(i) > c(i)/1000))
%      Bpr = quadgk(@(x) g(x,beta(i),c(i)),0,Inf);
%      Bpp = quadgk(@(x) f(x,beta(i),c(i)),0,Inf);
%      H = quadgk(@(x) h(x,beta(i),c(i)),0,Inf);
%      Dt(i) = (Bpr-H);
%      Xt(i) = (Bpp-Bpr);
%    elseif(beta(i) > c(i)/1000)   % much larger than cutoff, apply
%                                  % asymptotic result
%      Dt(i) = perp_asympt(beta(i));
%      Xt(i) = par_asympt(beta(i))-perp_asympt(beta(i));
%    else    % much smaller than cutoff,

      Bpr = quadgk(@(x) g(x,c/1000,c),0,Inf);
      H = quadgk(@(x) h(x,c/1000,c),0,Inf);
      Dt(beta <= bsmin) = Bpr-H;
      Xt(beta <= bsmin) = 0;
      
      %Dt(i) = (Bpr-H);
      %Xt(i) = 0;  % this is guaranteed
 %   end

  % if(rem(i,10)==0)
  %   i
  % end 
%end

