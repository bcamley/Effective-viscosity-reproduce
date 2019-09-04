function [xsp,ysp] = smooth_pair(spacing,R,addangle,nfactor)
[xs,ys] = smoother_circle(spacing,R);
% function [xs,ys] = hex_packed_circle(epsilon,R)
xsp = [];
ysp = [];
for j =1:nfactor
    xsp = [xsp xs+(2*R+spacing)*cos(addangle)*(j-1)]; % yes this is slow, but we don't run this very often
    ysp = [ysp ys+(2*R+spacing)*sin(addangle)*(j-1)];
end
xsp =xsp-mean(xsp);
ysp = ysp-mean(ysp);
    %ysp = ysp -mean(ysp);
%xs = XS(:);

