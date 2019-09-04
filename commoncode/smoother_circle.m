function [xs,ys] = smoother_circle(spacing,R)
% function [xs,ys] = hex_packed_circle(epsilon,R)
rs = linspace(spacing,R,round(R/spacing));
xs = []; ys = []; 
lastspace = 0;
for i = 1:length(rs)
    th = lastspace/2+linspace(0,2*pi-spacing/rs(i),floor((2*pi*rs(i)-spacing)/spacing)); 
    xs = [xs rs(i)*cos(th)]; 
    ys = [ys rs(i)*sin(th)];
    lastspace = floor((2*pi*rs(i)-spacing)/spacing);

end
    xs = [xs 0];
    ys = [ys 0];

%xs = XS(:);

