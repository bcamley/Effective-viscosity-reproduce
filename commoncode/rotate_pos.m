function [xsrot,ysrot] = rotate_pos(xs,ys,theta)

xs = xs - mean(xs);
ys = ys - mean(ys); 

Rm = [ [cos(theta)  -sin(theta)] ; [sin(theta) cos(theta)]];

%np = length(xs);

xsrot = NaN*ones(size(xs));
ysrot = NaN*ones(size(ys));

for k = 1:length(xs)
   pnew = Rm*[xs(k) ys(k)]';
   xsrot(k) = pnew(1); ysrot(k) = pnew(2);
end

xsrot = xsrot;
ysrot = ysrot;

end