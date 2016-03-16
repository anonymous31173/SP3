function area = voronoidens(k)
%
% function area = voronoidens(k);
%
% input:  k = kx + i ky is the  k-space trajectory
% output: area of cells for each point 
%           (if point doesn't have neighbors the area is NaN)

kx = real(k);
ky = imag(k);

[row,column] = size(kx);

kxy = [kx(:),ky(:)];
[V,C] = voronoin(kxy); 
area = [];
for j = 1:length(kxy)
  x = V(C{j},1); y = V(C{j},2); lxy = length(x);
  A = abs(sum( 0.5*(x([2:lxy 1]) - x(:)).*(y([2:lxy 1]) + y(:))));
  area = [area A];
end

area = reshape(area,row,column);

