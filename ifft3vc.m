function res = ifft3vc(x)
%
% res = ifft3vc(x)
%
% MATLAB function to perform 3D fast Fourier transform in the backward
% direction
%
% Input:
%   - x   : At least 3D; transform will be performed on the first 3 dimensions
%
% Output:
%   - res : Results after transformation
%
%   Y. Vince Chang, 20140314
%

assert(size(size(x),2) >= 3, 'image matrix cannot be less than 3D');

d1 = size(x,1); d2 = size(x,2); d3 = size(x,3);
rshx = reshape(x, d1, d2, d3, []);
res = zeros(size(rshx));
for ii = 1:size(rshx,4)
    res(:,:,:,ii) = sqrt(d1*d2*d3)*fftshift(ifftn(ifftshift(squeeze(rshx(:,:,:,ii)))));
end
res = reshape(res, size(x));

