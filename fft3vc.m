function res = fft3vc(x)

% res = fft3vc(x)
% 
% MATLAB function to perform 3D fast Fourier transform in the forward
% direction
%
% Input:
%   - x   : At least 3D; the first 3 dimensions must be the dimensions to be transformed
%
% Output:
%   - res : Resutls after transformation
%
%   Y. Vince Chang, 20140313
%


assert(size(size(x),2) >= 3, 'image matrix must be at least 3D');

d1 = size(x,1); d2 = size(x,2); d3 = size(x,3);
rshx = reshape(x, d1, d2, d3, []);
res = zeros(size(rshx));
for ii = 1:size(rshx,4)
    res(:,:,:,ii) = 1/sqrt(d1*d2*d3)*fftshift(fftn(ifftshift(squeeze(rshx(:,:,:,ii)))));
end
res = reshape(res,size(x));

