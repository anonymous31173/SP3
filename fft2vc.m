function res = fft2vc(x)

% res = fft2vc(x)
% MATLAB routine for orthonormal forward 2D FFT with single- or multi-slice
% input (x)

dim = size(x);
if size(dim,2) < 3  % x is 2D
    res = 1/sqrt(length(x(:)))*fftshift(fft2(ifftshift(x)));
else
    dim1 = size(x,1); dim2 = size(x,2);
    rshx = reshape(x,dim1,dim2,[]);
    res = zeros(size(rshx));
    for ii = 1:size(rshx,3)
        res(:,:,ii) = 1/sqrt(dim1*dim2)*fftshift(fft2(ifftshift(rshx(:,:,ii))));
    end
    res = reshape(res,dim);
end