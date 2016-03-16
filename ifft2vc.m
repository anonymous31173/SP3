function res = ifft2vc(x)

% res = ifft2vc(x)
% 

dim = size(x);
if size(dim,2) < 3
    res = sqrt(length(x(:)))*fftshift(ifft2(ifftshift(x)));
else
    dim1 = size(x,1); dim2 = size(x,2);
    rshx = reshape(x,dim1,dim2,[]);
    res = zeros(size(rshx));
    for ii = 1:size(rshx,3)
        res(:,:,ii) = sqrt(dim1*dim2)*fftshift(ifft2(ifftshift(rshx(:,:,ii))));
    end
    res = reshape(res,dim);
end