function ress = mtimes(a,bb)
% performs the normal nufft
% Modified by Yulin V. Chang, 20141023

sz = size(bb);
bb = reshape(bb,sz(1),sz(2),[]);
n = size(bb,3);

if a.adjoint
    ress = zeros([a.imSize n]);
    sz(1) = a.imSize(1); sz(2) = a.imSize(2);
    bb = reshape(bb, [], n).*repmat(a.w(:),[1 n]);
    for ii = 1:n
        ress(:,:,ii) = nufft_adj(bb(:,ii), a.st)/sqrt(prod(a.imSize));
    end
    ress = reshape(ress,sz);
else
    ress = zeros([prod(a.dataSize) n]);
    sz(1) = a.dataSize(1); sz(2) = a.dataSize(2);
    for ii = 1:n
        ress(:,ii) = nufft(bb(:,:,ii), a.st)/sqrt(prod(a.imSize)).*a.w(:);
    end
    ress = reshape(ress,sz);
end

