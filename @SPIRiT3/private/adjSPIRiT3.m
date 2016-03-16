function res = adjSPIRiT3(kData, kernel)
%
% Inputs:
%   kData  - k-space data in a single channel
%   kernel -
%

kernel = conj(kernel(end:-1:1,end:-1:1,end:-1:1,:));
res = zeros(size(kData));

Q = size(kernel,4);
for ii = 1:Q
    res(:,:,:,ii) = imfilter(kData, kernel(:,:,:,ii));
end
