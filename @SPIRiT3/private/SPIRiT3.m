function res = SPIRiT3(kData, kernel)
%
%
%

res = zeros(size(kData));
Q = size(kData,4);

for ii = 1:Q
    res(:,:,:,ii) = imfilter(kData(:,:,:,ii),kernel(:,:,:,ii));
end

res = sum(res,4);
