function res = mtimes(sp3,x)
%
% Implementation of SPIRiT3 operations
%

kn3 = sp3.kernel;
Q   = size(kn3,4);
ksz = [size(kn3,1),size(kn3,2),size(kn3,3)];

switch sp3.mth
    case 'conv'
        res = zeros(size(x));
        if sp3.adjoint
            for ii = 1:Q
                tmpk = kn3(:,:,:,:,ii);
                tmpk(floor(ksz(1)/2)+1,floor(ksz(2)/2)+1,floor(ksz(3)/2)+1,ii) = tmpk(floor(ksz(1)/2)+1,floor(ksz(2)/2)+1,floor(ksz(3)/2)+1,ii)-1;
                res = res + adjSPIRiT3(x(:,:,:,ii),tmpk);
            end
        else
            for ii = 1:Q
                tmpk = kn3(:,:,:,:,ii);
                tmpk(floor(ksz(1)/2)+1,floor(ksz(2)/2)+1,floor(ksz(3)/2)+1,ii) = tmpk(floor(ksz(1)/2)+1,floor(ksz(2)/2)+1,floor(ksz(3)/2)+1,ii)-1;
                res(:,:,:,ii) = SPIRiT3(x, tmpk);
            end
        end
    case 'fft'
        res = zeros(size(x));
        xx  = ifft3vc(x);
        if sp3.adjoint
            for ii = 1:Q
                tmpk = squeeze(conj(sp3.K3(:,:,:,ii,:))); % transpose
                res(:,:,:,ii) = sum(tmpk.*xx, 4);
            end
        else
            for ii = 1:Q
                tmpk = sp3.K3(:,:,:,:,ii);
                res(:,:,:,ii) = sum(tmpk.*xx, 4);
            end
        end
        res = fft3vc(res) - x;
    case 'image'
        res = zeros(size(x));
        if sp3.adjoint
            for ii = 1:Q
                tmpk = squeeze(conj(sp3.K3(:,:,:,ii,:)));
                res(:,:,:,ii) = sum(tmpk.*x,4);
            end
        else
            for ii = 1:Q
                tmpk = sp3.K3(:,:,:,:,ii);
                res(:,:,:,ii) = sum(tmpk.*x,4);
            end
        end
        res = res - x;
end

