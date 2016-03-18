% Example code for image reconstruction of 3D non-Cartesian data using
% SPIRiT3

load sngshot.mat

% plot k-space
figure; plot(k); axis image;

% Calibration
tik = 0.02; % Tikhnov regularization
mth = '\';  % \ is the fastest inversion method
ksh = ones(5,5,3); % kernel shape
tic, kn = compkn(d,k,w,mdh,ksh,tik,mth); toc; % kernel calibration
save('kernel553.mat','kn','tik','-v7.3'); % save kernel

% Reconstruction
lmd = 2;    % Balance between kernel consistency and data consistency
nit = 50;   % Number of iterations
mt0 = 'image'; % reconstruct in image space
tic, [im, RV] = SOSASLacc3(d,k,w,mdh,kn,nit,lmd,mt0); toc;

% Save results
save('kernel553.mat','kn','tik','-v7.3'); % save kernel
save('results.mat','im','lmd','nit','RV','-v7.3'); % save reconstructed image




