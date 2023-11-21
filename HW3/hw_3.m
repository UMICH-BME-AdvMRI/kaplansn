function hw_3()
close all;
q_1()
q_2()
end

%% Question 1 code
function q_1()
load('Data_Assignment3_Problem1.mat')

%%% Question 1a
% partial fourier undersampling
pffact = 5 / 8;
origSize = size(kspaceData_SingleCoil, 1);
numLines = origSize * pffact;
pf5_8 = zeros(origSize);
pf5_8(1:numLines, 1:numLines) = kspaceData_SingleCoil(1:numLines, 1:numLines);

% display partial fourier image
PFimg = ifftshift(ifft2(pf5_8));
figure;
subplot(1, 2, 1)
imagesc(abs(PFimg)); axis square; title('PF 5/8 Magnitude'); colormap('gray')
subplot(1, 2, 2)
imagesc(angle(PFimg)); axis square; title('PF 5/8 Phase'); colormap('gray')

% display original image
FSimg = ifftshift(ifft2(kspaceData_SingleCoil));
figure;
subplot(1, 2, 1)
imagesc(abs(FSimg)); axis square; title('Fully Sampled Magnitude'); colormap('gray')
subplot(1, 2, 2)
imagesc(angle(FSimg)); axis square; title('Fully Sampled Phase'); colormap('gray')

% display difference image
diffimg = FSimg - PFimg;
figure;
subplot(1, 2, 1)
imagesc(abs(diffimg)); axis square; title('Difference Zero-Filled Magnitude'); colormap('gray'); caxis([0 0.0005]); colorbar
subplot(1, 2, 2)
imagesc(angle(diffimg)); axis square; title('Difference Zero-Filled Phase'); colormap('gray'); caxis([-pi/2 pi/2]); colorbar

%%% Question 1b
pf5_8recon = POCS(pf5_8, pffact);
PFrecon = ifftshift(ifft2(pf5_8recon));
figure;
subplot(1, 2, 1)
imagesc(abs(PFrecon)); axis square; title('POCS Recon Magnitude'); colormap('gray')
subplot(1, 2, 2)
imagesc(angle(PFrecon)); axis square; title('POCS Recon Phase'); colormap('gray')

% display difference image
diffimg = FSimg - PFrecon;
figure;
subplot(1, 2, 1)
imagesc(abs(diffimg)); axis square; title('Difference POCS Magnitude'); colormap('gray'); caxis([0 0.0005]); colorbar
subplot(1, 2, 2)
imagesc(angle(diffimg)); axis square; title('Difference POCS Phase'); colormap('gray'); caxis([-pi/2 pi/2]); colorbar
end

% function to compute POCS partial fourier recon
function estk = POCS(pfk, pffact)
tol = 1e-6; % tolerance
maxIter = 200; % maximum number of iterations

% get center lines of k-space
imgSize = size(pfk, 1);
numLines = pffact * imgSize;
centLine = floor(imgSize / 2);
kcentLines = (centLine - (numLines - centLine)):numLines;

% estimate phase
phaseEst = zeros(imgSize);
phaseEst(kcentLines, :) = pfk(kcentLines, :);

iter = 0;
kprev = pfk; % current estimated k-space
diffk = ones(imgSize);
while(max(abs(diffk(:))) > tol && iter < maxIter)
    % compute new kspace
    mag = abs(ifftshift(ifft2(kprev)));
    phase = angle(ifftshift(ifft2(phaseEst)));
    estk = fft2(fftshift(mag .* exp(1i * phase)));
    estk(1:numLines, 1:numLines) = pfk(1:numLines, 1:numLines);
    diffk = estk - kprev;

    % update iterators
    iter = iter + 1;
    kprev = estk;
end
end

%% Question 2
function q_2()
load('Data_Assignment3_Problem2.mat')

%%% Question 2a
% coil-combined image
[imgSize, ~, numCoils] = size(kspaceData);
I = zeros(imgSize, imgSize, numCoils);
for k = 1:numCoils
    I(:, :, k) = conj(coilmaps(:, :, k)) .* ifftshift(ifft2(kspaceData(:, :, k))); % adjusted based on Jiayao's code
end
I = sum(I, 3);

figure;
imagesc(abs(I)); axis square; title('Coil Combined Image'); colormap('gray');

%%% Question 2b
R2k = kspaceData;
R2k(1:2:end, :, :) = 0;
IR2 = zeros(imgSize, imgSize, numCoils);
for k = 1:numCoils
    IR2(:, :, k) = conj(coilmaps(:, :, k)) .* ifftshift(ifft2(R2k(:, :, k)));
end
IR2 = sum(IR2, 3);

figure;
imagesc(abs(IR2)); axis square; title('Coil Combined Image R=2'); colormap('gray');

%%% Question 2c
IR2_sense = SENSE(R2k, coilmaps, 2);

figure;
imagesc(abs(IR2_sense)); axis square; title('SENSE Recon R=2'); colormap('gray');

% display difference image
diffImg = abs(I) - abs(IR2_sense);
figure;
imagesc(diffImg); axis square; title('Difference SENSE R=2'); colormap('gray');

%%% Question 2d
R4k = kspaceData;
R4k(1:4:end, :, :) = 0;
R4k(2:4:end, :, :) = 0;
R4k(3:4:end, :, :) = 0;
IR4 = zeros(imgSize, imgSize, numCoils);
for k = 1:numCoils
    IR4(:, :, k) = conj(coilmaps(:, :, k)) .* ifftshift(ifft2(R4k(:, :, k)));
end
IR4 = sum(IR4, 3);

figure;
imagesc(abs(IR4)); axis square; title('Coil Combined Image R=4'); colormap('gray');

% fix mismatched coilmap indices
IR4_sense = SENSE(R4k, coilmaps, 4);

figure;
imagesc(abs(IR4_sense)); axis square; title('SENSE Recon R=4'); colormap('gray');

% display difference image
diffImg = abs(I) - abs(IR4_sense);
figure;
imagesc(diffImg); axis square; title('Difference SENSE R=4'); colormap('gray');

end

function rho = SENSE(ksp, coilmaps, R)
[imgSize, ~, numCoils] = size(ksp);
FOVR = floor(imgSize / R);
rho = zeros(imgSize);

% get all coil images
imgs = zeros(imgSize, imgSize, numCoils);
for k = 1:numCoils
    imgs(:, :, k) = ifftshift(ifft2(ksp(:, :, k)));
end

for y = 1:FOVR
    yshift = y:FOVR:imgSize;
    for x = 1:imgSize
        I = squeeze(imgs(y, x, :));
        C = transpose(squeeze(coilmaps(yshift, x, :)));
        rho(yshift, x) = C \ I;
    end
end
end
