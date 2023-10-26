function hw_2()
%% clear workspace
clear;
close all;
clc

% run question 1
q1()

% run question 2
q2()
end

%% Question 1
function q1()
% Question 1a
% (see functions below)
T1s = [200 500 800 1100 1500];
T2s = [50 110 170 230 290];
LW = 2;
numPlots = 5;
ETL = 64;

% Question 1ai
alpha = 180;
legstr = cell(1, numPlots);
figure;
for i = 1:numPlots
    T1 = T1s(i);
    T2 = T2s(i);
    M = epg_fse(alpha * pi / 180, T1, T2, ETL);
    plot(abs(M), 'LineWidth', LW)
    hold on
    xlabel('Echo')
    ylabel('abs(S)')
    title(['\alpha=' num2str(alpha) char(176)])
    legstr{i} = sprintf('T1=%dms, T2=%dms', T1, T2);
end
legend(legstr)

% Question 1aii
alpha = 120;
legstr = cell(1, numPlots);
figure;
for i = 1:numPlots
    T1 = T1s(i);
    T2 = T2s(i);
    M = epg_fse(alpha * pi / 180, T1, T2, ETL);
    plot(abs(M), 'LineWidth', LW)
    hold on
    xlabel('Echo')
    ylabel('abs(S)')
    title(['\alpha=' num2str(alpha) char(176)])
    legstr{i} = sprintf('T1=%dms, T2=%dms', T1, T2);
end
legend(legstr)

% Question 1aiii
alpha = 60;
legstr = cell(1, numPlots);
figure;
for i = 1:numPlots
    T1 = T1s(i);
    T2 = T2s(i);
    M = epg_fse(alpha * pi / 180, T1, T2, ETL);
    plot(abs(M), 'LineWidth', LW)
    hold on
    xlabel('Echo')
    ylabel('abs(S)')
    title(['\alpha=' num2str(alpha) char(176)])
    legstr{i} = sprintf('T1=%dms, T2=%dms', T1, T2);
end
legend(legstr)

% Question 1b
T1s = 200:100:1500;
T2s = 50:30:300;
alphas = [180 120 60];
echonum = [6 16 32 48];
for i = 1:length(alphas)
    alpha = alphas(i);
    cntplt = zeros(length(T1s), length(T2s), length(echonum));
    for j = 1:length(T1s)
        T1 = T1s(j);
        for k = 1:length(T2s)
            T2 = T2s(k);
            M = epg_fse(alpha * pi / 180, T1, T2, ETL);
            M = abs(M);
            cntplt(j, k, :) = M(echonum);
        end
    end
    for j = 1:length(echonum)
        figure;
        imagesc(abs(cntplt(:, :, j)));
        caxis([0 0.8]);
        xlabel('T2 (ms)')
        ylabel('T1 (ms)')
        set(gca, 'XTickLabels', T2s)
        set(gca, 'YTickLabels', T1s)
        c = colorbar;
        c.Label.String = 'abs(S)';
        set(gca, 'YDir', 'normal')
        title(['\alpha = ' num2str(alpha) char(176) ', Echo = ' num2str(echonum(j))])
    end
end
end

%% Functions for question 1
% code in this section built off of epg simulations in Rad229
% simulate FSE
function S = epg_fse(alpha, T1, T2, ETL)
% initial parameters
alpha_init = pi / 2; % 90 deg
dt = 5; % ms
P = zeros(3, 2 * ETL);
P(3, 1) = 1;
alpha1 = alpha + ((pi / 2) - (alpha / 2));
S = zeros(1, ETL);

% initial flip
P = rf(P, alpha_init, alpha_init);
for i = 1:ETL
    P = relax(P, T1, T2, dt);
    if i == 1
        P = rf(P, abs(alpha1), angle(alpha1));
    else
        P = rf(P, abs(alpha), angle(alpha));
    end
    P = relax(P, T1, T2, dt);

    S(i) = P(1, 1);
end
end

% simulate rf pulse
function P = rf(P0, alpha, phi)
R = [((cos(alpha / 2)) ^ 2) (exp(2 * 1i * phi) * (sin(alpha / 2)) ^ 2) (-1i * exp(1i * phi) * sin(alpha));
    (exp(-2 * 1i * phi) * (sin(alpha / 2)) ^ 2) ((cos(alpha / 2)) ^ 2) (1i * exp(-1i * phi) * sin(alpha));
    (-1i / 2 *exp(-1i * phi) * sin(alpha)) (1i / 2 * exp(1i * phi) * sin(alpha)) (cos(alpha))];

P = R * P0;
end

% simulate relaxation
function P = relax(P0, T1, T2, dt)
E2 = exp(-dt / T2);
E1 = exp(-dt / T1);

A = [E2 0 0;
    0 E2 0;
    0 0 E1];

B = 1 - E1;
P = A * P0;
P(3, 1) = P(3, 1) + B;

% include gradients
P(1, :) = circshift(P(1, :), [0 1]);
P(2, :) = circshift(P(2, :), [0 -1]);
P(2, end) = 0;
P(1, 1) = conj(P(2, 1));
end

%% Question 2
function q2()
% (see functions below)
% load in data
load('brain_maps.mat')
LW = 2;

% Question 2a
% T1-weighting
TR = 500; % ms
TE = 15; % ms
t1w = zeros(size(T1map));
for i = 1:numel(t1w)
    t1w(i) = spin_echo(T1map(i), T2map(i), TE, TR, 0);
end
t1w = M0map .* t1w;
figure;
imagesc(abs(t1w)); colormap('gray');
title('T1-weighted')

% PD-weighting
TR = 4000; % ms
TE = 15; % ms
pdw = zeros(size(T1map));
for i = 1:numel(pdw)
    pdw(i) = spin_echo(T1map(i), T2map(i), TE, TR, 0);
end
pdw = M0map .* pdw;
figure;
imagesc(abs(pdw)); colormap('gray');
title('PD-weighted')

% T2-weighting
TR = 6000; % ms
TE = 100; % ms
t2w = zeros(size(T1map));
for i = 1:numel(t1w)
    t2w(i) = spin_echo(T1map(i), T2map(i), TE, TR, 0);
end
t2w = M0map .* t2w;
figure;
imagesc(abs(t2w)); colormap('gray');
title('T2-weighted')

% Question 2bi
T1s = [1000 1000 2000 2000];
T2s = [50 100 50 100];
ETL = 32;
TR = 3000; % ms
TE = 5; % ms
nTRs = 5;
figure;
legstr = cell(1, length(T1s));
for i = 1:length(T1s)
    T1 = T1s(i);
    T2 = T2s(i);
    S = fast_spin_echo(T1, T2, TE, TR, 0, ETL, nTRs);
    plot(abs(S((end-ETL+1):end)), 'LineWidth', LW)
    xlabel('Echo')
    ylabel('abs(S)')
    legstr{i} = sprintf('T1=%dms, T2=%dms', T1, T2);
    hold on
end
legend(legstr)

% Question 2bii
TEeff = 80; % ms
ESP = 5; % ms
ETL = 32;
TEs = 0:ESP:(ESP*ETL);
TEs(1) = [];
kcentecho = find(TEs == TEeff);
kspacetot = zeros([size(T1map), length(TEs)]);
imgs = zeros([size(T1map), length(TEs)]);
for j = 1:numel(t1w)
    [row, col] = ind2sub(size(T1map), j);
    M = fast_spin_echo(T1map(j), T2map(j), ESP, TR, 0, ETL, 1);
    imgs(row, col, :) = M;
end
imgs(isnan(imgs)) = 0;
for i = 1:length(TEs)
    imgs(:, :, i) = M0map .* imgs(:, :, i);
    kspacetot(:, :, i) = fftshift(fft2(imgs(:, :, i)));
end

kspace = zeros(size(T1map));
slabsize = size(T1map, 1) / ETL;
slabstart = 1:slabsize:size(T1map, 1);
slabend = slabstart + slabsize - 1;
kspacelines = kcentecho;
echonum = kcentecho;
for i = 1:(ETL/2)
    min1 = kcentecho - i;
    if min1 < 1
        min1 = ETL + kcentecho - i;
    end
    plu1 = kcentecho + i;
    if plu1 > ETL
        plu1 = kcentecho + i - ETL;
    end
    kspacelines = [kspacelines min1 plu1];
    echonum = [echonum plu1 plu1];
end
kspacelines(kspacelines < 1 | kspacelines > ETL) = [];
for i = 1:length(TEs)
    curslab = kspacelines(i);
    kspace(slabstart(curslab):slabend(curslab), :) = kspacetot(slabstart(curslab):slabend(curslab), :, echonum(i));
end
imgeff = ifft2(kspace);
figure;
imagesc(abs(imgeff)); colormap('gray');
title(sprintf('TE_{eff} = %dms', TEeff))

figure;
imagesc(abs(imgeff)); colormap('gray');
title(sprintf('ETL = %d', ETL))

% Question 2biii
TEeff = 40; % ms
kcentecho = find(TEs == TEeff);
kspace = zeros(size(T1map));
echonum = kcentecho;
for i = 1:(ETL/2)
    plu1 = kcentecho + i;
    if plu1 > ETL
        plu1 = kcentecho + i - ETL;
    end
    echonum = [echonum plu1 plu1];
end
for i = 1:length(TEs)
    curslab = kspacelines(i);
    kspace(slabstart(curslab):slabend(curslab), :) = kspacetot(slabstart(curslab):slabend(curslab), :, echonum(i));
end
imgeff = ifft2(kspace);
figure;
imagesc(abs(imgeff)); colormap('gray');
title(sprintf('TE_{eff} = %dms', TEeff))

TEeff = 120; % ms
kcentecho = find(TEs == TEeff);
kspace = zeros(size(T1map));
echonum = kcentecho;
for i = 1:(ETL/2)
    plu1 = kcentecho + i;
    if plu1 > ETL
        plu1 = kcentecho + i - ETL;
    end
    echonum = [echonum plu1 plu1];
end
for i = 1:length(TEs)
    curslab = kspacelines(i);
    kspace(slabstart(curslab):slabend(curslab), :) = kspacetot(slabstart(curslab):slabend(curslab), :, echonum(i));
end
imgeff = ifft2(kspace);
figure;
imagesc(abs(imgeff)); colormap('gray');
title(sprintf('TE_{eff} = %dms', TEeff))

% Question 2biv
TEeff = 80; % ms
ESP = 5; % ms
ETLs = [16 64 128];
for k = 1:length(ETLs)
    ETL = ETLs(k);
    TEs = 0:ESP:(ESP*ETL);
    TEs(1) = [];
    kcentecho = find(TEs == TEeff);
    kspacetot = zeros([size(T1map), length(TEs)]);
    imgs = zeros([size(T1map), length(TEs)]);
    for j = 1:numel(t1w)
        [row, col] = ind2sub(size(T1map), j);
        M = fast_spin_echo(T1map(j), T2map(j), ESP, TR, 0, ETL, 1);
        imgs(row, col, :) = M;
    end
    imgs(isnan(imgs)) = 0;
    for i = 1:length(TEs)
        imgs(:, :, i) = M0map .* imgs(:, :, i);
        kspacetot(:, :, i) = fftshift(fft2(imgs(:, :, i)));
    end

    kspace = zeros(size(T1map));
    slabsize = size(T1map, 1) / ETL;
    slabstart = 1:slabsize:size(T1map, 1);
    slabend = slabstart + slabsize - 1;
    kspacelines = kcentecho;
    for i = 1:(ETL/2)
        min1 = kcentecho - i;
        if min1 < 1
            min1 = ETL + kcentecho - i;
        end
        plu1 = kcentecho + i;
        if plu1 > ETL
            plu1 = kcentecho + i - ETL;
        end
        kspacelines = [kspacelines min1 plu1];
    end
    kspacelines(kspacelines < 1 | kspacelines > ETL) = [];
    for i = 1:length(TEs)
        curslab = kspacelines(i);
        kspace(slabstart(curslab):slabend(curslab), :) = kspacetot(slabstart(curslab):slabend(curslab), :, curslab);
    end
    imgeff = ifft2(kspace);
    figure;
    imagesc(abs(imgeff)); colormap('gray');
    title(sprintf('ETL = %d', ETL))
end

% single echo image
imgse = zeros(size(T1map));
TE = 80; % ms
for i = 1:numel(imgse)
    imgse(i) = spin_echo(T1map(i), T2map(i), TE, TR, 0);
end
imgse = M0map .* imgse;
figure;
imagesc(abs(imgse)); colormap('gray');
title('Single Echo TE=80ms')
end

%% Functions for question 2
% code in this section built off of bloch simulations in Rad229 
function M = spin_echo(T1, T2, TE, TR, dfreq)
% rf pulses
Flip90 = rotate_y(pi / 2);
Refoc180 = rotate_x(pi);

% precession
[R_tr, b_tr] = freeprecess(TR - TE, dfreq, T1, T2);
[R_te_2, b_te_2] = freeprecess(TE / 2, dfreq, T1, T2);

% apply RF pulses
R_tr = [0 0 0; 0 0 0; 0 0 1] * R_tr;
Mss = inv(eye(3) - R_te_2 * Refoc180 * R_te_2 * Flip90 * R_tr) * (b_te_2 + R_te_2 * Refoc180 * (b_te_2 + R_te_2 * Flip90 * b_tr));
M = Mss(1) + 1i * Mss(2);
end

function Rx = rotate_x(alpha)
% function to rotate alpha (in radians) around x-axis
Rx = [1 0 0;
    0 cos(alpha) -sin(alpha);
    0 sin(alpha) cos(alpha)];
end

function Ry = rotate_y(alpha)
% function to rotate alpha (in radians) around y-axis
Ry = [cos(alpha) 0 sin(alpha);
    0 1 0;
    -sin(alpha) 0 cos(alpha)];
end

function Rz = rotate_z(phi)
% function to rotate phi (in radians) around z-axis
Rz = [cos(phi) -sin(phi) 0;
    sin(phi) cos(phi) 0;
    0 0 1];
end

function [A, b] = freeprecess(T,df,T1,T2)
% function to simulate precession
phi = 2*pi*df*(T * 10 ^ -3);
E1 = exp(-T/T1);
E2 = exp(-T/T2);

A = [E2 0 0;0 E2 0;0 0 E1]*rotate_z(phi);
b = [0 0 1-E1]';
end

function Mtot = fast_spin_echo(T1, T2, TE, TR, dfreq, ETL, nTRs)
Mtot = [];
% rf pulses
Flip90 = rotate_y(pi / 2);
Refoc180 = rotate_x(pi);

% precession
[R_tr, b_tr] = freeprecess(TR - ETL * TE, dfreq, T1, T2);
[R_te_2, b_te_2] = freeprecess(TE / 2, dfreq, T1, T2);

% apply RF pulses
R_tr = [0 0 0; 0 0 0; 0 0 1] * R_tr;
R = eye(3);
b = [0; 0; 0];

for k = 1:nTRs
    % echos
    for j = 1:ETL
        R = R_te_2 * Refoc180 * R_te_2 * R;
    	b = b_te_2 + R_te_2 * Refoc180 * (b_te_2 + R_te_2 * b);
    end

    % steady state sig
    R = Flip90 * R_tr * R;
    b = Flip90 * (b_tr + R_tr * b);
    Mss = inv(eye(3) - R) * b;
    M = Mss;

    % total signal
    for j = 1:ETL
        M = R_te_2 * Refoc180 * R_te_2 * M + b_te_2 + R_te_2 * Refoc180 * b_te_2;
    	Mss(:, j)=M;
    	Mtot(end+1) = M(1) + 1i * M(2);
    end
end
end