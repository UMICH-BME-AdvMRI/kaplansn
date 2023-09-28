function hw_1()
%% clear workspace
clear;
close all;
clc

%% plotting properties
LW = 2;

%% Question 1
rBW = 750; % Hz/px
N = 256; % Nx = Ny
dx = 1.2; % mm
dy = 1.2; % mm
dz = 3; % mm
Gmax = 25; % mT/m
s = 180; % T / m / s
gammabar = 42.58; % MHz / T
FOV = dx * N; % mm

% Question 1a
fprintf('<strong>Question 1a</strong>\n')
TBW = 2;
tau_RF = 1; % ms

%%% compute slice select amplitude
Gss = (TBW / (tau_RF * 10 ^ -3)) / ((dz * 10 ^ -3) * (gammabar * 10 ^ 6)); % T / m
fprintf('The slice-select gradient amplitude is %.2fmT/m\n', Gss * 10 ^ 3)

%%% compute slice select duration
tau_ss_rise = Gss / s;
%%tau_ss = 1 / ((dz * 10 ^ -3) * (gammabar * 10 ^ 6) * Gss);
tau_ss = tau_RF * 10 ^ -3;
tau_ss_tot = 2 * tau_ss_rise + tau_ss;
fprintf('The total duration of the slice-select gradient is %.2fms\n', tau_ss_tot * 10 ^ 3)

% Question 1b
fprintf('\n<strong>Question 1b</strong>\n')
tau_pe = 1 / (2 * (dy * 10 ^ -3) * (gammabar * 10 ^ 6) * (Gmax * 10 ^ -3));
tau_pe_rise = (Gmax * 10 ^ -3) / s;
tau_pe_tot = 2 * tau_pe_rise + tau_pe;
fprintf('The minimum duration of the phase encoding gradient is %.2fms\n', tau_pe_tot * 10 ^ 3)

% Question 1c
fprintf('\n<strong>Question 1c</strong>\n')
%%% compute readout amplitude
Gro = (rBW * N) / ((gammabar * 10 ^ 6) * (FOV * 10 ^ -3));
fprintf('The frequency encoding gradient amplitude is %.2fmT/m\n', Gro * 10 ^ 3)

%%% compute readout duration
tau_ro = 1 / rBW;
tau_ro_rise = Gro / s;
tau_ro_tot = 2 * tau_ro_rise + tau_ro;
fprintf('The total duration of the frequency encoding gradient is %.2fms\n', tau_ro_tot * 10 ^ 3)

% Question 1d
fprintf('\n<strong>Question 1d</strong>\n')
tau_overlap = max([tau_pe_tot (tau_ss_tot / 2) (tau_ro_tot / 2)]);
TEs = (tau_ss_tot / 2) + tau_overlap + (tau_ro_tot / 2);
fprintf('The minimum TE is %.2fms\n', TEs * 10 ^ 3)
TRs = (tau_ss_tot / 2) + tau_overlap + tau_ro_tot + tau_overlap + (tau_ss_tot / 2);
fprintf('The minimum TR is %.2fms\n', TRs * 10 ^ 3)

% Question 1e
fprintf('\n<strong>Question 1e</strong>\n')
fprintf('1.\tIncrease receiver bandwidth\n')
fprintf('2.\tDecrease time bandwidth\n')
fprintf('3.\tDecrease RF pulse duration\n')
%%fprintf('4.\tIncrease slew rate\n')

%% Question 2
% code in this section built off of that located here:
% http://mrsrl.stanford.edu/~brian/bloch/
% Question 2a
TRs = [5 10 20]; % ms
TEs = [2.5 5 10]; % ms
FA = 60; % deg
freq = -100:100; % Hz
T1 = 1000; % ms
T2 = 100; % ms

% repeat for each TR/TE pair
for k = 1:length(TRs)
    TR = TRs(k);
    TE = TEs(k);
    S = zeros(1, length(freq)); % initalize signal

    % compute response
    for j = 1:length(freq)
        S(j) = bSSFP(FA, TE, TR, freq(j), 1, T1, T2);
    end

    % plot freq response
    figure;
    plot(freq, abs(S), 'k-', 'LineWidth', LW)
    title({'Frequency Response', sprintf('TR=%.2fms, TE=%.2fms', TR, TE)})
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    if k == 1
        ylimvals = get(gca, 'YLim');
    else
        ylim(ylimvals)
    end
end

% Question 2bi
fprintf('\n<strong>Question 2bi</strong>\n')
TR = 10; % ms
TE = 5; % ms
FA = 10; % deg
T1 = 1000; % ms
T2 = 100; % ms
perSpoil = 1; % assume perfect spoiler

% compute the steady-state signal across all isochromats
M_ss = FLASH(FA, TE, TR, 0, T1, T2, perSpoil);
fprintf('The steady-state magnetization for a perfect-spoiler is:\nM=[%.2f; %.2f; %.2f]\n', M_ss(1), M_ss(2), M_ss(3))
fprintf('The steady-state Mxy=%.2f\n', abs(M_ss(1) + 1i * M_ss(2)))

% Question 2bii
Nspins = 200; % number of spins to simulate
perSpoil = 0; % assume nonideal spoiler
numpi = (2 .^ (1:4));
pTwists = numpi * pi;
M_ss_tot = zeros(3, length(pTwists));

for j = 1:length(pTwists)
    M_ss = zeros(3, Nspins); % steady-state signal
    dfreq = linspace(0, pTwists(j), Nspins); % off-resonance for each spin
    for k = 1:Nspins
        % simulate with perfect spoiler and plot
        M_tot = FLASH(FA, TE, TR, dfreq(k), T1, T2, perSpoil, pTwists(j));
        M_ss(:, k) = M_tot;
    end

    % compute the steady-state signal across all isochromats
    M_ss = mean(M_ss');

    % add to total
    M_ss_tot(:, j) = M_ss;
end

figure;
plot(pTwists, abs(M_ss_tot(1, :) + 1i * M_ss_tot(2, :)), 'k-', 'LineWidth', LW)
xticks(pTwists);
xticklabels(strcat(strsplit(num2str(numpi)), '\pi'))
xlabel('Dephasing Moment')
ylabel('Steady-state Mxy')

% Question 2biii
Nex = 50; % number of excitations
inc = 2 * pi / (Nex * (Nex + 1)); % rf phase increment
perSpoil = 0;
pTwist = 2*pi;
dfreq = linspace(0, pTwist, Nspins); % off-resonance for each spin

S_tot = zeros(3, Nex+1, Nspins);
for j = 1:Nspins
    [~, S] = rf_plus_grad(FA, TE, TR, dfreq(j), T1, T2, Nex, inc, perSpoil, pTwist);
    S_tot(:, :, j) = S;
end
S_tot = mean(S_tot, 3);

% plot
figure;
rfphase = (0:Nex) .* (1:(Nex+1)) / 2 * inc;
plot(rfphase, abs(S_tot(1, :) + 1i * abs(S_tot(2, :))), 'k-', 'LineWidth', LW)
xticks(rfphase(1:5:end))
xlim([0 rfphase(end)])
xlabel('RF-phase (rad)')
ylabel('Steady-state Mxy')

% best rfphase
[~, idx] = max(abs(S_tot(1, :) + 1i * abs(S_tot(2, :))));
fprintf('The best RF phase is=%.2f\n', rfphase(idx))

%% Question 3
s = 180; % mT / (m * ms)
Gmax = 25; % mT / m
TBW = 8;
tau_RF = 2; % ms
dz = 5; % mm
gammabar = 42.58; % MHz / T

% Question 3.1
fprintf('\n<strong>Question 3.1</strong>\n')
%%% compute slice select gradient amplitude
Gss = (TBW / (tau_RF * 10 ^ -3)) / ((dz * 10 ^ -3) * (gammabar * 10 ^ 6)); % T / m
fprintf('The slice-select gradient amplitude is %.2fmT/m\n', Gss * 10 ^ 3)

%%% compute slice select duration
tau_ss_rise = Gss / s;
tau_ss = tau_RF * 10 ^ -3;
tau_ss_tot = 2 * tau_ss_rise + tau_ss;
fprintf('The total duration of the slice-select gradient is %.2fms\n', tau_ss_tot * 10 ^ 3)
fprintf('The RF pulse follows sinc(4t) but is truncated at 1ms\n')

% Question 3.2
FA = 90; % deg
dfreq = [0 200]; % off resonance
T1 = 1000; % ms
T2 = 100; % ms
zpos = -10:10;
BW = TBW / tau_RF;
t = -1:0.01:1; % rf pulse duration (ms)
t_shift = (0:0.01:2) .* (10 ^ -3); % shift to start from 0 and in (s)
rf = sinc(BW .* t);
rf_int = cumsum(t_shift .* rf);
rf_amp90 = FA / (gammabar * rf_int(end));
dt = t(2) - t(1);
position = -10:10; % position (mm)
grad = Gss * ones(1, length(t)) * 10 ^ 3; % mT/m
grad = [grad -grad(1:2:end)];
rf90 = zeros(1, length(grad));
rf90(1:length(t)) = rf_amp90 * sinc(BW .* t);

% simulation
for k = 1:length(dfreq)
    M = sliceprofile(T1, T2, rf90 / 10, dt, grad, dfreq(k), position);
    
    % plot slice profile
    figure;
    tt = tiledlayout(4, 1);
    nexttile;
    plot(position, M(1, :), '-k', 'LineWidth', LW)
    ylabel('Mx')
    xlabel('Position (mm)')
    nexttile;
    plot(position, M(2, :), '-k', 'LineWidth', LW)
    ylabel('My')
    xlabel('Position (mm)')
    nexttile;
    plot(position, M(3, :), '-k', 'LineWidth', LW)
    ylabel('Mz')
    xlabel('Position (mm)')
    if k == 1
        M90 = M;
    end
    nexttile;
    plot(position, abs(M(1, :) + 1i * M(2, :)), '-k', 'LineWidth', LW)
    ylabel('Mxy')
    xlabel('Position (mm)')
    title(tt, {['Off-Resonance = ' num2str(dfreq(k)) 'Hz'], ['FA = ' num2str(FA)]})
end

% plot pulse and gradient
figure;
tiledlayout(2, 1);

% rf pulse
nexttile;
plot((t_shift + tau_ss_rise) * 10 ^ 3, rf_amp90 * sinc(BW .* t), 'k-', 'LineWidth', LW)
xlim([0 ((tau_ss_tot + 2 * tau_ss_rise + tau_ss / 2) * 10 ^ 3)])
xlabel('Time (ms)')
ylabel('RF pulse (\muT)')

% gradient
nexttile;
plot([0 (tau_ss_rise * 10 ^ 3)], [0 Gss], 'k-', 'LineWidth', LW)
hold on
plot(([tau_ss_rise (tau_ss_rise + tau_ss)]) * 10 ^ 3, [Gss Gss], 'k-', 'LineWidth', LW)
plot(([(tau_ss_rise + tau_ss) tau_ss_tot]) * 10 ^ 3, [Gss 0], 'k-', 'LineWidth', LW)
plot(([tau_ss_tot (tau_ss_tot + tau_ss_rise)]) * 10 ^ 3, [0 -Gss], 'k-', 'LineWidth', LW)
plot(([(tau_ss_tot + tau_ss_rise) (tau_ss_tot + tau_ss_rise + tau_ss / 2)]) * 10 ^ 3, [-Gss -Gss], 'k-', 'LineWidth', LW)
plot(([(tau_ss_tot + tau_ss_rise + tau_ss / 2) (tau_ss_tot + 2 * tau_ss_rise + tau_ss / 2)]) * 10 ^ 3, [-Gss 0], 'k-', 'LineWidth', LW)
xlim([0 ((tau_ss_tot + 2 * tau_ss_rise + tau_ss / 2) * 10 ^ 3)])
xlabel('Time (ms)')
ylabel('Gss (mT/m)')

% Question 3.3
FA = 30; % deg
rf_amp30 = FA / (gammabar * rf_int(end));
rf30 = zeros(1, length(grad));
rf30(1:length(t)) = rf_amp30 * sinc(BW .* t);
% simulation
for k = 1:length(dfreq)
    M = sliceprofile(T1, T2, rf30 / 10, dt, grad, dfreq(k), position);
    
    % plot slice profile
    figure;
    tt = tiledlayout(4, 1);
    nexttile;
    plot(position, M(1, :), '-k', 'LineWidth', LW)
    ylabel('Mx')
    xlabel('Position (mm)')
    nexttile;
    plot(position, M(2, :), '-k', 'LineWidth', LW)
    ylabel('My')
    xlabel('Position (mm)')
    nexttile;
    plot(position, M(3, :), '-k', 'LineWidth', LW)
    ylabel('Mz')
    xlabel('Position (mm)')
    nexttile;
    plot(position, abs(M(1, :) + 1i * M(2, :)), '-k', 'LineWidth', LW)
    ylabel('Mxy')
    xlabel('Position (mm)')
    title(tt, {['Off-Resonance = ' num2str(dfreq(k)) 'Hz'], ['FA = ' num2str(FA)]})
end

FA = 10; % deg
rf_amp10 = FA / (gammabar * rf_int(end));
rf10 = zeros(1, length(grad));
rf10(1:length(t)) = rf_amp10 * sinc(BW .* t);
% simulation
for k = 1:length(dfreq)
    M = sliceprofile(T1, T2, rf10 / 10, dt, grad, dfreq(k), position);
    
    % plot slice profile
    figure;
    tt = tiledlayout(4, 1);
    nexttile;
    plot(position, M(1, :), '-k', 'LineWidth', LW)
    ylabel('Mx')
    xlabel('Position (mm)')
    nexttile;
    plot(position, M(2, :), '-k', 'LineWidth', LW)
    ylabel('My')
    xlabel('Position (mm)')
    nexttile;
    plot(position, M(3, :), '-k', 'LineWidth', LW)
    ylabel('Mz')
    xlabel('Position (mm)')
    if k == 1
        M10 = M;
    end
    nexttile;
    plot(position, abs(M(1, :) + 1i * M(2, :)), '-k', 'LineWidth', LW)
    ylabel('Mxy')
    xlabel('Position (mm)')
    title(tt, {['Off-Resonance = ' num2str(dfreq(k)) 'Hz'], ['FA = ' num2str(FA)]})
end

% overlay w/ FT
figure;
tt = tiledlayout(1, 2);
% 90
nexttile;
normM90 = -1 * M90(3, :);
normM90 = normM90 - min(normM90);
normM90 = normM90 / max(normM90);
normft = ft(rf90(1:length(t)));
normft = normft / max(normft);
plot(linspace(-10, 10, length(t)), abs(normft), 'r-', 'LineWidth', LW)
hold on
plot(position, normM90, 'b-', 'LineWidth', LW)
xlabel('Position (mm)')
ylabel('Normalized Slice Profile')
title('FA=90')
% 10
nexttile;
normM10 = -1 * M10(3, :);
normM10 = normM10 - min(normM10);
normM10 = normM10 / max(normM10);
normft = ft(rf10(1:length(t)));
normft = normft / max(normft);
plot(linspace(-10, 10, length(t)), abs(normft), 'r-', 'LineWidth', LW)
hold on
plot(position, normM10, 'b-', 'LineWidth', LW)
xlabel('Position (mm)')
ylabel('Normalized Slice Profile')
title('FA=10')
title(tt, 'T2 = 100ms')

% Set T2 = 2
T2 = 2; % ms
M90 = sliceprofile(T1, T2, rf90 / 10, dt, grad, 0, position);
M10 = sliceprofile(T1, T2, rf10 / 10, dt, grad, 0, position);
% overlay w/ FT
figure;
tt = tiledlayout(1, 2);
% 90
nexttile;
normM90 = -1 * M90(3, :);
normM90 = normM90 - min(normM90);
normM90 = normM90 / max(normM90);
normft = ft(rf90(1:length(t)));
normft = normft / max(normft);
plot(linspace(-10, 10, length(t)), abs(normft), 'r-', 'LineWidth', LW)
hold on
plot(position, normM90, 'b-', 'LineWidth', LW)
xlabel('Position (mm)')
ylabel('Normalized Slice Profile')
title('FA=90')
% 10
nexttile;
normM10 = -1 * M10(3, :);
normM10 = normM10 - min(normM10);
normM10 = normM10 / max(normM10);
normft = ft(rf10(1:length(t)));
normft = normft / max(normft);
plot(linspace(-10, 10, length(t)), abs(normft), 'r-', 'LineWidth', LW)
hold on
plot(position, normM10, 'b-', 'LineWidth', LW)
xlabel('Position (mm)')
ylabel('Normalized Slice Profile')
title('FA=10')
title(tt, 'T2 = 2ms')

% Question 3.4
FA = 90;
M = sliceprofile(T1, T2, rf90(1:length(t)) / 10, dt, grad(1:length(t)), 0, position);

% plot slice profile
figure;
tt = tiledlayout(4, 1);
nexttile;
plot(position, M(1, :), '-k', 'LineWidth', LW)
ylabel('Mx')
xlabel('Position (mm)')
nexttile;
plot(position, M(2, :), '-k', 'LineWidth', LW)
ylabel('My')
xlabel('Position (mm)')
nexttile;
plot(position, M(3, :), '-k', 'LineWidth', LW)
ylabel('Mz')
xlabel('Position (mm)')
nexttile;
plot(position, abs(M(1, :) + 1i * M(2, :)), '-k', 'LineWidth', LW)
ylabel('Mxy')
xlabel('Position (mm)')
title(tt, {['Off-Resonance = ' num2str(0) 'Hz'], ['FA = ' num2str(FA)]})

% Question 3.5
T2 = 100; % ms
numSlices = 5;
phaseshift = 0;
position = -100:100;
slicepos = linspace(-50, 50, numSlices) / 1000;
for k = 1:numSlices
    phaseshift = phaseshift + exp(1i * (gammabar * 10 ^ 6 * Gss * slicepos(k) * t));
end
rf_sms = zeros(1, length(grad));
rf_sms(1:length(t)) = rf_amp90 * sinc(BW * t) .* phaseshift;
M = sliceprofile(T1, T2, rf_sms / 10, dt, grad, 0, position);
figure;
tt = tiledlayout(4, 1);
nexttile;
plot(position, M(1, :), '-k', 'LineWidth', LW)
ylabel('Mx')
xlabel('Position (mm)')
nexttile;
plot(position, M(2, :), '-k', 'LineWidth', LW)
ylabel('My')
xlabel('Position (mm)')
nexttile;
plot(position, M(3, :), '-k', 'LineWidth', LW)
ylabel('Mz')
xlabel('Position (mm)')
nexttile;
plot(position, abs(M(1, :) + 1i * M(2, :)), '-k', 'LineWidth', LW)
ylabel('Mxy')
xlabel('Position (mm)')
end

function M_tot = sliceprofile(T1, T2, rf, dt, grad, dfreq, position)
% function to simulate slice profile of given RF pulse
%
% convert units
gammabar = 42.58 * 10 ^ 2; % unit conversion
dt = dt * 10 ^ -3; % sec
position = position / 10; % unit conversion
rf_radrot = 2*pi*rf .* gammabar * dt;
M_tot = zeros(3, length(position));
for k = 1:length(position)
    M = [0; 0; 1];
    [R_rf, b_rf] = freeprecess((dt / 2 * 10 ^ 3), dfreq, 1, T1, T2);
    for j = 1:length(rf)
        % initial precession
        M = R_rf * M + b_rf;

        % add in position
        R_grad = rotate_z(2 * pi * gammabar * position(k) * grad(j) * (dt / 2));
        M = R_grad * M;

        % rf rotation
        M = (inv(rotate_z(-angle(rf_radrot(j)))) * rotate_x_rad(abs(rf_radrot(j))) * rotate_z(-angle(rf_radrot(j)))) * M;
        M = R_rf * M + b_rf;

        % add in position
        M = R_grad * M;
    end
    M_tot(:, k) = M;
end

end

function [M, M_tot] = rf_plus_grad(FA, TE, TR, dfreq, T1, T2, Nex, inc, perSpoil, pTwist)
% function to simulate sequence with both RF and gradient spoiling
% perSpoil = 1/0 assume perfect spoiling
% pTwist (optional iff perSpoil = 0)
relax = 1; % only allowing option that includes relaxation
% initial magnetization
M = [0; 0; 1];

% magnetization at each excitation
M_tot = zeros(3, Nex + 1);
M_tot(:, 1) = M;

% inital phase
Rf_ph = 0;

% precession after TE
[R_te, b_te] = freeprecess(TE, dfreq, relax, T1, T2);

% precession at start of next TR
[R_tr, b_tr] = freeprecess(TR-TE, dfreq, relax, T1, T2);

% add gradient spoiling
if perSpoil
    % remove transverse magnetization
    R_tr = ([0 0 0; 0 0 0; 0 0 1]) * R_tr;
else
    R_tr = rotate_z(pTwist) * R_tr;
end

for k = 1:Nex
    % rf-spoiling
    R_rf = inv(rotate_z(-Rf_ph)) * rotate_x(FA) * rotate_z(-Rf_ph);

    % magnetization after TE
    M = R_te * R_rf * M + b_te;

    % magnetization after TR
    M = R_tr * M + b_tr;

    % tip magnetization
    M = rotate_z(dfreq) * M;

    % magnetization over time
    M_tot(:, k + 1) = M;

    % increment phase
    Rf_ph = Rf_ph * (k * (k + 1) / 2) * inc;
end
end

function M_tot = FLASH(FA, TE, TR, freq, T1, T2, perSpoil, pTwist)
% function to simulate FLASH sequence
% perSpoil = 1/0 assume perfect spoiling
% pTwist (optional iff perSpoil = 0)
relax = 1; % only allowing option that includes relaxation

% apply RF pulse
R_rf = rotate_y(FA);
% precession after TE
[R_te, b_te] = freeprecess(TE, freq, relax, T1, T2);

% precession at start of next TR
[R_tr, b_tr] = freeprecess(TR-TE, freq, relax, T1, T2);

if perSpoil
    % remove transverse magnetization
    R_tr = ([0 0 0; 0 0 0; 0 0 1]) * R_tr;
else
    R_tr = rotate_z(pTwist) * R_tr;
end

% solve for total magnetization
M_tot = inv(eye(3) - R_te * R_rf * R_tr) * (R_te * R_rf * b_tr + b_te);

end

function S = bSSFP(FA, TE, TR, freq, relax, T1, T2)
% function to simulate bSSFP
% relax = 0/1 if including relaxation
% T1 (optional iff relax = 0)
% T2 (optional iff relax = 0)
warning('off','all')

% apply RF pulse
R_rf = rotate_y(FA);

if relax
    % precession after TE
    [R_te, b_te] = freeprecess(TE, freq, relax, T1, T2);

    % precession at start of next TR
    [R_tr, b_tr] = freeprecess(TR-TE, freq, relax, T1, T2);
else
    % precession after TE
    [R_te, b_te] = freeprecess(TE, freq, relax);

    % precession at start of next TR
    [R_tr, b_tr] = freeprecess(TR-TE, freq, relax);
end

% solve for total magnetization
M_tot = inv(eye(3) - R_te * R_rf * R_tr) * (R_te * R_rf * b_tr + b_te);
S = M_tot(1) + 1i * M_tot(2);

warning('on','all')
end

function Rz = rotate_z(phi)
% function to rotate phi (in radians) around z-axis
Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
end

function Rx = rotate_x(alpha)
% function to rotate alpha (in degrees) around x-axis
Rx = [1 0 0; 0 cosd(alpha) sind(alpha); 0 -sind(alpha) cosd(alpha)];
end

function Rx = rotate_x_rad(alpha)
% function to rotate alpha (in radians) around x-axis
Rx = [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
end

function Ry = rotate_y(alpha)
% function to rotate alpha (in degrees) around y-axis
Ry = [cosd(alpha) 0 sind(alpha); 0 1 0; -sind(alpha) 0 cosd(alpha)];
end

function [A, b] = freeprecess(T,df,relax,T1,T2)
% function to simulate precession
% relax = 0/1 if including relaxation
% T1 (optional iff relax = 0)
% T2 (optional iff relax = 0)
phi = 2*pi*df*T*(10 ^ -3);
if relax
    E1 = exp(-T/T1);
    E2 = exp(-T/T2);
else
    E1 = 1;
    E2 = 1;
end

A = [E2 0 0;0 E2 0;0 0 E1]*rotate_z(phi);
b = [0 0 1-E1]';
end

function output=ft(in)
% Performs fftshift(fft(ifftshift(input)))
% 1-D forward FT
output = fftshift(fft(ifftshift(in)));
end
