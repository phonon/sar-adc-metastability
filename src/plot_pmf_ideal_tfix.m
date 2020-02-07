% plot ideal noiseless ASAR PMF with different TFIX, TSAR ratios
clear all; close all; clc; format compact
[FONTSIZE, LINEWIDTH, FIGSIZE, SCATTER ] = figure_settings(14, 1.1, [300 320], 30);

%% parameters

% true - rerun the model simulation
% false - use saved data
RUN_SIM = true;

% bits to simulate
B = 8;

% average regeneration time required (sets floor for total time)
TREG = (1 + log(2)/2) * B + log(2)/2 * B^2

% reset code
VRESET = -2^(B-1);    % 0000...
% VRESET = 2^(B-1) - 1; % 1111...
% VRESET = 0;           % 1000...

% =====================================
% TIME SET 1: Compare constant TADC with TFIX change
% time normalized by tau (sample timings, all time in ps)
Navg = 25;     % average TOTAL regeneration time
TS = 0;        % sampling time
TSAR = 3;     % SAR loop delay
TLATCH = 0;    % time for comparator bit to propogate output
TADC = Navg + TREG;
Navg_1 = TADC - TS - (B-1)*TSAR - TLATCH - TREG

% timing set 1: constant TADC
TIMING_a0 = [0, TADC, 0, 0, 0];
TIMING_a1 = [0, TADC, TS, TSAR, TLATCH];

% =====================================
% TIME SET 2: Compare fixed <Ntotal> with/without TFIX
% time normalized by tau (sample timings, all time in ps)
Navg = 25;     % average TOTAL regeneration time
TS = 0;        % sampling time
TSAR = 3;     % SAR loop delay
TLATCH = 0;    % time for comparator bit to propogate output
% TADC = TS + (B-1)*TSAR + TLATCH + Navg   % total ADC period
% TADC = 64;
% Navg = TADC - TS - (B-1)*TSAR - TLATCH - TREG

% timing set 1: constant TADC
TIMING_b0 = [0, Navg + TREG, 0, 0, 0];
TIMING_b1 = [0, Navg + TREG + (B-1)*TSAR, TS, TSAR, TLATCH];
% =====================================


%%
if RUN_SIM == false
    % load saved data
    load('./results/asar_pmf_ideal_8b_compare_TFIX.mat')
else % re run simulation
    %% do first order estimate of error
    VDD = 1;
    VFS = 1;
    Treg_avg = B * (1 + log(2)/2 + log(VDD/VFS)) + B^2 * log(2) / 2
    TADC_AVG = TS + (B-1) * TSAR + TLATCH + Treg_avg
    Pmeta_avg = 2*(2^B-1) .* exp(- (TADC - TADC_AVG) )

    %% run model + plot

    % TIMING SET 1
    % TADC = CONSTANT, TFIX INCREASE
    % pmf (no comparator thermal noise)
    [err_ideal_a0, pmf_ideal_a0, err_folded_ideal_a0, pmf_folded_ideal_a0] = asar_meta_pmf_ideal(B, VRESET, TIMING_a0);
    [err_ideal_a1, pmf_ideal_a1, err_folded_ideal_a1, pmf_folded_ideal_a1] = asar_meta_pmf_ideal(B, VRESET, TIMING_a1);

    % TIMING SET 2
    % Ntotal = CONSTANT, TFIX INCREASE
    % pmf (no comparator thermal noise)
    [err_ideal_b0, pmf_ideal_b0, err_folded_ideal_b0, pmf_folded_ideal_b0] = asar_meta_pmf_ideal(B, VRESET, TIMING_b0);
    [err_ideal_b1, pmf_ideal_b1, err_folded_ideal_b1, pmf_folded_ideal_b1] = asar_meta_pmf_ideal(B, VRESET, TIMING_b1);
    
    save './results/asar_pmf_ideal_8b_compare_TFIX.mat'
end

%% plot pmf

FIGSIZE = [640 280];
SIZE = 35;

format_error = @(err) [-1-log2(-err(err < 0)), 0, 1+log2(err(err > 0))];
xticks = [-8:1:8];

% plot full pmf (pos + neg side) (uncomment stem for different aesthetic)
figure; hold on
% stem(err_ideal_a0, pmf_ideal_a0, 'r', 'fill', 'MarkerSize', 6)
% stem(err_ideal_a0, pmf_ideal_a1, 'k', 'fill', 'MarkerSize', 6)
% stem(format_error(err_ideal_a0), pmf_ideal_a0, 'r', 'fill', 'MarkerSize', 6)
% stem(format_error(err_ideal_a0), pmf_ideal_a1, 'k', 'fill', 'MarkerSize', 6)
scatter(format_error(err_ideal_a0), pmf_ideal_a0, SIZE, 'filled', 'MarkerFaceColor', [1 0 0])
scatter(format_error(err_ideal_a0), pmf_ideal_a1, SIZE, 'filled', 'MarkerFaceColor', [0 0 0])

set(gcf, 'position', [100 100 FIGSIZE])
set(gca, 'yscale', 'log')
xlim([-2^(B-1) - 4, 2^(B-1) + 4])
% xlim([-150, 150])
xlim([-9, 9])
ylim([1e-40, 1e0])
xlabel('log_2(\epsilon)')
ylabel('Pr(\epsilon|v_{id})')
grid on 

set(gca, 'ytick', 10.^[-140:10:0])
% set(gca, 'xtick', [-150:50:150])
set(gca, 'xtick', xticks)
% set(gca, 'xticklabel', [])
% set(gca, 'yticklabel', [])
set(gca, 'TickLength', [0.015 0.005])


%% timing set 2 (uncomment stem for different aesthetic)
figure; hold on
% stem(err_ideal_b0, pmf_ideal_b0, 'r', 'fill', 'MarkerSize', 6)
% stem(err_ideal_b0, pmf_ideal_b1, 'k', 'fill', 'MarkerSize', 6)
% stem(format_error(err_ideal_b0), pmf_ideal_b0, 'r', 'fill', 'MarkerSize', 6)
% stem(format_error(err_ideal_b0), pmf_ideal_b1, 'k', 'fill', 'MarkerSize', 6)
scatter(format_error(err_ideal_b0), pmf_ideal_b0, SIZE, 'filled', 'MarkerFaceColor', [1 0 0])
scatter(format_error(err_ideal_b0), pmf_ideal_b1, SIZE, 'filled', 'MarkerFaceColor', [0 0 0])

set(gcf, 'position', [100 100 FIGSIZE])
set(gca, 'yscale', 'log')
xlim([-2^(B-1) - 4, 2^(B-1) + 4])
% xlim([-150, 150])
xlim([-9, 9])
ylim([1e-40, 1e0])
xlabel('log_2(\epsilon)')
ylabel('Pr(\epsilon|v_{id})')
grid on 

set(gca, 'ytick', 10.^[-140:10:0])
% set(gca, 'xtick', [-150:50:150])
set(gca, 'xtick', xticks)
% set(gca, 'xticklabel', [])
% set(gca, 'yticklabel', [])
set(gca, 'TickLength', [0.015 0.005])
