% plot ideal noiseless SAR PMF with different reset codes
clear all; close all; clc; format compact
[FONTSIZE, LINEWIDTH, FIGSIZE, SCATTER ] = figure_settings(14, 1.1, [600 400], 40);

%% parameters

% true - rerun the model simulation
% false - use saved data
RUN_SIM = true;

% bits to simulate
B = 8;

% time normalized by tau (sample timings, all time in ps)
TAU = 4;              % comparator regeneration tau
TS = 0 / TAU;       % sampling time
TSAR = 0 / TAU;     % SAR loop delay
TLATCH = 0 / TAU;    % time for comparator bit to propogate output
TREG = (1 + log(2)/2) * B + log(2)/2 * B^2

Navg = 25
TADC = Navg + (B-1)*TSAR + TLATCH + TREG    % total ADC period

% reset code
% VRESET = -2^(B-1);    % 0000...
% VRESET = 2^(B-1) - 1; % 1111...
% VRESET = 0;           % 1000...

% input format for functions (DO NOT CHANGE)
TIMING = [TAU, TADC, TS, TSAR, TLATCH];

%% run model + plot
if RUN_SIM == false
    % load saved data
    load('./results/asar_pmf_ideal_8b_reset.mat')
else % re run simulation

    % run err pmf model
    VRESET = -2^(B-1);    % 0000...
    [ err_ideal_1, pmf_ideal_1, err_folded_ideal_1, pmf_folded_ideal_1 ] = asar_meta_pmf_ideal(B, VRESET, TIMING);
    
    VRESET = 2^(B-1) - 1; % 1111...
    [ err_ideal_2, pmf_ideal_2, err_folded_ideal_2, pmf_folded_ideal_2 ] = asar_meta_pmf_ideal(B, VRESET, TIMING);
    
    VRESET = -2^(B-1) + 170;           % 10101010
    [ err_ideal_3, pmf_ideal_3, err_folded_ideal_3, pmf_folded_ideal_3 ] = asar_meta_pmf_ideal(B, VRESET, TIMING);

    save './results/asar_pmf_ideal_8b_reset.mat'
end

%% reset code 1 plot

format_error = @(err) [-1-log2(-err(err < 0)), 0, 1+log2(err(err > 0))];
xticks = [-8:1:8];

% plot full pmf (pos + neg side)
figure; hold on
stem(format_error(err_ideal_1), pmf_ideal_1, 'k', 'fill', 'MarkerSize', 3)

set(gcf, 'position', [200, 200, FIGSIZE]);
set(gca, 'yscale', 'log')
set(gca, 'ytick', 10.^[-150:10:0])
% set(gca, 'xtick', [-64, -32, 0, 32, 64])
set(gca, 'xtick', [-150:50:150])
set(gca, 'xtick', xticks)
% set(gca, 'xticklabel', [])
% set(gca, 'yticklabel', [])
set(gca, 'TickLength', [0.008 0.005])

% xlim([-2^(B-1) - 12, 2^(B-1) + 12])
% xlim([-8, 8])
xlim([-150, 150])
xlim([-9, 9])
ylim([1e-30, 10^10])
xlabel('log_2(\epsilon)')
ylabel('Pr(\epsilon|v_{id})')
grid on

%% reset code 2
% plot full pmf (pos + neg side)
figure; hold on
stem(format_error(err_ideal_2), pmf_ideal_2, 'k', 'fill', 'MarkerSize', 3)

set(gcf, 'position', [200, 200, FIGSIZE]);
set(gca, 'yscale', 'log')
set(gca, 'ytick', 10.^[-150:10:0])
% set(gca, 'xtick', [-64, -32, 0, 32, 64])
% set(gca, 'xtick', [-150:50:150])
set(gca, 'xtick', xticks)
% set(gca, 'xticklabel', [])
% set(gca, 'yticklabel', [])
set(gca, 'TickLength', [0.008 0.005])

% xlim([-2^(B-1) - 12, 2^(B-1) + 12])
% xlim([-150, 150])
xlim([-9, 9])
ylim([1e-30, 10^10])
xlabel('log_2(\epsilon)')
ylabel('Pr(\epsilon|v_{id})')
grid on


%% reset code 3
figure; hold on
stem(format_error(err_ideal_3), pmf_ideal_3, 'k', 'fill', 'MarkerSize', 3)

set(gcf, 'position', [200, 200, FIGSIZE]);
set(gca, 'yscale', 'log')
set(gca, 'ytick', 10.^[-150:10:0])
% set(gca, 'xtick', [-64, -32, 0, 32, 64])
% set(gca, 'xtick', [-150:50:150])
set(gca, 'xtick', xticks)
% set(gca, 'xticklabel', [])
% set(gca, 'yticklabel', [])
set(gca, 'TickLength', [0.008 0.005])

% xlim([-2^(B-1) - 12, 2^(B-1) + 12])
xlim([-150, 150])
xlim([-9, 9])
ylim([1e-30, 10^10])
xlabel('log_2(\epsilon)')
ylabel('Pr(\epsilon|v_{id})')
grid on
