%% plot average total SAR regeneration time
clear all; close all; clc; format compact;
[FONTSIZE, LINEWIDTH, FIGSIZE, SCATTER ] = figure_settings(14, 1.2, [400 300], 40);

% colors for plot
COLOR = {[1 0 0], [0 0.7 0], [0 0 1], [0.8 0 0.7], [0.0 0.9 0.9], [0.5 0.9 0.9]};

%% settings

% analytic functions (INPUT/OUTPUT SHOULD BE NORMALIZED BY VLSB)
vdac_func = @(B,n,vid) floor(heaviside(n)) .* -2.^(B-n) .* (floor(vid ./ 2.^(B-n)) + 0.5);
vout0_func = @(B,n,vid) 2.^(B-n) .* (floor(vid ./ 2.^(B-n)));
vreset_n_func = @(B,n,vreset) 2.^(B-n) .* (floor(vreset ./ 2.^(B-n)));
vout_n_func = @(B,n,vid,vreset) vreset + 2.^(B-n) .* (floor(vid ./ 2.^(B-n)) - floor(vreset ./ 2.^(B-n)));
vout_ideal_func = @(vid) floor(vid);
vdac_k_func = @(B,n,k) floor(heaviside(n)) .* -2.^(B-n) .* (floor(k ./ 2.^(B-n)) + 0.5);

% bits
B = 9;

% full scale and ref
VFS = 1;      % [V]
VDD = 1;
VLSB = VFS/2^B;

VID = linspace(-VFS/2, VFS/2, 10000);

% DAC function (digital string of bits -> analog output in LSBs)
DAC = @(s) bin2dec(s) - 2^(B-1) + 0.5;

VID_VLSB = VID/VLSB;
VOUT_IDEAL_VLSB = vout_ideal_func(VID_VLSB);

% comparator gain
A = 1;

%% plot timing "Teasy" figure and average

% full version of Tavg (for any VDD, VFS)
B_range = 1:1:B;
Tavg = B_range .* (1 + log(2) / 2 + log(VDD/VFS)) + B_range.^2 .* log(2) ./ 2;
Tavg_numerical = zeros(1, length(B_range));

% higher resolution vres
VID_HIGH_RES_LOG = logspace(-8, -4, 100);
VID_HIGH_RES_LIN = linspace(10^-3.9, 0.5, 100);
VID_HIGH_RES = [VID_HIGH_RES_LOG, VID_HIGH_RES_LIN];

figure; hold on;
for i = 1:length(B_range)
    b = B_range(i)
    
    VLSB = 1 / 2^b;
    
    Treg = zeros(size(VID_VLSB)) + b * log(VDD/VLSB/A);
    for k = 0:b-1
        Treg = Treg - log(VID_VLSB + vdac_func(B,k,VID_VLSB));
    end
    
    % color for plotting
    col = COLOR{ceil(i/2)};
    col_darker = 0.5 .* col;
    
    % running integral across k
    Treg_integral = 0;
    
    % create higher resolution vid and integrate in each region
    for k = -2^(b-1) + 0.5 : 1 : 2^(b-1) - 0.5
        
        Treg_high_res_p = zeros(size(VID_HIGH_RES)) + b * log(VDD/VLSB/A);
        Treg_high_res_m = zeros(size(VID_HIGH_RES)) + b * log(VDD/VLSB/A);
        
        for n = 0:b-1
            vdac = vdac_func(b,n,k);
            
            vres_p = k + vdac + 0.5 - VID_HIGH_RES;
            vres_m = k + vdac - 0.5 + VID_HIGH_RES;
            vres_p(~vres_p) = realmin;
            vres_m(~vres_m) = realmin;

        	Treg_high_res_p = Treg_high_res_p - log(abs(vres_p));
            Treg_high_res_m = Treg_high_res_m - log(abs(vres_m));
        end
        
        Treg_integral = Treg_integral + abs(trapz(- VID_HIGH_RES, Treg_high_res_p)) + abs(trapz(VID_HIGH_RES, Treg_high_res_m));
        
        if mod(i,2) == 1
            plot((k + 0.5 - VID_HIGH_RES) * VLSB, Treg_high_res_p, 'Color', col);
            plot((k - 0.5 + VID_HIGH_RES) * VLSB, Treg_high_res_m, 'Color', col);
        end
    end
    
    % store numerical result
    Tavg_numerical(i) = Treg_integral ./ 2^b;
    
    % plot analytic average
    if mod(i,2) == 1
        plot([-0.5, 0.5], [Tavg(i), Tavg(i)], ':', 'Color', col_darker);
    end
end

Tavg
Tavg_numerical

% set(gca, 'yscale', 'log')
set(gca, 'ytick', 0:8:64)
xlim([-0.5, 0.5])
ylim([0, 50])
set(gcf, 'position', [400 200 300 300])
xlabel('v_{id}/V_{FS}')
ylabel('T_{reg,n}/\tau')


%% compare simplified versions of Tavg vs. Teasy
% VDD = VFS

% Tavg analytic form
B_range_continuous = linspace(0, 12, 100);
T_avg = t_avg(B_range_continuous);

% Teasy analytic
T_easy = t_easy(B_range_continuous);

% plot
figure; hold on
h(1) = plot(B_range_continuous, T_avg ,'k-');
h(2) = plot(B_range_continuous, T_easy ,'r-');
scatter(B_range, Tavg_numerical , 45, 'o', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.2);

set(gca, 'xtick', 0:2:12)
set(gca, 'ytick', 0:10:60)
xlim([0, 12])
ylim([0, 60])
set(gcf, 'position', [600 200 300 300])
xlabel('B')
ylabel('T/\tau')
legend(h, {'Avg', 'Easy'}, 'Location', 'SouthEast')