%% used to plot ideal metastability timing bands visualization
clear all; close all; clc; format compact;
[FONTSIZE, LINEWIDTH, FIGSIZE, SCATTER ] = figure_settings(14, 1.2, [800 300], 40);

% colors for plot
COLOR = {[1 0 0], [0 0.8 0], [0 0 1], [0.8 0 0.7], [0.0 0.9 0.9], [0.5 0.9 0.9]};

%% analytic functions

% functions (INPUT/OUTPUT SHOULD BE NORMALIZED BY VLSB)
vdac_func = @(B,n,vid) floor(heaviside(n)) .* -2.^(B-n) .* (floor(vid ./ 2.^(B-n)) + 0.5);
vout0_func = @(B,n,vid) 2.^(B-n) .* (floor(vid ./ 2.^(B-n)));
vreset_n_func = @(B,n,vreset) 2.^(B-n) .* (floor(vreset ./ 2.^(B-n)));
vout_n_func = @(B,n,vid,vreset) vreset + 2.^(B-n) .* (floor(vid ./ 2.^(B-n)) - floor(vreset ./ 2.^(B-n)));
vout_ideal_func = @(vid) floor(vid);

% bits
B = 3;

% full scale and ref
VFS = 1;      % [V]
VDD = 1;
VLSB = VFS/2^B;

VID = linspace(-VFS/2, VFS/2, 1000);
VID_VLSB = VID ./ VLSB;

% create VID space
VID_VLSB = [];
for i = -2^(B-1):2^(B-1)-1
    VID_VLSB = [VID_VLSB, linspace(i + 1e-6, i+1-1e-6, 200)];
end
% manually fix endpoints
VID_VLSB = [-2^(B-1), VID_VLSB, 2^(B-1)];

VID = VID_VLSB .* VLSB;

% DAC function (digital string of bits -> analog output in LSBs)
DAC = @(s) bin2dec(s) - 2^(B-1) + 0.5;

%% plot residuals

figure; hold on;
for n = 1:B
    vres = VID_VLSB + vdac_func(B,n-1,VID_VLSB);
    plot(VID_VLSB * VLSB, vres, 'Color', COLOR{n});
end

% set(gca, 'yscale', 'log')
% set(gca, 'ytick', 10.^[-8:2:8])
xlim([-0.5, 0.5])
ylim([-2^(B-1), 2^(B-1)])
set(gca, 'ytick', -2^(B-1) : 2 : 2^(B-1))
set(gcf, 'position', [600 400 FIGSIZE])
xlabel('v_{id}/V_{FS}')
ylabel('v_{res,n}/V_{LSB}')

%% plot total time

tres = zeros(B, length(VID_VLSB));

TAU = 2;
TS = 4; % sample
TSAR = 4;
TLATCH = 0;

figure; hold on;
for n = 1:B
    vres = VID_VLSB + vdac_func(B,n-1,VID_VLSB);
    tres(n+1,:) = tres(n,:) + TAU .* log( VDD ./ abs(VLSB .* vres) );
end

for n = B:-1:1
    % bit time
    t_bit = TS + tres(n+1,:) + TSAR * (n-1) + TLATCH;
    plot([-0.5, VID, 0.5], [0, t_bit, 0], 'Color', COLOR{n})
    
    if n < B
        t_cmp = TS + tres(n+1,:) + TSAR * n;
        plot([-0.5, VID, 0.5], [0, t_cmp, 0], 'Color', COLOR{n} ./ 1.5)
    end
%     fill([-0.5, VID, 0.5], [0, t, 0], COLOR{n}, 'LineStyle', 'none')
end

% finally plot sample
plot([-0.5, -0.5, 0.5, 0.5], [0, TS, TS, 0], 'Color', [0.5, 0.5, 0.5])
% fill([-0.5, VID, 0.5], [0, TS .* ones(size(VID)), 0], [0.5, 0.5, 0.5], 'LineStyle', 'none')

% set(gca, 'yscale', 'log')
% set(gca, 'ytick', 10.^[-8:2:8])
set(gca,'Visible','off')
xlim([-0.5, 0.5])
ylim([0, 36])
set(gcf, 'Renderer', 'painters');
set(gca, 'ytick', -2^(B-1) : 2 : 2^(B-1))
set(gcf, 'position', [600 200 FIGSIZE])
xlabel('v_{id}/V_{FS}')

