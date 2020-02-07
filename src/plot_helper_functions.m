%% plot analytic helper functions for metastability models
% DAC, residual voltages and output bit functions
clear all; close all; clc; format compact;
[FONTSIZE, LINEWIDTH, FIGSIZE, SCATTER ] = figure_settings(14, 1.2, [400 300], 40);

% colors for plot
COLOR = {[1 0 0], [0 0.8 0], [0 0 1], [0.8 0 0.7], [0.0 0.9 0.9], [0.5 0.9 0.9]};

SIZE = 25

%% analytic functions

VFS = 1;
% functions (INPUT/OUTPUT SHOULD BE NORMALIZED BY VLSB)
vdac_func_LSB = @(B,n,vid) floor(heaviside(n)) .* 2.^(B-n) .* (floor(vid ./ 2.^(B-n)) + 0.5);
vdac_func_VFS = @(B,n,vid) floor(heaviside(n)) .* VFS/2^n .* (floor( 2^(n) * vid/ VFS) + 0.5);

vdac_k_func_VFS = @(B,n,k) floor(heaviside(n)) .* VFS/2^n .* (floor(k ./ 2.^(B-n)) + 0.5);

VID_VFS = linspace(-0.5, 0.5, 10000);

B = 4
VLSB = VFS/2^B

k = -2^(B-1) + 0.5 : 1 : 2^(B-1) - 0.5;

%% DAC voltage
figure; hold on
for n = 1:B-1
    vdac = vdac_func_VFS(B, n, VID_VFS);
    plot(VID_VFS, vdac, 'Color', COLOR{n});
    
    scatter(k ./ 2^B, vdac_k_func_VFS(B, n, k), SIZE, 'o', 'MarkerEdgeColor', COLOR{n})
end

ylim([-0.5, 0.5])
xlabel('v_{id}/V_{FS}')
ylabel('v_{DAC}/V_{FS}')

%% residual
figure; hold on
for n = 1:B
    vres = vdac_func_VFS(B, n-1, VID_VFS) - VID_VFS;
    vres_k = vdac_k_func_VFS(B, n-1, k) - k/2^B;
    plot(VID_VFS, vres, 'Color', COLOR{n});
    
    scatter(k ./ 2^B, vres_k, SIZE, 'o', 'MarkerEdgeColor', COLOR{n})
end

ylim([-0.5, 0.5])
xlabel('v_{id}/V_{FS}')
ylabel('v_{res}/V_{FS}')

%% absolute value residual
figure; hold on
for n = 1:B
    vres = abs(vdac_func_VFS(B, n-1, VID_VFS) - VID_VFS);
    vres_k = abs(vdac_k_func_VFS(B, n-1, k) - k/2^B);
    plot(VID_VFS, vres, 'Color', COLOR{n});
    
    scatter(k ./ 2^B, vres_k, SIZE, 'o', 'MarkerEdgeColor', COLOR{n})
end

ylim([0, 0.5])
xlabel('v_{id}/V_{FS}')
ylabel('|v_{res}|/V_{FS}')

%% output codes
DRESET = -2^(B-1) % 0000...
% DRESET = 0 % 1000...
% DRESET = -1 % 0111...
% DRESET = 2^(B-1) - 1 % 1111...

D_out0_n_k = @(B,n,k) 2.^(B-n) .* (floor(k ./ 2.^(B-n)));
D_reset_n = @(Dreset,B,n) 2.^(B-n) .* (floor(Dreset ./ 2.^(B-n)));
D_out_n_k = @(Dreset,B,n,k) Dreset + floor(heaviside(n)) .* 2.^(B-n) .* ( floor(k ./ 2.^(B-n)) - floor(Dreset ./ 2.^(B-n)) );

% output when reset = 0
figure; hold on
plot([-0.5, 0.5], [DRESET, DRESET], 'k--')
for n = 1:B
    Dout0 = D_out0_n_k(B,n,k);
    scatter(k ./ 2^B, Dout0, SIZE, 'o', 'MarkerEdgeColor', COLOR{n})
end
xlabel('v_{id}/V_{FS}')
ylabel('D_{out0}')

% leftover reset code voltage levels on each cycle
figure; hold on
plot([-0.5, 0.5], [DRESET, DRESET], 'k--')
for n = 1:B
    Dresetn = D_reset_n(DRESET,B,n)
    plot([-0.5, 0.5], [Dresetn, Dresetn], 'Color', COLOR{n})
end
xlabel('v_{id}/V_{FS}')
ylabel('D_{RESET}')

% actual output code
figure; hold on
plot([-0.5, 0.5], [DRESET, DRESET], 'k--')
for n = 1:B
    Dout = D_out_n_k(DRESET,B,n,k);
    scatter(k ./ 2^B, Dout, SIZE, 'o', 'MarkerEdgeColor', COLOR{n})
end
xlabel('v_{id}/V_{FS}')
ylabel('D_{out}')
