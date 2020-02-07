% plot FIRST ORDER Pmeta expressions for async and sync SAR
clear all; close all; clc; format compact
[FONTSIZE, LINEWIDTH, FIGSIZE, SCATTER ] = figure_settings(14, 1.1, [600 300], 40);

%%

% bits to simulate
B = 1:1:10;

Teasy = t_easy(B);

N = linspace(0, 40, 200);

for i = 1:length(B)
    Pmeta_sync(i,:) = 2 .* (2.^B(i) - 1) .* exp(-N);
    Pmeta_sync(i,Pmeta_sync(i,:) > 1) = 1;
    
    Pmeta_async(i,:) = 2.^(B(i)+1) .* exp( -(B(i).*N - Teasy(i)) );
    Pmeta_async(i,Pmeta_async(i,:) > 1) = 1;
end

%% plot

FIGSIZE = [420 400];

figure; hold on

h(1) = plot(N, Pmeta_sync(6,:), 'k');
plot(N, Pmeta_sync(7,:), 'k')
plot(N, Pmeta_sync(8,:), 'k')

h(2) = plot(N, Pmeta_async(6,:), 'r');
plot(N, Pmeta_async(7,:), 'r')
plot(N, Pmeta_async(8,:), 'r')

set(gcf, 'position', [200, 200, FIGSIZE]);
set(gca, 'yscale', 'log')

xlim([0, 40])
ylim([1e-20, 10^0])
xlabel('N')
ylabel('Pr(\epsilon|v_{id})')
grid on
legend(h, {'Sync (6b-8b)', 'Async (6b-8b)'}, 'Location', 'NorthEast')