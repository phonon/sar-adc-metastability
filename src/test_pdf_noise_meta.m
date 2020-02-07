% script used for testing/analyzing the derived analytic PDF for
% ASAR metastability with thermal noise
clear all; close all; clc; format compact;
[FONTSIZE, LINEWIDTH, FIGSIZE, SCATTER ] = figure_settings(14, 1.2, [400 280], 40);

% colors for plot
COLOR = {[1 0 0], [0 0.7 0], [0 0 1], [0.8 0 0.7], [0.0 0.9 0.9], [0.5 0.9 0.9]};

%% thermal noise voltage PDF, CDF functions
B = 8
TAU = 1;
VAVG = 0;
ENOB = 6
SIGMA = sqrt(0.5*(0.5^2)*(10.0^(-(6.02*ENOB + 1.76)/10)) - (1.0/(2.0^B))^2/12)
VDD = 1;

NPOINTS = 200;
v = [linspace(-0.5, 0, NPOINTS), 0, linspace(0, 0.5, NPOINTS)];

v_PDF = @(v) 1 ./ sqrt(2 * pi * SIGMA^2) .* exp( -(v - VAVG).^2 ./ 2 ./ SIGMA^2 );
v_PDF_half = @(v) floor(heaviside(v-VAVG)) .* ( sqrt(2) / SIGMA / sqrt(pi) .* exp( -(v - VAVG).^2 ./ 2 ./ SIGMA^2 ) );

v_CDF = @(v) 0.5 * (1 + erf( (v - VAVG) ./ SIGMA ./ sqrt(2) ) );
v_CDF_half = @(v) floor(heaviside(v-VAVG)) .* erf( (v - VAVG) ./ SIGMA ./ sqrt(2) );

figure; hold on
plot(v, v_PDF(v), 'r--');
plot(v, v_PDF_half(v), 'r');

figure; hold on
plot(v, v_CDF(v), 'r--');
plot(v, v_CDF_half(v), 'r');

%% thermal noise time PDF, CDF functions
XLIM = [1e-4, 1e2];
YLIM = [1e-100, 1e2];

t = logspace(0, 2, 20000);

% t_CDF_mu0 = @(t, vavg, sigma) 1 - erf( (VDD .* exp(-t./TAU) - vavg) / sigma / sqrt(2) );
% t_PDF_mu0 = @(t, vavg, sigma) 2 * VDD / sqrt(2*pi) / sigma / TAU .* exp( -( (VDD .* exp(-t./TAU) - vavg) ./ sigma ./ sqrt(2) ).^2 - t ./ TAU );

t_CDF = @(t, vavg, sigma) 1 - erf( - vavg / sigma / sqrt(2) ) + 0.5 * ( erf( (-VDD .* exp(-t./TAU) - vavg) / sigma / sqrt(2) ) - erf( (VDD .* exp(-t./TAU) - vavg) / sigma / sqrt(2) ) );
t_PDF = @(t, vavg, sigma) VDD / sqrt(2*pi) / sigma / TAU .* ( exp( -( (VDD .* exp(-t./TAU) - vavg) ./ sigma ./ sqrt(2) ).^2 - t ./ TAU ) + exp( -( (-VDD .* exp(-t./TAU) - vavg) ./ sigma ./ sqrt(2) ).^2 - t ./ TAU ));

t_PDF_L = @(t, tau, vavg, sigma) 2*VDD / sqrt(2*pi) / sigma / tau ./ ( 1 + erf(-vavg./sqrt(2)./sigma ) ) .* exp( -( (VDD .* exp(-t./tau) + vavg) ./ sigma ./ sqrt(2) ).^2 - t ./ tau );
% t_PDF_R = @(t, tau, vavg, sigma) 2*VDD / sqrt(2*pi) / sigma / tau ./ ( 1 + erf(-vavg./sqrt(2)./sigma ) .* ( exp( -( (VDD .* exp(-t./tau) + vavg) ./ sigma ./ sqrt(2) ).^2 - t ./ tau ) );

figure; hold on
plot(t, t_CDF(t, VAVG, SIGMA), 'b');
set(gca, 'xscale', 'log')
xlim(XLIM)

figure; hold on
plot(t, t_PDF(t, VAVG, SIGMA), 'b');
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
xlim(XLIM)
ylim(YLIM)

%% time CDF, PDF with different params

ENOB = 6
SIGMA = sqrt(0.5*(0.5^2)*(10.0^(-(6.02*ENOB + 1.76)/10)) - (1.0/(2.0^B))^2/12);
% SIGMA = 1e-5
% VAVG = [0, 0.1, 0.5, 0.999];
VAVG = [0.1, 0.05, 0.01, 1e-5, 1e-10];
T_expected = TAU .* log(VDD ./ VAVG);

figure; hold on
for i = 1:length(VAVG)
    plot(t, t_CDF(t, VAVG(i), SIGMA), 'Color', COLOR{i});
    plot([T_expected(i), T_expected(i)], [realmin, realmax], '--', 'Color', COLOR{i})
    set(gca, 'xscale', 'log')
    xlim(XLIM)
    ylim([0 1])
end

figure; hold on
for i = 1:length(VAVG)
    plot(t, t_PDF(t, VAVG(i), SIGMA), 'Color', COLOR{i});
    plot([T_expected(i), T_expected(i)], [realmin, realmax], '--', 'Color', COLOR{i})
    set(gca, 'xscale', 'log')
    set(gca, 'yscale', 'log')
    xlim([1e0, 1e2])
    ylim([1e-160, 1e2])
end

t_sat = -TAU .* log(sqrt(2) * SIGMA ./ VDD .* erfinv(0.5) )
plot([t_sat, t_sat], [realmin, realmax], 'k--')

figure; hold on
for i = 1:length(VAVG)
    plot(t, t_PDF(t, VAVG(i), SIGMA), 'Color', COLOR{i});
    plot([T_expected(i), T_expected(i)], [realmin, realmax], '--', 'Color', COLOR{i})
%     set(gca, 'xscale', 'log')
%     set(gca, 'yscale', 'log')
    xlim([1e0, 1e2])
    ylim([0 10])
end
plot([t_sat, t_sat], [realmin, realmax], 'k--')

%% single cycle timing spread from PDF

ENOB = [5, 6, 7, 7.5];
SIGMA = sqrt(0.5.*(0.5.^2).*(10.0.^(-(6.02.*ENOB + 1.76)./10)) - (1.0./(2.0.^B)).^2/12);
VAVG = 0.1;
TAVG = TAU .* log(VDD ./ VAVG);

figure; hold on
for i = 1:length(SIGMA)
    plot(t, t_PDF(t, VAVG, SIGMA(i)), 'Color', COLOR{i});
    plot([TAVG, TAVG], [realmin, realmax], 'k--');
    set(gca, 'xscale', 'log')
    set(gca, 'yscale', 'log')
    xlim([1e-1, 1e2])
    ylim([1e-200, 100])
end

figure; hold on
for i = 1:length(SIGMA)
    plot(t, t_PDF(t, VAVG, SIGMA(i)), 'Color', COLOR{i});
    plot([TAVG, TAVG], [realmin, realmax], 'k--');
%     set(gca, 'xscale', 'log')
%     set(gca, 'yscale', 'log')
    xlim([1e-1, 1e2])
    ylim([0 30])
end

%% convolving pdf test (to sum different cycles)

NPOINTS = 10000;

ENOB = 6;
SIGMA = sqrt(0.5.*(0.5.^2).*(10.0.^(-(6.02.*ENOB + 1.76)./10)) - (1.0./(2.0.^B)).^2/12);
t = linspace(0, 100, 10000);
dt = t(2) - t(1);
t_conv = dt .* [ 0 : 1 : 2*NPOINTS-2 ];

v1 = -0.1;
v2 = -0.5;
v3 = -0.7;

t_reg_1 = TAU .* log(VDD ./ abs(v1));
t_reg_2 = TAU .* log(VDD ./ abs(v2));
t_reg_3 = TAU .* log(VDD ./ abs(v3));
t_reg_1_2 = t_reg_1 + t_reg_2;
t_reg_1_2_3 = t_reg_1 + t_reg_2 + t_reg_3;

f1 = t_PDF_L(t, TAU, -abs(v1), SIGMA);
f2 = t_PDF_L(t, TAU, -abs(v2), SIGMA);
f3 = t_PDF_L(t, TAU, -abs(v3), SIGMA);
f12 = dt .* conv(f1, f2);
f123 = dt .* conv(f12, f3);

figure; hold on
plot(t, f1, 'r-')
plot(t, f2, 'b-')
plot(t_conv, f12, 'k-')
plot(t_conv, f123(1:length(t_conv)), 'g-')

% plot mean from adding expected regeneration times
% -> verify these coincide with PDF peak (e.g. mean)
plot([t_reg_1 t_reg_1], [realmin 1e2], 'r--')
plot([t_reg_2 t_reg_2], [realmin 1e2], 'b--')
plot([t_reg_1_2 t_reg_1_2], [realmin 1e2], 'k--')
plot([t_reg_1_2_3 t_reg_1_2_3], [realmin 1e2], 'g--')

% set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
ylim([1e-100, 1e2])