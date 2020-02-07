% test model code for metastability + thermal noise in 
% asar_meta_pmf_noise_combined, included here in non-function form
% - plots PDFs of t > TSAR during each cycle
% - plots error PMF for single set of parameters
clear all; close all; clc; format compact;
[FONTSIZE, LINEWIDTH, FIGSIZE, SCATTER ] = figure_settings(14, 1.2, [400 280], 40);

% colors for plot
COLOR0 = [0.5, 0.5, 0.5];
COLOR = {[1 0 0], [0 0.7 0], [0 0 1], [0.8 0 0.7], [0.0 0.9 0.9], [0.5 0.9 0.9]};

%%

% bits
B = 4;

% reset code
VRESET = -2^(B-1);    % 0000...
% VRESET = 2^(B-1) - 1; % 1111...
% VRESET = 0;           % 1000...

% time normalized by tau (sample timings)
N = 10;
TAU = 1;              % comparator regeneration tau
TLOGIC = 10;     % SAR loop delay
TADC = B * N + (B-1)*TLOGIC;   % total ADC period

% (for comparison with old model)
TIMING = [0, TADC, 0, TLOGIC, 0];

% 
% % get timing from inputs
% TAU = TIMING(1);
% TADC = TIMING(2);
% TS = TIMING(3);
% TSAR = TIMING(4);
% TLATCH = TIMING(5);

% full scale and ref
VDD = 1;
VFS = 1;         % [V]
VLSB = VFS/2^B;
VLSB_HALF = VLSB / 2; % half lsb

HALF_VLSB_RANGE_HALF = [logspace(-10, -2, 50), linspace(1.001e-2, 0.5, 50)];

VLSB_RANGE = [ HALF_VLSB_RANGE_HALF, -fliplr(HALF_VLSB_RANGE_HALF(2:end)) ];
VLSB_VFS_RANGE = VLSB .* VLSB_RANGE;

VOFFSET = VLSB_HALF .* [ -ones(size(HALF_VLSB_RANGE_HALF)), ones(size(HALF_VLSB_RANGE_HALF(2:end))) ];


% time range for t > T
% -> dt SETS ACCURACY OF PMF AT COMPUTATION SPEED TRADEOFF
dt = 0.1;

T_RANGE = dt .* [ 0 : 1 : ceil(1.5*TADC / dt) ];
NPOINTS_T = length(T_RANGE);

% std deviation of noise for comparator in terms of vlsb
% 0.5 -> VLSB/2
ENOB = 3.8
STDCOMP = sqrt(0.5*(0.5^2)*(10.0^(-(6.02*ENOB + 1.76)/10)) - (1.0/(2.0^B))^2/12)

%% probability distribution funcs

% t_CDF = @(t, tau, vavg, sigma) 1 + 0.5 * ( erf( (-VDD .* exp(-t./tau) - vavg) / sigma / sqrt(2) ) - erf( (VDD .* exp(-t./tau) - vavg) / sigma / sqrt(2) ) );
% t_CDF_inv = @(t, tau, vavg, sigma) - 0.5 * ( erf( (-VDD .* exp(-t./tau) - vavg) / sigma / sqrt(2) ) - erf( (VDD .* exp(-t./tau) - vavg) / sigma / sqrt(2) ) );
% t_CDF_inv = @(t, tau, vavg, sigma) 0.5 * ( erf( (VDD .* exp(-t./tau) + vavg) / sigma / sqrt(2) ) + erf( (VDD .* exp(-t./tau) - vavg) / sigma / sqrt(2) ) );
% t_CDF_inv = @(t, tau, vavg, sigma) 0.5 * ( erf( (VDD .* exp(-t./tau) + vavg) / sigma / sqrt(2) ) + erf( (VDD .* exp(-t./tau) - vavg) / sigma / sqrt(2) ) );

% vDAC - vid - vn < 0 : left conditional
% t_PDF_L = @(t, tau, vavg, sigma) 2*VDD / sqrt(2*pi) / sigma / tau ./ ( 1 + erf(-vavg./sqrt(2)./sigma ) ) .* exp( -( (VDD .* exp(-t./tau) + vavg) ./ sigma ./ sqrt(2) ).^2 - t ./ tau );
t_PDF_L_pr0 = @(t, tau, vavg, sigma) dt .* VDD / sqrt(2*pi) / sigma / tau .* exp( -( (VDD .* exp(-t./tau) + vavg) ./ sigma ./ sqrt(2) ).^2 - t ./ tau );

% vDAC - vid - vn > 0 : right conditional
% t_PDF_R = @(t, tau, vavg, sigma) 2*VDD / sqrt(2*pi) / sigma / tau ./ ( 1 - erf(-vavg./sqrt(2)./sigma ) ) .* exp( -( (VDD .* exp(-t./tau) - vavg) ./ sigma ./ sqrt(2) ).^2 - t ./ tau );
t_PDF_R_pr1 = @(t, tau, vavg, sigma) dt .* VDD / sqrt(2*pi) / sigma / tau .* exp( -( (VDD .* exp(-t./tau) - vavg) ./ sigma ./ sqrt(2) ).^2 - t ./ tau );

% combined PDF (not used by code at moment)
% t_PDF = @(t, tau, vavg, sigma) dt .* VDD / sqrt(2*pi) / sigma / tau .* ( exp( -( (VDD .* exp(-t./tau) - vavg) ./ sigma ./ sqrt(2) ).^2 - t ./ tau ) + exp( -( (-VDD .* exp(-t./tau) - vavg) ./ sigma ./ sqrt(2) ).^2 - t ./ tau ));

%% voltage/output/error functions

% discrete DAC for branch b
vdac_b_func = @(B,n,b) floor(heaviside(n)) .* ( 2.^(B-n) .* floor( b ) - 2.^(B-n-1) - 2.^(B-1) );

% discrete output and errors functions
verr_nbk_func = @(vreset,B,n,b,k) vreset - floor(heaviside(n)) .* ( 2.^(B-n) .* ( floor(b) + floor(vreset ./ 2.^(B-n)) - 2^(n-1) ) ) - floor(k);

% error voltage -> error probability vector index function
error_bin_func = @(k) 1 + 2^B + floor(k);

%% run

% error domain and output pmf (to calculate)
err_range = -2^B : 1 : 2^B;
pmf = zeros(1, 2^(B+1) + 1);

% VLSB index (half-integer k VLSB)
k_range = -2^(B-1) + 0.5 : 1 : 2^(B-1) - 0.5;
kVLSB = k_range .* VLSB;
len_k_range = length(k_range);

% prealloc arrays
vcenter = {};
err = {};
bits = {}; % tracks bit decisions as booleans, used for probability calc
ft_vid = {}; % pdf time distribution in branch

% t1 + t2 + ... > T - (n-1)Tlogic : metastable during regeneration time
pr_meta_during_reg = {};

% t1 + t2 + ... > T - nTlogic : metastable during logic between cycles
pr_meta_during_logic = {};

% prealloc arrays for "n = 0, b = 0"
vcenter{1}{1} = zeros(B+1, length(k_range));
err{1}{1} = zeros(B+1, length(k_range));
bits{1}{1} = [];
ft_vid{1}{1} = zeros(B+1, length(k_range));

pr_meta_during_reg{1}{1} = zeros(1, length(k_range));
pr_meta_during_logic{1}{1} = zeros(1, length(k_range));
pr_meta_during_logic{1}{2} = zeros(1, length(k_range));

% n = 0 initial values
vcenter{1}{1}(1,:) = VLSB .* k_range;
err{1}{1}(1,:) = verr_nbk_func(VRESET, B, 0, 1, k_range);

% metastability probability on first cycle
pr_MSB = 0;

% figure for plotting PDF distributions
figure; hold on
set(gca, 'yscale', 'log')
ylim([1e-40, 9.99e-9])
xlim([-0.5, 0.5])
xlabel('v_{id}/V_{FS}')
ylabel('f_{t > T}(v_{id})')

for n = 1:1:B
    
    % first generate next layer in tree
    % -> used for binning after the comparator resolves
    for b = 1 : 1 : 2^n
        % get parent branch index
        b_parent = round(b/2);

        % copy arrays from parent
        vcenter{n+1}{b} = vcenter{n}{b_parent};
        err{n+1}{b} = err{n}{b_parent};
        bits{n+1}{b} = logical([bits{n}{b_parent}, mod(b,2)==1]);

        % residual for this branch
        vcenter{n+1}{b}(n+1,:) = VLSB .* (k_range + vdac_b_func(B,n,b));

        % error for this cycle + branch
        err{n+1}{b}(n+1,:) = verr_nbk_func(VRESET, B, n, b, k_range);
        
        % next cycle prob 
        pr_meta_during_reg{n+1}{b} = zeros(1, length(k_range));
        pr_meta_during_logic{n+1}{2*b-1} = zeros(1, length(k_range));
        pr_meta_during_logic{n+1}{2*b} = zeros(1, length(k_range));
        
        % next time distributions
        ft_vid{n+1}{b} = {};
    end
    
    
    % go through current layer in tree, do metastable
    for b = 1 : 1 : 2^(n-1)
        
        % get parent branch index
        b_parent = round(b/2);
        
        % get branch center voltages
        vcenter_n_b = vcenter{n}{b}(n,:);
        
        % error pmf bin index during comparison phase
        error_bins = error_bin_func(err{n}{b}(n,:));
        
        % error bins when bit resolves to 1 or 0
        error_bins_1 = error_bin_func(err{n+1}{b*2-1}(n+1,:));
        error_bins_0 = error_bin_func(err{n+1}{b*2}(n+1,:));
        
        % iterate over k and calc
        for k = 1:len_k_range
            
            % voltage range in this k-index during cycle n / branch b
            v_range = (vcenter_n_b(k) + VOFFSET) + VLSB_VFS_RANGE;
            
            % for integrating and plotting, current chunk in input range
            vid_interval = kVLSB(k) + VOFFSET + VLSB_VFS_RANGE;
            
            % actual voltage will be defined as
            % v = (vres + voffset) + vid
            
            f0_greater_t = zeros(size(VLSB_VFS_RANGE));
            f1_greater_t = zeros(size(VLSB_VFS_RANGE));
            f_greater_t = zeros(size(VLSB_VFS_RANGE));
            f_resolves_0_meta = zeros(size(VLSB_VFS_RANGE));
            f_resolves_1_meta = zeros(size(VLSB_VFS_RANGE));
            
            % to save distributions for next cycle
%             ft_vid{n+1}{2*b-1}{k} = zeros(length(VLSB_VFS_RANGE), NPOINTS_T);
%             ft_vid{n+1}{2*b}{k} = zeros(length(VLSB_VFS_RANGE), NPOINTS_T);
            
            
            % get probability t1 + t2 + ... > T in this interval
            
            ft_vid_1_save = zeros(length(VLSB_VFS_RANGE), NPOINTS_T);
            ft_vid_0_save = zeros(length(VLSB_VFS_RANGE), NPOINTS_T);
            
            % first cycle, do not need convs or parallelization
            if n == 1

                for i = 1:length(VLSB_VFS_RANGE)

                    v = v_range(i);

                    f1 = t_PDF_R_pr1(T_RANGE, TAU, v, STDCOMP);
                    f0 = t_PDF_L_pr0(T_RANGE, TAU, v, STDCOMP);

                    % cache for next cycle
                    ft_vid{n+1}{2*b-1}{k}(i,:) = f1;
                    ft_vid{n+1}{2*b}{k}(i,:) = f0;
                    
                    % ===============================
                    % probability during conversion
                    % ===============================
                    
                    index = T_RANGE > TADC;
                    
                    f0_greater_t(i) = trapz(f0(index));
                    f1_greater_t(i) = trapz(f1(index));

                    f_greater_t(i) = f0_greater_t(i) + f1_greater_t(i);

                    % ===============================
                    % probability of resolving, then metastable
                    % ===============================
                    
                    if B == 1
                        index = T_RANGE < TADC;
                    else % B > 1
                        index = T_RANGE > (TADC - TLOGIC);
                    end
                    
                    % pr -> 0
                    f_resolves_0_meta(i) = trapz(f0(index));

                    % pr -> 1
                    f_resolves_1_meta(i) = trapz(f1(index));

                end
                
            % 2nd+ cycles, parallelize convolutions
            else
                f_prev_n_b = ft_vid{n}{b}{k};
                
                for i = 1:length(VLSB_VFS_RANGE)
                
                    % previous distribution (cached)
                    f_prev = f_prev_n_b(i,:);

                    % current vres
                    v = v_range(i);
                    f0n = t_PDF_L_pr0(T_RANGE, TAU, v, STDCOMP);
                    f1n = t_PDF_R_pr1(T_RANGE, TAU, v, STDCOMP);

                    f0 = conv(f_prev, f0n);
                    f1 = conv(f_prev, f1n);

                    % save for next cycle
                    % -> truncate the convolution at NPOINTS_T
                    %    for performance (most of the large T range
                    %    should not matter, adjust TADC range in
                    %    parameters)
%                     ft_vid{n+1}{2*b-1}{k}(i,:) = f1(1:NPOINTS_T);
%                     ft_vid{n+1}{2*b}{k}(i,:) = f0(1:NPOINTS_T);

                    ft_vid_1_save(i,:) = f1(1:NPOINTS_T);
                    ft_vid_0_save(i,:) = f0(1:NPOINTS_T);

                    
                    % ===============================
                    % probability during conversion
                    % ===============================
                    index = T_RANGE > (TADC - (n-1) * TLOGIC);
                    f0_greater_t(i) = trapz(f0(index));
                    f1_greater_t(i) = trapz(f1(index));

                    f_greater_t(i) = f0_greater_t(i) + f1_greater_t(i);

                    % ===============================
                    % probability of resolving, then metastable
                    % ===============================

                    if n < B
                        index = T_RANGE > (TADC - n * TLOGIC);

                        % pr -> 0
                        f_resolves_0_meta(i) = trapz(f0(index));

                        % pr -> 1
                        f_resolves_1_meta(i) = trapz(f1(index));

                    else % LAST CYCLE n == B

                        % bin remaining time into non-metastable cycles
                        index = T_RANGE < ( TADC - (B-1) * TLOGIC );

                        % pr -> 0
                        f_resolves_0_meta(i) = trapz(f0(index));

                        % pr -> 1
                        f_resolves_1_meta(i) = trapz(f1(index));

                    end
                end

                % save results from for loop
                ft_vid{n+1}{2*b-1}{k} = ft_vid_1_save;
                ft_vid{n+1}{2*b}{k} = ft_vid_0_save;
            end
            
            % ============
            % PLOT
            % ============
            fig = plot(vid_interval, f_greater_t, 'Color', COLOR{n});
            uistack(fig, 'bottom');
            
            % ================================
            % Part I: t1 + t2 + t3 + ... + tn > T - (n-1)*Tlogic
            %         ends during regeneration
            % ===============================
            
            % probability of ending during the cycle
            pr_ends_during_reg_0 = abs(trapz(vid_interval, f0_greater_t));
            pr_ends_during_reg_1 = abs(trapz(vid_interval, f1_greater_t));
            
            % save TOTAL probability of metastability:
            % t1 + t2 + ... > T - (n-1) * TLOGIC
            pr_meta_during_reg{n}{b}(k) = pr_ends_during_reg_0 + pr_ends_during_reg_1;
            
            % first cycle MSB unique, do not need subtract previous 
            if n == 1
                % bin error
                err_index = error_bins(k);
                pmf(err_index) = pmf(err_index) + pr_ends_during_reg_0 + pr_ends_during_reg_1;
                
                % first cycle sanity checking for error/accuracy
                pr_MSB = pr_MSB + pr_ends_during_reg_0 + pr_ends_during_reg_1;
            else
                % bin error, subtract Pr of metastability from in logic delay from previous cycle
                err_index = error_bins(k);
                pmf(err_index) = pmf(err_index) + pr_ends_during_reg_0 + pr_ends_during_reg_1 - pr_meta_during_logic{n-1}{b}(k);
            end
            
            % ================================
            % Part II: t1 + t2 + t3 + ... + tn > T - n*Tlogic
            %         ends before next cycle
            %         (for n < B)
            % ===============================
            if n < B
                % branch pr output 1
                pr_resolves_1 = abs(trapz(vid_interval, f_resolves_1_meta));
                err_index = error_bins_1(k);
                pmf(err_index) = pmf(err_index) + pr_resolves_1 - pr_ends_during_reg_1;
                
                % branch pr output 0
                pr_resolves_0 = abs(trapz(vid_interval, f_resolves_0_meta));
                err_index = error_bins_0(k);
                pmf(err_index) = pmf(err_index) + pr_resolves_0 - pr_ends_during_reg_0;
                
                % save results (for next cycles)
                pr_meta_during_logic{n}{2*b-1}(k) = pr_resolves_1;
                pr_meta_during_logic{n}{2*b}(k) = pr_resolves_0;
                
            else % bin remaining cycles
                
                % branch pr output 1
                pr_resolves_1 = abs(trapz(vid_interval, f_resolves_1_meta));
                err_index = error_bins_1(k);
                pmf(err_index) = pmf(err_index) + pr_resolves_1;
                
                % branch pr output 0
                pr_resolves_0 = abs(trapz(vid_interval, f_resolves_0_meta));
                err_index = error_bins_0(k);
                pmf(err_index) = pmf(err_index) + pr_resolves_0;
                
            end
        end
        
        
    end
    
end

% get "folded" (absolute value) code error
index0 = error_bin_func(0);
err_folded = err_range(index0:end);
pmf_folded = pmf(index0:end) + [0, fliplr(pmf(1:index0-1))];

% sanity checks:
% probability of MSB error
pr_MSB
pr_conventional_MSB = 2.* (2^1 - 1) * exp(-TADC/TAU)

% probability of exactly 1 LSB code error
P1LSB_analytic = sqrt(2/pi) .* (1 - 2^(-B)) .* STDCOMP ./ VLSB
P1LSB = pmf_folded(2)

%% get ideal noiseless SAR pmf

[ err_ideal, pmf_ideal, err_folded_ideal, pmf_folded_ideal ] = asar_meta_pmf_ideal(B, VRESET, TIMING);

%% plot PMF

FIGSIZE = [460 360];

% plot full pmf (pos + neg side)
figure; hold on
h(1) = stem(err_range, pmf, 'r', 'fill', 'MarkerSize', 4);
h(2) = stem(err_ideal, pmf_ideal, 'k', 'fill', 'MarkerSize', 4);

set(gcf, 'position', [200, 200, FIGSIZE]);
set(gca, 'yscale', 'log')
set(gca, 'ytick', 10.^[-150:10:0])
% set(gca, 'xtick', [-64, -32, 0, 32, 64])
set(gca, 'xtick', [-150:10:150])
% set(gca, 'xticklabel', [])
% set(gca, 'yticklabel', [])
set(gca, 'TickLength', [0.008 0.005])

% xlim([-2^(B-1) - 12, 2^(B-1) + 12])
xlim([-20, 20])
ylim([1e-100, 10^0])

% set(gca, 'ytick', 10.^[-40:4:0])
% set(gca, 'xtick', [-16:4:16])

xlabel('\epsilon')
ylabel('Pr(\epsilon|v_{id})')
legend(h, {'Noise', 'Ideal'}, 'Location', 'NorthEast')
grid on
