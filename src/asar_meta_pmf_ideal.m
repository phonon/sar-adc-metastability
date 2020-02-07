function [ err, err_pmf, err_folded, err_pmf_folded ] = asar_meta_pmf_ideal(B, VRESET, TIMING)
% generates error PMF for ideal noiseless async SAR

% fzero optimization options, see:
% https://www.mathworks.com/help/matlab/ref/fzero.html
% https://www.mathworks.com/help/matlab/math/setting-options.html
OPTIMIZER = optimset(@fzero);
OPTIMIZER.('TolX') = realmin;

%% model parameters

% bits
% B = 6;

% reset code
% VRESET = -2^(B-1);    % 0000...
% VRESET = 2^(B-1) - 1; % 1111...
% VRESET = 0;           % 1000...

% input probability distribution function
fvid = @(x) 1/2^B;

% time normalized by tau (sample timings)
% TAU = 4;              % comparator regeneration tau
% TADC = 1000 / TAU;   % total ADC period
% TS = 100 / TAU;       % sampling time
% TSAR = 100 / TAU;     % SAR loop delay
% TLATCH = 40 / TAU;    % time for comparator bit to propogate output

% get timing from inputs
TAU = TIMING(1);
TADC = TIMING(2);
TS = TIMING(3);
TSAR = TIMING(4);
TLATCH = TIMING(5);

% full scale and ref
VDD = 1;
VFS = 1;         % [V]
VLSB = VFS/2^B;

% probability X state fails to resolve (stays in reset code)
PR_X = 1;

%%

% VLSB index (half-integer k VLSB)
k_range = -2^(B-1) + 0.5 : 1 : 2^(B-1) - 0.5;
len_k_range = length(k_range);

% error domain and output pmf (to calculate)
err = -2^B : 1 : 2^B;
err_pmf = zeros(1, 2^(B+1) + 1);

% error voltage -> error probability vector index function
error_bin = @(k) 1 + 2^B + floor(k);

% discrete DAC and residual functions
vdac_k_func = @(B,n,k) floor(heaviside(n)) .* -2.^(B-n) .* (floor(k ./ 2.^(B-n)) + 0.5);
vres_k_func = @(B,n,k) k + floor(heaviside(n-1)) .* -2.^(B-n+1) .* (floor(k ./ 2.^(B-n+1)) + 0.5);

% discrete output functions
vout_ideal_k_func = @(k) floor(k);
vout_k_func = @(vreset,B,n,k) vreset + 2.^(B-n) .* (floor(k ./ 2.^(B-n)) - floor(vreset ./ 2.^(B-n)));

% discrete output and errors functions
verr_k_func = @(vreset,B,n,k) vreset + floor(heaviside(n)) * ( 2.^(B-n) .* (floor(k ./ 2.^(B-n)) - floor(vreset ./ 2.^(B-n))) ) - floor(k);

% prealloc arrays
vmeta_bit_p = zeros(B+1, length(k_range));   % comparison phase (X state)
vmeta_bit_m = zeros(B+1, length(k_range));
vmeta_cmp_p = zeros(B+1, length(k_range));   % bit out delay phase (0/1, no X)
vmeta_cmp_m = zeros(B+1, length(k_range));
fvals = zeros(B+1, length(k_range));
vres = zeros(B, length(k_range));
verr = zeros(B+1, length(k_range));
T_bit = zeros(1, B);
T_cmp = zeros(1, B);

% initialize vmeta edge cases at boundaries of full scale range
vmeta_bit_m(:,1) = 0.5;
vmeta_bit_p(:,end) = 0.5;
vmeta_cmp_m(:,1) = 0.5;
vmeta_cmp_p(:,end) = 0.5;

% n = 0 initial values
verr(1,:) = verr_k_func(VRESET, B, 0, k_range);

for n = 1:1:B
    
    % residual and error for this cycle
    vres(n,:) = vres_k_func(B, n, k_range);
    verr(n+1,:) = verr_k_func(VRESET, B, n, k_range);
    
    % error bins for this cycle
    error_bins_cycle_start = error_bin(verr(n,:));
    error_bins_cycle_end = error_bin(verr(n+1,:));

    % T_bit: X STATE: ADC ends during comparison
    % T_cmp: 0/1 STATE: finished comparison, waiting for next cycle
    T_bit_tau = TS + (n-1) * TSAR + n * log(VDD/VLSB) + TLATCH - TADC;
    T_cmp_tau = TS + n * TSAR + n * log(VDD/VLSB) - TADC;
    
    % save (for debugging purposes)
    T_bit(n) = T_bit_tau + TADC;
    T_cmp(n) = T_cmp_tau + TADC;
    
    for k_index = 1:1:length(k_range)

        % get residual for bits 1...n
        vres_n = vres(1:n, k_index);

        % ============================================================
        % phase 1: X STATE, ADC ends during comparison
        % ============================================================

        % ==============================
        % positive side
        % ==============================
        f = @(x) sum(log(abs(vres_n + 0.5 - x))) - T_bit_tau;
        
        % solve
        % first check if endpoint and middle are greater than zero
        if f(0.5) < 0 && f(realmin) < 0
            if k_index == len_k_range
                vmeta_bit_p(n+1, k_index) = 0;
            else
                vmeta_bit_p(n+1, k_index) = 0.5;
            end
        else
            try
                [vid, fval, exitflag] = fzero(f, [realmin, 0.5], OPTIMIZER);
                vmeta_bit_p(n+1, k_index) = vid;
                fvals(n+1, k_index) = fval;
            catch flag_err
                % failed because endpoints
            end
        end
        
        % ==============================
        % negative side
        % ==============================
        f = @(x) sum(log(abs(vres_n - 0.5 + x))) - T_bit_tau;
        
        % solve
        % first check if endpoint greater than zero
        if f(0.5) < 0 && f(realmin) < 0
            if k_index == 1
                vmeta_bit_m(n+1, k_index) = 0;
            else
                vmeta_bit_m(n+1, k_index) = 0.5;
            end
        else
            try
                [vid, fval, exitflag] = fzero(f, [realmin, 0.5], OPTIMIZER);
                vmeta_bit_m(n+1, k_index) = vid;
                fvals(n+1, k_index) = fval;
            catch flag_err
                % failed because endpoints
            end
        end
        
        % ==============================
        % convert edges to error probability
        % ==============================
        pr_p = abs(integral(fvid, vmeta_cmp_p(n, k_index), vmeta_bit_p(n+1, k_index), 'ArrayValued', true, 'AbsTol', realmin, 'RelTol', 1e-100));
        pr_m = abs(integral(fvid, vmeta_cmp_m(n, k_index), vmeta_bit_m(n+1, k_index), 'ArrayValued', true, 'AbsTol', realmin, 'RelTol', 1e-100));

        err_bin_x_fails = error_bins_cycle_start(k_index);
        err_pmf(err_bin_x_fails) = err_pmf(err_bin_x_fails) + PR_X * (pr_p + pr_m);
        
        err_bin_x_resolves = error_bins_cycle_end(k_index);
        err_pmf(err_bin_x_resolves) = err_pmf(err_bin_x_resolves) + (1 - PR_X) * (pr_p + pr_m);
        
        % ============================================================
        % phase 2: 0,1 ADC ends during delay to next stage
        %          -> only occur for n = 1..B-1
        % ============================================================
        if n < B
            % ==============================
            % positive side
            % ==============================
            f = @(x) sum(log(abs(vres_n + 0.5 - x))) - T_cmp_tau;

            % solve
            % first check if endpoint and middle are greater than zero
            if f(0.5) < 0 && f(realmin) < 0
                if k_index == len_k_range
                    vmeta_cmp_p(n+1, k_index) = 0;
                else
                    vmeta_cmp_p(n+1, k_index) = 0.5;
                end
            else
                try
                    [vid, fval, exitflag] = fzero(f, [realmin, 0.5], OPTIMIZER);
                    vmeta_cmp_p(n+1, k_index) = vid;
                    fvals(n+1, k_index) = fval;
                catch flag_err
                    % failed because endpoints
                end
            end

            % ==============================
            % negative side
            % ==============================
            f = @(x) sum(log(abs(vres_n - 0.5 + x))) - T_cmp_tau;

            % solve
            % first check if endpoint greater than zero
            if f(0.5) < 0 && f(realmin) < 0
                if k_index == 1
                    vmeta_cmp_m(n+1, k_index) = 0;
                else
                    vmeta_cmp_m(n+1, k_index) = 0.5;
                end
            else
                try
                    [vid, fval, exitflag] = fzero(f, [realmin, 0.5], OPTIMIZER);
                    vmeta_cmp_m(n+1, k_index) = vid;
                    fvals(n+1, k_index) = fval;
                catch flag_err
                    % failed because endpoints
                end
            end

            % ==============================
            % convert edges to error probability
            % ==============================
            pr_p = abs(integral(fvid, vmeta_bit_p(n+1, k_index), vmeta_cmp_p(n+1, k_index), 'ArrayValued', true, 'AbsTol', realmin, 'RelTol', 1e-100));
            pr_m = abs(integral(fvid, vmeta_bit_m(n+1, k_index), vmeta_cmp_m(n+1, k_index), 'ArrayValued', true, 'AbsTol', realmin, 'RelTol', 1e-100));

            err_bin = error_bins_cycle_end(k_index);
            err_pmf(err_bin) = err_pmf(err_bin) + pr_p + pr_m;
        end
    end
    
end

% give 0 error the rest
err_bin = error_bin(0);
err_pmf(err_bin) = 1.0 - (sum(err_pmf) - err_pmf(err_bin));

% get "folded" (absolute value) code error
index0 = error_bin(0);
err_folded = err(index0:end);
err_pmf_folded = err_pmf(index0:end) + [0, fliplr(err_pmf(1:index0-1))];

% vmeta_cmp_m
% vmeta_cmp_p

% optional: plot
% figure; hold on
% stem(err, err_pmf)
% set(gca, 'yscale', 'log')
% xlim([err(1), err(end)])
% xlabel('\epsilon')
% ylabel('Pr(\epsilon|v_{id})')
% set(gca, 'YTick', 10.^[-100:20:0])

end