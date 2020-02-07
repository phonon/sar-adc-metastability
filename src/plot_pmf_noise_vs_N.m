% Run and plot ASAR combined thermal noise model vs. timing margin N
% NOTE: this code is compute intensive and will take a while to run at >6 bits
% 
% By default uses parallelized version of combined noise PMF model
% which uses MATLAB parfor function.
% If MATLAB version does not support parfor, replace the function call
%     asar_meta_pmf_noise_combined_parallel
% with
%     asar_meta_pmf_noise_combined
clear all; close all; clc; format compact
[FONTSIZE, LINEWIDTH, FIGSIZE, SCATTER ] = figure_settings(14, 1.1, [600 300], 40);

%% parameters

% true - rerun the model simulation
% false - use saved data
RUN_SIM = true;

% results save folder
SAVE_FOLDER = './results'

% bits to simulate
B = 7;

% time normalized by tau (sample timings, all time in ps)
TAU = 0;              % comparator regeneration tau
TS = 0;         % sampling time
TLOGIC = [8, 16];     % SAR loop delay

% B*N = TADC - TFIXED -> TADC = TFIXED + B*N
N = 1:1:10;

% reset code
VRESET = -2^(B-1);    % 0000...
% VRESET = 2^(B-1) - 1; % 1111...
% VRESET = 0;           % 1000...

% std deviation of noise for comparator in terms of VFS
% 0.5 -> VFS/2
ENOB = 6
STDCOMP = sqrt(0.5*(0.5^2)*(10.0^(-(6.02*ENOB + 1.76)/10)) - (1.0/(2.0^B))^2/12)

% size of PMFs (for initializing arrays)
FOLDED_PMF_SIZE = 2^B + 1;
FULL_PMF_SIZE = 2^(B+1) + 1;
NUM_TIME_POINTS = length(N);

for k = 1:length(TLOGIC)

    TADC = B*N + TS + (B-1)*TLOGIC(k)

    % input format for functions (DO NOT CHANGE)
    TIMING = zeros(length(TADC), 5);
    TIMING(:,1) = TAU;
    TIMING(:,2) = TADC;
    TIMING(:,3) = TS;
    TIMING(:,4) = TLOGIC(k);
    TIMING(:,5) = 0;
    
    if RUN_SIM == false
        
        % errors/pmfs for current Tlogic (so we can save them separately
        % and re-generate as needed)
        all_err_ideal{k} = zeros(NUM_TIME_POINTS, FULL_PMF_SIZE);
        all_pmf_ideal{k} = zeros(NUM_TIME_POINTS, FULL_PMF_SIZE);
        all_err_noise{k} = zeros(NUM_TIME_POINTS, FULL_PMF_SIZE);
        all_pmf_noise{k} = zeros(NUM_TIME_POINTS, FULL_PMF_SIZE);
        
        all_err_folded_ideal{k} = zeros(NUM_TIME_POINTS, FOLDED_PMF_SIZE);
        all_pmf_folded_ideal{k} = zeros(NUM_TIME_POINTS, FOLDED_PMF_SIZE);
        all_err_folded_noise{k} = zeros(NUM_TIME_POINTS, FOLDED_PMF_SIZE);
        all_pmf_folded_noise{k} = zeros(NUM_TIME_POINTS, FOLDED_PMF_SIZE);
        
        for i = 1:length(TADC)

            % load saved data
            fname = [SAVE_FOLDER, '/asar_pmf_', num2str(B, 2) ,'b_ENOB_', num2str(ENOB, 3) ,'_Tlogic_', num2str(TLOGIC(k), 3), '_N_', num2str(N(i), 4), '.mat'];
            load(fname)

            % save pmfs
            all_err_ideal{k}(i,:) = err_ideal;
            all_pmf_ideal{k}(i,:) = pmf_ideal;
            all_err_noise{k}(i,:) = err_noise;
            all_pmf_noise{k}(i,:) = pmf_noise;

            all_err_folded_ideal{k}(i,:) = err_folded_ideal;
            all_pmf_folded_ideal{k}(i,:) = pmf_folded_ideal;
            all_err_folded_noise{k}(i,:) = err_folded_noise;
            all_pmf_folded_noise{k}(i,:) = pmf_folded_noise;
            
        end
        
    else % run simulation
        
        % errors/pmfs for current Tlogic (so we can save them separately
        % and re-generate as needed)
        all_err_ideal{k} = zeros(NUM_TIME_POINTS, FULL_PMF_SIZE);
        all_pmf_ideal{k} = zeros(NUM_TIME_POINTS, FULL_PMF_SIZE);
        all_err_noise{k} = zeros(NUM_TIME_POINTS, FULL_PMF_SIZE);
        all_pmf_noise{k} = zeros(NUM_TIME_POINTS, FULL_PMF_SIZE);
        
        all_err_folded_ideal{k} = zeros(NUM_TIME_POINTS, FOLDED_PMF_SIZE);
        all_pmf_folded_ideal{k} = zeros(NUM_TIME_POINTS, FOLDED_PMF_SIZE);
        all_err_folded_noise{k} = zeros(NUM_TIME_POINTS, FOLDED_PMF_SIZE);
        all_pmf_folded_noise{k} = zeros(NUM_TIME_POINTS, FOLDED_PMF_SIZE);

        % loop and calculate full pmf for all timings
        for i = 1:length(TADC)
            disp(['TSAR = ', num2str(TLOGIC(k), 3), ' | TADC = ', num2str(TADC(i), 3), ' tau'])

            % pmf ideal
            [err_ideal, pmf_ideal, err_folded_ideal, pmf_folded_ideal] = asar_meta_pmf_ideal(B, VRESET, TIMING(i,:));

            % pmf with comparator thermal noise) use parallelized version
            % by default
            [err_noise, pmf_noise, err_folded_noise, pmf_folded_noise] = asar_meta_pmf_noise_combined_parallel(B, VRESET, TIMING(i,:), STDCOMP);
            
            % non-parallelized version (uncomment if parfor does not work)
%             [err_noise, pmf_noise, err_folded_noise, pmf_folded_noise] = asar_meta_pmf_noise_combined(B, VRESET, TIMING(i,:), STDCOMP);

            fname = [SAVE_FOLDER, '/asar_pmf_', num2str(B, 2) ,'b_ENOB_', num2str(ENOB, 3) ,'_Tlogic_', num2str(TLOGIC(k), 3), '_N_', num2str(N(i), 4), '.mat'];
            save(fname, 'err_ideal', 'pmf_ideal', 'err_noise', 'pmf_noise', 'err_folded_ideal', 'pmf_folded_ideal', 'err_folded_noise', 'pmf_folded_noise');
            
            all_err_ideal{k}(i,:) = err_ideal;
            all_pmf_ideal{k}(i,:) = pmf_ideal;
            all_err_noise{k}(i,:) = err_noise;
            all_pmf_noise{k}(i,:) = pmf_noise;

            all_err_folded_ideal{k}(i,:) = err_folded_ideal;
            all_pmf_folded_ideal{k}(i,:) = pmf_folded_ideal;
            all_err_folded_noise{k}(i,:) = err_folded_noise;
            all_pmf_folded_noise{k}(i,:) = pmf_folded_noise;
        
        end
        
    end
    
end

%% plot err > X vs  N

for k = 1:length(TLOGIC)
    
TADC = N + TS + (B-1)*TLOGIC(k);
index = TADC > 0;

pmf_err_greater_than_ideal = zeros(1, FOLDED_PMF_SIZE);
pmf_err_greater_than_noise = zeros(1, FOLDED_PMF_SIZE);

for i = 1:length(N);
    for n = 1:FOLDED_PMF_SIZE
        pmf_err_greater_than_ideal(i,n) = sum(all_pmf_folded_ideal{k}(i,n:end));
        pmf_err_greater_than_noise(i,n) = sum(all_pmf_folded_noise{k}(i,n:end));
    end
end


figure; hold on
for i = 1:8
    plot(N(index), pmf_err_greater_than_noise(index,i), 'r')
    plot(N(index), pmf_err_greater_than_ideal(index,i), 'k-')
end

grid on

set(gcf, 'Position', [800 100 FIGSIZE])
set(gca, 'TickLength', [0.015 0.005])
set(gca, 'YScale', 'log')
set(gca, 'YTick', 10.^[-80:4:0])
set(gca, 'XTick', -10:1:10)
xlim( [-1 8] )
ylim( [1e-20, 1e0])

xlabel('N')
ylabel('Pr(\epsilon \geq X)')
% set(gca, 'xticklabel', [])
% set(gca, 'yticklabel', [])

end