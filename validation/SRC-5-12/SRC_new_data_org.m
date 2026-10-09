clear all
close all

%% Load data
% FIX: "NA" strings in the CSV are treated as missing so isnan()/filters work.
total_MS_data = readtable("Original_matrix.csv", "TreatAsMissing", "NA");

%% Step 1: log2 transform
data_cols = [2,3,4,11:19];                               % NT + SRC columns
log_data  = log2(table2array(total_MS_data(:, data_cols)));

%% Step 2: Quantile normalisation
total_MS_data_norm = total_MS_data;
total_MS_data_norm(:, data_cols) = array2table(quantileNormalization(log_data));

%% Replicate definitions
% Each row: [ProteinID, NT(t0), SRC30(0.5h), SRC5(5h), SRC24(24h)]
replicates = {[1,2,17,11,14], [1,3,18,12,15], [1,4,19,13,16]};
times      = [0 0.5 5 24];

rebound_names = cell(1, numel(replicates));   % qualifying ProteinIDs per replicate
valid_names   = cell(1, numel(replicates));   % all QC-passing ProteinIDs per replicate

for j = 1:numel(replicates)

    replicate_raw = total_MS_data_norm(:, replicates{j});
    replicate_num = 2.^(table2array(replicate_raw(:, 2:5)));    % back to linear scale

    % Remove rows where baseline (t0) is NaN or Inf
    bad_t0 = isnan(replicate_num(:,1)) | isinf(replicate_num(:,1));

    % Remove rows with >=2 NaNs across the three post-treatment timepoints
    many_nan = sum(isnan(replicate_num(:,2:end)), 2) >= 2;

    keep          = ~(bad_t0 | many_nan);
    replicate_raw = replicate_raw(keep, :);
    replicate     = replicate_num(keep, :);

    valid_names{j} = replicate_raw{:,1};

    %% Detect rebound ("bounce back")
    is_rebound = false(size(replicate,1), 1);
    for i = 1:size(replicate,1)
        reduction_phase = [
            0.7*replicate(i,1) > replicate(i,2)      % t0  -> 0.5h drop
            0.7*replicate(i,1) > replicate(i,3)      % t0  -> 5h  drop
            0.7*replicate(i,2) > replicate(i,3) ];   % 0.5h-> 5h  drop

        if any(reduction_phase)
            if reduction_phase(1) == 1
                if 1.2*replicate(i,2) <= replicate(i,4) || ...
                   1.2*replicate(i,2) <= replicate(i,3) || ...
                   1.2*replicate(i,3) <= replicate(i,4)
                    is_rebound(i) = true;
                end
            elseif reduction_phase(2) == 1 || reduction_phase(3) == 1
                if 1.3*replicate(i,3) <= replicate(i,4)
                    is_rebound(i) = true;
                end
            end
        end
    end

    rebound_names{j} = replicate_raw{is_rebound, 1};   % ProteinIDs only
end

%% Define shared sets by ProteinID
% Rebound in all three replicates
shared_three = intersect(intersect(rebound_names{1}, rebound_names{2}), rebound_names{3});

% Rebound in at least two of three replicates (union of pairwise intersections)
pair_12    = intersect(rebound_names{1}, rebound_names{2});
pair_13    = intersect(rebound_names{1}, rebound_names{3});
pair_23    = intersect(rebound_names{2}, rebound_names{3});
shared_two = unique([pair_12; pair_13; pair_23]);

% All sites that passed QC in at least one replicate (the measured universe)
total_valid = unique([valid_names{1}; valid_names{2}; valid_names{3}]);

%% Save results
% FIX: select from the full normalised table by name, so a site is NOT dropped
% just because it was filtered out of replicate 3, and ALL replicate columns
% are written (not only replicate 3's values).
out_cols  = [1, 2,3,4, 11,12,13, 14,15,16, 17,18,19];   % ProteinID + NT + SRC (normalised log2)
all_names = total_MS_data_norm{:,1};

final_shared_three = total_MS_data_norm(ismember(all_names, shared_three), out_cols);
final_shared_two   = total_MS_data_norm(ismember(all_names, shared_two),   out_cols);
final_total_valid  = total_MS_data_norm(ismember(all_names, total_valid),  out_cols);

writetable(final_shared_three, "bounceback_src_shared_three.csv")
writetable(final_shared_two,   "bounceback_src_shared_two.csv")
writetable(final_total_valid,  "bounceback_src_total_valid.csv")

%% Sanity check: how many rebound genes are in the full gene list
base_name = @(id) extractBefore([id '_'], '_');            % gene before first "_"
genes_all    = unique(cellfun(base_name, total_MS_data{:,1},        'UniformOutput', false));
genes_bounce = unique(cellfun(base_name, final_shared_two{:,1},     'UniformOutput', false));
fprintf('Rebound genes in full list: %d of %d\n', ...
        sum(ismember(genes_bounce, genes_all)), numel(genes_bounce));

%% Plot trajectories of the "shared two" sites
% Mean across the three replicates at each timepoint (geometric, done in log2 space).
sel = ismember(all_names, shared_two);
t0  = mean(table2array(total_MS_data_norm(sel, [2 3 4])),    2, 'omitnan');
t05 = mean(table2array(total_MS_data_norm(sel, [17 18 19])), 2, 'omitnan');
t5  = mean(table2array(total_MS_data_norm(sel, [11 12 13])), 2, 'omitnan');
t24 = mean(table2array(total_MS_data_norm(sel, [14 15 16])), 2, 'omitnan');

traj = 2.^[t0 t05 t5 t24];        % undo log2
traj = traj ./ traj(:,1);         % normalise to baseline

rng_ = max(traj,[],2) - min(traj,[],2);
m01  = 2*(traj - min(traj,[],2)) ./ rng_;
m01(rng_ == 0, :) = 1;            % guard flat rows (avoid 0/0)

figure('Position',[1175 658 859 289]);
plot(times, m01', 'linewidth', 1);
xlabel('Time (h)'); ylabel('Scaled intensity');
set(gca, 'linewidth', 1.5, 'fontsize', 12);