clear all
close all
[total_MS_data0] = readtable("CDK46_clean.csv");

%% Data preparation

vars = 5:16;   % column indices

total_MS_data = total_MS_data0;
total_MS_data{:, vars}(total_MS_data{:, vars} == 0) = NaN;
idx = cellfun(@isempty, total_MS_data.GeneNames);
total_MS_data(idx, :) = [];
% Select columns 4 to 18 (15 columns)
cols = total_MS_data(:, vars);

% Convert table to numeric array
X = table2array(cols);

% Count how many values per row are NaN or zero
numBad = sum( isnan(X) | X == 0, 2 );

% Drop row if more than 12 out of 15 are NaN or zero
rowsToDrop = numBad > 12;

% Remove those rows
total_MS_data(rowsToDrop, :) = [];
total_MS_data.GeneNames = cellfun( ...
    @(x) strtok(x, ';'), ...
    total_MS_data.GeneNames, ...
    'UniformOutput', false);

gene  = string(total_MS_data.GeneNames);
aa    = string(total_MS_data.AminoAcid);
pos   = string(total_MS_data.PositionsWithinProteins);

% Identify rows where AA or Position is missing
isMissing = ismissing(pos);

% Preallocate
total_MS_data.Gene_phospho = strings(height(total_MS_data), 1);

% Full format
total_MS_data.Gene_phospho(~isMissing) = ...
    gene(~isMissing) + "_" + aa(~isMissing) + "_" + pos(~isMissing);

% Gene only format
total_MS_data = total_MS_data(~isMissing,:);


%% Step 2: Quantile Normalization
log_data = log2(table2array(total_MS_data(:,vars)));

total_MS_data_norm = total_MS_data;
% Apply Quantile Normalization
total_MS_data_norm(:,[5:16]) = array2table(quantileNormalization(log_data));


%%
replicates = {[17,5,8,14,11],[17, 6,9,15, 12],[17,7,10,16,13]};

times = [0 1 8 24];


for j = 1:size(replicates,2)

    replicate_raw1 = total_MS_data_norm(:,replicates{j});
    replicate_num = [2.^(table2array(replicate_raw1(1:end,[2:5])))];
    nan_indices = isnan(abs(replicate_num(:,1)));
    inf_indices = isinf(abs(replicate_num(:,1)));


    % first column
    replicate_raw_col_1 = replicate_raw1 (~(nan_indices | inf_indices),:);
    replicate = replicate_num (~(nan_indices | inf_indices),:);

    % rest of columns
    nan_count_per_row = sum(isnan(replicate_num(:,2:end)), 2);

    % Find rows with at least two NaN values
    rows_to_remove = nan_count_per_row >= 2;
    replicate = replicate_num(~(rows_to_remove|inf_indices|nan_indices), :);
    replicate_raw = replicate_raw1 (~(rows_to_remove|inf_indices|nan_indices),:);

    %% finding bounce back in MS data

    bounceback10nm = [];
    for i = 1:size(replicate,1)
        reduction_phase = [
        0.7*replicate(i,1)> replicate(i,2)
        0.7*replicate(i,1)> replicate(i,3)
        0.7*replicate(i,2)> replicate(i,3)
        ];

        if sum(reduction_phase)
            if reduction_phase(1)==1

                if      1.2*replicate(i,2) <= replicate(i,4) |...
                        1.2*replicate(i,2) <= replicate(i,3) |...
                        1.2*replicate(i,3) <= replicate(i,4)


                    bounceback10nm = [ bounceback10nm; [replicate_raw(i,1), num2cell(replicate(i,:))]];

                end
            elseif reduction_phase(2)==1 | reduction_phase(3)==1
                if  1.3*replicate(i,3) <= replicate(i,4)

                    bounceback10nm = [ bounceback10nm; [replicate_raw(i,1), num2cell(replicate(i,:))]];

                end
            end

        end
        replicate_rebound{j} = bounceback10nm;
    end
end
% Common among three replicates
common_replicate_gene_names=intersect(intersect(replicate_rebound{1,1}(:,1),replicate_rebound{1,2}(:,1)),replicate_rebound{1,3}(:,1));

final_data = total_MS_data_norm(ismember(total_MS_data_norm(:,17),common_replicate_gene_names),:);
writetable(final_data,"bounceback_CDK_shared_three.csv")

% shared two

common_replicate_gene_names_1 = intersect(replicate_rebound{1,1}(:,1),replicate_rebound{1,2}(:,1));
common_replicate_gene_names_2 = intersect(replicate_rebound{1,1}(:,1),replicate_rebound{1,3}(:,1));
common_replicate_gene_names_3 = intersect(replicate_rebound{1,2}(:,1),replicate_rebound{1,3}(:,1));

share_two = unique([common_replicate_gene_names_1; common_replicate_gene_names_2; common_replicate_gene_names_3]);
final_data_shared_two = total_MS_data_norm(ismember(total_MS_data_norm(:,17),share_two),:);
writetable(final_data_shared_two,"bounceback_CDK_shared_two.csv")
writetable(total_MS_data_norm,"bounceback_CDK_total_valid.csv")



unique_genes = (total_MS_data(:,17));
% Gene_names_no_site = {};
% for i = 1:size(unique_genes,1)
% 
%     test = split(unique_genes{i,1},"_");
%     Gene_names_no_site{i} =test{1};
% 
% end

unique_genes_1 = unique(unique_genes);

unique_genes = final_data_shared_two(:,1);
% Gene_names_no_site = {};
% for i = 1:size(unique_genes,1)
% 
%     test = split(unique_genes{i,1},"_");
%     Gene_names_no_site{i} =test{1};
% 
% end

unique_genes_bounce = unique(unique_genes);

size(unique_genes_bounce)

%
plot_BB_10nm = table2array(final_data_shared_two(:,2:end));
plot_BB_10nm = 2.^plot_BB_10nm;
plot_BB_10nm = plot_BB_10nm./plot_BB_10nm(:,1);

range = max(plot_BB_10nm') - min(plot_BB_10nm');
a = (plot_BB_10nm' - min(plot_BB_10nm'))';
m01 = 2*a./ range' ;

% BB_protein = bounceback10nm(:,1);

figure('Position',[1175         658         859         289]);

hold on
% for c = 1:size(plot_BB_10nm,1)
plot(m01','linewidth',1);
% end
hold off
set(gca,'linewidth',1.5)
set(gca,'fontsize',12);







