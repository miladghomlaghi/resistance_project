clear all
close all
[total_MS_data] = readtable("Results.csv");

replicates = {[3,4,5,7,8]};

times = [0 0.5 5 24];


for j = 1:size(replicates,2)

    replicate_raw = total_MS_data(:,replicates{j});
    replicate_raw = replicate_raw(:,[4,5,1,2,2]);

    replicate = [repmat(1,size(replicate_raw,1),1),1./2.^(table2array(replicate_raw(1:end,[3:4 4])))];

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


                    bounceback10nm = [ bounceback10nm; [replicate_raw(i,1:2), num2cell(replicate(i,:))]];

                end
            elseif reduction_phase(2)==1 | reduction_phase(3)==1
                if  1.3*replicate(i,3) <= replicate(i,4)

                    bounceback10nm = [ bounceback10nm; [replicate_raw(i,1:2), num2cell(replicate(i,:))]];

                end
            end

        end
        replicate_rebound{j} = bounceback10nm;
    end
end
% Common among three replicates
common_replicate_gene_names=intersect(intersect(replicate_rebound{1,1}(:,1),replicate_rebound{1,2}(:,1)),replicate_rebound{1,3}(:,1));

common_replicate_gene_names = replicate_rebound{1,1}(:,1);
length(unique(common_replicate_gene_names))
final_data = total_MS_data(ismember(total_MS_data(:,1),common_replicate_gene_names),:);
writecell(final_data,"bounceback_src_shared_three.csv")

% shared two

common_replicate_gene_names_1 = intersect(replicate_rebound{1,1}(:,1),replicate_rebound{1,2}(:,1));
common_replicate_gene_names_2 = intersect(replicate_rebound{1,1}(:,1),replicate_rebound{1,3}(:,1));
common_replicate_gene_names_3 = intersect(replicate_rebound{1,2}(:,1),replicate_rebound{1,3}(:,1));

share_two = [common_replicate_gene_names_1; common_replicate_gene_names_2; common_replicate_gene_names_3];
final_data_shared_two = total_MS_data(ismember(total_MS_data(:,1),share_two),:);
writecell(final_data_shared_two,"bounceback_src_shared_two.csv")



unique_genes = total_MS_data(2:end,1);
Gene_names_no_site = {};
for i = 1:size(unique_genes,1)

    test = split(unique_genes{i,1},"_");
    Gene_names_no_site{i} =test{2};

end

unique_genes_1 = unique(Gene_names_no_site);

%
plot_BB_10nm = table2array(bounceback10nm(:,3:end));
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
% legend(BB_protein{:,:},'Location','north','Orientation','horizontal');
% lgd = legend;
% lgd.NumColumns = 4;
% lgd.FontSize = 8;
% legend boxoff

















