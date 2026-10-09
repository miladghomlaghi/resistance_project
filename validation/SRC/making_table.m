clear all
close all
[MS_numbers MS_text total_MS_data] = xlsread('Lapatinib.xlsx');


table1= {'Gene','Phosphosite','Lapatinib 1uM  1min/Ctr','Lapatinib 1uM  5min/Ctr','Lapatinib 1uM  10min/Ctr','Lapatinib 1uM  60min/Ctr','Lapatinib 10uM  1min/Ctr','Lapatinib 10uM  5min/Ctr','Lapatinib 10uM  10min/Ctr','Lapatinib 10uM  60min/Ctr'};


%%%%%% Making the table containing MS data of YAP1: phosphorylationsites vs
%%%%%% treatments

for i = 2:size(total_MS_data,1)
    
    if ismember((total_MS_data(i,2)),(table1(2:end,1))) % checking protein
        
        if ismember(cell2mat(total_MS_data(i,3)),cell2mat(table1(2:end,2))) % checking phosphosite
            
            protein_index = find(strcmp(total_MS_data(i,2),table1(2:end,1)))+1;
            phosphosite_index = find(cell2mat(total_MS_data(i,3))==cell2mat(table1(2:end,2)))+1;
            targeted_index = intersect(phosphosite_index , protein_index);
            table1(targeted_index,find(strcmp(table1(1,2:end),total_MS_data(i,7)))+1) = mat2cell(2^cell2mat(total_MS_data(i,9)),1,1);
            
        else
            table1(end+1,1) = total_MS_data(i,2);
            table1(end,2) = total_MS_data(i,3);
            table1(end,find(strcmp(table1(1,2:end),total_MS_data(i,7)))+1) = mat2cell(2^cell2mat(total_MS_data(i,9)),1,1);
        end
        
    else % protein is not in the list
        table1(end+1,1) = total_MS_data(i,2);
        table1(end,2) = total_MS_data(i,3);
        table1(end,find(strcmp(table1(1,2:end),total_MS_data(i,7)))+1) = mat2cell(2^cell2mat(total_MS_data(i,9)),1,1);
        
        
    end
end


table1nm = table1(:,1:6);

tf1nm = cellfun('isempty',table1nm); % true for empty cells
table1nm(tf1nm) = {NaN};

table_number = cell2mat(table1nm(2:end,3:end));
nanIndices = any(isnan(table_number),2);
nanIndices = [0; nanIndices];
table1nm(logical(nanIndices),:) = [];
table_number_final_1nm = cell2mat(table1nm(2:end,3:end));


corrDist = pdist(table_number_final_1nm,'corr');
clusterTree = linkage(corrDist,'average');

clusters = cluster(clusterTree,'maxclust',4);
times = [1 5 10 60];
figure
for c = 1:4
    subplot(2,2,c);
    plot(times,table_number_final_1nm((clusters == c),:)');
    axis tight
end
subtitle('Hierarchical Clustering of Profiles');

%% clustering 1nm lapatinib
table10nm = table1(:,[1 2 7:10]);
tf10nm = cellfun('isempty',table10nm); % true for empty cells
table10nm(tf10nm) = {NaN};
table_number = [];
nanIndices = [];
table_number = cell2mat(table10nm(2:end,3:end));
nanIndices = any(isnan(table_number),2);
nanIndices = [0; nanIndices];
table10nm(logical(nanIndices),:) = [];
table_number_final_10nm = cell2mat(table10nm(2:end,3:end));
table_number_final_10nm_norm = normalize(table_number_final_10nm')';
corrDist = pdist(table_number_final_10nm_norm,'corr');
clusterTree = linkage(corrDist,'average');

clusters = cluster(clusterTree,'maxclust',10);
figure
for c = 1:10
    subplot(2,5,c);
    plot(times,table_number_final_10nm_norm((clusters == c),:)');
    axis tight
end
subtitle('Hierarchical Clustering of Profiles');


h = clustergram (table_number_final_10nm_norm,'Colormap',redbluecmap,'Standardize','Row');
% saveas(h,'cluster10nm.png')


%% finding bounce back in MS data

%-----10nm
bounceback10nm = [];
for i = 1:size(table_number_final_10nm,1)
    reduction_phase = [0.7*table_number_final_10nm(i,1)> table_number_final_10nm(i,2) 0.7*table_number_final_10nm(i,1)> table_number_final_10nm(i,3) 0.7*table_number_final_10nm(i,2)> table_number_final_10nm(i,3)];
    
    if sum(reduction_phase)
        if reduction_phase(1)==1
            if 1.2*table_number_final_10nm(i,3) <= table_number_final_10nm(i,4) | 1.2*table_number_final_10nm(i,2) <= table_number_final_10nm(i,4) | 1.2*table_number_final_10nm(i,2) <= table_number_final_10nm(i,3)
                bounceback10nm = [ bounceback10nm; table10nm(i+1,:)];
                
            end
            
        elseif reduction_phase(2)==1 | reduction_phase(3)==1
            if  1.2*table_number_final_10nm(i,3) <= table_number_final_10nm(i,4)
                bounceback10nm = [ bounceback10nm; table10nm(i+1,:)];
                
            end
        end
    end
end
plot_BB_10nm = cell2mat(bounceback10nm(:,3:end));
plot_BB_10nm = plot_BB_10nm./plot_BB_10nm(:,1);

BB_protein = bounceback10nm(:,1);

figure('Position',[1175         658         859         289]);

hold on
for c = 1:size(plot_BB_10nm,1)
    plot(times,plot_BB_10nm,'linewidth',1);
end
hold off
set(gca,'linewidth',1.5)
set(gca,'fontsize',12);
legend(BB_protein{:,:},'Location','north','Orientation','horizontal');
lgd = legend;
lgd.NumColumns = 4;
lgd.FontSize = 8;
legend boxoff
saveas(gcf,'timecourse_BB.png')

%------1nm

bounceback1nm = [];
for i = 2:size(table_number_final_1nm,1)
    reduction_phase = [table_number_final_1nm(i,1)> 2*table_number_final_1nm(i,2) table_number_final_1nm(i,1)> 2*table_number_final_1nm(i,3) table_number_final_1nm(i,2)> 2*table_number_final_1nm(i,3)];
    
    if sum(reduction_phase)
        if reduction_phase(1)==1
            if 1.3*table_number_final_1nm(i,3) <= table_number_final_1nm(i,4) | 1.3*table_number_final_1nm(i,2) <= table_number_final_1nm(i,4) | 1.3*table_number_final_1nm(i,2) <= table_number_final_1nm(i,3)
                bounceback1nm = [ bounceback1nm; table1nm(i+1,:)];
                
            end
            
        elseif reduction_phase(2)==1 | reduction_phase(3)==1
            if  1.3*table_number_final_1nm(i,3) <= table_number_final_1nm(i,4)
                bounceback1nm = [ bounceback1nm; table1nm(i+1,:)];
                
            end
        end
    end
end






