function [] = Asymmetric_heatmap(~,~,...
    Matrix_high,Matrix_low,Unique_high_all,Unique_low_all,pheno_name,...
    pixelexpansion,permutations,custom_gatesfolder,Extra_information,pVal_sig,cut_off_percent)

% ASYMMETRIC_HEATMAP Generates "asymmetric" heatmap for neighborhood results
%
% Input:
% Matrix_high, Matrix_low --> Matrix from left/right tailed permutation test
% 1 = significant; 0 = present not significant; -1 = not present
% Unique_high_all, Unique_low_all --> Unique interactions from right/left
% tailed permutation test
% pheno_name --> PhenoGraph name used
% pixelexpansion --> pixel expansion used
% permutations --> amount permutations used
% custom_gatesfolder --> location custom gates folder
% Extra_information --> name for the special cluster selected (annotation)
% pVal_sig --> p-value considered to be significant
% cut_off_percent --> how many images (ratio: 0-1) need to be represented
%
% Histology Topography Cytometry Analysis Toolbox (histoCAT)
% Denis Schapiro - Bodenmiller Group - UZH

set(gca,'TickLabelInterpreter','none');

% Calculate sum or mean of all present clusters using the significants
% output. -1 should be transformed to NaN
Matrix_high(Matrix_high == -1) = NaN;
Matrix_low(Matrix_low == -1) = NaN;



% Exclude clusters which are present in less than 2 patients and
% cut-off images
% We are using (cut-off) images at least
Clusters_to_exclude_high_nrImages = sum(isnan(Matrix_high))>(size(Matrix_high,1)-(size(Matrix_high,1)*0.000001));
Clusters_to_exclude_high = Clusters_to_exclude_high_nrImages;
% Clusters_to_exclude_high = or(Clusters_to_exclude_high_nrImages,exclude_based_on_nrPatients_high);
% Get mean using only present clusters
Sum_all_matrix_high = nanmean(Matrix_high,1);
% Remove clusters not present in (cut-off) images
Sum_all_matrix_high(Clusters_to_exclude_high) = 0;
Sum_all_for_transformation_high = [Unique_high_all Sum_all_matrix_high'];

n = max(Sum_all_for_transformation_high(:, 1));
m = max(Sum_all_for_transformation_high(:, 2));
SymmetricHeatMap_high = nan(n, m);
for i = 1:length(Sum_all_for_transformation_high)
    SymmetricHeatMap_high(Sum_all_for_transformation_high(i, 1), Sum_all_for_transformation_high(i, 2)) ...
        =Sum_all_for_transformation_high(i, 3);
end
% 
% 
% %Exclude clusters if not in a certain amount of patients
% all_nr_low = [];
% for cols = 1:size(Matrix_low,2)
%     cur_col = Matrix_low(:,cols);
%     non_nan = find(~isnan(cur_col));
%     found_patient = cellfun(@(x) find(ismember(non_nan,x)),patients_rows,'UniformOutput',false);
%     patient_not_empty = ~cellfun(@isempty, found_patient);
%     nr_patients = sum(patient_not_empty);
%     all_nr_low = [all_nr_low, nr_patients];  
% end
% exclude_based_on_nrPatients_low = all_nr_low < 2;

%Matrix_low = Matrix_low(grade_order == 3,:);

% Exclude clusters which are only present in a few images
% We are using (cut-off) images at least
Clusters_to_exclude_low_nrImages = sum(isnan(Matrix_low))>(size(Matrix_low,1)-(size(Matrix_low,1)*0.000001));
Clusters_to_exclude_low = Clusters_to_exclude_low_nrImages;
% Clusters_to_exclude_low = or(Clusters_to_exclude_low_nrImages,exclude_based_on_nrPatients_low);

Sum_all_matrix_low = nanmean(Matrix_low,1);
% Remove clusters not present in (cut-off) images
Sum_all_matrix_low(Clusters_to_exclude_low) = 0;



Sum_all_for_transformation_low = [Unique_low_all Sum_all_matrix_low'];
n = max(Sum_all_for_transformation_low(:, 1));
m = max(Sum_all_for_transformation_low(:, 2));
SymmetricHeatMap_low = nan(n, m);
for i = 1:length(Sum_all_for_transformation_low)
    SymmetricHeatMap_low(Sum_all_for_transformation_low(i, 1), Sum_all_for_transformation_low(i, 2)) ...
        =Sum_all_for_transformation_low(i, 3);
end

% For correlation, neighborhood overlay -> dot size
%Tot number images in which interaction is significant low
temp = Matrix_low;
temp(isnan(Matrix_low)) = 0;
tot_number_images_low = [Unique_low_all sum(temp)'];
tot_number_images_low(Clusters_to_exclude_low_nrImages,3) = 0;

temp2 = Matrix_high;
temp2(isnan(Matrix_high)) = 0;
tot_number_images_high = [Unique_high_all sum(temp2)'];
tot_number_images_high(Clusters_to_exclude_high_nrImages,3) = 0;

neg_or_pos = [Sum_all_for_transformation_high(:,1:2),Sum_all_for_transformation_high(:,3) - Sum_all_for_transformation_low(:,3)];
neg = neg_or_pos(:,3) < 0;
pos = neg_or_pos(:,3) > 0;

diff_nr_images = [Unique_high_all zeros(size(Unique_high_all,1),1)];
diff_nr_images(neg,3) = tot_number_images_low(neg,3);
diff_nr_images(pos,3) = tot_number_images_high(pos,3);

% cond1 = ~ismember(diff_nr_images(:,1),32:40);
% cond2 = ~ismember(diff_nr_images(:,2),32:40);
% cond = cond1 & cond2;
nr_images = diff_nr_images %(cond,:);
% all_possible = nchoosek(1:22,2);
% not_there = all_possible(~ismember(all_possible,nr_images(:,1:2),'rows'),:);
% nr_images = [nr_images; [not_there, zeros(size(not_there,1),1)]];
% nr_images = [nr_images; [[not_there(:,2),not_there(:,1)], zeros(size(not_there,1),1)]];
% nr_images = sortrows(nr_images);
writetable(array2table(nr_images),'nr_images.csv');

figure()

SymmetricHeatMap_high(isnan(SymmetricHeatMap_high)) = 0;
SymmetricHeatMap_low(isnan(SymmetricHeatMap_low)) = 0;

Delta_allvsall = SymmetricHeatMap_high-SymmetricHeatMap_low;
x = 1:size(Delta_allvsall,2);
y = (1:size(Delta_allvsall,1))';

% %get rid of between 13 and 100
% x( :,14:99) = [];
% y( 14:99, : ) = [];
% Delta_allvsall( 14:99, : ) = [];
% Delta_allvsall( :, 14:99) = [];

%get rid of rows and columns that are all empty
x( :, ~any(Delta_allvsall,1) ) = [];
y( ~any(Delta_allvsall,2), : ) = [];
Delta_allvsall( ~any(Delta_allvsall,2), : ) = [];
Delta_allvsall( :, ~any(Delta_allvsall,1) ) = [];

imagesc(Delta_allvsall);
colormap(b2r(min(min(Delta_allvsall)),max(max(Delta_allvsall))));
writetable(array2table(Delta_allvsall),'neighborhood_heatmap.csv');
title(['Heatmap_Pixel',num2str(pixelexpansion),'_',Extra_information,'_Perm_',permutations,'_',pheno_name,'_p-value',num2str(pVal_sig)]);
% x = unique(Unique_all_string(:,1),'stable');
set(gca,'Xtick',1:length(x),'Ytick',1:length(y),'XtickLabel',x,'YtickLabel',y);


%custom_gatesfolder =  retr('custom_gatesfolder');
saveas(gca,[custom_gatesfolder,'/','Heatmap_Pixel',num2str(pixelexpansion),'_Perm_',num2str(permutations),'_',pheno_name,'_',Extra_information,'_p-value',num2str(pVal_sig),'.fig']);


end

