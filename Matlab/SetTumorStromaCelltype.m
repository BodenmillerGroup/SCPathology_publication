%Set all tumor or all stromal cells to same cell type for neighbrohood
%analysis with specific question (for the analysis to get the TMEs
%irrespective of the tumor celltypes, set all Tumor cell types to a
%separate cluster number, that is not already used, here 100)

%Retrieve session data
ses = retr('sessionData');
gates = retr('gates');

%Hard coded!! Enter column of PhenoGraph clusters of interest/ celltype
%numbers to be adapted
tumor_pheno = ses(:,74);%Basel
tumor_pheno = ses(:,73);%Zuri

%Set all tumor cell types to label 100 (here cluster numbers 14:26)
tumor_pheno(ismember(tumor_pheno,14:27)) = 100; %1:13 for stroma instead of tumor celltypes Basel
tumor_pheno(ismember(tumor_pheno,[4,6,8:15,17,19:20,22:27,29:38,40])) = 100;%Zuri

%Add to session as new Clustering column
ses(:,end+1) = tumor_pheno;
et_gates = cellfun(@(x) [x(1:(end-1)), 'PhenographTumor100'],gates(:,3),'UniformOutput',false);
gates(:,3) = et_gates;
put('gates',gates);
put('sessionData',ses);