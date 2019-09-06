%Retrieve sessionData
ses = retr('sessionData');
gates = retr('gates');
handles = gethand;

Tiff_all = retr('Tiff_all');
Tiff_name = retr('Tiff_name');

%Get full Image marker means
all_mat = [];
for i = 1:size(Tiff_all,1)

    cur_tiff = Tiff_all(i,:);
    all_mean_intensities = cellfun(@(x) sum(sum(x))/(size(x,1)*size(x,2)),cur_tiff,'UniformOutput',false);
    all_mat = [all_mat; cell2mat(all_mean_intensities)];
    
end

names = gates{1,3}(3:44);
intensity_table = array2table(all_mat,'VariableNames',names);
intensity_table = [intensity_table,gates(1:381,1)];
writetable(intensity_table,'ImageMeanMarkerIntensity.csv');




%Get tumor region marker means


%Retrieve necessary variables
global Mask_all;
masks = Mask_all;
handles = gethand;
Tiff_all = retr('Tiff_all');
Tiff_name = retr('Tiff_name');
Tiff_all = Tiff_all;
Tiff_name = Tiff_name;
gates = retr('gates');
sessionData = retr('sessionData');

%Select folder containing the tumor masks and get all files corresponding
%to the masks in question
CNNmaskfolder = uipickfiles('Prompt','Select folder containging CNN masks'); 
all_files = dir(char(CNNmaskfolder));
all_file_names = {all_files(:).name};
all_file_names = all_file_names(~cellfun(@(x) contains(x,'._'),all_file_names));
all_file_names = all_file_names(3:end);
all_mask_files = all_file_names(~cellfun('isempty',regexpi(all_file_names,'_mask_AllTumorFilled.tiff')));
all_mask_files = all_mask_files (contains(all_mask_files ,'BaselTMA'));

%Store masks in temporary variables
temp_Mask_all = struct('Image',[]);
image_names = cellfun(@(x) extractBefore(x,'a0') ,all_mask_files,'UniformOutput',false);
unique_image_name = unique(image_names);
store_preds = [];

%Loop through all images and find the corresponding masks based on the
%image name
for i=1:length(unique_image_name)
    
    cur_image = unique_image_name(i);
    all_masks_fo_curr_image_names = all_mask_files(contains(all_mask_files,cur_image));
    
    %Read in masks
    MultiMasks = cellfun(@(x) imread(fullfile(char(CNNmaskfolder),x)), all_masks_fo_curr_image_names,'UniformOutput',false);
    names_keep = all_masks_fo_curr_image_names;
    
    %Don't read in single-cell maak again
    kick_out_normal_seg = cellfun(@(x) contains(x,'full_mask.tif'),names_keep);
    kick_out_test = cellfun(@(x) contains(x,'test'),names_keep);
    MultiMasks = MultiMasks(~(kick_out_normal_seg | kick_out_test));
    
    %Set 65535 to 1 incase necessary
    for j=1:length(MultiMasks)
        [a,b] = ismember(MultiMasks{j},65535);
        MultiMasks{j}(a) = 1;
    end
    
    %Find image in session corresponding to mask and add new masks to
    %single-cell mask
    corresp_mask = find(contains(gates(:,1),cur_image));
    MultiMasks_ext = [{masks(corresp_mask).Image},MultiMasks];
    temp_Mask_all(i,1,:).Image = MultiMasks_ext;
    
end

%Find all affected images in session (for which additional masks exist)
found_name_in_session = sum(cell2mat(cellfun(@(x) contains(gates(:,1),x) ,unique_image_name,'UniformOutput',false)),2);

%Loop through all new masks and extract tumor region mean marker expression
image_names = {};
mean_mat = [];
not_working = {};
for k=1:length(temp_Mask_all)
    
    %Get other mask layers (additional to single-cell mask)
    other_Masks = temp_Mask_all(k).Image(2:end);
    
    %Find index in tiffs
    cur_image_name = unique_image_name(k);
    tiff_idx = contains(gates(:,1),cur_image_name);
    
    %Get Tiffs corresponding to mask
    cur_marker_names = Tiff_name(tiff_idx,:);
    cur_tiffs = Tiff_all(tiff_idx,:);

    try
        %Extract marker mean of tumor regions for every channel
        marker_means = cellfun(@(x) cell2mat(struct2cell(regionprops(other_Masks{1},x, 'MeanIntensity'))),cur_tiffs,'UniformOutput',false);
        em = cellfun(@isempty, marker_means);
        marker_means(em) = {0};

        image_names = [image_names, cur_image_name];
        mean_mat = [mean_mat; cell2mat(marker_means)];
    catch
        
        not_working = [not_working,cur_image_name];
    end

end

names = gates{1,3}(3:54);%44 in Basel,59 uncleaned Zurich
intensity_table = array2table(mean_mat,'VariableNames',names);
intensity_table = [image_names',intensity_table];
writetable(intensity_table,'ZuriTumorHollowRegionMeanMarkerIntensity.csv');


%Loop through all new masks and extract stromal region area
image_names = {};
areas = [];
not_working = {};
for k=1:length(temp_Mask_all)
    
    %Get other mask layers (additional to single-cell mask)
    other_Masks = temp_Mask_all(k).Image(2:end);
    
    %Find index in tiffs
    cur_image_name = unique_image_name(k);


%     try
        %Extract marker mean of tumor regions for every channel
        stroma_area = sum(cell2mat(struct2cell(regionprops(~other_Masks{1},'Area'))));

        image_names = [image_names, cur_image_name];
        areas = [areas; stroma_area];
%     catch
%         
%         not_working = [not_working,cur_image_name];
%     end

end

area_table = array2table(image_names','VariableNames',{'core'});
area_table(:,'area') = array2table(areas);
writetable(area_table,'Basel_stroma_area.csv');

