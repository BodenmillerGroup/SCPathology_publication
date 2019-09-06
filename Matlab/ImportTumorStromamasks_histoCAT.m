%Imports the Tumor-Stroma or other masks (in addition to single-cell masks)
%into the histoCAT session and calculates distance of every single cell to 
%the tumor-stroma boundary, to display distances or masks in images or
%extract them for plotting in R (Mask logic vector and single cell 
%distances vector will be added as additional channels in the session)

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
all_mask_files = all_mask_files (contains(all_mask_files ,'Basel'));

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

%Loop through all new masks and extract single cell distances to rim
row_store = [];
full_session = [];
not_working = [];
name_not_working  = {};
for k=1:length(temp_Mask_all)
    
    %Get other mask layers (additional to single-cell mask)
    other_Masks = temp_Mask_all(k).Image(2:end);
    singlecellmask = temp_Mask_all(k).Image{1};
    
    %Find single-cell rows in session corresponding to current mask image
    cur_image_name = unique_image_name(k);
    rows = [gates{contains(gates(:,1),cur_image_name),2}];
    row_store = [row_store,rows];
    cur_session = sessionData(rows,:);
    

     try   
        layer_names = {};
        sc_names = {};
        for layer = 1:length(other_Masks)

            store_dist = {};

            %Get single cells within each object of other layer
            get_single_cells_within = cellfun(@unique, struct2cell(regionprops(other_Masks{layer}, singlecellmask, 'PixelValues')),'UniformOutput',false);
            new_matrix = zeros(size(cur_session,1),1);

            for obj = 1:length(get_single_cells_within)
                curr_obj = get_single_cells_within(obj);
                rowfound_in = ismember(double(cur_session(:,2)), curr_obj{1});
                new_matrix(rowfound_in,1) = obj;
            end

            %Get all pixel distances first and then extract for regions of
            %interest
            all_pos = singlecellmask;
            all_pos(:) = 1;

            %Get distances of each pixel to edge of mask of current layer object from
            %inside
            pixel_dists = bwdistgeodesic(logical(all_pos),logical(other_Masks{layer}));

            %Get smallest pixeldist for each cell
            cell_dist = cellfun(@min, struct2cell(regionprops(singlecellmask, pixel_dists, 'PixelValues')),'UniformOutput',false);
            cell_dist = cell2mat(cell_dist)';
            cell_dist(isnan(cell_dist)) = 0;
            store_dist = cell_dist;

            %Get distances of each pixel to edge of closest mask object of current layer from outside mask
            pixel_dists_outside = bwdistgeodesic(logical(all_pos),logical(~other_Masks{layer}));
            cell_dist_out = cellfun(@min, struct2cell(regionprops(singlecellmask, pixel_dists_outside, 'PixelValues')),'UniformOutput',false);
            cell_dist_out = cell2mat(cell_dist_out)';
            cell_dist_out(isnan(cell_dist_out)) = 0;
            store_dist_out = cell_dist_out;

            %Overlay the distances for in and outside of mask into one column
            distances =sum(horzcat(store_dist,store_dist_out),2);                   
            if isempty(distances)
                distances = zeros(length(new_matrix),1);
            end

            %Add to fcs matrix
            cur_session = [cur_session, new_matrix,distances];

            %Add variable names
            current_name_dist = {'Mask_AllTumorFilled',strcat('Mask_AllTumorFilled','_distance_to_edge')};
            layer_names = [layer_names,current_name_dist];
        end

     catch
         %Record failed images/masks
         not_working = [not_working,k];
         name_not_working = [name_not_working,cur_image_name];
         new_matrix = zeros(size(cur_session,1),1);
         distances = zeros(length(new_matrix),1);
         cur_session = [cur_session, new_matrix,distances];

     end
    %Add current mask layer distances to session
    full_session = [full_session; cur_session];

end

%Add Full_session to sessionData in rows that correspond
sessionData(row_store,1:size(full_session,2)) = full_session;

%Same for gate names
et_gates = cellfun(@(x) [x, layer_names],gates(:,3),'UniformOutput',false);
gates(:,3) = et_gates;

%Update variables
Mask_all(logical(found_name_in_session)) = temp_Mask_all;
put('Tiff_all',Tiff_all);
put('Tiff_name',Tiff_name);
put ('sessionData',sessionData);
put('gates',gates);

list_samples_Callback;



