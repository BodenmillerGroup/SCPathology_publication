function plot_additionalLayer_mask_Callback(hObject, eventdata, handles)
% PLOT_ADDITONALLAYER_MASK_CALLBACK: This function is executed upon checking/unchecking of 
% the 'Additional layer mask on/off' checkbox. It displayes the segmentation mask outlines
% (outlines of the additional masks such as tumor region) on top of the currently selected image tab.
% This function or the checkbox are not in  the publicly available version
% of histoCAT but can be added to it.

%Retrieve GUI and global variables
tabmaster_histonetiff = retr('tabmaster_histonetiff');
global Mask_all
global Sample_Set_arranged

%Get the current tab number
tabnum = find(tabmaster_histonetiff.Children == tabmaster_histonetiff.SelectedTab);

maskoutline_additional = retr('maskoutline_additional');
    
%If the checkbox is checked, show mask
if handles.additional_layer_mask_onoff.Value == 1
    
    if size(Mask_all(1).Image,2) > 2
        numMasks = 1:size(Mask_all(1).Image,2)-2;

        choices = cellfun(@(x) strcat('Layer_',num2str(x)), num2cell(numMasks),'UniformOutput',false);
        Selection = listdlg('PromptString','Select a Layer:','SelectionMode','single','ListString',choices);


    elseif size(Mask_all(1).Image,2) == 2
        Selection = 1;

    else
        errordlg('There is no additional mask in dataset');
        return;
    end

    %Split the filepaths and extract the sample name of all samples
    splitSamplename = cellfun(@(x) strsplit(x,fullfile('/')),Sample_Set_arranged,'UniformOutput',false);
    allcutnames = cellfun(@(x) x(end),splitSamplename);

    %Find the index of the sample that corresponds to the currently
    %visualized image
    idxfound_name = find(~cellfun('isempty',regexpi(allcutnames,tabmaster_histonetiff.SelectedTab.Title)));

    %Store the corresponding single-cell mask (each pixel of a
    %cell is marked with the corresponding cell number)
    lblImg_filled = Mask_all(idxfound_name).Image{1+Selection};

    %If there is no mask and hence no single-cell data, return
    if isempty(lblImg_filled) == 1
        return;
    end

    %Get only the outlines of the individual cells (not all the pixels of a
    %cell, but only the edges)
    lblImg=conv2(single(lblImg_filled),[0 -1 0; -1 4 -1;0 -1 0],'same')>0;

    %Set focus on current axes and hold on to it
    axes(tabmaster_histonetiff.SelectedTab.Children.findobj('Type','axes'));
    hold on;

    %Display the mask outline image on top of the current image
    %axes, and set the transparancy of the mask to a level such
    %that both the cell outlines and the background image are visible
    cmap = colormap;
    lblImg = gray2ind(lblImg,200);
    maskoutline_additional{tabnum} = imshow(lblImg);
    freezeColors;
    set(maskoutline_additional{tabnum},'AlphaData',0.4);
    put('maskoutline_additional',maskoutline_additional);

else
    maskoutline_additional{tabnum}.Visible = 'off';
end

try
    %Store the colormap if it has been generated
    if isempty(cmap) ~= 1
        put('cmap',cmap);
        tabmaster_histonetiff.SelectedTab.Children.findobj('type','colorbar');
    else
        colorbar(tabmaster_histonetiff.SelectedTab.Children.findobj('Type','axes'),'off');
    end
catch
    return;
end

%Set axes position
tabmaster_histonetiff.SelectedTab.Children.findobj('Type','axes').Position = [0 0 1.1 1];


end