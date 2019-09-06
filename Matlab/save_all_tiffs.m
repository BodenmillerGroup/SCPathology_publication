
handles = gethand;
global Sample_Set_arranged;
global Fcs_Interest_all;
global Mask_all;
global HashID;

%Set gating variable to zero (previous gates)
gatedontiff = 0;
put('gatedontiff',gatedontiff);


for i= 1:381    

    set(handles.list_visual,'Value',i+1);
    handles = gethand;

    %Function call to get the imageIDs of the selected samples
    [imageids,gate_names, SGsof_imageids_open,sample_orderIDX ] = choosetiffs_overlay_Callback;
    
    %Function call to apply the channel colors to Image (RGBCMY) based on
    %user selection
    overlay_maskandchannels( Mask_all,Fcs_Interest_all,imageids,gate_names,SGsof_imageids_open,sample_orderIDX );
    
    %Function call to highlight the cells of interest(from selected gate)
    show_selected_area_onTiff( Sample_Set_arranged,HashID,Fcs_Interest_all,Mask_all );
    
    %Set the channels list selection (between 1 and 6 channels can be
    %visualized simultaneouslynin RGBCMY)
    set(handles.list_channels,'Min',1,'Max',6);
    
    %Reset channels list incase RGBCMY image was used before
    channels = retr('list_channels');
    set(handles.list_channels,'String',channels);
    put('valchannel',{});
    
    
    for i=1:size(all_tabs,1)

        %Open new figure
        currentfig=figure;

        %Get only the image in the current tab of iteraction
        tab = get(all_tabs(i,:).Parent);


        %Copy the tiff image to the figure
        copyobj(tab.Children,currentfig);
% 
%         if ~exist('filename','var')  
%             %Ask user for the folder to save it in
%             [filename,path]=uiputfile('*.jpeg','Save Tiff');
%         end

%         If no path found
%         if path==0
%             return;
%         end

        %Save the figure as png
        %saveas(currentfig,fullfile(path,strcat(tab.Title,'_',filename)));
        path = retr('custom_gatesfolder');
        filename = tab.Title;
        saveas(currentfig,fullfile(path,strcat(tab.Title,'_',filename,'.png')))

        %Close the current figure
        close(currentfig);
    end
    
    
end

