function [ tiff_matrix ] = fuse_images(tabchild,imh)
% FUSE_IMAGES: Slightly adapted to our case to not display legend on top of
% image because we batch saved all images and didn't want parts to be
% hidden.

%Get GUI variables
handles = gethand;
tab_axes = retr('tab_axes1');

%Delete javawrapper classes
delete(tabchild.Children.findobj('Units','pixels'));

%Retrieve global variables
global Mask_all

%Function call to get the index and the tiff name of the selected channels
[sel_channels,tiff_matrix] = Comparetiffnames_tolistchannels(Mask_all);

%Store the colormap based on the number of selected channels
%in order: 'r','g','b','c','m','y'
colorstouse = [[1 0 0];[0 1 0];[0 0 1];[0 1 1];[1 0 1];[1 1 0]];

%If no axes found, create one
if isempty(tab_axes) == 1
    handles.panel_tiff_images;
    tab_axes = subplot(1,1,1,'Parent',tabchild);
    put('tab_axes1',tab_axes);
end

%Loop through the selected tiffs (channels)
for k=1:length(tiff_matrix{1,imh})
       
    %Scale image
    tiffimage_read = mat2gray(tiff_matrix{1,imh}{k});
    
    %Focus on axes
    handles.panel_tiff_images;
    axes(tab_axes);
    hold on;
    
    %If it is the first image, set background as the BWimage
    if k == 1
        blackim = imshow(tiffimage_read);
        set(blackim,'Tag','firstgrayimage');
        hold on;     
    end
    
    %Function call to convert image to RGB
    [rgb_Image] = make_rgb( tiffimage_read,colorstouse,k);
    hold on;
    
    %Display RGB image
    imagesh = imshow(rgb_Image);freezeColors;
    hold on;
    
    %Tag image
    set(imagesh,'Tag',strcat('rgbimage',int2str(k)));
    hold off;
    
    %Freeze colors
    freezeColors;
    
    %Adjust the intensity of the cell colors if multiple channels are
    %selected
    if length(tiff_matrix{1,imh}) ~= 1
        disp('Applying contrast to image to display all markers')
        intensemask =  imadjust(tiffimage_read);
    else
        intensemask =  tiffimage_read;
    end
 
    %Set the alphadata of the RGB image to the adjusted grayimage
    set(imagesh,'AlphaData',intensemask);
    freezeColors;
    
    hold off;
    
end

%If multiple channels are selected
if numel(sel_channels) > 1
    
    %Set up legend of which color correspond to which channel
%     string_channels = retr('list_channels');
%     La = line(ones(numel(sel_channels)),ones(numel(sel_channels)),'LineWidth',2,'Parent',tabchild.Children.findobj('Type','axes'));
%     set(La,{'color'},mat2cell(colorstouse(1:numel(sel_channels),:),ones(1,numel(sel_channels)),3)); freezeColors;
%     hla=legend(La,cellfun(@(n)(num2str(n)), string_channels(sel_channels), 'UniformOutput', false));
%     
%     %Define the location of the legend
%     set(hla, 'Location','South');
%     set(hla,'FontSize',8,'Interpreter','none');
    
end

end

