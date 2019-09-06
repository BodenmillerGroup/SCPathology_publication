function heatmap_images_overlay( labelimg, labelvec, valuevec, axis ,handles )
% HEATMAP_IMAGES_OVERLAY: Slightly adapted to our case to not display color
%bar on top of image (inconvenient when batch savin).

%Initialize amount of colors
ncols =100;

%Retrieve current percentile cut-off slider value
global currentsliderValue;
perc = currentsliderValue;

%If a percentile cut-off has been set
if ~(isempty(perc) == 1)
    
    %Function call to cut off the values above a given percentile
    %(outliers)
    valuevec = percentile_cutoff(labelvec, valuevec, handles, perc);
    
end

%Sort the label vector and the value vector according to the same order
[labelvec, ord] = sort(labelvec);
valuevec = valuevec(ord);

%Normalize the value vector
maxval = max(valuevec);
minval = min(valuevec);
res_valuevec = (valuevec-minval)/(maxval-minval);

%If there are no values found, return
if isnan(res_valuevec) == 1
    disp('No Data found to visualize');
    return;
end

%Make a 'full vector' in case some cell labels were missing
full_valuevec = zeros(max(labelimg(:)),1);
full_valuevec(labelvec) = res_valuevec;
full_valuevec = (full_valuevec-min(full_valuevec(:))) ./ (max(full_valuevec(:)-min(full_valuevec(:))));

%Define the color map
colmap = jet(ncols+1);

%Assign the colors to the values
full_valuevec = round(full_valuevec*ncols)+1;
colmap_lab = colmap(full_valuevec,:);

%Remove labels that are not in the labelvector from the image mask
labelimg(~ismember(labelimg, labelvec)) = 0;

%Apply the colormap
rgb_img = label2rgb(labelimg, colmap_lab, [0,0,0]);

%Set focus on axis and hold on to it
axes(axis);
hold on;

%Display image
intenseim = imshow(rgb_img);
hold on;

%Remove colorbar because when batch saving we didn't want it to hide parts
%of the image.

% %Set colorbar
% colormap(axis,colmap);
% cbr = colorbar(axis);
% cbr.Location = 'SouthOutside';
% hold on;

%Set labels, lims and ticks
drawnow;
% lims = get(cbr,'Limits');
% yval = linspace(lims(1), lims(2),11);
% set(cbr,'ytick',yval);
ylab=linspace(minval,maxval,11);
ylab =round(ylab, 2, 'significant');
%set(cbr,'YTickLabel',ylab);
freezeColors;

%Set position of colorbar
%cbr.Position = [0.1542 0.022 0.7274 0.0200];

%Tag the image
set(intenseim,'Tag','rgbimage1');

end

