%Script to run tsne without GUI and adapted to subsample cells

%fast tsne is only implemented for 2 dims.
ndims = 2;

%retrieve gui data
handles = gethand;
selected_gates = get(handles.list_samples,'Value');
selected_channels = get(handles.list_channels,'Value');
gates = retr('gates');
sessionData = retr('sessionData');
gate_context = retr('gateContext');
custom_gatesfolder = retr('custom_gatesfolder');

%HARDCODED which channels to use, adapt if anything changes
rows = [gates{selected_gates,2}];
sel_channels = [9:20,23:25,28,31:33,42:44,46,49:51,53:54,60:62,65,72];
selectedset = ses(rows,sel_channels);

%Subsample 20% of cells from each image
comb = [sessionData(:,1:2), selectedset];
un_im = unique(comb(:,1),'stable');
store_sampled_id =[];
data = [];
for i = 1:length(un_im)
    cur_im = un_im(i);
    cur_dat = comb(comb(:,1) == cur_im,:);
    [sampled,idx] = datasample(cur_dat, ceil((size(cur_dat,1) * 0.2)),'Replace',false);
    store_sampled_id = [store_sampled_id;cur_dat(idx,1:2)];
    data = [data; sampled];
    
end

%Run tsne on subsampled data
map = fast_tsne(data, 110);

%Write out results
id_map = [array2table(store_sampled_id), map];
id_map = table2array(id_map);
names = gates(:,1);
store_names = {};
for i = 1:length(un_im)
    cur_rows = id_map(:,1) == un_im(i);
    store_names(cur_rows,1) = names(i);
    
end
t = [store_names,array2table(id_map(:,2:end))];
t.Properties.VariableNames = {'core','CellId','tsne1','tsne2'};

writetable(t,'tsne_combined_20perc.csv');
