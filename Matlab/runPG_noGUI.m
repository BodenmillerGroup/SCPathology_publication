%Run PG without using the GUI

%%%%%%%%%%%%%Basel

%Retrieve sessionData
ses = retr('sessionData');
gates = retr('gates');
handles = gethand;

%Select samples to use for PhenoGraph in GUI or enter the index of
%samples to be used
sel = get(handles.list_samples,'Value');
rows = [gates{sel,2}];

%Hardcoded to specific session
sel_channels = [3:18,20:21,23:25,33,35,37,39,40,42:46,48,57];
data_cut = ses(rows,sel_channels);

%Normalize data to the 99th percentile
percentile = 99;  
data_cut = mynormalize(data_cut, percentile);

%Run Phenograph
near_neighbors = 20;
[labels,~,~] = phenograph(data_cut, near_neighbors,'random_seed','Yes');

%Write out data as csv to save and use in R, first column are sample names,
%second cell IDs and third the PG cluster labels
names = gates(1:381,1);
un = unique(ses(rows,1),'stable');
for i=1:length(un)
    cur_un = un(i);
    cur_rows = ses(rows,1) == cur_un;
    out_table(cur_rows,1) = names(i);
end

out_table(:,2) = num2cell(ses(rows,2));
out_table(:,3) = num2cell(labels);
t = array2table(out_table,'VariableNames',{'core','CellId','PhenoGraph'});
writetable(t,'PhenoGraph_test.csv');




%%%%%%%%%%%%%%%%%Zurich


ses = retr('sessionData');
gates = retr('gates');
handles = gethand;

%Select samples to use for PhenoGraph in GUI or enter the index of
%samples to be used
sel = get(handles.list_samples,'Value');
rows = [gates{sel,2}];


sel_channels = [3:18,20:21,23:25,33,35,37,39,40,42:46,48,57];
data_cut = ses(rows,sel_channels);

%Normalize data to the 99th percentile
percentile = 99;  
data_cut = mynormalize(data_cut, percentile);

%Run Phenograph
near_neighbors = 30;
[labels,~,~] = phenograph(data_cut, near_neighbors,'random_seed','Yes');

%Write out data as csv to save and use in R, first column are sample names,
%second cell IDs and third the PG cluster labels
names = gates(sel,1);
un = unique(ses(rows,1),'stable');
for i=1:length(un)
    cur_un = un(i);
    cur_rows = ses(rows,1) == cur_un;
    out_table(cur_rows,1) = names(i);
end

out_table(:,2) = num2cell(ses(rows,2));
out_table(:,3) = num2cell(labels);
t = array2table(out_table,'VariableNames',{'core','CellId','PhenoGraph'});
writetable(t,'PhenoGraph_zurich_x1_k30_nomets.csv');



