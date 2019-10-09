%Script to convert tissues to topological neighborhood graph representation
%and then run Louvain community detection on either all cells to detect
%microenvironment communities or only tumor cells for tumor communities.

%Retrieve session Data from histoCA
gates = retr('gates');
sessionData = retr('sessionData');

%Fragmentation score & communities on tumor cells only
%Loop through images
mean_size = [];
failed = [];
label_out = [];
for image_num = 1:381

    %Get single-cell data rows of current image in session Data
    rows = sessionData(gates{image_num,2},2);
    selectedall_gates = image_num;
    
    %Get neighbor column names
    pixelexpansion = 4;
    expansion_name = ['neighbour_',num2str(pixelexpansion),'_CellId'];
    
    %Get index of neighbor columns
    neigb_index = cellfun(@(x) find(~cellfun('isempty',regexp(x,expansion_name))),...
    gates(selectedall_gates,3),'UniformOutput',false);
    
    %Get neighbor matrix containing info of neighbors of each single-cell
    Neighbor_Matrix = sessionData(gates{selectedall_gates,2},[neigb_index{1}]);
    
    %Get column containing the PG of interest (relevant to distinguish
    %tumor from stromal cells)
    PG_idx = 74; 
    Phenograph_Vector_orig = sessionData(gates{selectedall_gates,2},PG_idx);
    Phenograph_Vector = Phenograph_Vector_orig;
    
    %Save for plotting
    rows_orig = rows;
    curIDs = sessionData(gates{selectedall_gates,2},2);
    
    %Use only tumor cell types
    idx_tumor_pheno = ismember(Phenograph_Vector_orig, 14:27);
    curIDs_tumor = curIDs(idx_tumor_pheno);
    Phenograph_Vector = Phenograph_Vector_orig(idx_tumor_pheno);
    Neighbor_Matrix = Neighbor_Matrix(idx_tumor_pheno,:);
    Neighbor_Matrix(~ismember(Neighbor_Matrix,curIDs_tumor)) = 0;
    rows = rows(idx_tumor_pheno);

    %Convert to topological neighborhood graph
    S = sparse(repmat(rows,1,size(Neighbor_Matrix,2)),(Neighbor_Matrix+1),1);
    fs = full(S);
    cutS = fs(:,2:size(fs,2));
    size1 = size(cutS,1);
    size2 = size(cutS,2);
    if size1 > size2
        cutS = [cutS, zeros(size1,size1-size2)];
    elseif size2 > size1
        cutS = [cutS; zeros(size1-size2, size2)];
    end
    
    %%%%%%%%%Plot graph, uncomment this part when looping over all images
    handles = gethand;
    idx_gates = selectedall_gates;
    global Mask_all
    mask = Mask_all(idx_gates).Image;
    mask(~ismember(mask,curIDs_tumor)) = 0;
    
    if ~isempty(mask)

        %Get single-cell centroids to plot graph nodes
        centroids_cell = struct2cell(regionprops(mask,'Centroid'));
        centroids = cell2mat(centroids_cell');
        centroids(isnan(centroids)) = 0;
        
        %Plot
        G = graph(cutS);
        figure;
        p = plot(G,'XData',centroids(:,1),'YData',centroids(:,2),'MarkerSize',3,'LineWidth',2);
        set(gca,'Ydir','reverse');
    end
    %%%%%%%%%
try
    %Run louvain community detection
    custom_gates = retr('custom_gatesfolder');	
    Graph2Binary(sparse(cutS),fullfile(custom_gates,'S'));	
    niter = 20;
    if ispc == 1	
        [c,Q,labels,communities] = LouvainfromBin_Windows(fullfile(custom_gates,'S.bin'),niter);	
    else	
        [c,Q,labels,communities] = LouvainfromBin_ubuntu(fullfile(custom_gates,'S.bin'),niter,'Yes');	
    end	
    llim = max([ceil(length(labels)./1e4) 1]);	

    %Get all detected community labels
    unLabel = unique(labels);	
    
    %Get sizes of the communities and find all communities that include
    %more than 1 cells and calculate average community size (for
    %fragmentation/cohesiveness score)
    sizes = arrayfun(@(x) length(find(labels == x)) ,unLabel);
    idxgr1 = sizes > 1;
    curr_mean = mean(arrayfun(@(x) length(find(labels == x)) ,unLabel(idxgr1)));

    fragm_out = array2table(gates(image_num,1),'VariableNames',{'core'});
    fragm_out(:,'frag') = {curr_mean};
    mean_size = [mean_size; fragm_out];
    
    %Return all communities above certain size for analysis in R pipeline
    larger = arrayfun(@(x) length(find(labels == x)) > 10 ,unLabel);
    labels_use = unLabel(larger);
    
    %Save Communities of each image and involved single-cells to dataset
    currPheno = {};	
    ClusteringCoeff_perCommunity = {};	
    comm = {};
    core = {};
    id = {};
    for c=1:length(labels_use)
        currlab = labels == labels_use(c);	
        tempS = cutS;	
        tempS(~currlab,~currlab) = 0;		
        [acc, ~ ] = avgClusteringCoefficient(tempS);	
        ClusteringCoeff_perCommunity{c} = repmat(acc,sum(currlab),1);	
        currPheno{c} = Phenograph_Vector_orig(currlab);	
        comm{c} = repmat(labels_use(c),sum(currlab),1);
        core{c} = repmat(gates(image_num,1),sum(currlab),1);
        id{c} = rows_orig(currlab);

    end
    out_put = array2table([cell2mat(comm'),cell2mat(currPheno'),cell2mat(ClusteringCoeff_perCommunity'),cell2mat(id')],'VariableNames',{'Community','Pheno','ClusteringCoef','CellId'});
    out_put(:,'core') = vertcat(core{:});
    label_out = [label_out;out_put];
    
catch
    failed = [failed, image_num];
    continue
end

end

%Write out mean size of communities as indicator of tumor fragmentation
writetable(mean_size,'fragmentation.csv');
%Write out community data for analysis in R
writetable(label_out,'nodules_stromal_basel.csv');


%Highlight tumor communities in different colors on graph
custom_gates = retr('custom_gatesfolder');	
Graph2Binary(sparse(cutS),fullfile(custom_gates,'S'));	
niter = 20;
if ispc == 1	
    [c,Q,labels,communities] = LouvainfromBin_Windows(fullfile(custom_gates,'S.bin'),niter);	
else	
    [c,Q,labels,communities] = LouvainfromBin_ubuntu(fullfile(custom_gates,'S.bin'),niter,'Yes');	
end	
llim = max([ceil(length(labels)./1e4) 1]);	
    	
%Visualize communities of at least a certain amount of nodes
unLabel = unique(labels);	
larger = arrayfun(@(x) length(find(labels == x)) > 10,unLabel);
labels_use = unLabel(larger);

%For each community highlight it on graph
currPheno = {};	
ClusteringCoeff_perCommunity = {};	
comm = {};
for c=1:length(labels_use)
    currlab = labels == labels_use(c);	
    tempS = cutS;	
    tempS(~currlab,~currlab) = 0;	
    tG = graph(tempS);
    tp = plot(tG,'XData',centroids(:,1),'YData',centroids(:,2),'MarkerSize',3,'LineWidth',3);
        set(gca,'Ydir','reverse');           
    hold on;

end



%Visualize example images picked in R pipeline
highlight = readtable('/home/jana/Desktop/R_dat/community_examples.csv');
cellid = table2array(highlight(:,'CellId'));
show = ismember(curIDs,cellid);
cutcutS = cutS;
cutcutS(~show,~show) = 0;

%Visualize individual communities of a given type from R pipeline
labels_use = unique(table2array(highlight(:,'Community')));
labels = table2array(highlight(:,'Community'));
currPheno = {};	
ClusteringCoeff_perCommunity = {};	
comm = {};
core = {};
id = {};
for c=1:length(labels_use)
    currlab = ismember(curIDs,cellid(labels == labels_use(c)));	
    tempS = cutcutS;	
    tempS(~currlab,~currlab) = 0;	

    tG = graph(tempS);
    tp = plot(tG,'XData',centroids(:,1),'YData',centroids(:,2),'MarkerSize',3,'LineWidth',3);
    set(gca,'Ydir','reverse');

    hold on;

    [acc, ~ ] = avgClusteringCoefficient(tempS);	
    ClusteringCoeff_perCommunity{c} = repmat(acc,sum(currlab),1);	
    currPheno{c} = Phenograph_Vector_orig(currlab);	
    comm{c} = repmat(labels_use(c),sum(currlab),1);
    core{c} = repmat(gates(image_num,1),sum(currlab),1);
    id{c} = rows_orig(currlab);

end

%Export visualized community graph to overlay onto label image
export_fig('/home/jana/Desktop/R_dat/filename.pdf','-transparent','-dpdf')



%Microenvironment communities on all cell types
%Loop through all images
mean_size = [];
failed = [];
label_out = [];
for image_num = 1:381

    %Get single-cell data rows of current image in session Data
    rows = sessionData(gates{image_num,2},2);
    selectedall_gates = image_num;
    
    %Get neighbor column names
    pixelexpansion = 4;
    expansion_name = ['neighbour_',num2str(pixelexpansion),'_CellId'];
    
    %Get index of neighbor columns
    neigb_index = cellfun(@(x) find(~cellfun('isempty',regexp(x,expansion_name))),...
    gates(selectedall_gates,3),'UniformOutput',false);
    
    %Get neighbor matrix containing info of neighbors of each single-cell
    Neighbor_Matrix = sessionData(gates{selectedall_gates,2},[neigb_index{1}]);
    
    %Get the cell type label column that has the same label for all tumor
    %cells and separate labels for stroma and immune cells    
    PG_idx = 75; 
    Phenograph_Vector_orig = sessionData(gates{selectedall_gates,2},PG_idx);
    Phenograph_Vector = Phenograph_Vector_orig; 
    rows_orig = rows;
    curIDs = sessionData(gates{selectedall_gates,2},2);
    
    %Convert to topological neighborhood graph
    S = sparse(repmat(rows,1,size(Neighbor_Matrix,2)),(Neighbor_Matrix+1),1);
    fs = full(S);
    cutS = fs(:,2:size(fs,2));
    size1 = size(cutS,1);
    size2 = size(cutS,2);
    if size1 > size2
        cutS = [cutS, zeros(size1,size1-size2)];
    elseif size2 > size1
        cutS = [cutS; zeros(size1-size2, size2)];
    end
    
    %%%%%%%%%Plot graph, uncomment this part when looping over all images
    handles = gethand;
    idx_gates = selectedall_gates;
    mask = Mask_all(idx_gates).Image;

    if ~isempty(mask)
        
        centroids_cell = struct2cell(regionprops(mask,'Centroid'));
        centroids = cell2mat(centroids_cell');
        centroids(isnan(centroids)) = 0;

        G = graph(cutS);
        figure;
        pg_ranks = centrality(G,'pagerank');
        closeness = centrality(G,'closeness');
        betweenness = centrality(G,'betweenness');
        eigenvector = centrality(G,'eigenvector');
        degree = centrality(G,'degree');

        p = plot(G,'XData',centroids(:,1),'YData',centroids(:,2),'MarkerSize',3,'LineWidth',2);
        set(gca,'Ydir','reverse');
    end
    %%%%%%%%%%%
try
    %Run louvain community detection
    custom_gates = retr('custom_gatesfolder');	
    Graph2Binary(sparse(cutS),fullfile(custom_gates,'S'));	
    niter = 20;
    if ispc == 1	
        [c,Q,labels,communities] = LouvainfromBin_Windows(fullfile(custom_gates,'S.bin'),niter);	
    else	
        [c,Q,labels,communities] = LouvainfromBin_ubuntu(fullfile(custom_gates,'S.bin'),niter,'Yes');	
    end	
    llim = max([ceil(length(labels)./1e4) 1]);	

    %Return communities including at least a minimum number of nodes
    unLabel = unique(labels);	
    larger = arrayfun(@(x) length(find(labels == x)) > 10 ,unLabel);
    labels_use = unLabel(larger);
    
    %Save Communities of each image and involved single-cells to dataset
    currPheno = {};	
    ClusteringCoeff_perCommunity = {};	
    comm = {};
    core = {};
    id = {};
    for c=1:length(labels_use)
        currlab = labels == labels_use(c);	
        tempS = cutS;	
        tempS(~currlab,~currlab) = 0;	
        [acc, ~ ] = avgClusteringCoefficient(tempS);	
        ClusteringCoeff_perCommunity{c} = repmat(acc,sum(currlab),1);	
        currPheno{c} = Phenograph_Vector_orig(currlab);	
        comm{c} = repmat(labels_use(c),sum(currlab),1);
        core{c} = repmat(gates(image_num,1),sum(currlab),1);
        id{c} = rows_orig(currlab);

    end
    out_put = array2table([cell2mat(comm'),cell2mat(currPheno'),cell2mat(ClusteringCoeff_perCommunity'),cell2mat(id')],'VariableNames',{'Community','Pheno','ClusteringCoef','CellId'});
    out_put(:,'core') = vertcat(core{:});
    label_out = [label_out;out_put];
    
catch
    failed = [failed, image_num];
    continue
end

end

%Write out community data for further analysis in R
writetable(label_out,'nodules_stromal_basel.csv');


%Highlight tumor communities in different colors on graph
custom_gates = retr('custom_gatesfolder');	
Graph2Binary(sparse(cutS),fullfile(custom_gates,'S'));	
niter = 20;
if ispc == 1	
    [c,Q,labels,communities] = LouvainfromBin_Windows(fullfile(custom_gates,'S.bin'),niter);	
else	
    [c,Q,labels,communities] = LouvainfromBin_ubuntu(fullfile(custom_gates,'S.bin'),niter,'Yes');	
end	
llim = max([ceil(length(labels)./1e4) 1]);	
    	
unLabel = unique(labels);	
larger = arrayfun(@(x) length(find(labels == x)) > 10,unLabel);
labels_use = unLabel(larger);

%Highlight each community on graph	
currPheno = {};	
ClusteringCoeff_perCommunity = {};	
comm = {};
for c=1:length(labels_use)
    currlab = labels == labels_use(c);	
    tempS = cutS;	
    tempS(~currlab,~currlab) = 0;	

    tG = graph(tempS);
    tp = plot(tG,'XData',centroids(:,1),'YData',centroids(:,2),'MarkerSize',3,'LineWidth',3);
    set(gca,'Ydir','reverse');
            
    hold on;

end

    
