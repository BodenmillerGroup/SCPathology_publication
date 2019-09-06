
gates = retr('gates');
sessionData = retr('sessionData');


idx = find(contains(gates(:,1),'178_X4Y9'));%_X16Y7_130 %_X11Y6_78 %_X11Y1_6
image_num = idx

%Fragmentation score & communities, tumor cells only
mean_size = [];
failed = [];
label_out = [];
for image_num = 1:381



%     disp(image_num);    
    rows = sessionData(gates{image_num,2},2);
    selectedall_gates = image_num;
    
    pixelexpansion = 4;
    expansion_name = ['neighbour_',num2str(pixelexpansion),'_CellId'];
    
    neigb_index = cellfun(@(x) find(~cellfun('isempty',regexp(x,expansion_name))),...
    gates(selectedall_gates,3),'UniformOutput',false);
    
    Neighbor_Matrix = sessionData(gates{selectedall_gates,2},[neigb_index{1}]);
    
     PG_idx = 74; %For tumor only
    %PG_idx = 75; %For stromal - tumor 100
    Phenograph_Vector_orig = sessionData(gates{selectedall_gates,2},PG_idx);
    Phenograph_Vector = Phenograph_Vector_orig; %For stromal - tumor 100
    
    rows_orig = rows;
    curIDs = sessionData(gates{selectedall_gates,2},2);
    
    %Only tumor cell types
    idx_tumor_pheno = ismember(Phenograph_Vector_orig, 14:27);
    
    curIDs_tumor = curIDs(idx_tumor_pheno);
    
    Phenograph_Vector = Phenograph_Vector_orig(idx_tumor_pheno);
    Neighbor_Matrix = Neighbor_Matrix(idx_tumor_pheno,:);
    Neighbor_Matrix(~ismember(Neighbor_Matrix,curIDs_tumor)) = 0;
    
    rows = rows(idx_tumor_pheno);

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
    
        %Plotting out
    handles = gethand;
    idx_gates = selectedall_gates;
    global Mask_all
    mask = Mask_all(idx_gates).Image;
    mask(~ismember(mask,curIDs_tumor)) = 0; %curIDs_tumor
    

    if ~isempty(mask)


        centroids_cell = struct2cell(regionprops(mask,'Centroid'));
        centroids = cell2mat(centroids_cell');
        centroids(isnan(centroids)) = 0;
        

        G = graph(cutS);
%         figure;
        pg_ranks = centrality(G,'pagerank');
        closeness = centrality(G,'closeness');
        betweenness = centrality(G,'betweenness');
        eigenvector = centrality(G,'eigenvector');
        degree = centrality(G,'degree');

        p = plot(G,'XData',centroids(:,1),'YData',centroids(:,2),'MarkerSize',3,'LineWidth',2);
        set(gca,'Ydir','reverse');

%         catch
%             continue
% %         end
%         p.NodeCData = pg_ranks;
%          p.NodeCData = Phenograph_Vector;
%         colorbar

    end
        
try
    custom_gates = retr('custom_gatesfolder');	
    Graph2Binary(sparse(cutS),fullfile(custom_gates,'S'));	
    % Run Louvain on file for multiple iterations	
    niter = 20;
    if ispc == 1	
        [c,Q,labels,communities] = LouvainfromBin_Windows(fullfile(custom_gates,'S.bin'),niter);	
    else	
        [c,Q,labels,communities] = LouvainfromBin_ubuntu(fullfile(custom_gates,'S.bin'),niter,'Yes');	
    end	
    llim = max([ceil(length(labels)./1e4) 1]);	

    unLabel = unique(labels);	
    sizes = arrayfun(@(x) length(find(labels == x)) ,unLabel);
    idxgr1 = sizes > 1;
    curr_mean = mean(arrayfun(@(x) length(find(labels == x)) ,unLabel(idxgr1)));

    fragm_out = array2table(gates(image_num,1),'VariableNames',{'core'});
    fragm_out(:,'frag') = {curr_mean};
    mean_size = [mean_size; fragm_out];
    
    
    
    %return community labels above certain size	
    larger = arrayfun(@(x) length(find(labels == x)) > 200 ,unLabel);
    labels_use = unLabel(larger);
    
    currPheno = {};	
    ClusteringCoeff_perCommunity = {};	
    comm = {};
    core = {};
    id = {};
    for c=1:length(labels_use)
        currlab = labels == labels_use(c);	
        tempS = cutS;	
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
    out_put = array2table([cell2mat(comm'),cell2mat(currPheno'),cell2mat(ClusteringCoeff_perCommunity'),cell2mat(id')],'VariableNames',{'Community','Pheno','ClusteringCoef','CellId'});
    out_put(:,'core') = vertcat(core{:});
    label_out = [label_out;out_put];
    
catch
    failed = [failed, image_num];
    continue
end

end


%write out mean size of communities as indicator of tumor fragmentation
writetable(mean_size,'fragmentation.csv');

%Write out nodule data
writetable(label_out,'nodules_stromal_basel.csv');




custom_gates = retr('custom_gatesfolder');	
Graph2Binary(sparse(cutS),fullfile(custom_gates,'S'));	
% Run Louvain on file for multiple iterations	
niter = 20;
if ispc == 1	
    [c,Q,labels,communities] = LouvainfromBin_Windows(fullfile(custom_gates,'S.bin'),niter);	
else	
    [c,Q,labels,communities] = LouvainfromBin_ubuntu(fullfile(custom_gates,'S.bin'),niter,'Yes');	
end	
llim = max([ceil(length(labels)./1e4) 1]);	
    	
unLabel = unique(labels);	
larger = arrayfun(@(x) length(find(labels == x)) > 1,unLabel);
labels_use = unLabel(larger);

%for each community highlight it on graph, calc clustering coeff and	
%find involved celltypes	
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

    %     gplot(tempS,centroids,'-*');
    %     set(gca,'Ydir','reverse');	

        [acc, ~ ] = avgClusteringCoefficient(tempS);	
        ClusteringCoeff_perCommunity{c} = repmat(acc,sum(currlab),1);	
        currPheno{c} = Phenograph_Vector(currlab);	
        comm{c} = repmat(labels_use(c),sum(currlab),1);

end

out_put = [cell2mat(comm'),cell2mat(currPheno'),cell2mat(ClusteringCoeff_perCommunity')];
writetable(array2table(out_put,'VariableNames',{'Community','Pheno','ClusteringCoef'}),'communities.csv');






%Visualize example images
highlight = readtable('/home/jana/Desktop/R_dat/community_examples.csv');
cellid = table2array(highlight(:,'CellId'));
show = ismember(curIDs,cellid);
cutcutS = cutS;
cutcutS(~show,~show) = 0;


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




export_fig('filename.pdf','-transparent','-dpdf')
    
    
