%Script to automatically select specific images based on different criteria
%in list boxes (so one doesn't have to manually search through the image
%names)

%Retrieve session data
ses = retr('sessionData');
gates = retr('gates');

%Import patient groups made in R
patient_groups = readtable('/home/jana/Desktop/R_dat/patient_groups_tumors_meta.csv');

%Import and clean metadata
basel_meta = readtable('/home/jana/Desktop/bb_volume_jana/Data/2018/Clinical_paper/SpilloverCorrected/R_dat/Metadata_final.csv');
names = gates(:,1);

split_name = cellfun(@(x) strsplit(x,'_'),names,'UniformOutput',false);
split_name = cellfun(@(x) strcat(x{1},'_',x{2},'_',x{8},'_',x{9},'_',x{10}), split_name,'UniformOutput',false);
split_name = regexprep(split_name,'000','');

%Match SC data to metadata
found_names = cellfun(@(x) contains(split_name,x),table2array(basel_meta(:,1)),'UniformOutput',false);
idx_store = {};
for i=1:length(found_names)
    curr = found_names{i};
    idx = find(curr);
    if ~isempty(idx)
        idx_store{idx} = i;
    else
        continue
    end
end
em = cellfun(@isempty ,idx_store);

%Label controls and nomral samples accordingly for our case
table_ordered_according_gates_order(~em,:) = basel_meta(cell2mat(idx_store'),:);
table_ordered_according_gates_order(em,'core') = {'control'};
table_ordered_according_gates_order(em,'patientcode') = {'control'};
table_ordered_according_gates_order(em,'grade') = {4};
table_ordered_according_gates_order(em,'diseasestatus') = {'control'};

normals = ismember(table2array(table_ordered_according_gates_order(:,'diseasestatus')),'non-tumor');

%Select certain images based on some criteria to be selected in gates box
idx = ismember(table2array(patient_groups(:,'patient_pheno')),[12]); %For example select patients from patient group 12
idx_patients = ismember(table_ordered_according_gates_order(:,'patientcode'),patient_groups(idx,'patientcode'));
idx_patients(normals) = 0;
set(handles.list_samples,'Value',find(idx_patients));

%Select certain images based on some criteria to be selected in visualize samples box
names = get(handles.list_visual,'String');
split_names = cellfun(@(x) strsplit(x,'_'),table2array(patient_groups(idx,'patientcode')),'UniformOutput',false);
patients = cellfun(@(x) x{3}, split_names,'UniformOutput',false);
idx_patients = cellfun(@(x) contains(names,x),patients,'UniformOutput',false)';
idx_patients = sum(cell2mat(idx_patients),2);
set(handles.list_visual,'Value',find(idx_patients));

%Select specific list of images if e.g. neighborhood analysis should only
%be run on a subset
names = get(handles.list_visual,'String');
names = names(2:end);
split_name = cellfun(@(x) strsplit(x,'_'),names,'UniformOutput',false);
split_name = cellfun(@(x) strcat(x{1},'_',x{2},'_',x{8},'_',x{9},'_',x{10}), split_name,'UniformOutput',false);
split_name = regexprep(split_name,'000','');

to_find = {'BaselTMA_SP41_086_X12Y4','BaselTMA_SP41_074_X14Y6','BaselTMA_SP42_366_X4Y9','BaselTMA_SP42_360_X2Y9','BaselTMA_SP43_387_X14Y7','BaselTMA_SP42_230_X15Y5_235'...
    ,'BaselTMA_SP42_357_X1Y9','BaselTMA_SP43_500_X1Y9','BaselTMA_SP42_267_X3Y7','BaselTMA_SP42_361_X3Y9','BaselTMA_SP41_094_X6Y7','BaselTMA_SP42_245_X14Y4','BaselTMA_SP43_490_X11Y5'...
    ,'BaselTMA_SP43_494_X3Y4','BaselTMA_SP43_431_X8Y6','BaselTMA_SP41_082_X5Y4','BaselTMA_SP42_237_X11Y6','BaselTMA_SP41_045_X15Y5','BaselTMA_SP42_193_X3Y5','BaselTMA_SP41_132_X13Y7'...
    ,'BaselTMA_SP41_042_X11Y2','BaselTMA_SP42_227_X9Y4','BaselTMA_SP41_057_X10Y3','BaselTMA_SP41_027_X10Y5','BaselTMA_SP42_215_X5Y4','BaselTMA_SP41_061_X5Y6','BaselTMA_SP42_207_X9Y5'...
    ,'BaselTMA_SP41_083_X3Y7','BaselTMA_SP42_249_X5Y6','BaselTMA_SP41_065_X8Y6','BaselTMA_SP41_043_X13Y2','BaselTMA_SP41_085_X10Y4_244','BaselTMA_SP43_450_X11Y6','BaselTMA_SP41_054_X5Y3'...
    ,'BaselTMA_SP43_410_X2Y5','BaselTMA_SP41_138_X16Y7','BaselTMA_SP41_171_X2Y9','BaselTMA_SP41_069_X11Y6','BaselTMA_SP42_181_X1Y3','BaselTMA_SP42_214_X4Y5','BaselTMA_SP42_343_X11Y8'...
    ,'BaselTMA_SP42_278_X8Y7','BaselTMA_SP42_325_X6Y8'};

found = sum(cell2mat(cellfun(@(x) contains(split_name,x),to_find,'UniformOutput',false)),2);
set(handles.list_visual,'Value',(find(found)+1));

%Get met and control to select everything but them for neighborhood
%analysis if necessary
control_idx = [contains(table2array(table_ordered_according_gates_order(:,'core')),'control'); [1;1;1]]; %Last three controls not in meta

mets = {'Ay12x4','Ay15x6','By1x1','By1x7','By3x7','By4x5','By10x3','By12x3','Cy1x4','Cy11x5','Cy13x5','Cy14x8'};
met_idx = [sum(cell2mat(cellfun(@(x) contains(table2array(table_ordered_according_gates_order(:,'core')),x),mets,'UniformOutput',false)),2); [0;0;0]];

exclude = or(control_idx,met_idx);
set(handles.list_samples,'Value',find(~exclude));

