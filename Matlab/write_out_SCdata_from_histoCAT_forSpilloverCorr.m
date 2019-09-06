%Script to write out uncompensated single-cell data from histoCAT session 
%for spillover compensation with R scripts and then import the resulting 
%compensated single cell data again to use in histoCAT this only works for 
%our case where there were no other channels/gates in the histoCAT session
%but the ones used for compensation, otherwise adapt script.

%Export
ses = retr('sessionData');
gates = retr('gates');
%Amount of images hardcoded in case there are already custom channels
%bellow that shouldn't be used, otherwise just take full session
lens = cellfun(@(x) sum(contains(x,'Cell_')),gates(:,3),'UniformOutput',false);
%Should be same for all
max = unique(cell2mat(lens));
writetable(array2table(ses(:,1:(max+2)),'VariableNames',gates{1,3}(1:(max+2))),'curr_single_cell.csv');

%Import

%Read compensated data back in and replace uncompensated session data with
%compensated data (in order to have the compensated single-cell data for
%histoCAT analyses) 
comp = readtable('compensated_correct_basel.csv');
ses = retr('sessionData');
gates = retr('gates');
ses(:,1:(max+2)) = table2array(comp);
%Replace sessionData with compensated data
put('sessionData',ses);
%Replace table Fcs_Interest_all with compensated data
global Fcs_Interest_all;
temp = Fcs_Interest_all;
un = unique(table2array(comp(:,1)),'stable');
for i=1:length(temp)
    cur_un = un(i);
    temp{i}(:,1:(max+2)) = comp(table2array(comp(:,1)) == cur_un,:);  
end
Fcs_Interest_all = temp;