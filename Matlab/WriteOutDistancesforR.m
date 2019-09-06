%Write  out single cell distance information from histoCAT to use in R for 
%plotting

%Get session info
ses = retr('sessionData');
gates = retr('gates');

%Convert image hash ID to original file name and put into first column of
%output table
names = gates(:,1);
names = names(:);
rows = [gates{:,2}];
un = unique(ses(rows,1),'stable');
for i=1:length(un)
    cur_un = un(i);
    cur_rows = ses(rows,1) == cur_un;
    out_table(cur_rows,1) = names(i);
end

%Add the cell IDs from the session to the second column of the output table
out_table(:,2) = num2cell(ses(rows,2));

%Hardcoded!! Adapt to column numbers of distances and tumor region mask and
%add those to third and fourth column of output table
out_table(:,3) = num2cell((ses(rows,74)));
out_table(:,4) = num2cell((ses(rows,73)));

%Write out as csv to read in R for plotting
t = array2table(out_table,'VariableNames',{'core','CellId','Distances','Mask'});
writetable(t,'Tumor_mask.csv');