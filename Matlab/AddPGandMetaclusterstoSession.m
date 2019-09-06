%Import phenograph (or any other clustering) info into session if not 
%already there and add metaclusters for our case too (only necessary when 
%a clustering was run or save somewhere else and needs to be added to 
%current histoCAT session). Columns are added to session.

%Retrieve session info
ses = retr('sessionData');
gates = retr('gates');

%Import PG info (this script assumes that the images and single-cells are
%already in the same order for both the session and the PG fcs file, if
%this is not the case match cells first based on Image and CellID);
phenoBasel = readtable('/home/jana/Desktop/R_dat/PhenoGraphBasel_k20_2.csv','Delimiter',',');
phenoBasel.Properties.VariableNames = {'core','CellId','cluster'};
ses(:,end+1) = table2array(phenoBasel(:,'cluster'));
et_gates = cellfun(@(x) [x, 'PhenographBasel'],gates(:,3),'UniformOutput',false);
gates(:,3) = et_gates;
put('gates',gates);
put('sessionData',ses);

%Metaclusters, according to the way we grouped the small clusters,
%HARDCODED!!!
pheno_col = ses(:,end);
meta_pheno(pheno_col == 25) = 1;
meta_pheno(pheno_col == 19) = 2;
meta_pheno(pheno_col == 2) = 3;
meta_pheno(pheno_col == 6) = 4;
meta_pheno(pheno_col == 38) = 5;
meta_pheno(pheno_col == 70) = 6;
meta_pheno(pheno_col == 10) = 7;
meta_pheno(pheno_col == 3) = 8;
meta_pheno(pheno_col == 4) = 9;
meta_pheno(pheno_col == 1) = 10;
meta_pheno(pheno_col == 15) = 11;
meta_pheno(pheno_col == 71) = 12;
meta_pheno(pheno_col == 36) = 13;
meta_pheno(pheno_col == 11) = 14;
meta_pheno(pheno_col == 66 | pheno_col == 43 | pheno_col == 32 | pheno_col == 49 | pheno_col == 56) = 15;
meta_pheno(pheno_col == 35 | pheno_col == 68 | pheno_col == 8 | pheno_col == 60) = 16;
meta_pheno(pheno_col == 23 | pheno_col == 33 | pheno_col == 26) = 17;
meta_pheno(pheno_col == 47 | pheno_col == 57 | pheno_col == 17 | pheno_col == 50) = 18;
meta_pheno(pheno_col == 16 | pheno_col == 69 | pheno_col == 67) = 27;
meta_pheno(pheno_col == 52 | pheno_col == 28 | pheno_col == 64 | pheno_col == 45) = 19;
meta_pheno(pheno_col == 55 | pheno_col == 13 | pheno_col == 40 | pheno_col == 51 | pheno_col == 42) = 20;
meta_pheno(pheno_col == 21 | pheno_col == 65 | pheno_col == 7) = 21;
meta_pheno(pheno_col == 5 | pheno_col == 41 | pheno_col == 30 | pheno_col == 27 | pheno_col == 58) = 22;
meta_pheno(pheno_col == 59 | pheno_col == 61 | pheno_col == 14 | pheno_col == 48 | pheno_col == 62 | pheno_col == 20 | pheno_col == 24) = 23;
meta_pheno(pheno_col == 18 | pheno_col == 46 | pheno_col == 53 | pheno_col == 37 | pheno_col == 31) = 24;
meta_pheno(pheno_col == 39 | pheno_col == 9 | pheno_col == 12 | pheno_col == 29 | pheno_col == 34 | pheno_col == 22) = 25;
meta_pheno(pheno_col == 44 | pheno_col == 54 | pheno_col == 63) = 26;

%Store in session
ses(:,end+1) = meta_pheno';
et_gates = cellfun(@(x) [x, 'PhenographMetaclusters'],gates(:,3),'UniformOutput',false);
gates(:,3) = et_gates;
put('gates',gates);
put('sessionData',ses);


