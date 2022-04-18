libpath = getenv('G2S3LIB')
datapath = getenv('MYDATA')
savepath = getenv('MYSAVEPATH')

parts = strsplit(libpath, filesep);
parent_path = strjoin(parts(1:end-1), filesep);
addpath(genpath(parent_path))

init_unlocbox();
gsp_start;addpath(".");
M = readtable(datapath,'Delimiter',',','ReadRowNames',true,'ReadVariableNames',true);% M is raw data in table form, rows are genes and columns are cells.
M0 = table2array(M); % data matrix

[X,network] = G2S3(M0);

writetable(cell2table(M.Properties.VariableNames','VariableNames',{'cellID'}),strcat(savepath,'/cellids.csv'));
writematrix(X,strcat(savepath,'/imputed.csv'));
genetable = cell2table(M.Properties.RowNames,'VariableNames',{'geneID'});
writetable(genetable,strcat(savepath,'/geneids.csv'));