libpath = getenv('CMFLIB')
datapath = getenv('MYDATA')
savepath = getenv('MYSAVEPATH')

parts = strsplit(libpath, filesep);
parent_path = strjoin(parts(1:end-1), filesep);
addpath(genpath(parent_path))

// fprintf('path is:\n')
// fprintf(path)

%load data
M = readtable(datapath,'Delimiter',',','ReadRowNames',true,'ReadVariableNames',true);% M is raw data in table form, rows are genes and columns are cells.
M0 = table2array(M); % data matrix
processed_data=process(M0');
[score1]=CMF(processed_data',1,1,0.0001,0.0001);% return the imputated data matrix
X = max(10.^score1-1,0);

writetable(cell2table(M.Properties.VariableNames','VariableNames',{'cellID'}),strcat(savepath,'/cellids.csv'));
writematrix(X,strcat(savepath,'/imputed.csv'));
genetable = cell2table(M.Properties.RowNames,'VariableNames',{'geneID'});
writetable(genetable,strcat(savepath,'/geneids.csv'));