origPath = pwd;
cd /pool1/unshared/Active_Project_Backup/DURIAN/DURIAN_paper_clean/G2S3/unlocbox/;
init_unlocbox();
cd /pool1/unshared/Active_Project_Backup/DURIAN/DURIAN_paper_clean/G2S3/gspbox/;
gsp_start;addpath(".");
cd(origPath);
x = csvread('x.csv');
[x_impute,network] = G2S3(x);
csvwrite('x_impute.csv', x_impute);