%% Data processing function
function [Xdat,Ydat,Dat] = loadData(fname)
fileID = fopen(fname);
FDATA=readmatrix(fname);
Xdat = FDATA(:,2);
Ydat = FDATA(:,3);
Dat = FDATA(:,4);
fclose all;
end