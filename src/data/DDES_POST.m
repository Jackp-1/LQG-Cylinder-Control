%% Fig_20_DDES_Validate

addpath('data\');

DNS_Cp = readmatrix('Re1000_3D_DNS_Cp');
EXP_Cp = readmatrix('Re1000_Experimental_Cp');
URANS_Cp = readmatrix('Theta_Cp_base');
DDES_Cp = readmatrix('Re1000_DDES_Cp_TEMP');
DDES_Cp = sortrows(DDES_Cp,1);
%%

scatter(EXP_Cp(:,1),EXP_Cp(:,2));
hold on;
plot(DNS_Cp(:,1),DNS_Cp(:,2));
plot(URANS_Cp(:,1)./180,URANS_Cp(:,2));
plot(DDES_Cp(:,1),DDES_Cp(:,2));
grid on;