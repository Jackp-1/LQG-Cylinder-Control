%% Cylinder Post Processing

%% Initialise Data
clc; close all; clear all;
addpath utils\;
addpath simulink\;
addpath data\;
load('train.mat');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
% Parameters;
nx = 999; %Number of grid points in x direction
ny = 499; %Number of grid points in y direction
N = size(Data,1); %Size of data
L = size(Data,2); %Length of training data
r_max = 18; %Max truncation rank size
o_max = 11; %Max sensor size
truncation_range = 1:18;
U_inf = 0.01460734693878; %Freestream velocity
q = 0.5*1.225*U_inf^2; %Dynamic Head
Ts = 0.4;  %Sampling time
fs=1/Ts;
t=0:Ts:Ts*200000-Ts;
f_shed=0.00305176;
Tau = t*f_shed;
ReadData
U = Train_input(:,2);
Ut=U';
meanData = mean(Data(:,1:1553),2);
meanLift = mean(Train_lift(1:1553,2));
Data_submu = Data-meanData;
lift_submu = Train_lift(:,2)-meanLift;
%% DMD of Cylinder
mod_ss = input('Use modified state space (1=yes)?'); %1 for yes 0 for no.
thresh=0;
Upr = cell(1,r_max);
x_train = cell(1,r_max);
sys = cell(1,r_max);
Phi = cell(1,r_max);
CKF = cell(1,r_max);
CLQR = cell(1,r_max);
Sys_mod = cell(1,r_max);

for r=truncation_range
    [A,B,U_1,Vr,Sd,Upr{r},Sr,Phi{r}] = DMDc(Data_submu,U',thresh,r,Ts);
    x_train{r}=Upr{r}'*(Data_submu); %PINV Method
    x0=x_train{r}(:,1);
    CKF{r}=Data_submu(4107,:)*pinv(x_train{r}); %Sensing
    CLQR{r}=(lift_submu)'*pinv(x_train{r}); %Direct lift measurement 
    sys_plant_force=ss(A,B,CLQR{r},0,Ts); %Force Sensing
    %sys_plant_vel=ss(A,B,CKF{r},0,Ts); %Variable (Velocity or Force) Sensing
    sys{r}=sys_plant_force;
    controllability_check(sys{r}) %Only returns if uncontrollable or unobservable
    if r==11    % Print Singular Value Spectrum, lift estimate
        FIG_4_Spectrum
        FIG_3_Lift
    end
    %Modify SS models enforcing stability
    [v,d] = eig(sys{r}.A);
    d = diag(d);
    if mod_ss == 1 %adjusts DMD models to fix poles on unit circle without modifying arg(Pole)
        clear scale Mod_A
        scale = array(1,r_max);
        Mod_A = array(1,r_max);
        for i = 1:r
            scale(i) = abs(d(i))-1;
            if scale(i)>0
                Mod_A(i) = dot((1-scale(i)/abs(d(i))),d(i)/abs(d(i)));
            else
                Mod_A(i)=d(i);
            end
        end
        Mod_A = diag(Mod_A);
        Mod_A = real(v*Mod_A/v);
        Sys_mod{r}=ss(Mod_A, sys{r}.B,sys{r}.C,0,Ts);
    end
    fprintf('Rank: %.0f\n',r) %Printing progress to console
end




%% Tuning gains
close all
n='SISO';%input('SISO or MIMO?','s');
KLQRi = cell(1,r_max);
Kf =  cell(1,r_max);
L_lqgi = cell(1,r_max);
L_lqri = cell(1,r_max);
H_lqgi = cell(1,r_max);
H2_lqgi = cell(1,r_max);
H3_lqgi = cell(1,r_max);
C_lqg = cell(1,r_max);
C_lqgi = cell(1,r_max);
KLQRi_MIMO = cell(1,r_max);
Qi_R = cell(1,r_max);
Kf_MIMO = cell(1,o_max);
switch n
    case 'SISO'
        for r = truncation_range
            close all
            Q=(CLQR{r}'*CLQR{r}); %State penalty Regulator
            Qi = 1e-6;
            Qi_R{r} = blkdiag(Q,Qi); %Integral action penalty
            R = 2e-5; %Actuation penalty Regulator
            W = 1; % process covariance Kalman filter
            V = 3e-4;% noise covariance Kalman filter
            [KLQRi{r}, Kf{r}, L_lqgi{r}, L_lqri{r}, H_lqgi{r}, H2_lqgi{r}, H3_lqgi{r}, C_lqg{r}, C_lqgi{r}] = Sys_response(sys{r}.A,sys{r}.B,CLQR{r},CKF{r},Ts,Q,Qi_R{r},R,W,V,r);
            %Write_UDF_Files; %Exporting files to UDF text
        end
    case 'MIMO'
        for r = truncation_range
            for o = truncation_range
                close all
                Q = (sys{r}.C'*sys{r}.C);
                Qi = blkdiag(Q,0);
                R = 3e-5;
                W = 1; 
                V = 3e-4;
                [KLQRi_MIMO{r}, Kf_MIMO{o}] = Sys_response_MIMO(sys{r}.A,sys{r}.B,sys{r}.C,Theta{o},Ts,Q,Qi,R,W,V,r);
                %Write_UDF_Files_MIMO; %Exporting files to UDF text
            end
        end
end

%% Sine Forcing
Sin_force_sim=sim(Sinforce);

%% Calculate Efficiencies
close all
Results = cell(1,r_max);
R_Sweep = {R5e_6 R75e_7 R1e_5 R2e_5 Mode_11 R5e_5};
R = [5e-6 75e-7 1e-5 2e-5 3e-5 5e-5];
Mean_Drag = cell(1,r_max);
Min_drag = mean(Half_drag(7500:end,2))*2;
P_actuation = cell(1,r_max);
Delta_P_Drag = cell(1,r_max);
P_Saved = cell(1,r_max);
Cmu = cell(1,r_max);
Cmu_RSweep = cell(1,r_max);
Mean_Cmu = cell(1,r_max);
Max_Cmu = cell(1,r_max);
Eta = cell(1,r_max);
Mean_Drag_RSweep = cell(1,6);
P_actuation_RSweep = cell(1,6);
Delta_P_Drag_RSweep = cell(1,6);
P_Saved_RSweep = cell(1,6);
Eta_RSweep = cell(1,6);
Mean_Cmu_RSweep = cell(1,6);
Max_Cmu_RSweep = cell(1,6);
Mean_Base_Drag = mean(Base_drag(:,2));
Constant_mean_drag = cell(1,4);
Eta_constant = cell(1,4);
Delta_P_drag_constant = cell(1,4);
P_Saved_constant = cell(1,4);

%SRS Calcs
Eta_SRS = Eta_calc(mean(SRS_Base_Drag(10000:end,2))/(q*2),mean(SRS_Control(10000:end,end))/(q*2),Min_drag); %Need to find min drag for 3d
P_actuation_SRS = 1.225*mean(abs(SRS_Control(:,end-1)))^3*0.05*2; %accounting for span area of control 
Delta_P_Drag_SRS = (mean(SRS_Base_Drag(10000:end,2))-mean(SRS_Control(10000:end,end)))*U_inf;
P_Saved_SRS = (Delta_P_Drag_SRS-P_actuation_SRS)./(mean(SRS_Base_Drag(10000:end,2))*U_inf)*100;

%rank sweep calcs
for i = 4:16
    Results{i} = eval(['Mode_',num2str(i)]);
    Cmu{i} =   (2*Results{i}(:,end-1).^2*0.05)/U_inf^2;
    Mean_Drag{i} = mean(Results{i}(15477:23216,end)); % 20 <= Tau <=30 == 15477:23216
    P_actuation{i} = 1.225*mean(abs(Results{i}(15477:23216,end-1)))^3*0.05; 
    Delta_P_Drag{i} = (Mean_Base_Drag-mean(Results{i}(15477:23216,end)))*U_inf;
    P_Saved{i} = (Delta_P_Drag{i}-P_actuation{i})./(Mean_Base_Drag*U_inf)*100;
    Mean_Cmu{i} = mean(Cmu{i});
    Max_Cmu{i} = max(Cmu{i});
    Eta{i} = Eta_calc(Mean_Base_Drag,Mean_Drag{i},Min_drag);
end

%R sweep calcs
for i = 1:6
    Cmu_RSweep{i} = (2*R_Sweep{i}(:,end-1).^2*0.05)/U_inf^2;
    Mean_Cmu_RSweep{i} = mean(Cmu_RSweep{i}(15477:23216));
    Max_Cmu_RSweep{i} = max(Cmu_RSweep{i});
    Mean_Drag_RSweep{i} = mean(R_Sweep{i}(15477:23216,end));
    P_actuation_RSweep{i} = 1.225*mean(abs(R_Sweep{i}(15477:23216,end-1))).^3*0.05;
    Delta_P_Drag_RSweep{i} = (Mean_Base_Drag-mean(R_Sweep{i}(:,end)))*U_inf;
    P_Saved_RSweep{i} = (Delta_P_Drag_RSweep{i}-P_actuation_RSweep{i})./(Mean_Base_Drag*U_inf)*100;
    Eta_RSweep{i} = Eta_calc(Mean_Base_Drag,Mean_Drag_RSweep{i},Min_drag);
end

%Constant input calcs
Constant_Drag = {Drag_constant_quarter(15477:23216,2), Drag_constant_half(15477:23216,2), Drag_constant_average(15477:23216,2), Drag_constant_double(15477:23216,2)};
Constant_Lift = {Lift_constant_quarter(15477:23216,2), Drag_constant_half(15477:23216,2), Lift_constant_average(15477:23216,2), Drag_constant_double(15477:23216,2)};

P_actuation_constant = 1.225.*[0.002,0.0028,0.0039,0.0055].^3.*0.1;

for i = 1:4
    Constant_mean_drag{i}= mean(Constant_Drag{i});
    Eta_constant{i} = Eta_calc(Mean_Base_Drag,Constant_mean_drag{i},Min_drag);
    Delta_P_drag_constant{i} = (Mean_Base_Drag-Constant_mean_drag{i})*U_inf;
    P_Saved_constant{i} = (Delta_P_drag_constant{i}-P_actuation_constant(i))./(Mean_Base_Drag*U_inf)*100;
end

P_actuation_noise = mean(1.225*0.05*(abs(Disturbance_Noise(12000:end,end-1).^3)));
Delta_P_drag_noise = mean((Mean_Base_Drag-Disturbance_Noise(12000:end,end))*U_inf);
P_Saved_noise = ((Delta_P_drag_noise-P_actuation_noise)./(Mean_Base_Drag*U_inf))*100;

%% Plotting
close all
%{
for i=truncation_range
    plotContourPhi(Phi{18}(:,i),Xdat,Ydat,1) %DMD Modes
end

for i=truncation_range
    plotContourPhi(Upr{18}(:,i),Xdat,Ydat,1)
end
%}

FIG_3_Lift
FIG_4_Spectrum
FIG_5_Bode
Fig_6_ControllerPerformance
FIG_7_Sensitivity2
Fig_8_Sensitivity
Fig_9_OpenLoop
FIG_10_Constant
FIG_11_rsweep
FIG_12_Cp_theta_wake
FIG_13_Cumsum
Fig_14_snapshots
FIG_15_Integral
FIG_16_DisturbanceNoise
Fig_17_Reynolds
FIG_18_SRS
FIG_19_PSD
