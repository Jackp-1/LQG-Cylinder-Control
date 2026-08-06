function [KLQRi,Kf] = Sys_response_MIMO(A,B,C,Theta,Ts,Q,Qi,R,W,V,r)

clc; close all;

nx = size(A,1); nu = size(B,2); ny = size(C,1);
Aint = [A zeros(size(A,1),1); C,0.99999999.*eye(ny,nu)];
Bint = [B;zeros(ny,nu)];
CKFi = [C,zeros(ny,ny)]; %mapping to sensors

KLQR = lqr(ss(A,B,C,0,Ts),Q,R);

KLQRi =lqr(ss(Aint,Bint,[C 0],0,Ts),Qi,R);

[Kf,~,~] =  dlqr(A',Theta',W,V);
Kf=Kf';
KFi = dlqr(Aint',CKFi',W,V);
KFi = KFi';
%L_lqr = ss(A,B,KLQR,0,Ts); %GC
%L_lqg = (ss(A-B*KLQR-Kf*CKF,Kf,KLQR,0,Ts)*ss(A,B,CKF,0,Ts));

%L_lqri = ss(Aint,Bint,KLQRi,0,Ts);
L_lqgi = (ss(Aint-Bint*KLQRi-KFi*CKFi,KFi,KLQRi,0,Ts)*ss(Aint,Bint,CKFi,0,Ts));

%C_lqgi = ss(Aint-Bint*KLQRi-KFi*CKFi,KFi,KLQRi,0,Ts); % Kalman measure y_tilde to LQR input u transfer function C(s)

%H_lqr = minreal((ss(A,B,CKF,0,Ts))/((L_lqr)+1)); %complementary sensitivity
%H_lqg = minreal((ss(A,B,CKF,0,Ts))/((L_lqg)+1)); %disturbance to output (G/1+GC)

%H_lqri = minreal((ss(Aint,Bint,CKFi,0,Ts))/((L_lqri)+1));
H_lqgi = minreal((ss(Aint,Bint,CKFi,0,Ts))/((L_lqgi)+1)); %disturbance to output (G/1+GC), PS Astrom

%H2_lqr = minreal((L_lqr)/((L_lqr)+1)); 
%H2_lqg = minreal((L_lqg)/((L_lqg)+1)); %reference to output (GC/(1+GC))

%H2_lqri = minreal((L_lqri)/((L_lqri)+1));
H2_lqgi = minreal((L_lqgi)/((L_lqgi)+1)); %reference to output (GC/(1+GC)) , T Astrom

%H3_lqri = minreal(1/((L_lqri)+1)); %sensitivity
H3_lqgi = minreal(1/((L_lqgi)+1)); %reference to output (1/(1+GC))

%H3_lqr = minreal(1/((L_lqr)+1));
%H3_lqg = minreal(1/((L_lqg)+1)); %reference to output (1/(1+GC)), S Astrom

opts = bodeoptions;
opts.FreqUnits = 'Hz';
opts.Title.FontSize = 15;
opts.XLabel.FontSize = 13;
opts.YLabel.FontSize = 13;
opts.TickLabel.FontSize = 11;
opts.Grid = 'on';
if r == 11
    figure
    margin(L_lqgi); annotation('textbox',[.2 .25 .3 .3],'String','LQRi Margin','FitBoxToText','on')
    set(gcf,'Position',[0 0 500 500])
    
    % figure
    % margin(L_lqgi); annotation('textbox',[.2 .25 .3 .3],'String','LQGi Margin','FitBoxToText','on')
    % set(gcf,'Position',[0 0 500 500])
    % 
    figure
    bode(ss(A,B,C,0,Ts), L_lqgi,opts); 
    ax = findall(gcf,'type','axes');
    legend(ax(3),'Open Loop','L_{lqri}','Location','northeast'); 
    title('Open vs I/O TF')
    set(gcf,'Position',[0 0 500 500])
    
    % figure
    % bode(ss(A,B,CLQR,0,Ts),H_lqr,H_lqg,opts);
    % ax = findall(gcf,'type','axes');
    % legend(ax(3),'Open Loop','H_{lqr}','H_{lqg}'); 
    % set(gcf,'Position',[0 0 500 500])
    % 
    % figure
    % bode(ss(A,B,CLQR,0,Ts),H2_lqr,H2_lqg,opts);
    % ax = findall(gcf,'type','axes');
    % legend(ax(3),'Open Loop','H2_{lqr}','H2_{lqg}'); 
    % set(gcf,'Position',[0 0 500 500])
    % 
    % figure
    % bode(H3_lqr,H3_lqg,opts);
    % ax = findall(gcf,'type','axes');
    % legend(ax(3),'H3_{lqr}','H3_{lqg}'); 
    % set(gcf,'Position',[0 0 500 500])
    
    figure
    bode(ss(A,B,C,0,Ts),H_lqgi,opts);
    ax = findall(gcf,'type','axes');
    legend(ax(3),'Open Loop','H_{lqri}','Location','northwest'); 
    title('Open vs Complimentary TF')
    set(gcf,'Position',[0 0 500 500])
    
    figure
    bode(ss(A,B,C,0,Ts), H2_lqgi,opts);
    ax = findall(gcf,'type','axes');
    legend(ax(3),'Open Loop','H2_{lqri}','Location','northeast'); 
    title('Open vs Closed-loop TF')
    set(gcf,'Position',[0 0 500 500])
    
    figure
    bode(H3_lqgi,opts);
    ax = findall(gcf,'type','axes');
    legend(ax(3),'H3_{lqri}','Location','northwest');
    title('Open vs Sensitivity TF')
    set(gcf,'Position',[0 0 500 500])
end