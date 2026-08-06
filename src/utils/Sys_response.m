function [KLQRi, Kf, L_lqgi, L_lqri, H_lqgi, H2_lqgi, H3_lqgi, C_lqg, C_lqgi] = Sys_response(A,B,CLQR,CKF,Ts,Q,Qi,R,W,V,r)

clc;

nx = size(A,1); nu = size(B,2); ny = size(CLQR,1);
Aint = [A zeros(nx,1); CLQR,0.99999999.*eye(ny,ny)];
Bint = [B;zeros(ny,nu)];
CKFi = [CLQR,zeros(ny,ny)];

KLQR = lqr(ss(A,B,CLQR,0,Ts),Q,R);

[KLQRi, ~, ~] =lqr(ss(Aint,Bint,[CLQR 0],0,Ts),Qi,R);

zeros(size(A,1),1);

[Kf,~,~] =  dlqr(A',CKF',W,V);
Kf=Kf';
KFi = dlqr(Aint',CKFi',W,V);
KFi = KFi';

L_lqr = ss(A,B,KLQR,0,Ts); %GC
L_lqg = (ss(A-B*KLQR-Kf*CKF,Kf,KLQR,0,Ts)*ss(A,B,CKF,0,Ts));

L_lqri = ss(Aint,Bint,KLQRi,0,Ts);
L_lqgi = (ss(Aint-Bint*KLQRi-KFi*CKFi,KFi,KLQRi,0,Ts)*ss(Aint,Bint,CKFi,0,Ts));

C_lqg = ss(A-B*KLQR-Kf*CKF,Kf,KLQR,0,Ts); % Kalman measure y_tilde to LQR input u transfer function C(s)
C_lqgi = ss(Aint-Bint*KLQRi-KFi*CKFi,KFi,KLQRi,0,Ts); % Kalman measure y_tilde to LQR input u transfer function C(s)

H_lqr = minreal((ss(A,B,CKF,0,Ts))/((L_lqr)+1)); %complementary sensitivity
H_lqg = minreal((ss(A,B,CKF,0,Ts))/((L_lqg)+1)); %disturbance to output (G/1+GC)

H_lqri = minreal((ss(Aint,Bint,CKFi,0,Ts))/((L_lqri)+1));
H_lqgi = minreal((ss(Aint,Bint,CKFi,0,Ts))/((L_lqgi)+1)); %disturbance to output (G/1+GC), PS Astrom

H2_lqr = minreal((L_lqr)/((L_lqr)+1)); 
H2_lqg = minreal((L_lqg)/((L_lqg)+1)); %reference to output (GC/(1+GC))

H2_lqri = minreal((L_lqri)/((L_lqri)+1));
H2_lqgi = minreal((L_lqgi)/((L_lqgi)+1)); %reference to output (GC/(1+GC)) , T Astrom

H3_lqri = minreal(1/((L_lqri)+1)); %sensitivity
H3_lqgi = minreal(1/((L_lqgi)+1)); %reference to output (1/(1+GC))

H3_lqr = minreal(1/((L_lqr)+1));
H3_lqg = minreal(1/((L_lqg)+1)); %reference to output (1/(1+GC)), S Astrom
opts = bodeoptions;
opts.FreqUnits = 'Hz';
opts.Title.FontSize = 15;
opts.XLabel.FontSize = 13;
opts.YLabel.FontSize = 13;
opts.TickLabel.FontSize = 11;
opts.Grid = 'on';

if r == 12
    % figure
    % margin(L_lqri); annotation('textbox',[.2 .25 .3 .3],'String','LQRi Margin','FitBoxToText','on')
    % set(gcf,'Position',[0 0 500 500])
    
    figure
    margin(L_lqgi); annotation('textbox',[.2 .25 .3 .3],'String','LQGi Margin','FitBoxToText','on')
    set(gcf,'Position',[0 0 500 500])
    
    figure
    bode(ss(A,B,CLQR,0,Ts), L_lqri,L_lqgi,opts); 
    ax = findall(gcf,'type','axes');
    legend(ax(3),'Open Loop','L_{lqri}','L_{lqgi}','Location','northeast'); 
    title('Open vs I/O TF')
    set(gcf,'Position',[0 0 500 500])
   
    %figure
    %bode(KalmanTF,opts); 
    %ax = findall(gcf,'type','axes');
    %legend(ax(3),'Open Loop','L_{lqri}','L_{lqgi}','Location','northeast'); 
    %title('Open vs I/O TF')
    %set(gcf,'Position',[0 0 500 500])
   
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
    bode(ss(A,B,CLQR,0,Ts),H_lqri,H_lqgi,opts);
    ax = findall(gcf,'type','axes');
    legend(ax(3),'Open Loop','H_{lqri}','H_{lqgi}','Location','northwest'); 
    title('Open vs Complimentary TF')
    set(gcf,'Position',[0 0 500 500])
    
    figure
    bode(ss(A,B,CLQR,0,Ts), H2_lqri,H2_lqgi,opts);
    ax = findall(gcf,'type','axes');
    legend(ax(3),'Open Loop','H2_{lqri}','H2_{lqgi}','Location','northeast'); 
    title('Open vs Closed-loop TF')
    set(gcf,'Position',[0 0 500 500])
    
    figure
    bode(H3_lqri,H3_lqgi,opts);
    ax = findall(gcf,'type','axes');
    legend(ax(3),'H3_{lqri}','H3_{lqgi}','Location','northwest');
    title('Open vs Sensitivity TF')
    set(gcf,'Position',[0 0 500 500])
end
