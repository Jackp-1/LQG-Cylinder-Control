function []=controllability_check(sys)
A = sys.A;
B = sys.B;
C = sys.C;
% Split into controllability form
[~,~,~,~,k]=ctrbf(A,B,C);
% Then the noncontrollable modes are contained in first size(A,1)-sum(k)
% rows of abar, bbar e.g.
nc=size(A,1)-sum(k);

% Split into observability form
[~,~,~,~,k]=obsvf(A,B,C);
% Then the noncontrollable modes are contained in first size(A,1)-sum(k)
% rows of abar, bbar e.g.
no=size(A,1)-sum(k);

if nc~=0 || no~=0
    disp(['Model of size ',num2str(size(A,1)) ,' has some uncontrollable or unobservable modes'])
end
