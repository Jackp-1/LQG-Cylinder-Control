function [A,B,U_1,Vr,Sd,Upr,Sr,Phi] = DMDc(X,U,thresh,r,dt)
%DMDc  Dynamic Mode Decomposition with control

% Split snapshots
X1 = X(:,1:end-1);
X2 = X(:,2:end);
Ups = U(:,1:end-1);

% Build augmented input matrix
Omega = [X1; Ups];

% SVD of augmented matrix
[Uo,S,V] = svd(Omega,'econ');

% SVD of X2
[U2s,~,~] = svd(X2,'econ');

% Singular values
Sd = diag(S);

% Rank selection
if thresh ~= 0
    energy = cumsum(Sd) / sum(Sd);
    r = find(energy >= thresh, 1, 'first');
end

% Truncate
Ur  = Uo(:,1:r);
Sr  = S(1:r,1:r);
Vr  = V(:,1:r);
Upr = U2s(:,1:r);

% Dimensions
n = size(X1,1);
q = size(Ups,1);

% Partition left singular vectors of Omega
U_1 = Ur(1:n,:);
U_2 = Ur(n+1:n+q,:);

% DMDc operators
A = Upr' * X2 * Vr * (Sr \ U_1') * Upr;
B = Upr' * X2 * Vr * (Sr \ U_2');


[W,D] = eig(A);

Phi = real(X2*Vr*Sr^(-1)*U_1'*Upr*W);
lambda = diag (D);
omega = log ( lambda )/dt;
% Compute DMD Solution
x1 = X(:, 1);
b = Phi \x1;
t = 200;
time_dynamics = zeros(r, length (t));
for iter = 1: length (t)
    time_dynamics (:, iter ) = (b.*exp ( omega*t(iter )));
end
X_dmd = Phi * time_dynamics ;


end