function [ Uhat ] = func_Algeb_TensorCPD( Y, L, AlgebOpt )

[M,T,K] = size(Y);
K3 = AlgebOpt.K3;
L3 = AlgebOpt.L3;

SelectionMatrix3 = zeros(K3*T, K*T, L3);
for lll = 1:L3
    JJ_temp = [zeros(K3,lll-1) eye(K3) zeros(K3,L3-lll)];
    SelectionMatrix3(:,:,lll) = kron(JJ_temp, eye(T));
end

YY_Sel = zeros(K3*T, L3*M);
Y_mode1 = tens2mat(Y,1);
for lll = 1:L3
    YY_Sel(:,(lll-1)*M+1:lll*M) = SelectionMatrix3(:,:,lll) * transpose(Y_mode1);
end

%% Factor matrix 3
[UU, Sigma, VV] = svd(YY_Sel);
UU1 = UU(1:(K3-1)*T,1:L);
UU2 = UU(T+1:K3*T,1:L);

[MM, ZZ] = eig(pinv(UU1)*UU2);

% est tau
z_est = diag(ZZ);

BB3_hat = zeros(K,L);
for ll = 1:L
    z_est(ll) = z_est(ll)/norm(z_est(ll));
    BB3_hat(:,ll) = z_est(ll).^([1:1:K]);
end

%% Factor matrix 2
BB2_hat = zeros(T,L);
for ll = 1:L
    bb3_K3_hat = BB3_hat(1:K3,ll);
    BB2_hat(:,ll) = kron(bb3_K3_hat', eye(T)) * UU(:,1:L) * MM(:,ll);
end

%% Factor matrix 1
BB1_hat = zeros(M,L);
TT = inv(MM).';
for ll = 1:L
    bb3_L3_hat = BB3_hat(1:L3,ll);
    BB1_hat(:,ll) = kron(bb3_L3_hat', eye(M)) * conj(VV(:,1:L)) * Sigma(1:L,1:L) * TT(:,ll);
end

Uhat.BB1_hat = BB1_hat;
Uhat.BB2_hat = BB2_hat;
Uhat.BB3_hat = BB3_hat;
Uhat.BB3_hatMM = MM;
Uhat.BB3_hatZZ = ZZ;

end

