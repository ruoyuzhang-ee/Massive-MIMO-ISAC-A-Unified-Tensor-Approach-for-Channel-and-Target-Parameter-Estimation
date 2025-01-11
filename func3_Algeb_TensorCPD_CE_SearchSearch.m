function [ sin_AoA_est, sin_AoD_est, tau_est, alpha_est, H_est, Doppler_est ] = func3_Algeb_TensorCPD_CE_SearchSearch( Uhat, CEOptions, optionsAltOpt  )


F = CEOptions.F;
N_BS = CEOptions.N_BS;
W = CEOptions.W;
N_MS = CEOptions.N_MS;
search_dg = CEOptions.search_dg;
L = CEOptions.L;
fs = CEOptions.fs;
K = CEOptions.K;
K_bar = CEOptions.K_bar;
tau_max = CEOptions.tau_max;
T = CEOptions.T;
M = size(W,2);

% AAhat = Uhat{1,1};
% BBhat = Uhat{1,2};
% CChat = Uhat{1,3};

BB1_hat = Uhat.BB1_hat;
BB2_hat = Uhat.BB2_hat;
BB3_hat = Uhat.BB3_hat;
ZZ = Uhat.BB3_hatZZ;   %第三个因子矩阵的特征值

AAhat = BB1_hat;
BBhat = BB2_hat;
CChat = BB3_hat;

% est tau
z_est = diag(ZZ);
z_phase_est = atan2(imag(z_est),real(z_est)); %z_phase_est = atan(imag(z_est)./real(z_est));
% tau_est = -K_bar/2/pi/fs * z_phase_est.';
[ tau_est, f_corr2 ] = func_1D_searchTau( BB3_hat, tau_max, K, K_bar, fs, search_dg);

% est AoA
[ sin_AoA_est, f_corr2 ] = func_1D_searchAngle( BB1_hat, W, N_MS, search_dg);

% est AoD
%% off Grid AoD and Doppler Est
[sin_AoD_est, Doppler_est] = func_DopplerAngle_AltOpt2_Exhaust( BB2_hat, conj(F), N_BS, optionsAltOpt );


% est alpha
a_MS_matrix_est = 1/sqrt(N_MS)*exp(1j*pi*(0:N_MS-1).'*sin_AoA_est);
a_BS_matrix_est = 1/sqrt(N_BS)*exp(1j*pi*(0:N_BS-1).'*sin_AoD_est);
AAest = W' * a_MS_matrix_est;
BBest = zeros(T, L);
for l = 1:L
    w_list_est = Doppler_est(l) * [1:1:T];
    BBest(:,l) = diag(exp(1i*w_list_est)) * F.' * a_BS_matrix_est(:,l);
end

g_tau_est = exp(-1i*2*pi* [1:K]' * tau_est *fs/K_bar);
Y = tens2mat(CEOptions.Y,3);
Y_transpose = Y.'; 
Phi = zeros(T*M*K, L);
for ll = 1:L
    BB_AA_est_ll = kron(BBest(:,ll), AAest(:,ll));
    Phi(:,ll) = kron(g_tau_est(:,ll) , BB_AA_est_ll);
end
alpha_est = pinv(Phi) * reshape(Y_transpose,[],1);
alpha_est = alpha_est.';


% est H
H_est = zeros(N_MS, N_BS, K, T);
for tt = 1:T
    for kk = 1:K
        H_kk_est = zeros(N_MS,N_BS);
        for ll = 1:L
            H_kk_est = H_kk_est + alpha_est(ll) * exp(1i*Doppler_est(ll)*tt) * exp(-1i*2*pi*tau_est(ll)*fs*kk/K_bar) * a_MS_matrix_est(:,ll) * a_BS_matrix_est(:,ll).';
        end
        H_est(:,:,kk,tt) = H_kk_est;
    end
end


end

