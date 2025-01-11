% The demo file for MassiveMIMO ISAC with Unified Tensor Framework 


% Setting
c_light = 3*1e8;
fc = 28*1e9;
K_bar = 128; 
N_BS = 64;
N_MS = 8;

L = 4;
fs = 100 * 1e6;
delta_f = fs/K_bar;
tau_max = 0.5 * 1/delta_f; 
N_CP = ceil(fs * tau_max); 

%Training parameters
T = 16; 
M = N_MS; 
K = 16; 
 
SNR_dB = 10; 

% Doppler generate
T_sampling = 1/(fs); 
v_max = 2*108; %Target velocity, double due to the round trip
w_max = 2*pi * fc/c_light * (v_max)/3.6 * T_sampling * (N_CP+K_bar); 



%% 2. Generate sensing channel G
sin_AoD_true = sin(random('unif',-pi/3,pi/3,1,L));
sin_AoA_true = sin_AoD_true;
tau_true = sort(tau_max * random('unif',0,1,1,L));
alpha_true = 1/sqrt(2) * (randn(1,L)+ 1i*randn(1,L));
velocity_ll = sort(random('unif',-v_max,v_max,1,L));
Doppler_true = 2*pi * fc/c_light * velocity_ll/3.6 * T_sampling * (N_CP+K_bar);
posSort_true = [1:L];

a_MS_matrix = 1/sqrt(N_MS)*exp(1j*pi*(0:N_MS-1).'*sin_AoA_true);
a_BS_matrix = 1/sqrt(N_BS)*exp(1j*pi*(0:N_BS-1).'*sin_AoD_true);

H = zeros(N_MS, N_BS, K_bar, T);
for tt = 1:T
    for kk = 1:K_bar
        H_kk = zeros(N_MS,N_BS);
        for ll = 1:L
            H_kk = H_kk + alpha_true(ll) * exp(1i*Doppler_true(ll)*tt) * exp(-1i*2*pi*tau_true(ll)*fs*kk/K_bar) * a_MS_matrix(:,ll) * a_BS_matrix(:,ll).';
        end
        H(:,:,kk,tt) = H_kk;
    end
end


% Training
phase = 2*pi * rand(N_BS, T);
F_temp = exp(1j*phase);
F = sqrt(T)*F_temp/norm(F_temp, 'fro');
W = eye(N_MS);

search_dg = 1e-5;
CEOptions.F = F;
CEOptions.N_BS = N_BS;
CEOptions.W = W;
CEOptions.N_MS = N_MS;
CEOptions.search_dg = search_dg;
CEOptions.L = L;
CEOptions.fs = fs;
CEOptions.K = K;
CEOptions.K_bar = K_bar;
CEOptions.tau_max = tau_max;
CEOptions.T = T;

% Generate the tensor
Y = zeros(N_MS,T,K);
for tt = 1:T
    for kk = 1:K
        H_kk = H(:,:,kk,tt);
        Y_kk_tt = H_kk * F(:,tt);
        sigma2 = norm(Y_kk_tt,'fro')^2 * 10^(-SNR_dB/10);
        noise_vector = 1/sqrt(2) * (randn(N_MS,1)+ 1i*randn(N_MS,1));
        noise_vector = sqrt(sigma2) * noise_vector/norm(noise_vector,'fro');
        Y(:,tt,kk) = Y_kk_tt + noise_vector;
    end
end



%% 5. Proposed Tensor Method
% Factor matrix estimation
K3 = 8;
L3 = K+1-K3;
AlgebOpt.K3 = K3;
AlgebOpt.L3 = L3;
[Uhat] = func_Algeb_TensorCPD( Y, L, AlgebOpt );

% Parameter estimation
CEOptions.Y = Y;
optionsAltOpt.T = T;
optionsAltOpt.a_BS_matrix = a_BS_matrix;
optionsAltOpt.search_dg = search_dg;
optionsAltOpt.w_max = w_max;
optionsAltOpt.Num_Iter = 30;
[sin_AoA_est, sin_AoD_est, tau_est, alpha_est, H_est, Doppler_est] = func3_Algeb_TensorCPD_CE_SearchSearch( Uhat, CEOptions, optionsAltOpt );


%display true and estimated parameters
sin_AoA_true
sin_AoA_est

tau_true
tau_est

Doppler_true
Doppler_est

H_est_T = H_est(:,:,:,T);
H_true_T = H(:,:,1:K,T);
mse_H = norm(H_est_T(:) - H_true_T(:))^2 / norm(H_true_T(:))^2

