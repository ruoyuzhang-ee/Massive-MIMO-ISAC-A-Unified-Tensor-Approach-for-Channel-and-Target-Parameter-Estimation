function [ tau_est, f_corr2 ] = func_1D_searchTau( CChat, tau_max, K, K_bar, fs, search_dg)

[~,L] = size(CChat);
tau_est = zeros(1,L);
resol = tau_max/K_bar/4;
candiGroupMax = 1;
% 1 Coarse grid search
tau_Dict = [0+resol: resol: tau_max-resol];
g_tau_coarse = exp(-1i*2*pi* [1:K]' * tau_Dict *fs/K_bar);

tau_coarse = zeros(candiGroupMax,L);
for l = 1:L
    f_corr1 = zeros(length(tau_Dict),1);
    cc_hat = CChat(:,l);
    for ii = 1:length(tau_Dict)
        f_corr1(ii) = abs(cc_hat' * g_tau_coarse(:,ii))^2/(norm(cc_hat)^2 * norm(g_tau_coarse(:,ii))^2);
    end
    [val1, pos1] = sort(f_corr1, 'descend');
    tau_coarse(:,l) = tau_Dict(pos1(1:candiGroupMax)); % figure; plot(tau_Dict,f_corr1);
end

% 2 Refining the search in the vicinity of possible grid points.
f_corr_candiGroup = zeros(candiGroupMax,L);
tauest_candiGroup = zeros(candiGroupMax,L);
for l = 1:L
    cc_hat = CChat(:,l);
    for jj = 1:candiGroupMax
        flag = 0;
        search_dg_temp = resol/10; % The inital stepsize
        search_tau_temp = [tau_coarse(jj,l)-resol; tau_coarse(jj,l)+resol]; % The initial range
        while flag == 0
            tau_grid = [search_tau_temp(1,:) : search_dg_temp : search_tau_temp(2,:)];
            Num_Taucandi = length(tau_grid);
            f_corr2 = zeros(Num_Taucandi,1);
            g_tau_refine = exp(-1i*2*pi* [1:K]' * tau_grid *fs/K_bar);
            for ii = 1:Num_Taucandi
                f_corr2(ii) = abs(cc_hat' * g_tau_refine(:,ii))^2/(norm(cc_hat)^2 * norm(g_tau_refine(:,ii))^2);
            end
            [val, pos] = sort(f_corr2,'descend'); % Find the maxmum one
            hat_tau = tau_grid(pos(1));
            if search_dg_temp <= search_dg * tau_max * 1e-2;
                break;
            end
            search_tau_temp(1,:) = hat_tau - search_dg_temp;
            search_tau_temp(2,:) = hat_tau + search_dg_temp;
            search_dg_temp = search_dg_temp/10;
        end
        f_corr_candiGroup(jj,l) = f_corr2(pos(1));
        tauest_candiGroup(jj,l) = tau_grid(pos(1));
    end
end
for l = 1:L
    f_corrFinal = f_corr_candiGroup(:,l);
    [valFinal, posFinal] = sort(f_corrFinal,'descend');
    posFinalmax = posFinal(1);
    tau_est(l) = tauest_candiGroup(posFinalmax,l);
end
    
end