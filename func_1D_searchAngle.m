function [ sin_AoA_est, f_corr2 ] = func_1D_searchAngle( AAhat, W, N_MS, search_dg)

[~,L] = size(AAhat);
sin_AoA_est = zeros(1,L);
resol = 2/N_MS /4;
candiGroupMax = 2;
% 1 Coarse grid search
sin_DictAoA = [-1+resol: resol: 1-resol];
a_MS_coarse = W' *  1/sqrt(N_MS)*exp(1j*pi*(0:N_MS-1).'*sin_DictAoA);
sin_AoA_coarse = zeros(candiGroupMax,L);
for l = 1:L
    f_corr1 = zeros(length(sin_DictAoA),1);
    aa_hat = AAhat(:,l);
    for ii = 1:length(sin_DictAoA)
        f_corr1(ii) = abs(aa_hat' * a_MS_coarse(:,ii))^2/(norm(aa_hat)^2 * norm(a_MS_coarse(:,ii))^2);
    end
    [val1, pos1] = sort(f_corr1, 'descend');
    sin_AoA_coarse(:,l) = sin_DictAoA(pos1(1:candiGroupMax)); % figure; plot(sin_DictAoA,f_corr1);
end

% 2 Refining the search in the vicinity of possible grid points.
f_corr_candiGroup = zeros(candiGroupMax,L);
sinAoAest_candiGroup = zeros(candiGroupMax,L);
for l = 1:L
    aa_hat = AAhat(:,l);
    for jj = 1:candiGroupMax
        flag = 0;
        search_dg_temp = resol/10; % The inital stepsize
        search_DOA_temp = [sin_AoA_coarse(jj,l)-resol; sin_AoA_coarse(jj,l)+resol]; % The initial range
        while flag == 0
            degree_grid = [search_DOA_temp(1,:) : search_dg_temp : search_DOA_temp(2,:)];
            Num_DOAcandi = length(degree_grid);
            f_corr2 = zeros(Num_DOAcandi,1);
            a_MS_refine = W' *  1/sqrt(N_MS)*exp(1j*pi*(0:N_MS-1).'*degree_grid);
            for ii = 1:Num_DOAcandi
                f_corr2(ii) = abs(aa_hat' * a_MS_refine(:,ii))^2/(norm(aa_hat)^2 * norm(a_MS_refine(:,ii))^2);
            end
            [val, pos] = sort(f_corr2,'descend'); % Find the maxmum one
            hat_DOA = degree_grid(pos(1));
            if search_dg_temp <= search_dg
                break;
            end
            search_DOA_temp(1,:) = hat_DOA - search_dg_temp;
            search_DOA_temp(2,:) = hat_DOA + search_dg_temp;
            search_dg_temp = search_dg_temp/10;
        end
        f_corr_candiGroup(jj,l) = f_corr2(pos(1));
        sinAoAest_candiGroup(jj,l) = degree_grid(pos(1));
    end
end
for l = 1:L
    f_corrFinal = f_corr_candiGroup(:,l);
    [valFinal, posFinal] = sort(f_corrFinal,'descend');
    posFinalmax = posFinal(1);
    sin_AoA_est(l) = sinAoAest_candiGroup(posFinalmax,l);
end
    
end