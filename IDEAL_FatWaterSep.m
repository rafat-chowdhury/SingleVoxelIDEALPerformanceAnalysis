function [water_s_norm fat_s_norm res_psi_norm] = ...
    IDEAL_FatWaterSep(W, F, df_W, df_F, TE, psi)

%% Description
% Script to model the fitting efficiency and distance between potential
% solutions identified by an IDEAL model described in literature (Yu 2005)
% as a function of varying field inhomogeneity (psi_range). See reference
% below for more technical details

%% Author
    % Rafat Chowdhury
    % email: rafat.u.a.chowdhury@gmail.com
    % mob: 07469889441
    
%% Structure
% The script works in three parts
% The first facilitates the generation of echo data 
% The second stores the fat and water signals the IDEAL model 
% estimates as a function of different starting field inhomogeneities
% (psi_range)
% The third involves calculating the residue between the original echo
% data and the estimated values - the lower the value the more likely we
% have approached our true solution

%% Reference
% Reeder, Scott B., et al. "Iterative decomposition of water and fat with 
% echo asymmetry and least?squares estimation (IDEAL): application with 
% fast spin?echo imaging." Magnetic Resonance in Medicine: An Official 
% Journal of the International Society for Magnetic Resonance in Medicine 
% 54.3 (2005): 636-644.

%% Inputs:
% W - water concentration; df_W - water resonant frequency
% F - fat concentration; df_F - fat resonant frequency
% TE - echo times in ms
% psi - the field inhomogeneity for the initial data

%% Outputs:
% water_s_norm - range of water signals estimated by the IDEAL model
% fat_s_norm - range of fat signals estimated by the IDEAL model
% ress_psi_norm - range of residues estimated by the IDEAL model

%% Extract values and define constants
    df(1) = df_W; 

    df(2) = df_F;

    const = i*2*pi*TE;

%% Generate echo data
    raw = (W.*exp(const*df(1)) + F.*exp(const*df(2)))...
        .*exp(const*psi);

%% Estimate water, fat andd residue as a function of psi_range
    S = (raw');
    psi_range = -1000:1000;

    A = exp(-1i*2*pi*TE'*df);
    s = S.*exp(1i*2*pi*TE'*psi_range);
    met_norm = A\s;
    A_ls = inv(A'*A)*A';
    rho = A_ls*s;

    water_s_norm = abs(rho(1,:));
    fat_s_norm = abs(rho(2,:));

%% Calculate residue
    estim = zeros(length(psi_range),length(TE));
    for p = 1:length(psi_range)
    estim(p,:) = (((water_s_norm(:,p).*exp(const*df(1))) + (fat_s_norm(:,p).*exp(const*df(2)))).*exp(const*psi_range(:,p)));
    end

    for p = 1:length(psi_range)
    res_ec(p,:) = raw-estim(p,:);
    res_psi_norm(p,:) = squeeze(norm(res_ec(p,:),2));
    end
    
%% Normalise data
    res_psi_norm = res_psi_temp./max(res_psi_temp);
    water_s_norm = water_s_temp./max(water_s_temp);
    fat_s_norm = fat_s_temp./max(fat_s_temp);

%% Visualise data
figure()
plot(psi_range, res_psi_norm,'linewidth',3.0)
hold on
    plot(psi_range, water_s_norm,'linewidth',3.0)
hold on
plot(psi_range, fat_s_norm,'linewidth',3.0)
legend('Res','Water','Fat')
xlabel('\psi (Hz)')
ylabel('||R_i||_2 (i = 3)')
title('Residual (S_{true} - S_{est}) versus \psi')
set(gca, 'fontsize', 16, 'fontweight', 'bold','linewidth',3.0);
set(gcf,'color','w');