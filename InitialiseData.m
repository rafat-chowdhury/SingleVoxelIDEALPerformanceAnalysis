%% Initialise data for Fat water separation via IDEAL analysis

% Define metabolite concentrations (a.u.)
    W = 2;
    F = 1;
    
% Define metabolite resonances (Hz)
    df_W = 0;
    df_F = 224;

% Define echo times
    TE = [-1.8e-3 0 2.5e-3];
    
% Define true field inhomogeneity
    psi = [100];

% Run analysis script
[water_s_norm fat_s_norm res_psi_norm] = ...
    IDEAL_FatWaterSep(W, F, df_W, df_F, TE,...
    psi);