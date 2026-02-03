%% 1. Start a UQLab session
uqlab 

%% 2. Surrogates for wind-off tip deflections (uncertain variables: EI_1, EI_2, G)
% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [60, 80]*1e9; % (uncertain Young's modulus E for OOP bending) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [60, 80]*1e9; % (uncertain Young's modulus E for IP bending) lower and upper uncertainty bound
InputOpts.Marginals(3).Type = 'Uniform';
InputOpts.Marginals(3).Parameters = [26, 32]*1e9; % (uncertain Young's modulus G for torsion) lower and upper uncertainty bound
myInput = uq_createInput(InputOpts); 
alp = linspace(pi/2,0,15);  %ranges of pitch angles...

set(0, 'DefaultFigureVisible', 'off');
surrogates_tip_deflection_windoff = cell(length(alp), 1);

for ii = 1:length(alp)
    % Description of the physical model for UQLab
    ModelOpts.mFile = 'model_tip_deflection_windoff';
    ModelOpts.isVectorized = false;
    ModelOpts.Parameters = [alp(ii)];
    myModel_tip_deflection_windoff = uq_createModel(ModelOpts);

    N_train = 5;                                            % initial training set size (the set will be updated until the surrogate validation error is low enough)
    MetaOpts.Type = 'Metamodel';                            % 'metamodel': another word for 'surrogate'    
    MetaOpts.MetaType = 'PCE';
    MetaOpts.Input = myInput;                               % probability distribution for the uncertain variables
    MetaOpts.FullModel = myModel_tip_deflection_windoff;    % the physical model as a UQLab object
    MetaOpts.ExpDesign.NSamples = N_train;                  % 'experimental design' (ExpDesign): another word for 'training set'
    if strcmp(MetaOpts.MetaType, 'Kriging')
        MetaOpts.ExpDesign.Sampling = 'User';
    end

    flag_parfor = true;             % can we run the physical model in parallel to build the training set? (True/False)
    seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
    N_train_increment = 8;          % we will increment the training set size until we reach convergence
    N_train_max = 500;              % training budget (i.e., maximum number of training points allowed)
    % run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
    % recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
    flag_test_for_mean_and_sigma = false;
    flag_test_set = true;          % will a test set be generated for further surrogate validation?

    % Plots generator for parameter sweeps for the uncertain variables
    inputs_name = ["E for OOP (Pa)", "E for IP (Pa)", "Torsion modulus (Pa)"];     % list of the names of the uncertain variables
    outputs_name = ["Tip deflection (m)"];      % list of the names of the QIs
    N_outputs = length(outputs_name);           % number of quantities of interest (QIs)
    descriptive_title_for_plots = sprintf('%s surrogate (pitch angle (rad):%.2e)', MetaOpts.MetaType, alp(ii));
    N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
    plotsfolderName = 'tip_deflection_windoff_uncertain_E'; 
    subfolder_plotsfolderName = sprintf('pitch_case_%u', ii); % each 'case' refers to one fixed combination of design variables
    fullPath = fullfile(plotsfolderName, subfolder_plotsfolderName);
    mkdir(fullPath);
    mkdir(fullPath, 'plots_uq');
    try
        surrogates_tip_deflection_windoff{ii, 1} =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, fullPath, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
        elementToSave = surrogates_tip_deflection_windoff{ii, 1}; 
        save(fullfile(fullPath, sprintf('surrogate_tip_deflection_windoff_pitch_case_%u.mat', ii)), 'elementToSave'); % save the surrogate
        uncertain_variables_exploration(elementToSave, inputs_name, outputs_name, descriptive_title_for_plots, N_eval, seed, fullPath); % plots generator using the surrogates 
        % add readme to explain each 'case'
        fileID = fopen(fullfile(fullPath, 'readme.txt'), 'w');
        fprintf(fileID, 'Surrogates for each quantity of interest (QI) as a function of the uncertain variables (the deterministic variables are fixed).\n\n');
        fprintf(fileID, 'Legend:\n');
        N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
        N_variables = length(inputs_name);            % number of uncertain variables
        for ll = 1:N_outputs
            fprintf(fileID, 'QI %d: %s\n', ll, outputs_name(ll));
        end   
        for ll = 1:N_variables
            fprintf(fileID, 'Uncertain variable %d: %s\n', ll, inputs_name(ll));
        end    
        fprintf(fileID, 'Deterministic variable 1: Pitch angle; Case: %d out of %d; Value: %.5f (rad)\n', ii, length(alp), alp(ii));
        fprintf(fileID, 'Methodology: %s\n\n', descriptive_title_for_plots);
        fprintf(fileID, 'The trained surrogates and most of the design exploration figures are stored externally due to size limits.\n'); 
        fclose(fileID);
        disp('Readme file for the surrogates has been created successfully.');
    catch ME
        fprintf('Error in sample (%d): %s\n', ii, ME.message);
    end
end    
set(0, 'DefaultFigureVisible', 'on');

%% 3. Surrogates for the modal frequencies (uncertain variables: E, EI_2, G)
% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [60, 80]*1e9; % (uncertain Young's modulus E for OOP bending) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [60, 80]*1e9; % (uncertain Young's modulus E for IP bending) lower and upper uncertainty bound
InputOpts.Marginals(3).Type = 'Uniform';
InputOpts.Marginals(3).Parameters = [26, 32]*1e9; % (uncertain Young's modulus G for torsion) lower and upper uncertainty bound
myInput = uq_createInput(InputOpts); 
alp = linspace(pi/2,0,15);  %ranges of pitch angles...

set(0, 'DefaultFigureVisible', 'off');
surrogates_modal_frequencies = cell(length(alp), 1);

for ii = 1:length(alp)
    % Description of the physical model for UQLab
    ModelOpts.mFile = 'model_evals_windoff';
    ModelOpts.isVectorized = false;
    ModelOpts.Parameters = [alp(ii)];
    myModel_evals_windoff = uq_createModel(ModelOpts);

    N_train = 5;                                            % initial training set size (the set will be updated until the surrogate validation error is low enough)
    MetaOpts.Type = 'Metamodel';                            % 'metamodel': another word for 'surrogate'    
    MetaOpts.MetaType = 'PCE';
    MetaOpts.Input = myInput;                               % probability distribution for the uncertain variables
    MetaOpts.FullModel = myModel_evals_windoff;    % the physical model as a UQLab object
    MetaOpts.ExpDesign.NSamples = N_train;                  % 'experimental design' (ExpDesign): another word for 'training set'
    if strcmp(MetaOpts.MetaType, 'Kriging')
        MetaOpts.ExpDesign.Sampling = 'User';
    end

    flag_parfor = true;             % can we run the physical model in parallel to build the training set? (True/False)
    seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
    N_train_increment = 8;          % we will increment the training set size until we reach convergence
    N_train_max = 500;              % training budget (i.e., maximum number of training points allowed)
    % run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
    % recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
    flag_test_for_mean_and_sigma = false;
    flag_test_set = false;          % will a test set be generated for further surrogate validation?

    % Plots generator for parameter sweeps for the uncertain variables
    inputs_name = ["E for OOP (Pa)", "E for IP (Pa)", "Torsion modulus (Pa)"];     % list of the names of the uncertain variables
    outputs_name = ["Mode 1 (Hz)", "Mode 2 (Hz)", "Mode 3 (Hz)", "Mode 4 (Hz)"];      % list of the names of the QIs
    N_outputs = length(outputs_name);           % number of quantities of interest (QIs)
    descriptive_title_for_plots = sprintf('%s surrogate (pitch angle (rad):%.2e)', MetaOpts.MetaType, alp(ii));
    N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
    plotsfolderName = 'modal_frequencies_uncertain_E'; 
    subfolder_plotsfolderName = sprintf('pitch_case_%u', ii); % each 'case' refers to one fixed combination of design variables
    fullPath = fullfile(plotsfolderName, subfolder_plotsfolderName);
    mkdir(fullPath);
    mkdir(fullPath, 'plots_uq');
    try
        surrogates_modal_frequencies{ii, 1} =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, fullPath, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
        elementToSave = surrogates_modal_frequencies{ii, 1}; 
        save(fullfile(fullPath, sprintf('surrogate_modal_frequencies_pitch_case_%u.mat', ii)), 'elementToSave'); % save the surrogate
        uncertain_variables_exploration(elementToSave, inputs_name, outputs_name, descriptive_title_for_plots, N_eval, seed, fullPath); % plots generator using the surrogates 
        % add readme to explain each 'case'
        fileID = fopen(fullfile(fullPath, 'readme.txt'), 'w');
        fprintf(fileID, 'Surrogates for each quantity of interest (QI) as a function of the uncertain variables (the deterministic variables are fixed).\n\n');
        fprintf(fileID, 'Legend:\n');
        N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
        N_variables = length(inputs_name);            % number of uncertain variables
        for ll = 1:N_outputs
            fprintf(fileID, 'QI %d: %s\n', ll, outputs_name(ll));
        end   
        for ll = 1:N_variables
            fprintf(fileID, 'Uncertain variable %d: %s\n', ll, inputs_name(ll));
        end    
        fprintf(fileID, 'Deterministic variable 1: Pitch angle; Case: %d out of %d; Value: %.5f (rad)\n', ii, length(alp), alp(ii));
        fprintf(fileID, 'Methodology: %s\n\n', descriptive_title_for_plots);
        fprintf(fileID, 'The trained surrogates and most of the design exploration figures are stored externally due to size limits.\n'); 
        fclose(fileID);
        disp('Readme file for the surrogates has been created successfully.');
    catch ME
        fprintf('Error in sample (%d): %s\n', ii, ME.message);
    end
end    
set(0, 'DefaultFigureVisible', 'on');

%% 4. Surrogates for the modal frequencies (uncertain variable: pitch angle)
% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0, pi/2]; % (pitch angle) lower and upper uncertainty bound
myInput = uq_createInput(InputOpts); 
alp = linspace(60*1e9,80*1e9,3);  %ranges of Young's modulus E

set(0, 'DefaultFigureVisible', 'off');
surrogates_modal_frequencies = cell(length(alp), 1);

for ii = 1:length(alp)
    % Description of the physical model for UQLab
    ModelOpts.mFile = 'model_evals_windoff_uncertain_pitch';
    ModelOpts.isVectorized = false;
    ModelOpts.Parameters = [alp(ii)];
    myModel_evals_windoff = uq_createModel(ModelOpts);

    N_train = 5;                                            % initial training set size (the set will be updated until the surrogate validation error is low enough)
    MetaOpts.Type = 'Metamodel';                            % 'metamodel': another word for 'surrogate'    
    MetaOpts.MetaType = 'PCE';
    MetaOpts.Input = myInput;                               % probability distribution for the uncertain variables
    MetaOpts.FullModel = myModel_evals_windoff;    % the physical model as a UQLab object
    MetaOpts.ExpDesign.NSamples = N_train;                  % 'experimental design' (ExpDesign): another word for 'training set'
    if strcmp(MetaOpts.MetaType, 'Kriging')
        MetaOpts.ExpDesign.Sampling = 'User';
    end

    flag_parfor = true;             % can we run the physical model in parallel to build the training set? (True/False)
    seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
    N_train_increment = 8;          % we will increment the training set size until we reach convergence
    N_train_max = 100;              % training budget (i.e., maximum number of training points allowed)
    % run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
    % recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
    flag_test_for_mean_and_sigma = false;
    flag_test_set = false;          % will a test set be generated for further surrogate validation?

    % Plots generator for parameter sweeps for the uncertain variables
    inputs_name = ["Pitch angle (rad)"];     % list of the names of the uncertain variables
    outputs_name = ["Mode 1 (Hz)", "Mode 2 (Hz)", "Mode 3 (Hz)", "Mode 4 (Hz)"];      % list of the names of the QIs
    N_outputs = length(outputs_name);           % number of quantities of interest (QIs)
    descriptive_title_for_plots = sprintf('%s surrogate (Young modulus E (Pa):%.2e)', MetaOpts.MetaType, alp(ii));
    N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
    plotsfolderName = 'modal_frequencies_uncertain_pitch_angle'; 
    subfolder_plotsfolderName = sprintf('Young_modulus_case_%u', ii); % each 'case' refers to one fixed combination of design variables
    fullPath = fullfile(plotsfolderName, subfolder_plotsfolderName);
    mkdir(fullPath);
    mkdir(fullPath, 'plots_uq');
    try
        surrogates_modal_frequencies{ii, 1} =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, fullPath, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
        elementToSave = surrogates_modal_frequencies{ii, 1}; 
        save(fullfile(fullPath, sprintf('surrogate_modal_frequencies_Young_modulus_case_%u.mat', ii)), 'elementToSave'); % save the surrogate
        uncertain_variables_exploration(elementToSave, inputs_name, outputs_name, descriptive_title_for_plots, N_eval, seed, fullPath); % plots generator using the surrogates 
        % add readme to explain each 'case'
        fileID = fopen(fullfile(fullPath, 'readme.txt'), 'w');
        fprintf(fileID, 'Surrogates for each quantity of interest (QI) as a function of the uncertain variables (the deterministic variables are fixed).\n\n');
        fprintf(fileID, 'Legend:\n');
        N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
        N_variables = length(inputs_name);            % number of uncertain variables
        for ll = 1:N_outputs
            fprintf(fileID, 'QI %d: %s\n', ll, outputs_name(ll));
        end   
        for ll = 1:N_variables
            fprintf(fileID, 'Uncertain variable %d: %s\n', ll, inputs_name(ll));
        end    
        fprintf(fileID, 'Deterministic variable 1: Young modulus; Case: %d out of %d; Value: %.5f (Pa)\n', ii, length(alp), alp(ii));
        fprintf(fileID, 'Methodology: %s\n\n', descriptive_title_for_plots);
        fprintf(fileID, 'The trained surrogates and most of the design exploration figures are stored externally due to size limits.\n'); 
        fclose(fileID);
        disp('Readme file for the surrogates has been created successfully.');
    catch ME
        fprintf('Error in sample (%d): %s\n', ii, ME.message);
    end
end    
set(0, 'DefaultFigureVisible', 'on');

%% 5. Surrogates for wind-on tip deflections (uncertain variables: EI_1, EI_2, G)
% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [60, 80]*1e9; % (uncertain Young's modulus E for OOP bending) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [60, 80]*1e9; % (uncertain Young's modulus E for IP bending) lower and upper uncertainty bound
InputOpts.Marginals(3).Type = 'Uniform';
InputOpts.Marginals(3).Parameters = [26, 32]*1e9; % (uncertain Young's modulus G for torsion) lower and upper uncertainty bound
myInput = uq_createInput(InputOpts); 
asp = linspace(50,175,10);  %ranges of airspeeds...

set(0, 'DefaultFigureVisible', 'off');
surrogates_tip_deflection_windon = cell(length(asp), 1);

for ii = 1:length(asp)
    % Description of the physical model for UQLab
    ModelOpts.mFile = 'model_tip_deflection_windon';
    ModelOpts.isVectorized = false;
    ModelOpts.Parameters = [asp(ii)];
    myModel_tip_deflection_windon = uq_createModel(ModelOpts);

    N_train = 5;                                            % initial training set size (the set will be updated until the surrogate validation error is low enough)
    MetaOpts.Type = 'Metamodel';                            % 'metamodel': another word for 'surrogate'    
    MetaOpts.MetaType = 'PCE';
    MetaOpts.Input = myInput;                               % probability distribution for the uncertain variables
    MetaOpts.FullModel = myModel_tip_deflection_windon;     % the physical model as a UQLab object
    MetaOpts.ExpDesign.NSamples = N_train;                  % 'experimental design' (ExpDesign): another word for 'training set'
    if strcmp(MetaOpts.MetaType, 'Kriging')
        MetaOpts.ExpDesign.Sampling = 'User';
    end

    flag_parfor = true;             % can we run the physical model in parallel to build the training set? (True/False)
    seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
    N_train_increment = 8;          % we will increment the training set size until we reach convergence
    N_train_max = 500;              % training budget (i.e., maximum number of training points allowed)
    % run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
    % recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
    flag_test_for_mean_and_sigma = false;
    flag_test_set = true;          % will a test set be generated for further surrogate validation?

    % Plots generator for parameter sweeps for the uncertain variables
    inputs_name = ["E for OOP (Pa)", "E for IP (Pa)", "Torsion modulus (Pa)"];     % list of the names of the uncertain variables
    outputs_name = ["Tip deflection (m)"];      % list of the names of the QIs
    N_outputs = length(outputs_name);           % number of quantities of interest (QIs)
    descriptive_title_for_plots = sprintf('%s surrogate (airspeed (m/s):%.2e)', MetaOpts.MetaType, asp(ii));
    N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
    plotsfolderName = 'tip_deflection_windon_uncertain_E'; 
    subfolder_plotsfolderName = sprintf('airspeed_case_%u', ii); % each 'case' refers to one fixed combination of design variables
    fullPath = fullfile(plotsfolderName, subfolder_plotsfolderName);
    mkdir(fullPath);
    mkdir(fullPath, 'plots_uq');
    try
        surrogates_tip_deflection_windon{ii, 1} =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, fullPath, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
        elementToSave = surrogates_tip_deflection_windon{ii, 1}; 
        save(fullfile(fullPath, sprintf('surrogate_tip_deflection_windon_airspeed_case_%u.mat', ii)), 'elementToSave'); % save the surrogate
        uncertain_variables_exploration(elementToSave, inputs_name, outputs_name, descriptive_title_for_plots, N_eval, seed, fullPath); % plots generator using the surrogates 
        % add readme to explain each 'case'
        fileID = fopen(fullfile(fullPath, 'readme.txt'), 'w');
        fprintf(fileID, 'Surrogates for each quantity of interest (QI) as a function of the uncertain variables (the deterministic variables are fixed).\n\n');
        fprintf(fileID, 'Legend:\n');
        N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
        N_variables = length(inputs_name);            % number of uncertain variables
        for ll = 1:N_outputs
            fprintf(fileID, 'QI %d: %s\n', ll, outputs_name(ll));
        end   
        for ll = 1:N_variables
            fprintf(fileID, 'Uncertain variable %d: %s\n', ll, inputs_name(ll));
        end    
        fprintf(fileID, 'Deterministic variable 1: Airspeed; Case: %d out of %d; Value: %.5f (m/s)\n', ii, length(asp), asp(ii));
        fprintf(fileID, 'Methodology: %s\n\n', descriptive_title_for_plots);
        fprintf(fileID, 'The trained surrogates and most of the design exploration figures are stored externally due to size limits.\n'); 
        fclose(fileID);
        disp('Readme file for the surrogates has been created successfully.');
    catch ME
        fprintf('Error in sample (%d): %s\n', ii, ME.message);
    end
end    
set(0, 'DefaultFigureVisible', 'on');

%% 6. Surrogates for wind-on tip deflections with different angles of attack (uncertain variables: EI_1, EI_2, G)
% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [60, 80]*1e9; % (uncertain Young's modulus E for OOP bending) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [60, 80]*1e9; % (uncertain Young's modulus E for IP bending) lower and upper uncertainty bound
InputOpts.Marginals(3).Type = 'Uniform';
InputOpts.Marginals(3).Parameters = [26, 32]*1e9; % (uncertain Young's modulus G for torsion) lower and upper uncertainty bound
myInput = uq_createInput(InputOpts); 
asp = linspace(50,175,10);  %ranges of airspeeds...
aoa = linspace(1,10,10);    %ranges of angles of attack

set(0, 'DefaultFigureVisible', 'off');
surrogates_tip_deflection_windon = cell(length(asp), length(aoa));

for ii = 1:length(asp)
    for jj = 1:length(aoa)
        % Description of the physical model for UQLab
        ModelOpts.mFile = 'model_tip_deflection_windon_different_aoa';
        ModelOpts.isVectorized = false;
        ModelOpts.Parameters = [asp(ii) aoa(jj)];
        myModel_tip_deflection_windon = uq_createModel(ModelOpts);

        N_train = 5;                                            % initial training set size (the set will be updated until the surrogate validation error is low enough)
        MetaOpts.Type = 'Metamodel';                            % 'metamodel': another word for 'surrogate'    
        MetaOpts.MetaType = 'PCE';
        MetaOpts.Input = myInput;                               % probability distribution for the uncertain variables
        MetaOpts.FullModel = myModel_tip_deflection_windon;     % the physical model as a UQLab object
        MetaOpts.ExpDesign.NSamples = N_train;                  % 'experimental design' (ExpDesign): another word for 'training set'
        if strcmp(MetaOpts.MetaType, 'Kriging')
            MetaOpts.ExpDesign.Sampling = 'User';
        end

        flag_parfor = true;             % can we run the physical model in parallel to build the training set? (True/False)
        seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
        N_train_increment = 8;          % we will increment the training set size until we reach convergence
        N_train_max = 150;              % training budget (i.e., maximum number of training points allowed)
        % run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
        % recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
        flag_test_for_mean_and_sigma = false;
        flag_test_set = true;          % will a test set be generated for further surrogate validation?

        % Plots generator for parameter sweeps for the uncertain variables
        inputs_name = ["E for OOP (Pa)", "E for IP (Pa)", "Torsion modulus (Pa)"];     % list of the names of the uncertain variables
        outputs_name = ["Tip deflection (m)"];      % list of the names of the QIs
        N_outputs = length(outputs_name);           % number of quantities of interest (QIs)
        descriptive_title_for_plots = sprintf('%s surrogate (airspeed (m/s):%.2e, aoa (deg):%.2e)', MetaOpts.MetaType, asp(ii), aoa(jj));
        N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
        plotsfolderName = 'tip_deflection_windon_different_aoa_uncertain_E'; 
        subfolder_plotsfolderName = sprintf('airspeed_%u_aoa_%u', ii, jj); % each 'case' refers to one fixed combination of design variables
        fullPath = fullfile(plotsfolderName, subfolder_plotsfolderName);
        mkdir(fullPath);
        mkdir(fullPath, 'plots_uq');
        try
            surrogates_tip_deflection_windon{ii, jj} =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, fullPath, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
            elementToSave = surrogates_tip_deflection_windon{ii, jj}; 
            save(fullfile(fullPath, sprintf('surrogate_tip_deflection_windon_airspeed_case_%u_aoa_case_%u.mat', ii, jj)), 'elementToSave'); % save the surrogate
            uncertain_variables_exploration(elementToSave, inputs_name, outputs_name, descriptive_title_for_plots, N_eval, seed, fullPath); % plots generator using the surrogates 
             % add readme to explain each 'case'
            fileID = fopen(fullfile(fullPath, 'readme.txt'), 'w');
            fprintf(fileID, 'Surrogates for each quantity of interest (QI) as a function of the uncertain variables (the deterministic variables are fixed).\n\n');
            fprintf(fileID, 'Legend:\n');
            N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
            N_variables = length(inputs_name);            % number of uncertain variables
            for ll = 1:N_outputs
                fprintf(fileID, 'QI %d: %s\n', ll, outputs_name(ll));
            end   
            for ll = 1:N_variables
                fprintf(fileID, 'Uncertain variable %d: %s\n', ll, inputs_name(ll));
            end    
            fprintf(fileID, 'Deterministic variable 1: Airspeed; Case: %d out of %d; Value: %.5f (m/s)\n', ii, length(asp), asp(ii));
            fprintf(fileID, 'Deterministic variable 2: AoA; Case: %d out of %d; Value: %.5f (deg)\n', jj, length(aoa), aoa(jj));
            fprintf(fileID, 'Methodology: %s\n\n', descriptive_title_for_plots);
            fprintf(fileID, 'The trained surrogates and most of the design exploration figures are stored externally due to size limits.\n'); 
            fclose(fileID);
            disp('Readme file for the surrogates has been created successfully.');
        catch ME
            fprintf('Error in sample (%d, %d): %s\n', ii, jj, ME.message);            
        end
    end
end    
set(0, 'DefaultFigureVisible', 'on');

%% 7. Surrogates for wind-off spanwise deflections (uncertain variables: EI_1, EI_2, G)
% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [60, 80]*1e9; % (uncertain Young's modulus E for OOP bending) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [60, 80]*1e9; % (uncertain Young's modulus E for IP bending) lower and upper uncertainty bound
InputOpts.Marginals(3).Type = 'Uniform';
InputOpts.Marginals(3).Parameters = [26, 32]*1e9; % (uncertain Young's modulus G for torsion) lower and upper uncertainty bound
myInput = uq_createInput(InputOpts); 
alp = linspace(pi/2,0,15);  %ranges of pitch angles...

for ii = 1:length(alp)
    % generate training set, reduce dimension for outputs (PCA), and train GP surrogates
    N_MC = 5e2;
    X_MC = uq_getSample(myInput, N_MC, 'MC');
    Y_MC = nan(N_MC, 50);
    fun_tip_deflection_windoff_spanwise = @(x) model_tip_deflection_windoff_spanwise([x(1) x(2) x(3)], alp(ii));
    parfor kk=1:N_MC
        try
            Y_MC(kk, :) = fun_tip_deflection_windoff_spanwise(X_MC(kk, :));
        catch ME
            fprintf('Error in sample %d: %s\n', kk, ME.message);
        end
    end
    [coeff, score, latent, ~, explained, mu] = pca(Y_MC);
    cumExplained = cumsum(explained);
    numComponents = find(cumExplained >= 95, 1);  % e.g., retain 95% variance
    Y_reduced = score(:, 1:numComponents);
    for jj = 1:numComponents
        gpModel{jj} = fitrgp(X_MC, Y_reduced(:, jj));
    end

    % generate test set
    N_MC = 1e2;
    X_MC_test = uq_getSample(myInput, N_MC, 'MC');
    Y_MC_test = nan(N_MC, 50);
    parfor kk=1:N_MC
        try
            Y_MC_test(kk, :) = fun_tip_deflection_windoff_spanwise(X_MC_test(kk, :));
        catch ME
            fprintf('Error in sample %d: %s\n', kk, ME.message);
        end
    end
    % Predict in reduced space
    Y_pred_reduced = zeros(N_MC, numComponents);
    for jj = 1:numComponents
        Y_pred_reduced(:, jj) = predict(gpModel{jj}, X_MC_test);
    end
    % Reconstruct full output
    Y_pred_full = Y_pred_reduced * coeff(:, 1:numComponents)' + mu;
    mae = mean(abs(Y_MC_test - Y_pred_full), 'all'); % mean absolute error
    plotsfolderName = 'spanwise_deflection_windoff_uncertain_E'; 
    subfolder_plotsfolderName = sprintf('pitch_case_%u', ii); % each 'case' refers to one fixed combination of design variables
    fullPath = fullfile(plotsfolderName, subfolder_plotsfolderName);
    mkdir(fullPath);
    save(fullfile(fullPath, sprintf('surrogate_spanwise_deflection_windoff_pitch_case_%u.mat', ii)), 'gpModel'); % save the surrogate
    save(fullfile(fullPath, 'training_set.mat'), 'X_MC', 'Y_MC');             % Save both to a .mat file
    save(fullfile(fullPath, 'validation_set.mat'), 'X_MC_test', 'Y_MC_test');          % Save both to a .mat file
    save(fullfile(fullPath, 'mean_absolute_error_validation.mat'), 'mae');             % Save test performance to a .mat file
    save(fullfile(fullPath, sprintf('pca_spanwise_deflection_windoff_pitch_case_%u.mat', ii)), 'coeff', 'score', 'latent', 'explained', 'mu');            
    % add readme to explain each 'case'
    fileID = fopen(fullfile(fullPath, 'readme.txt'), 'w');
    fprintf(fileID, 'PCA dimension reduction and GP surrogates for spanwise tip deflection as a function of the uncertain variables (the deterministic variables are fixed).\n\n');
    fprintf(fileID, 'Legend:\n');
    inputs_name = ["E for OOP (Pa)", "E for IP (Pa)", "Torsion modulus (Pa)"];     % list of the names of the uncertain variables
    N_variables = length(inputs_name);            % number of uncertain variables
    for ll = 1:N_variables
        fprintf(fileID, 'Uncertain variable %d: %s\n', ll, inputs_name(ll));
    end    
    fprintf(fileID, 'Deterministic variable 1: Pitch angle; Case: %d out of %d; Value: %.5f (rad)\n', ii, length(alp), alp(ii));
    fclose(fileID);
    disp('Readme file for the surrogates has been created successfully.');
end    

%% 8. Surrogates for wind-on spanwise deflections (uncertain variables: EI_1, EI_2, G)
% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [60, 80]*1e9; % (uncertain Young's modulus E for OOP bending) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [60, 80]*1e9; % (uncertain Young's modulus E for IP bending) lower and upper uncertainty bound
InputOpts.Marginals(3).Type = 'Uniform';
InputOpts.Marginals(3).Parameters = [26, 32]*1e9; % (uncertain Young's modulus G for torsion) lower and upper uncertainty bound
myInput = uq_createInput(InputOpts); 
asp = linspace(50,175,10);  %ranges of airspeeds...

for ii = 1:length(asp)
    % generate training set, reduce dimension for outputs (PCA), and train GP surrogates
    N_MC = 5e2;
    X_MC = uq_getSample(myInput, N_MC, 'MC');
    Y_MC = nan(N_MC, 50);
    fun_tip_deflection_windon_spanwise = @(x) model_tip_deflection_windon_spanwise([x(1) x(2) x(3)], asp(ii));
    parfor kk=1:N_MC
        try
            Y_MC(kk, :) = fun_tip_deflection_windon_spanwise(X_MC(kk, :));
        catch ME
            fprintf('Error in sample %d: %s\n', kk, ME.message);
        end
    end
    [coeff, score, latent, ~, explained, mu] = pca(Y_MC);
    cumExplained = cumsum(explained);
    numComponents = find(cumExplained >= 95, 1);  % e.g., retain 95% variance
    Y_reduced = score(:, 1:numComponents);
    for jj = 1:numComponents
        gpModel{jj} = fitrgp(X_MC, Y_reduced(:, jj));
    end

    % generate test set
    N_MC = 1e2;
    X_MC_test = uq_getSample(myInput, N_MC, 'MC');
    Y_MC_test = nan(N_MC, 50);
    parfor kk=1:N_MC
        try
            Y_MC_test(kk, :) = fun_tip_deflection_windon_spanwise(X_MC_test(kk, :));
        catch ME
            fprintf('Error in sample %d: %s\n', kk, ME.message);
        end
    end
    % Predict in reduced space
    Y_pred_reduced = zeros(N_MC, numComponents);
    for jj = 1:numComponents
        Y_pred_reduced(:, jj) = predict(gpModel{jj}, X_MC_test);
    end
    % Reconstruct full output
    Y_pred_full = Y_pred_reduced * coeff(:, 1:numComponents)' + mu;
    mae = mean(abs(Y_MC_test - Y_pred_full), 'all'); % mean absolute error
    plotsfolderName = 'spanwise_deflection_windon_uncertain_E'; 
    subfolder_plotsfolderName = sprintf('airspeed_case_%u', ii); % each 'case' refers to one fixed combination of design variables
    fullPath = fullfile(plotsfolderName, subfolder_plotsfolderName);
    mkdir(fullPath);
    save(fullfile(fullPath, sprintf('surrogate_spanwise_deflection_windon_airspeed_case_%u.mat', ii)), 'gpModel'); % save the surrogate
    save(fullfile(fullPath, 'training_set.mat'), 'X_MC', 'Y_MC');             % Save both to a .mat file
    save(fullfile(fullPath, 'validation_set.mat'), 'X_MC_test', 'Y_MC_test');          % Save both to a .mat file
    save(fullfile(fullPath, 'mean_absolute_error_validation.mat'), 'mae');             % Save test performance to a .mat file
    save(fullfile(fullPath, sprintf('pca_spanwise_deflection_windon_pitch_case_%u.mat', ii)), 'coeff', 'score', 'latent', 'explained', 'mu');            
    % add readme to explain each 'case'
    fileID = fopen(fullfile(fullPath, 'readme.txt'), 'w');
    fprintf(fileID, 'PCA dimension reduction and GP surrogates for spanwise tip deflection as a function of the uncertain variables (the deterministic variables are fixed).\n\n');
    fprintf(fileID, 'Legend:\n');
    inputs_name = ["E for OOP (Pa)", "E for IP (Pa)", "Torsion modulus (Pa)"];     % list of the names of the uncertain variables
    N_variables = length(inputs_name);            % number of uncertain variables
    for ll = 1:N_variables
        fprintf(fileID, 'Uncertain variable %d: %s\n', ll, inputs_name(ll));
    end    
    fprintf(fileID, 'Deterministic variable 1: Airspeed; Case: %d out of %d; Value: %.5f (m/s)\n', ii, length(asp), asp(ii));
    fclose(fileID);
    disp('Readme file for the surrogates has been created successfully.');
end    

%% 9. Surrogates for tip deflections due to 1MC gust (uncertain variables: EI_1, EI_2, G)
% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [60, 80]*1e9; % (uncertain Young's modulus E for OOP bending) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [60, 80]*1e9; % (uncertain Young's modulus E for IP bending) lower and upper uncertainty bound
InputOpts.Marginals(3).Type = 'Uniform';
InputOpts.Marginals(3).Parameters = [26, 32]*1e9; % (uncertain Young's modulus G for torsion) lower and upper uncertainty bound
myInput = uq_createInput(InputOpts); 
asp = linspace(50,175,10);  %ranges of airspeeds...

set(0, 'DefaultFigureVisible', 'off');
surrogates_tip_deflection_gust = cell(length(asp), 1);

for ii = 1:length(asp)
    % Description of the physical model for UQLab
    ModelOpts.mFile = 'model_tip_deflection_gust';
    ModelOpts.isVectorized = false;
    ModelOpts.Parameters = [asp(ii)];
    myModel_tip_deflection_gust = uq_createModel(ModelOpts);

    N_train = 5;                                            % initial training set size (the set will be updated until the surrogate validation error is low enough)
    MetaOpts.Type = 'Metamodel';                            % 'metamodel': another word for 'surrogate'    
    MetaOpts.MetaType = 'PCE';
    MetaOpts.Input = myInput;                               % probability distribution for the uncertain variables
    MetaOpts.FullModel = myModel_tip_deflection_gust;     % the physical model as a UQLab object
    MetaOpts.ExpDesign.NSamples = N_train;                  % 'experimental design' (ExpDesign): another word for 'training set'
    if strcmp(MetaOpts.MetaType, 'Kriging')
        MetaOpts.ExpDesign.Sampling = 'User';
    end

    flag_parfor = true;             % can we run the physical model in parallel to build the training set? (True/False)
    seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
    N_train_increment = 8;          % we will increment the training set size until we reach convergence
    N_train_max = 500;              % training budget (i.e., maximum number of training points allowed)
    % run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
    % recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
    flag_test_for_mean_and_sigma = false;
    flag_test_set = true;          % will a test set be generated for further surrogate validation?

    % Plots generator for parameter sweeps for the uncertain variables
    inputs_name = ["E for OOP (Pa)", "E for IP (Pa)", "Torsion modulus (Pa)"];     % list of the names of the uncertain variables
    outputs_name = ["Max tip deflection (m)"];      % list of the names of the QIs
    N_outputs = length(outputs_name);           % number of quantities of interest (QIs)
    descriptive_title_for_plots = sprintf('%s surrogate (airspeed (m/s):%.2e)', MetaOpts.MetaType, asp(ii));
    N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
    plotsfolderName = 'tip_deflection_gust_uncertain_E'; 
    subfolder_plotsfolderName = sprintf('airspeed_case_%u', ii); % each 'case' refers to one fixed combination of design variables
    fullPath = fullfile(plotsfolderName, subfolder_plotsfolderName);
    mkdir(fullPath);
    mkdir(fullPath, 'plots_uq');
    try
        surrogates_tip_deflection_gust{ii, 1} =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, fullPath, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
        elementToSave = surrogates_tip_deflection_gust{ii, 1}; 
        save(fullfile(fullPath, sprintf('surrogate_tip_deflection_gust_airspeed_case_%u.mat', ii)), 'elementToSave'); % save the surrogate
        uncertain_variables_exploration(elementToSave, inputs_name, outputs_name, descriptive_title_for_plots, N_eval, seed, fullPath); % plots generator using the surrogates 
        % add readme to explain each 'case'
        fileID = fopen(fullfile(fullPath, 'readme.txt'), 'w');
        fprintf(fileID, 'Surrogates for each quantity of interest (QI) as a function of the uncertain variables (the deterministic variables are fixed).\n\n');
        fprintf(fileID, 'Legend:\n');
        N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
        N_variables = length(inputs_name);            % number of uncertain variables
        for ll = 1:N_outputs
            fprintf(fileID, 'QI %d: %s\n', ll, outputs_name(ll));
        end   
        for ll = 1:N_variables
            fprintf(fileID, 'Uncertain variable %d: %s\n', ll, inputs_name(ll));
        end    
        fprintf(fileID, 'Deterministic variable 1: Airspeed; Case: %d out of %d; Value: %.5f (m/s)\n', ii, length(asp), asp(ii));
        fprintf(fileID, 'Methodology: %s\n\n', descriptive_title_for_plots);
        fprintf(fileID, 'The trained surrogates and most of the design exploration figures are stored externally due to size limits.\n'); 
        fclose(fileID);
        disp('Readme file for the surrogates has been created successfully.');
    catch ME
        fprintf('Error in sample (%d): %s\n', ii, ME.message);
    end
end    
set(0, 'DefaultFigureVisible', 'on');
