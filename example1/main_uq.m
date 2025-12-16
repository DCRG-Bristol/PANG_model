%% 1. Start a UQLab session
uqlab 

%% 2. Surrogates for wind-off tip deflections (uncertain variable: E)
% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [60, 80]*1e9; % (uncertain Young's modulus E) lower and upper uncertainty bound
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
    N_train_max = 40;               % training budget (i.e., maximum number of training points allowed)
    % run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
    % recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
    flag_test_for_mean_and_sigma = false;
    flag_test_set = false;          % will a test set be generated for further surrogate validation?

    % Plots generator for parameter sweeps for the uncertain variables
    inputs_name = ["Young's modulus (Pa)"];     % list of the names of the uncertain variables
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

%% 3. Surrogates for the modal frequencies (uncertain variable: E)
% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [60, 80]*1e9; % (uncertain Young's modulus E) lower and upper uncertainty bound
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
    N_train_max = 40;               % training budget (i.e., maximum number of training points allowed)
    % run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
    % recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
    flag_test_for_mean_and_sigma = false;
    flag_test_set = false;          % will a test set be generated for further surrogate validation?

    % Plots generator for parameter sweeps for the uncertain variables
    inputs_name = ["Young's modulus (Pa)"];     % list of the names of the uncertain variables
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