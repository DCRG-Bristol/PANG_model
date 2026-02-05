%% 1. Start a UQLab session
uqlab 

%% 2. Surrogates for tip deflections at the first Hopf point (uncertain variables: EI_1, EI_2, G)
% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.9, 1.1]; % (uncertain scaling factor Young's modulus E for OOP bending) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.9, 1.1]; % (uncertain scaling factor Young's modulus E for IP bending) lower and upper uncertainty bound
InputOpts.Marginals(3).Type = 'Uniform';
InputOpts.Marginals(3).Parameters = [0.9, 1.1]; % (uncertain scaling factor modulus G for torsion) lower and upper uncertainty bound
myInput = uq_createInput(InputOpts); 
alp = ([0.3, 1.0, 1.2])*pi/180; %angles of attack to run...

set(0, 'DefaultFigureVisible', 'off');
surrogates_tip_deflection_hopf = cell(length(alp), 1);

for ii = 1:length(alp)
    % Description of the physical model for UQLab
    ModelOpts.mFile = 'model_tip_deflection_first_hopf';
    ModelOpts.isVectorized = false;
    ModelOpts.Parameters = [alp(ii)];
    myModel_tip_deflection_hopf = uq_createModel(ModelOpts);

    N_train = 5;                                            % initial training set size (the set will be updated until the surrogate validation error is low enough)
    MetaOpts.Type = 'Metamodel';                            % 'metamodel': another word for 'surrogate'    
    MetaOpts.MetaType = 'Kriging';
    MetaOpts.Input = myInput;                               % probability distribution for the uncertain variables
    MetaOpts.FullModel = myModel_tip_deflection_hopf;       % the physical model as a UQLab object
    MetaOpts.ExpDesign.NSamples = N_train;                  % 'experimental design' (ExpDesign): another word for 'training set'
    if strcmp(MetaOpts.MetaType, 'Kriging')
        MetaOpts.ExpDesign.Sampling = 'User';
    end

    flag_parfor = true;             % can we run the physical model in parallel to build the training set? (True/False)
    seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
    N_train_increment = 8;          % we will increment the training set size until we reach convergence
    N_train_max = 600;              % training budget (i.e., maximum number of training points allowed)
    % run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
    % recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
    flag_test_for_mean_and_sigma = false;
    flag_test_set = true;          % will a test set be generated for further surrogate validation?

    % Plots generator for parameter sweeps for the uncertain variables
    inputs_name = ["E for OOP (Pa)", "E for IP (Pa)", "Torsion modulus (Pa)"];     % list of the names of the uncertain variables
    outputs_name = ["Tip deflection first Hopf (m)"];      % list of the names of the QIs
    N_outputs = length(outputs_name);           % number of quantities of interest (QIs)
    descriptive_title_for_plots = sprintf('%s surrogate (angle of attack(rad):%.2e)', MetaOpts.MetaType, alp(ii));
    N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
    plotsfolderName = 'tip_deflection_first_hopf'; 
    subfolder_plotsfolderName = sprintf('aoa_case_%u', ii); % each 'case' refers to one fixed combination of design variables
    fullPath = fullfile(plotsfolderName, subfolder_plotsfolderName);
    mkdir(fullPath);
    mkdir(fullPath, 'plots_uq');
    try
        surrogates_tip_deflection_hopf{ii, 1} =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, fullPath, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
        elementToSave = surrogates_tip_deflection_hopf{ii, 1}; 
        save(fullfile(fullPath, sprintf('surrogate_tip_deflection_first_hopf_aoa_case_%u.mat', ii)), 'elementToSave'); % save the surrogate
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
        fprintf(fileID, 'Deterministic variable 1: angle of attack; Case: %d out of %d; Value: %.5f (rad)\n', ii, length(alp), alp(ii));
        fprintf(fileID, 'Methodology: %s\n\n', descriptive_title_for_plots);
        fprintf(fileID, 'The trained surrogates and most of the design exploration figures are stored externally due to size limits.\n'); 
        fclose(fileID);
        disp('Readme file for the surrogates has been created successfully.');
    catch ME
        fprintf('Error in sample (%d): %s\n', ii, ME.message);
    end
end   
set(0, 'DefaultFigureVisible', 'on');

%% 3. PCE surrogates for tip deflections at the first Hopf point (uncertain variables: EI_1, EI_2, G)

set(0, 'DefaultFigureVisible', 'off');
load('C:\Users\ua24606\git_2\PANG_model\SJD_wing\tip_deflection_first_hopf\aoa_case_1\training_set.mat');
% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.9, 1.1]; % (uncertain scaling factor Young's modulus E for OOP bending) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.9, 1.1]; % (uncertain scaling factor Young's modulus E for IP bending) lower and upper uncertainty bound
InputOpts.Marginals(3).Type = 'Uniform';
InputOpts.Marginals(3).Parameters = [0.9, 1.1]; % (uncertain scaling factor modulus G for torsion) lower and upper uncertainty bound
myInput = uq_createInput(InputOpts); 

MetaOpts.ExpDesign.X = X;                             % 'experimental design' (ExpDesign): another word for 'training set'
MetaOpts.ExpDesign.Y = Y;
MetaOpts.ExpDesign.NSamples = size(X, 1);
MetaOpts.Type = 'Metamodel';                            % 'metamodel': another word for 'surrogate'    
MetaOpts.MetaType = 'PCE';
MetaOpts.Input = myInput;                               % probability distribution for the uncertain variables
MetaOpts.Degree = 3:15;
if strcmp(MetaOpts.MetaType, 'Kriging')
    MetaOpts.ExpDesign.Sampling = 'User';
end
% MetaOpts.TruncOptions.qNorm = 0.75;

elementToSave = uq_createModel(MetaOpts);                % train the surrogate

% Plots generator for parameter sweeps for the uncertain variables
inputs_name = ["E for OOP (Pa)", "E for IP (Pa)", "Torsion modulus (Pa)"];     % list of the names of the uncertain variables
outputs_name = ["Tip deflection first Hopf (m)"];      % list of the names of the QIs
N_outputs = length(outputs_name);           % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate (angle of attack(rad):%.2e)', MetaOpts.MetaType, 0.3*pi/180);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
plotsfolderName = 'tip_deflection_first_hopf'; 
subfolder_plotsfolderName = 'aoa_case_1_pce'; % each 'case' refers to one fixed combination of design variables
fullPath = fullfile(plotsfolderName, subfolder_plotsfolderName);
mkdir(fullPath);
mkdir(fullPath, 'plots_uq');

save(fullfile(fullPath, 'surrogate_tip_deflection_first_hopf_aoa_case_1_pce.mat'), 'elementToSave'); % save the surrogate
uncertain_variables_exploration(elementToSave, inputs_name, outputs_name, descriptive_title_for_plots, N_eval, seed, fullPath); % plots generator using the surrogates 
set(0, 'DefaultFigureVisible', 'on');

%% 4. Surrogates for flutter speed (first Hopf point) (uncertain variables: E, G)
% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.9, 1.1]; % (uncertain scaling factor Young's modulus E) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.9, 1.1]; % (uncertain scaling factor modulus G for torsion) lower and upper uncertainty bound
myInput = uq_createInput(InputOpts); 
alp = ([0.3, 1.0, 1.2])*pi/180; %angles of attack to run...

set(0, 'DefaultFigureVisible', 'off');
surrogates_flutter_hopf = cell(length(alp), 1);

for ii = 3:length(alp)
    % Description of the physical model for UQLab
    ModelOpts.mFile = 'model_flutter_speed_first_hopf';
    ModelOpts.isVectorized = false;
    ModelOpts.Parameters = [alp(ii)];
    myModel_flutter_hopf = uq_createModel(ModelOpts);

    N_train = 5;                                            % initial training set size (the set will be updated until the surrogate validation error is low enough)
    MetaOpts.Type = 'Metamodel';                            % 'metamodel': another word for 'surrogate'    
    MetaOpts.MetaType = 'Kriging';
    MetaOpts.Input = myInput;                               % probability distribution for the uncertain variables
    MetaOpts.FullModel = myModel_flutter_hopf;       % the physical model as a UQLab object
    MetaOpts.ExpDesign.NSamples = N_train;                  % 'experimental design' (ExpDesign): another word for 'training set'
    if strcmp(MetaOpts.MetaType, 'Kriging')
        MetaOpts.ExpDesign.Sampling = 'User';
    end

    flag_parfor = true;             % can we run the physical model in parallel to build the training set? (True/False)
    seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
    N_train_increment = 8;          % we will increment the training set size until we reach convergence
    N_train_max = 1000;             % training budget (i.e., maximum number of training points allowed)
    % run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
    % recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
    flag_test_for_mean_and_sigma = false;
    flag_test_set = true;          % will a test set be generated for further surrogate validation?

    % Plots generator for parameter sweeps for the uncertain variables
    inputs_name = ["Young's modulus (Pa)", "Torsion modulus (Pa)"];     % list of the names of the uncertain variables
    outputs_name = ["Flutter speed (m/s)"];      % list of the names of the QIs
    N_outputs = length(outputs_name);           % number of quantities of interest (QIs)
    descriptive_title_for_plots = sprintf('%s surrogate (angle of attack(rad):%.2e)', MetaOpts.MetaType, alp(ii));
    N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
    plotsfolderName = 'flutter_first_hopf'; 
    subfolder_plotsfolderName = sprintf('aoa_case_%u', ii); % each 'case' refers to one fixed combination of design variables
    fullPath = fullfile(plotsfolderName, subfolder_plotsfolderName);
    mkdir(fullPath);
    mkdir(fullPath, 'plots_uq');
    try
        surrogates_flutter_hopf{ii, 1} =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, fullPath, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
        elementToSave = surrogates_flutter_hopf{ii, 1}; 
        save(fullfile(fullPath, sprintf('surrogate_flutter_first_hopf_aoa_case_%u.mat', ii)), 'elementToSave'); % save the surrogate
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
        fprintf(fileID, 'Deterministic variable 1: angle of attack; Case: %d out of %d; Value: %.5f (rad)\n', ii, length(alp), alp(ii));
        fprintf(fileID, 'Methodology: %s\n\n', descriptive_title_for_plots);
        fprintf(fileID, 'The trained surrogates and most of the design exploration figures are stored externally due to size limits.\n'); 
        fclose(fileID);
        disp('Readme file for the surrogates has been created successfully.');
    catch ME
        fprintf('Error in sample (%d): %s\n', ii, ME.message);
    end
end   
set(0, 'DefaultFigureVisible', 'on');

