%% Title section - UQ analysis and parameter sweeps for the numerical models corresponding to the experimental tests
%{
--------------------------------------------------------
Comments:
* For each model, fast sensitivity analysis is performed using Polynomial Chaos Expansion (PCE) surrogates
* The sensitivity analysis informs which variables should be included in the optimisation formulation for model update
* The resulting surrogates are very accurate, and hence they are also used for fast model update via least-squares error minimisation
--------------------------------------------------------
%}

%% 1. Start a UQLab session
uqlab 

%% 2. Surrogates and parameter sweeps for the Hopf bifurcations speeds as a function of lambda_L, alpha_L, EI, GJ, Sxx, Szz (W1.2 test)

% This section generates many plots, which are saved rather than displayed on the screen
set(0, 'DefaultFigureVisible', 'off'); 

% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for lambda_L) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for alpha_L) lower and upper uncertainty bound
% InputOpts.Marginals(3).Type = 'Uniform';
% InputOpts.Marginals(3).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for EI) lower and upper uncertainty bound
% InputOpts.Marginals(4).Type = 'Uniform';
% InputOpts.Marginals(4).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for GJ) lower and upper uncertainty bound
% InputOpts.Marginals(5).Type = 'Uniform';
% InputOpts.Marginals(5).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for Sxx) lower and upper uncertainty bound
% The uncertain variables are inputs for physical maps that output QIs
myInput = uq_createInput(InputOpts); 

plotsfolderName = 'Hopf_bifurcations_speeds_w1_2_uq'; 
mkdir(plotsfolderName)

% Description of the physical model for UQLab
ModelOpts.mFile = 'model_W1_2_Uf';
ModelOpts.isVectorized = false;
myModel = uq_createModel(ModelOpts);

N_train = 5;                                        % initial training set size (the set will be updated until the surrogate validation error is low enough)
MetaOpts.Type = 'Metamodel';                        % 'metamodel': another word for 'surrogate'
MetaOpts.MetaType = 'PCE';
MetaOpts.Input = myInput;                           % probability distribution for the uncertain variables
MetaOpts.FullModel = myModel;                       % the physical model as a UQLab object
MetaOpts.ExpDesign.NSamples = N_train;              % 'experimental design' (ExpDesign): another word for 'training set'
if strcmp(MetaOpts.MetaType, 'Kriging')
    MetaOpts.ExpDesign.Sampling = 'User';
end

flag_parfor = true;             % can we run the physical model in parallel to build the training set? (True/False)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
N_train_increment = 8;          % we will increment the training set size until we reach convergence
N_train_max = 50;               % training budget (i.e., maximum number of training points allowed)
% run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
% recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
flag_test_for_mean_and_sigma = false;
flag_test_set = true;           % will a test set be generated for further surrogate validation?
experimental_data_set = [];   

% Plots generator for parameter sweeps for the uncertain variables
% inputs_name = ["lambdaL scaling factor", "alphaL scaling factor", "EI scaling factor", "GJ scaling factor", "Sxx scaling factor"];  % list of the names of the uncertain variables
inputs_name = ["lambdaL scaling factor", "alphaL scaling factor"];        % list of the names of the uncertain variables
outputs_name = ["First Hopf speed(m/s)", "Second Hopf speed(m/s)"];       % list of the names of the QIs
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate', MetaOpts.MetaType);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
plotsfolderName = 'Hopf_bifurcations_speeds_w1_2_uq';  
mkdir(plotsfolderName, 'plots_uq');
tic;
surrogates =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, plotsfolderName, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
totalTime = toc;
fprintf('Total surrogate building time: %.4f seconds\n', totalTime);
elementToSave = surrogates;
save(fullfile(plotsfolderName, 'surrogates_Hopf_bifurcations_speeds_w1_2.mat'), 'elementToSave'); % save the surrogate
tic;
uncertain_variables_exploration(elementToSave, inputs_name, outputs_name, descriptive_title_for_plots, N_eval, seed, plotsfolderName, experimental_data_set); % plots generator using the surrogates 
totalTime = toc;
fprintf('Total design space exploration time: %.4f seconds\n', totalTime);

% add readme to explain each 'case' (i.e., each fixed combination (ii, jj) of deterministic variables)
fileID = fopen(fullfile(plotsfolderName, 'readme.txt'), 'w');
fprintf(fileID, 'Surrogates for each quantity of interest (QI) as a function of the uncertain variables.\n\n');
fprintf(fileID, 'Legend:\n');
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
N_variables = length(inputs_name);            % number of uncertain variables
for kk = 1:N_outputs
    fprintf(fileID, 'QI %d: %s\n', kk, outputs_name(kk));
end   
for kk = 1:N_variables
    fprintf(fileID, 'Uncertain variable %d: %s\n', kk, inputs_name(kk));
end    
fprintf(fileID, 'Methodology: %s\n\n', descriptive_title_for_plots); 
fprintf(fileID, 'The trained surrogates and most of the design exploration figures are stored externally due to size limits.\n'); 
fclose(fileID);
disp('Readme file for the surrogates has been created successfully.'); 

set(0, 'DefaultFigureVisible', 'on');

%% 2.1. Surrogates and parameter sweeps for the Hopf bifurcations speeds as a function of lambda_L, alpha_L, EI, GJ, Sxx, Szz (W1.2 test)

% This section generates many plots, which are saved rather than displayed on the screen
set(0, 'DefaultFigureVisible', 'off'); 

% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.95, 1.05]; % (multiplicative scaling factor for lambda_L) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.95, 1.05]; % (multiplicative scaling factor for alpha_L) lower and upper uncertainty bound
% InputOpts.Marginals(3).Type = 'Uniform';
% InputOpts.Marginals(3).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for EI) lower and upper uncertainty bound
% InputOpts.Marginals(4).Type = 'Uniform';
% InputOpts.Marginals(4).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for GJ) lower and upper uncertainty bound
% InputOpts.Marginals(5).Type = 'Uniform';
% InputOpts.Marginals(5).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for Sxx) lower and upper uncertainty bound
% The uncertain variables are inputs for physical maps that output QIs
myInput = uq_createInput(InputOpts); 

plotsfolderName = 'Hopf_bifurcations_speeds_5_perc_unc_w1_2_uq'; 
mkdir(plotsfolderName)

% Description of the physical model for UQLab
ModelOpts.mFile = 'model_W1_2_Uf';
ModelOpts.isVectorized = false;
myModel = uq_createModel(ModelOpts);

N_train = 5;                                        % initial training set size (the set will be updated until the surrogate validation error is low enough)
MetaOpts.Type = 'Metamodel';                        % 'metamodel': another word for 'surrogate'
MetaOpts.MetaType = 'PCE';
MetaOpts.Input = myInput;                           % probability distribution for the uncertain variables
MetaOpts.FullModel = myModel;                       % the physical model as a UQLab object
MetaOpts.ExpDesign.NSamples = N_train;              % 'experimental design' (ExpDesign): another word for 'training set'
if strcmp(MetaOpts.MetaType, 'Kriging')
    MetaOpts.ExpDesign.Sampling = 'User';
end

flag_parfor = true;             % can we run the physical model in parallel to build the training set? (True/False)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
N_train_increment = 8;          % we will increment the training set size until we reach convergence
N_train_max = 50;               % training budget (i.e., maximum number of training points allowed)
% run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
% recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
flag_test_for_mean_and_sigma = false;
flag_test_set = true;           % will a test set be generated for further surrogate validation?
experimental_data_set = [];   

% Plots generator for parameter sweeps for the uncertain variables
% inputs_name = ["lambdaL scaling factor", "alphaL scaling factor", "EI scaling factor", "GJ scaling factor", "Sxx scaling factor"];  % list of the names of the uncertain variables
inputs_name = ["lambdaL scaling factor", "alphaL scaling factor"];        % list of the names of the uncertain variables
outputs_name = ["First Hopf speed(m/s)", "Second Hopf speed(m/s)"];       % list of the names of the QIs
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate', MetaOpts.MetaType);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
plotsfolderName = 'Hopf_bifurcations_speeds_5_perc_unc_w1_2_uq';  
mkdir(plotsfolderName, 'plots_uq');
tic;
surrogates =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, plotsfolderName, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
totalTime = toc;
fprintf('Total surrogate building time: %.4f seconds\n', totalTime);
elementToSave = surrogates;
save(fullfile(plotsfolderName, 'surrogates_Hopf_bifurcations_speeds_w1_2.mat'), 'elementToSave'); % save the surrogate
tic;
uncertain_variables_exploration(elementToSave, inputs_name, outputs_name, descriptive_title_for_plots, N_eval, seed, plotsfolderName, experimental_data_set); % plots generator using the surrogates 
totalTime = toc;
fprintf('Total design space exploration time: %.4f seconds\n', totalTime);

% add readme to explain each 'case' (i.e., each fixed combination (ii, jj) of deterministic variables)
fileID = fopen(fullfile(plotsfolderName, 'readme.txt'), 'w');
fprintf(fileID, 'Surrogates for each quantity of interest (QI) as a function of the uncertain variables.\n\n');
fprintf(fileID, 'Legend:\n');
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
N_variables = length(inputs_name);            % number of uncertain variables
for kk = 1:N_outputs
    fprintf(fileID, 'QI %d: %s\n', kk, outputs_name(kk));
end   
for kk = 1:N_variables
    fprintf(fileID, 'Uncertain variable %d: %s\n', kk, inputs_name(kk));
end    
fprintf(fileID, 'Methodology: %s\n\n', descriptive_title_for_plots); 
fprintf(fileID, 'The trained surrogates and most of the design exploration figures are stored externally due to size limits.\n'); 
fclose(fileID);
disp('Readme file for the surrogates has been created successfully.'); 

set(0, 'DefaultFigureVisible', 'on');

%% 2.2. Surrogates and parameter sweeps for the Hopf bifurcations speeds as a function of lambda_L, alpha_L, EI, GJ, Sxx, Szz (W1.2 test)

% This section generates many plots, which are saved rather than displayed on the screen
set(0, 'DefaultFigureVisible', 'off'); 

% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.95, 1.05]; % (multiplicative scaling factor for lambda_L) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.95, 1.05]; % (multiplicative scaling factor for alpha_L) lower and upper uncertainty bound
% InputOpts.Marginals(3).Type = 'Uniform';
% InputOpts.Marginals(3).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for EI) lower and upper uncertainty bound
% InputOpts.Marginals(4).Type = 'Uniform';
% InputOpts.Marginals(4).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for GJ) lower and upper uncertainty bound
% InputOpts.Marginals(5).Type = 'Uniform';
% InputOpts.Marginals(5).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for Sxx) lower and upper uncertainty bound
% The uncertain variables are inputs for physical maps that output QIs
myInput = uq_createInput(InputOpts); 

plotsfolderName = 'Hopf_bifurcations_speeds_5_perc_unc_lse_optim_w1_2_uq'; 
mkdir(plotsfolderName)

% Description of the physical model for UQLab
ModelOpts.mFile = 'model_W1_2_Uf_lse_optim';
ModelOpts.isVectorized = false;
myModel = uq_createModel(ModelOpts);

N_train = 5;                                        % initial training set size (the set will be updated until the surrogate validation error is low enough)
MetaOpts.Type = 'Metamodel';                        % 'metamodel': another word for 'surrogate'
MetaOpts.MetaType = 'PCE';
MetaOpts.Input = myInput;                           % probability distribution for the uncertain variables
MetaOpts.FullModel = myModel;                       % the physical model as a UQLab object
MetaOpts.ExpDesign.NSamples = N_train;              % 'experimental design' (ExpDesign): another word for 'training set'
if strcmp(MetaOpts.MetaType, 'Kriging')
    MetaOpts.ExpDesign.Sampling = 'User';
end

flag_parfor = true;             % can we run the physical model in parallel to build the training set? (True/False)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
N_train_increment = 8;          % we will increment the training set size until we reach convergence
N_train_max = 50;               % training budget (i.e., maximum number of training points allowed)
% run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
% recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
flag_test_for_mean_and_sigma = false;
flag_test_set = true;           % will a test set be generated for further surrogate validation?
experimental_data_set = [];   

% Plots generator for parameter sweeps for the uncertain variables
% inputs_name = ["lambdaL scaling factor", "alphaL scaling factor", "EI scaling factor", "GJ scaling factor", "Sxx scaling factor"];  % list of the names of the uncertain variables
inputs_name = ["lambdaL scaling factor", "alphaL scaling factor"];        % list of the names of the uncertain variables
outputs_name = ["First Hopf speed(m/s)", "Second Hopf speed(m/s)"];       % list of the names of the QIs
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate', MetaOpts.MetaType);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
plotsfolderName = 'Hopf_bifurcations_speeds_5_perc_unc_lse_optim_w1_2_uq'; 
mkdir(plotsfolderName, 'plots_uq');
tic;
surrogates =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, plotsfolderName, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
totalTime = toc;
fprintf('Total surrogate building time: %.4f seconds\n', totalTime);
elementToSave = surrogates;
save(fullfile(plotsfolderName, 'surrogates_Hopf_bifurcations_speeds_w1_2.mat'), 'elementToSave'); % save the surrogate
tic;
uncertain_variables_exploration(elementToSave, inputs_name, outputs_name, descriptive_title_for_plots, N_eval, seed, plotsfolderName, experimental_data_set); % plots generator using the surrogates 
totalTime = toc;
fprintf('Total design space exploration time: %.4f seconds\n', totalTime);

% add readme to explain each 'case' (i.e., each fixed combination (ii, jj) of deterministic variables)
fileID = fopen(fullfile(plotsfolderName, 'readme.txt'), 'w');
fprintf(fileID, 'Surrogates for each quantity of interest (QI) as a function of the uncertain variables.\n\n');
fprintf(fileID, 'Legend:\n');
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
N_variables = length(inputs_name);            % number of uncertain variables
for kk = 1:N_outputs
    fprintf(fileID, 'QI %d: %s\n', kk, outputs_name(kk));
end   
for kk = 1:N_variables
    fprintf(fileID, 'Uncertain variable %d: %s\n', kk, inputs_name(kk));
end    
fprintf(fileID, 'Methodology: %s\n\n', descriptive_title_for_plots); 
fprintf(fileID, 'The trained surrogates and most of the design exploration figures are stored externally due to size limits.\n'); 
fclose(fileID);
disp('Readme file for the surrogates has been created successfully.'); 

set(0, 'DefaultFigureVisible', 'on');

%% 2.4. Surrogates and parameter sweeps for the Hopf bifurcations speeds as a function of lambda_L, alpha_L, EI, GJ, Sxx, Szz (W1.2 test)

% This section generates many plots, which are saved rather than displayed on the screen
set(0, 'DefaultFigureVisible', 'off'); 

% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.95, 1.05]; % (multiplicative scaling factor for lambda_L) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.95, 1.05]; % (multiplicative scaling factor for alpha_L) lower and upper uncertainty bound
InputOpts.Marginals(3).Type = 'Uniform';
InputOpts.Marginals(3).Parameters = [0.95, 1.05]; % (multiplicative scaling factor for EI) lower and upper uncertainty bound
InputOpts.Marginals(4).Type = 'Uniform';
InputOpts.Marginals(4).Parameters = [0.95, 1.05]; % (multiplicative scaling factor for GJ) lower and upper uncertainty bound
InputOpts.Marginals(5).Type = 'Uniform';
InputOpts.Marginals(5).Parameters = [0.95, 1.05]; % (multiplicative scaling factor for Sxx) lower and upper uncertainty bound
% The uncertain variables are inputs for physical maps that output QIs
myInput = uq_createInput(InputOpts); 

plotsfolderName = 'Hopf_bifurcations_speeds_all_unc_w1_2_uq'; 
mkdir(plotsfolderName)

% Description of the physical model for UQLab
ModelOpts.mFile = 'model_W1_2_Uf_all_unc';
ModelOpts.isVectorized = false;
myModel = uq_createModel(ModelOpts);

N_train = 5;                                        % initial training set size (the set will be updated until the surrogate validation error is low enough)
MetaOpts.Type = 'Metamodel';                        % 'metamodel': another word for 'surrogate'
MetaOpts.MetaType = 'PCE';
MetaOpts.Input = myInput;                           % probability distribution for the uncertain variables
MetaOpts.FullModel = myModel;                       % the physical model as a UQLab object
MetaOpts.ExpDesign.NSamples = N_train;              % 'experimental design' (ExpDesign): another word for 'training set'
if strcmp(MetaOpts.MetaType, 'Kriging')
    MetaOpts.ExpDesign.Sampling = 'User';
end

flag_parfor = true;             % can we run the physical model in parallel to build the training set? (True/False)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
N_train_increment = 8;          % we will increment the training set size until we reach convergence
N_train_max = 50;               % training budget (i.e., maximum number of training points allowed)
% run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
% recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
flag_test_for_mean_and_sigma = false;
flag_test_set = true;           % will a test set be generated for further surrogate validation?
experimental_data_set = [];   

% Plots generator for parameter sweeps for the uncertain variables
inputs_name = ["lambdaL scaling factor", "alphaL scaling factor", "EI scaling factor", "GJ scaling factor", "Sxx scaling factor"];  % list of the names of the uncertain variables
% inputs_name = ["lambdaL scaling factor", "alphaL scaling factor"];        % list of the names of the uncertain variables
outputs_name = ["First Hopf speed(m/s)", "Second Hopf speed(m/s)"];       % list of the names of the QIs
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate', MetaOpts.MetaType);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
plotsfolderName = 'Hopf_bifurcations_speeds_all_unc_w1_2_uq';  
mkdir(plotsfolderName, 'plots_uq');
tic;
surrogates =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, plotsfolderName, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
totalTime = toc;
fprintf('Total surrogate building time: %.4f seconds\n', totalTime);
elementToSave = surrogates;
save(fullfile(plotsfolderName, 'surrogates_Hopf_bifurcations_speeds_w1_2.mat'), 'elementToSave'); % save the surrogate
tic;
uncertain_variables_exploration(elementToSave, inputs_name, outputs_name, descriptive_title_for_plots, N_eval, seed, plotsfolderName, experimental_data_set); % plots generator using the surrogates 
totalTime = toc;
fprintf('Total design space exploration time: %.4f seconds\n', totalTime);

% add readme to explain each 'case' (i.e., each fixed combination (ii, jj) of deterministic variables)
fileID = fopen(fullfile(plotsfolderName, 'readme.txt'), 'w');
fprintf(fileID, 'Surrogates for each quantity of interest (QI) as a function of the uncertain variables.\n\n');
fprintf(fileID, 'Legend:\n');
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
N_variables = length(inputs_name);            % number of uncertain variables
for kk = 1:N_outputs
    fprintf(fileID, 'QI %d: %s\n', kk, outputs_name(kk));
end   
for kk = 1:N_variables
    fprintf(fileID, 'Uncertain variable %d: %s\n', kk, inputs_name(kk));
end    
fprintf(fileID, 'Methodology: %s\n\n', descriptive_title_for_plots); 
fprintf(fileID, 'The trained surrogates and most of the design exploration figures are stored externally due to size limits.\n'); 
fclose(fileID);
disp('Readme file for the surrogates has been created successfully.'); 

set(0, 'DefaultFigureVisible', 'on');

%% 2.5. Surrogates and parameter sweeps for the Hopf bifurcations speeds as a function of lambda_L, alpha_L, EI, GJ, Sxx, Szz (W1.2 test)

% This section generates many plots, which are saved rather than displayed on the screen
set(0, 'DefaultFigureVisible', 'off'); 

% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.15, 0.33]; % (lambda_L) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.35, 0.55]; % (alpha_L) lower and upper uncertainty bound
InputOpts.Marginals(3).Type = 'Uniform';
InputOpts.Marginals(3).Parameters = [0.95, 1.05]; % (multiplicative scaling factor for EI) lower and upper uncertainty bound
InputOpts.Marginals(4).Type = 'Uniform';
InputOpts.Marginals(4).Parameters = [0.95, 1.05]; % (multiplicative scaling factor for GJ) lower and upper uncertainty bound
InputOpts.Marginals(5).Type = 'Uniform';
InputOpts.Marginals(5).Parameters = [0.95, 1.05]; % (multiplicative scaling factor for Sxx) lower and upper uncertainty bound
% The uncertain variables are inputs for physical maps that output QIs
myInput = uq_createInput(InputOpts); 

plotsfolderName = 'Hopf_bifurcations_speeds_all_unc_large_alpha_lambda_w1_2_uq'; 
mkdir(plotsfolderName)

% Description of the physical model for UQLab
ModelOpts.mFile = 'model_W1_2_Uf_all_unc';
ModelOpts.isVectorized = false;
myModel = uq_createModel(ModelOpts);

N_train = 5;                                        % initial training set size (the set will be updated until the surrogate validation error is low enough)
MetaOpts.Type = 'Metamodel';                        % 'metamodel': another word for 'surrogate'
MetaOpts.MetaType = 'PCE';
MetaOpts.Input = myInput;                           % probability distribution for the uncertain variables
MetaOpts.FullModel = myModel;                       % the physical model as a UQLab object
MetaOpts.ExpDesign.NSamples = N_train;              % 'experimental design' (ExpDesign): another word for 'training set'
if strcmp(MetaOpts.MetaType, 'Kriging')
    MetaOpts.ExpDesign.Sampling = 'User';
end

flag_parfor = true;             % can we run the physical model in parallel to build the training set? (True/False)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
N_train_increment = 8;          % we will increment the training set size until we reach convergence
N_train_max = 50;               % training budget (i.e., maximum number of training points allowed)
% run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
% recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
flag_test_for_mean_and_sigma = false;
flag_test_set = true;           % will a test set be generated for further surrogate validation?
experimental_data_set = [];   

% Plots generator for parameter sweeps for the uncertain variables
inputs_name = ["lambdaL scaling factor", "alphaL scaling factor", "EI scaling factor", "GJ scaling factor", "Sxx scaling factor"];  % list of the names of the uncertain variables
% inputs_name = ["lambdaL scaling factor", "alphaL scaling factor"];        % list of the names of the uncertain variables
outputs_name = ["First Hopf speed(m/s)", "Second Hopf speed(m/s)"];       % list of the names of the QIs
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate', MetaOpts.MetaType);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
plotsfolderName = 'Hopf_bifurcations_speeds_all_unc_large_alpha_lambda_w1_2_uq';  
mkdir(plotsfolderName, 'plots_uq');
tic;
surrogates =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, plotsfolderName, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
totalTime = toc;
fprintf('Total surrogate building time: %.4f seconds\n', totalTime);
elementToSave = surrogates;
save(fullfile(plotsfolderName, 'surrogates_Hopf_bifurcations_speeds_w1_2.mat'), 'elementToSave'); % save the surrogate
tic;
uncertain_variables_exploration(elementToSave, inputs_name, outputs_name, descriptive_title_for_plots, N_eval, seed, plotsfolderName, experimental_data_set); % plots generator using the surrogates 
totalTime = toc;
fprintf('Total design space exploration time: %.4f seconds\n', totalTime);

% add readme to explain each 'case' (i.e., each fixed combination (ii, jj) of deterministic variables)
fileID = fopen(fullfile(plotsfolderName, 'readme.txt'), 'w');
fprintf(fileID, 'Surrogates for each quantity of interest (QI) as a function of the uncertain variables.\n\n');
fprintf(fileID, 'Legend:\n');
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
N_variables = length(inputs_name);            % number of uncertain variables
for kk = 1:N_outputs
    fprintf(fileID, 'QI %d: %s\n', kk, outputs_name(kk));
end   
for kk = 1:N_variables
    fprintf(fileID, 'Uncertain variable %d: %s\n', kk, inputs_name(kk));
end    
fprintf(fileID, 'Methodology: %s\n\n', descriptive_title_for_plots); 
fprintf(fileID, 'The trained surrogates and most of the design exploration figures are stored externally due to size limits.\n'); 
fclose(fileID);
disp('Readme file for the surrogates has been created successfully.'); 

set(0, 'DefaultFigureVisible', 'on');

%% 3. Surrogates and parameter sweeps for the damping ratio as a function of lambda_L, alpha_L, EI, GJ, Sxx, Szz (W1.2 test)

% This section generates many plots, which are saved rather than displayed on the screen
set(0, 'DefaultFigureVisible', 'off'); 

% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.95, 1.05]; % (multiplicative scaling factor for lambda_L) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.95, 1.05]; % (multiplicative scaling factor for alpha_L) lower and upper uncertainty bound
% InputOpts.Marginals(3).Type = 'Uniform';
% InputOpts.Marginals(3).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for EI) lower and upper uncertainty bound
% InputOpts.Marginals(4).Type = 'Uniform';
% InputOpts.Marginals(4).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for GJ) lower and upper uncertainty bound
% InputOpts.Marginals(5).Type = 'Uniform';
% InputOpts.Marginals(5).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for Sxx) lower and upper uncertainty bound
% The uncertain variables are inputs for physical maps that output QIs
myInput = uq_createInput(InputOpts); 

plotsfolderName = 'tor1_damping_w2_uq'; 
mkdir(plotsfolderName)

% Description of the physical model for UQLab
ModelOpts.mFile = 'model_W2_tor1_damping';
ModelOpts.isVectorized = false;
myModel = uq_createModel(ModelOpts);

N_train = 5;                                        % initial training set size (the set will be updated until the surrogate validation error is low enough)
MetaOpts.Type = 'Metamodel';                        % 'metamodel': another word for 'surrogate'
MetaOpts.MetaType = 'PCE';
MetaOpts.Input = myInput;                           % probability distribution for the uncertain variables
MetaOpts.FullModel = myModel;                       % the physical model as a UQLab object
MetaOpts.ExpDesign.NSamples = N_train;              % 'experimental design' (ExpDesign): another word for 'training set'
if strcmp(MetaOpts.MetaType, 'Kriging')
    MetaOpts.ExpDesign.Sampling = 'User';
end

flag_parfor = true;             % can we run the physical model in parallel to build the training set? (True/False)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
N_train_increment = 8;          % we will increment the training set size until we reach convergence
N_train_max = 50;               % training budget (i.e., maximum number of training points allowed)
% run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
% recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
flag_test_for_mean_and_sigma = false;
flag_test_set = true;           % will a test set be generated for further surrogate validation?
experimental_data_set = [];   

% Plots generator for parameter sweeps for the uncertain variables
% inputs_name = ["lambdaL scaling factor", "alphaL scaling factor", "EI scaling factor", "GJ scaling factor", "Sxx scaling factor"];  % list of the names of the uncertain variables
inputs_name = ["lambdaL scaling factor", "alphaL scaling factor"];        % list of the names of the uncertain variables
outputs_name = ["Tor1 damping ratio (speed1)", "Tor1 damping ratio (speed2)", "Tor1 damping ratio (speed3)", "Tor1 damping ratio (speed4)", "Tor1 damping ratio (speed5)", "Tor1 damping ratio (speed6)", "Tor1 damping ratio (speed7)", "Tor1 damping ratio (speed8)", "Tor1 damping ratio (speed9)"];       % list of the names of the QIs
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate', MetaOpts.MetaType);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
plotsfolderName = 'tor1_damping_w2_uq'; 
mkdir(plotsfolderName, 'plots_uq');
tic;
surrogates =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, plotsfolderName, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
totalTime = toc;
fprintf('Total surrogate building time: %.4f seconds\n', totalTime);
elementToSave = surrogates;
save(fullfile(plotsfolderName, 'surrogates_tor1_damping.mat'), 'elementToSave'); % save the surrogate
tic;
uncertain_variables_exploration(elementToSave, inputs_name, outputs_name, descriptive_title_for_plots, N_eval, seed, plotsfolderName, experimental_data_set); % plots generator using the surrogates 
totalTime = toc;
fprintf('Total design space exploration time: %.4f seconds\n', totalTime);

% add readme to explain each 'case' (i.e., each fixed combination (ii, jj) of deterministic variables)
fileID = fopen(fullfile(plotsfolderName, 'readme.txt'), 'w');
fprintf(fileID, 'Surrogates for each quantity of interest (QI) as a function of the uncertain variables.\n\n');
fprintf(fileID, 'Legend:\n');
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
N_variables = length(inputs_name);            % number of uncertain variables
for kk = 1:N_outputs
    fprintf(fileID, 'QI %d: %s\n', kk, outputs_name(kk));
end   
for kk = 1:N_variables
    fprintf(fileID, 'Uncertain variable %d: %s\n', kk, inputs_name(kk));
end    
fprintf(fileID, 'Methodology: %s\n\n', descriptive_title_for_plots); 
fprintf(fileID, 'The trained surrogates and most of the design exploration figures are stored externally due to size limits.\n'); 
fclose(fileID);
disp('Readme file for the surrogates has been created successfully.'); 

set(0, 'DefaultFigureVisible', 'on');

%% 3.1. Surrogates and parameter sweeps for the damping ratio as a function of lambda_L, alpha_L, EI, GJ, Sxx, Szz (W1.2 test)

% This section generates many plots, which are saved rather than displayed on the screen
set(0, 'DefaultFigureVisible', 'off'); 

% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for lambda_L) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for alpha_L) lower and upper uncertainty bound
% InputOpts.Marginals(3).Type = 'Uniform';
% InputOpts.Marginals(3).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for EI) lower and upper uncertainty bound
% InputOpts.Marginals(4).Type = 'Uniform';
% InputOpts.Marginals(4).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for GJ) lower and upper uncertainty bound
% InputOpts.Marginals(5).Type = 'Uniform';
% InputOpts.Marginals(5).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for Sxx) lower and upper uncertainty bound
% The uncertain variables are inputs for physical maps that output QIs
myInput = uq_createInput(InputOpts); 

plotsfolderName = 'tor1_damping_20_perc_unc_w2_uq'; 
mkdir(plotsfolderName)

% Description of the physical model for UQLab
ModelOpts.mFile = 'model_W2_tor1_damping';
ModelOpts.isVectorized = false;
myModel = uq_createModel(ModelOpts);

N_train = 5;                                        % initial training set size (the set will be updated until the surrogate validation error is low enough)
MetaOpts.Type = 'Metamodel';                        % 'metamodel': another word for 'surrogate'
MetaOpts.MetaType = 'PCE';
MetaOpts.Input = myInput;                           % probability distribution for the uncertain variables
MetaOpts.FullModel = myModel;                       % the physical model as a UQLab object
MetaOpts.ExpDesign.NSamples = N_train;              % 'experimental design' (ExpDesign): another word for 'training set'
if strcmp(MetaOpts.MetaType, 'Kriging')
    MetaOpts.ExpDesign.Sampling = 'User';
end

flag_parfor = true;             % can we run the physical model in parallel to build the training set? (True/False)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
N_train_increment = 8;          % we will increment the training set size until we reach convergence
N_train_max = 50;               % training budget (i.e., maximum number of training points allowed)
% run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
% recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
flag_test_for_mean_and_sigma = false;
flag_test_set = true;           % will a test set be generated for further surrogate validation?
experimental_data_set = [];   

% Plots generator for parameter sweeps for the uncertain variables
% inputs_name = ["lambdaL scaling factor", "alphaL scaling factor", "EI scaling factor", "GJ scaling factor", "Sxx scaling factor"];  % list of the names of the uncertain variables
inputs_name = ["lambdaL scaling factor", "alphaL scaling factor"];        % list of the names of the uncertain variables
outputs_name = ["Tor1 damping ratio (speed1)", "Tor1 damping ratio (speed2)", "Tor1 damping ratio (speed3)", "Tor1 damping ratio (speed4)", "Tor1 damping ratio (speed5)", "Tor1 damping ratio (speed6)", "Tor1 damping ratio (speed7)", "Tor1 damping ratio (speed8)"];       % list of the names of the QIs
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate', MetaOpts.MetaType);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
plotsfolderName = 'tor1_damping_20_perc_unc_w2_uq'; 
mkdir(plotsfolderName, 'plots_uq');
tic;
surrogates =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, plotsfolderName, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
totalTime = toc;
fprintf('Total surrogate building time: %.4f seconds\n', totalTime);
elementToSave = surrogates;
save(fullfile(plotsfolderName, 'surrogates_tor1_damping.mat'), 'elementToSave'); % save the surrogate
tic;
uncertain_variables_exploration(elementToSave, inputs_name, outputs_name, descriptive_title_for_plots, N_eval, seed, plotsfolderName, experimental_data_set); % plots generator using the surrogates 
totalTime = toc;
fprintf('Total design space exploration time: %.4f seconds\n', totalTime);

% add readme to explain each 'case' (i.e., each fixed combination (ii, jj) of deterministic variables)
fileID = fopen(fullfile(plotsfolderName, 'readme.txt'), 'w');
fprintf(fileID, 'Surrogates for each quantity of interest (QI) as a function of the uncertain variables.\n\n');
fprintf(fileID, 'Legend:\n');
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
N_variables = length(inputs_name);            % number of uncertain variables
for kk = 1:N_outputs
    fprintf(fileID, 'QI %d: %s\n', kk, outputs_name(kk));
end   
for kk = 1:N_variables
    fprintf(fileID, 'Uncertain variable %d: %s\n', kk, inputs_name(kk));
end    
fprintf(fileID, 'Methodology: %s\n\n', descriptive_title_for_plots); 
fprintf(fileID, 'The trained surrogates and most of the design exploration figures are stored externally due to size limits.\n'); 
fclose(fileID);
disp('Readme file for the surrogates has been created successfully.'); 

set(0, 'DefaultFigureVisible', 'on');