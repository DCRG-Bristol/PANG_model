%% Title section - UQ analysis and parameter sweeps ...
%{
--------------------------------------------------------
Comments:
* ...
--------------------------------------------------------
%}

%% 1. Start a UQLab session
uqlab 

%% 2. Surrogates and parameter sweeps for the tip deflection (LE) as a function of EI, GJ (G1.1 test)

% This section generates many plots, which are saved rather than displayed on the screen
set(0, 'DefaultFigureVisible', 'off'); 

% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for EI) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for GJ) lower and upper uncertainty bound
% The uncertain variables are inputs for physical maps that output QIs
myInput = uq_createInput(InputOpts); 

plotsfolderName = 'tip_deflection_LE_g1_1_uq'; 
mkdir(plotsfolderName)

% Description of the physical model for UQLab
ModelOpts.mFile = 'model_G1_1_tip_deflection_LE';
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

flag_parfor = false;            % can we run the physical model in parallel to build the training set? (True/False)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
N_train_increment = 8;          % we will increment the training set size until we reach convergence
N_train_max = 50;               % training budget (i.e., maximum number of training points allowed)
% run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
% recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
flag_test_for_mean_and_sigma = false;
flag_test_set = true;           % will a test set be generated for further surrogate validation?
load('groundTests\testData\SJD_groundTestData.mat', 'exprData');
experimental_data_set = exprData{1, 1}.delta_LE;   

% Plots generator for parameter sweeps for the uncertain variables
inputs_name = ["EI scaling factor", "GJ scaling factor"];  % list of the names of the uncertain variables
outputs_name = ["Tip deflection(m) at 0deg", "Tip deflection(m) at 10deg", "Tip deflection(m) at 20deg", "Tip deflection(m) at 30deg", "Tip deflection(m) at 60deg"];       % list of the names of the QIs
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate', MetaOpts.MetaType);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
plotsfolderName = 'tip_deflection_LE_g1_1_uq';  
mkdir(plotsfolderName, 'plots_uq');
tic;
surrogates =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, plotsfolderName, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
totalTime = toc;
fprintf('Total surrogate building time: %.4f seconds\n', totalTime);
elementToSave = surrogates;
save(fullfile(plotsfolderName, 'surrogates_tip_deflection_LE_g1_1.mat'), 'elementToSave'); % save the surrogate
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

%% 3. Surrogates and parameter sweeps for the tip deflection (LE) as a function of EI, GJ (G1.2 test)

% This section generates many plots, which are saved rather than displayed on the screen
set(0, 'DefaultFigureVisible', 'off'); 

% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for EI) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for GJ) lower and upper uncertainty bound
% The uncertain variables are inputs for physical maps that output QIs
myInput = uq_createInput(InputOpts); 

plotsfolderName = 'tip_deflection_LE_g1_2_uq'; 
mkdir(plotsfolderName)

% Description of the physical model for UQLab
ModelOpts.mFile = 'model_G1_2_tip_deflection_LE';
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

flag_parfor = false;            % can we run the physical model in parallel to build the training set? (True/False)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
N_train_increment = 8;          % we will increment the training set size until we reach convergence
N_train_max = 50;               % training budget (i.e., maximum number of training points allowed)
% run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
% recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
flag_test_for_mean_and_sigma = false;
flag_test_set = true;           % will a test set be generated for further surrogate validation?
load('groundTests\testData\SJD_groundTestData.mat', 'exprData');
experimental_data_set = exprData{1, 2}.delta_LE;   

% Plots generator for parameter sweeps for the uncertain variables
inputs_name = ["EI scaling factor", "GJ scaling factor"];  % list of the names of the uncertain variables
outputs_name = ["Tip deflection(m) at 0deg", "Tip deflection(m) at 10deg", "Tip deflection(m) at 20deg", "Tip deflection(m) at 30deg", "Tip deflection(m) at 60deg"];       % list of the names of the QIs
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate', MetaOpts.MetaType);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
plotsfolderName = 'tip_deflection_LE_g1_2_uq';  
mkdir(plotsfolderName, 'plots_uq');
tic;
surrogates =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, plotsfolderName, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
totalTime = toc;
fprintf('Total surrogate building time: %.4f seconds\n', totalTime);
elementToSave = surrogates;
save(fullfile(plotsfolderName, 'surrogates_tip_deflection_LE_g1_2.mat'), 'elementToSave'); % save the surrogate
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

%% 4. Surrogates and parameter sweeps for the bending strain (beta_y) as a function of EI, GJ (G1.1 test)

% This section generates many plots, which are saved rather than displayed on the screen
set(0, 'DefaultFigureVisible', 'off'); 

% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for EI) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for GJ) lower and upper uncertainty bound
% The uncertain variables are inputs for physical maps that output QIs
myInput = uq_createInput(InputOpts); 

plotsfolderName = 'bending_strain_g1_1_uq'; 
mkdir(plotsfolderName)

% Description of the physical model for UQLab
ModelOpts.mFile = 'model_G1_1_beta_y';
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

flag_parfor = false;            % can we run the physical model in parallel to build the training set? (True/False)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
N_train_increment = 8;          % we will increment the training set size until we reach convergence
N_train_max = 50;               % training budget (i.e., maximum number of training points allowed)
% run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
% recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
flag_test_for_mean_and_sigma = false;
flag_test_set = true;           % will a test set be generated for further surrogate validation?
load('groundTests\testData\SJD_groundTestData.mat', 'exprData');
experimental_data_set = exprData{1, 1}.beta_y;   

% Plots generator for parameter sweeps for the uncertain variables
inputs_name = ["EI scaling factor", "GJ scaling factor"];  % list of the names of the uncertain variables
outputs_name = ["Bend.strain(1/m) at 0deg", "Bend.strain(1/m) at 10deg", "Bend.strain(1/m) at 20deg", "Bend.strain(1/m) at 30deg", "Bend.strain(1/m) at 60deg"];       % list of the names of the QIs
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate', MetaOpts.MetaType);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
plotsfolderName = 'bending_strain_g1_1_uq';   
mkdir(plotsfolderName, 'plots_uq');
tic;
surrogates =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, plotsfolderName, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
totalTime = toc;
fprintf('Total surrogate building time: %.4f seconds\n', totalTime);
elementToSave = surrogates;
save(fullfile(plotsfolderName, 'surrogates_bending_strain_g1_1.mat'), 'elementToSave'); % save the surrogate
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

%% 5. Surrogates and parameter sweeps for the torsional strain (beta_x) as a function of EI, GJ (G1.1 test)

% This section generates many plots, which are saved rather than displayed on the screen
set(0, 'DefaultFigureVisible', 'off'); 

% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for EI) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for GJ) lower and upper uncertainty bound
% The uncertain variables are inputs for physical maps that output QIs
myInput = uq_createInput(InputOpts); 

plotsfolderName = 'torsional_strain_g1_1_uq'; 
mkdir(plotsfolderName)

% Description of the physical model for UQLab
ModelOpts.mFile = 'model_G1_1_beta_x';
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

flag_parfor = false;            % can we run the physical model in parallel to build the training set? (True/False)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
N_train_increment = 8;          % we will increment the training set size until we reach convergence
N_train_max = 50;               % training budget (i.e., maximum number of training points allowed)
% run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
% recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
flag_test_for_mean_and_sigma = false;
flag_test_set = true;           % will a test set be generated for further surrogate validation?
load('groundTests\testData\SJD_groundTestData.mat', 'exprData');
experimental_data_set = exprData{1, 1}.beta_x;   

% Plots generator for parameter sweeps for the uncertain variables
inputs_name = ["EI scaling factor", "GJ scaling factor"];  % list of the names of the uncertain variables
outputs_name = ["Tors.strain(1/m) at 0deg", "Tors.strain(1/m) at 10deg", "Tors.strain(1/m) at 20deg", "Tors.strain(1/m) at 30deg", "Tors.strain(1/m) at 60deg"];       % list of the names of the QIs
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate', MetaOpts.MetaType);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
plotsfolderName = 'torsional_strain_g1_1_uq';   
mkdir(plotsfolderName, 'plots_uq');
tic;
surrogates =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, plotsfolderName, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
totalTime = toc;
fprintf('Total surrogate building time: %.4f seconds\n', totalTime);
elementToSave = surrogates;
save(fullfile(plotsfolderName, 'surrogates_torsional_strain_g1_1.mat'), 'elementToSave'); % save the surrogate
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

%% 6. Surrogates and parameter sweeps for the bending strain (beta_y) as a function of EI, GJ (G1.2 test)

% This section generates many plots, which are saved rather than displayed on the screen
set(0, 'DefaultFigureVisible', 'off'); 

% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for EI) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for GJ) lower and upper uncertainty bound
% The uncertain variables are inputs for physical maps that output QIs
myInput = uq_createInput(InputOpts); 

plotsfolderName = 'bending_strain_g1_2_uq'; 
mkdir(plotsfolderName)

% Description of the physical model for UQLab
ModelOpts.mFile = 'model_G1_2_beta_y';
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

flag_parfor = false;            % can we run the physical model in parallel to build the training set? (True/False)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
N_train_increment = 8;          % we will increment the training set size until we reach convergence
N_train_max = 50;               % training budget (i.e., maximum number of training points allowed)
% run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
% recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
flag_test_for_mean_and_sigma = false;
flag_test_set = true;           % will a test set be generated for further surrogate validation?
load('groundTests\testData\SJD_groundTestData.mat', 'exprData');
experimental_data_set = exprData{1, 2}.beta_y;   

% Plots generator for parameter sweeps for the uncertain variables
inputs_name = ["EI scaling factor", "GJ scaling factor"];  % list of the names of the uncertain variables
outputs_name = ["Bend.strain(1/m) at 0deg", "Bend.strain(1/m) at 10deg", "Bend.strain(1/m) at 20deg", "Bend.strain(1/m) at 30deg", "Bend.strain(1/m) at 60deg"];       % list of the names of the QIs
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate', MetaOpts.MetaType);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
plotsfolderName = 'bending_strain_g1_2_uq';   
mkdir(plotsfolderName, 'plots_uq');
tic;
surrogates =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, plotsfolderName, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
totalTime = toc;
fprintf('Total surrogate building time: %.4f seconds\n', totalTime);
elementToSave = surrogates;
save(fullfile(plotsfolderName, 'surrogates_bending_strain_g1_2.mat'), 'elementToSave'); % save the surrogate
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

%% 7. Surrogates and parameter sweeps for the torsional strain (beta_x) as a function of EI, GJ (G1.2 test)

% This section generates many plots, which are saved rather than displayed on the screen
set(0, 'DefaultFigureVisible', 'off'); 

% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for EI) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for GJ) lower and upper uncertainty bound
% The uncertain variables are inputs for physical maps that output QIs
myInput = uq_createInput(InputOpts); 

plotsfolderName = 'torsional_strain_g1_2_uq'; 
mkdir(plotsfolderName)

% Description of the physical model for UQLab
ModelOpts.mFile = 'model_G1_2_beta_x';
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

flag_parfor = false;            % can we run the physical model in parallel to build the training set? (True/False)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
N_train_increment = 8;          % we will increment the training set size until we reach convergence
N_train_max = 50;               % training budget (i.e., maximum number of training points allowed)
% run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
% recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
flag_test_for_mean_and_sigma = false;
flag_test_set = true;           % will a test set be generated for further surrogate validation?
load('groundTests\testData\SJD_groundTestData.mat', 'exprData');
experimental_data_set = exprData{1, 2}.beta_x;   

% Plots generator for parameter sweeps for the uncertain variables
inputs_name = ["EI scaling factor", "GJ scaling factor"];  % list of the names of the uncertain variables
outputs_name = ["Tors.strain(1/m) at 0deg", "Tors.strain(1/m) at 10deg", "Tors.strain(1/m) at 20deg", "Tors.strain(1/m) at 30deg", "Tors.strain(1/m) at 60deg"];       % list of the names of the QIs
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate', MetaOpts.MetaType);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
plotsfolderName = 'torsional_strain_g1_2_uq';   
mkdir(plotsfolderName, 'plots_uq');
tic;
surrogates =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, plotsfolderName, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
totalTime = toc;
fprintf('Total surrogate building time: %.4f seconds\n', totalTime);
elementToSave = surrogates;
save(fullfile(plotsfolderName, 'surrogates_torsional_strain_g1_2.mat'), 'elementToSave'); % save the surrogate
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

%% 8. Surrogates and parameter sweeps for the 1st modal frequency (OOP bending) as a function of EI, GJ, Sxx, Szz (G2 test)

% This section generates many plots, which are saved rather than displayed on the screen
set(0, 'DefaultFigureVisible', 'off'); 

% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for EI) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for GJ) lower and upper uncertainty bound
InputOpts.Marginals(3).Type = 'Uniform';
InputOpts.Marginals(3).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for Sxx) lower and upper uncertainty bound
InputOpts.Marginals(4).Type = 'Uniform';
InputOpts.Marginals(4).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for Szz) lower and upper uncertainty bound
% The uncertain variables are inputs for physical maps that output QIs
myInput = uq_createInput(InputOpts); 

plotsfolderName = 'frq_1_g2_uq'; 
mkdir(plotsfolderName)

% Description of the physical model for UQLab
ModelOpts.mFile = 'model_G2_frq_1';
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

flag_parfor = false;            % can we run the physical model in parallel to build the training set? (True/False)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
N_train_increment = 8;          % we will increment the training set size until we reach convergence
N_train_max = 50;               % training budget (i.e., maximum number of training points allowed)
% run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
% recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
flag_test_for_mean_and_sigma = false;
flag_test_set = true;           % will a test set be generated for further surrogate validation?
load('groundTests\testData\SJD_groundTestData.mat', 'exprData');
experimental_data_set = exprData{1, 3}.frequencies(1, :);   

% Plots generator for parameter sweeps for the uncertain variables
inputs_name = ["EI scaling factor", "GJ scaling factor", "Sxx scaling factor", "Szz scaling factor"];  % list of the names of the uncertain variables
outputs_name = ["1st freq(Hz) at 0deg", "1st freq(Hz) at 10deg", "1st freq(Hz) at 20deg", "1st freq(Hz) at 30deg", "1st freq(Hz) at 60deg", "1st freq(Hz) at 90deg"]; % list of the names of the QIs
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate', MetaOpts.MetaType);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
plotsfolderName = 'frq_1_g2_uq';   
mkdir(plotsfolderName, 'plots_uq');
tic;
surrogates =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, plotsfolderName, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
totalTime = toc;
fprintf('Total surrogate building time: %.4f seconds\n', totalTime);
elementToSave = surrogates;
save(fullfile(plotsfolderName, 'surrogates_frq_1_g2.mat'), 'elementToSave'); % save the surrogate
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

%% 9. Surrogates and parameter sweeps for the 2nd modal frequency (IP bending) as a function of EI, GJ, Sxx, Szz (G2 test)

% This section generates many plots, which are saved rather than displayed on the screen
set(0, 'DefaultFigureVisible', 'off'); 

% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for EI) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for GJ) lower and upper uncertainty bound
InputOpts.Marginals(3).Type = 'Uniform';
InputOpts.Marginals(3).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for Sxx) lower and upper uncertainty bound
InputOpts.Marginals(4).Type = 'Uniform';
InputOpts.Marginals(4).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for Szz) lower and upper uncertainty bound
% The uncertain variables are inputs for physical maps that output QIs
myInput = uq_createInput(InputOpts); 

plotsfolderName = 'frq_2_g2_uq'; 
mkdir(plotsfolderName)

% Description of the physical model for UQLab
ModelOpts.mFile = 'model_G2_frq_2';
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

flag_parfor = false;            % can we run the physical model in parallel to build the training set? (True/False)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
N_train_increment = 8;          % we will increment the training set size until we reach convergence
N_train_max = 50;               % training budget (i.e., maximum number of training points allowed)
% run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
% recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
flag_test_for_mean_and_sigma = false;
flag_test_set = true;           % will a test set be generated for further surrogate validation?
load('groundTests\testData\SJD_groundTestData.mat', 'exprData');
experimental_data_set = exprData{1, 3}.frequencies(2, :);   

% Plots generator for parameter sweeps for the uncertain variables
inputs_name = ["EI scaling factor", "GJ scaling factor", "Sxx scaling factor", "Szz scaling factor"];  % list of the names of the uncertain variables
outputs_name = ["2nd freq(Hz) at 0deg", "2nd freq(Hz) at 10deg", "2nd freq(Hz) at 20deg", "2nd freq(Hz) at 30deg", "2nd freq(Hz) at 60deg", "2nd freq(Hz) at 90deg"]; % list of the names of the QIs
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate', MetaOpts.MetaType);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
plotsfolderName = 'frq_2_g2_uq';   
mkdir(plotsfolderName, 'plots_uq');
tic;
surrogates =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, plotsfolderName, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
totalTime = toc;
fprintf('Total surrogate building time: %.4f seconds\n', totalTime);
elementToSave = surrogates;
save(fullfile(plotsfolderName, 'surrogates_frq_2_g2.mat'), 'elementToSave'); % save the surrogate
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

%% 10. Surrogates and parameter sweeps for the 3rd modal frequency (2nd OOP bending) as a function of EI, GJ, Sxx, Szz (G2 test)

% This section generates many plots, which are saved rather than displayed on the screen
set(0, 'DefaultFigureVisible', 'off'); 

% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for EI) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for GJ) lower and upper uncertainty bound
InputOpts.Marginals(3).Type = 'Uniform';
InputOpts.Marginals(3).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for Sxx) lower and upper uncertainty bound
InputOpts.Marginals(4).Type = 'Uniform';
InputOpts.Marginals(4).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for Szz) lower and upper uncertainty bound
% The uncertain variables are inputs for physical maps that output QIs
myInput = uq_createInput(InputOpts); 

plotsfolderName = 'frq_3_g2_uq'; 
mkdir(plotsfolderName)

% Description of the physical model for UQLab
ModelOpts.mFile = 'model_G2_frq_3';
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

flag_parfor = false;            % can we run the physical model in parallel to build the training set? (True/False)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
N_train_increment = 8;          % we will increment the training set size until we reach convergence
N_train_max = 50;               % training budget (i.e., maximum number of training points allowed)
% run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
% recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
flag_test_for_mean_and_sigma = false;
flag_test_set = true;           % will a test set be generated for further surrogate validation?
load('groundTests\testData\SJD_groundTestData.mat', 'exprData');
experimental_data_set = exprData{1, 3}.frequencies(3, :);   

% Plots generator for parameter sweeps for the uncertain variables
inputs_name = ["EI scaling factor", "GJ scaling factor", "Sxx scaling factor", "Szz scaling factor"];  % list of the names of the uncertain variables
outputs_name = ["3rd freq(Hz) at 0deg", "3rd freq(Hz) at 10deg", "3rd freq(Hz) at 20deg", "3rd freq(Hz) at 30deg", "3rd freq(Hz) at 60deg", "3rd freq(Hz) at 90deg"]; % list of the names of the QIs
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate', MetaOpts.MetaType);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
plotsfolderName = 'frq_3_g2_uq';   
mkdir(plotsfolderName, 'plots_uq');
tic;
surrogates =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, plotsfolderName, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
totalTime = toc;
fprintf('Total surrogate building time: %.4f seconds\n', totalTime);
elementToSave = surrogates;
save(fullfile(plotsfolderName, 'surrogates_frq_3_g2.mat'), 'elementToSave'); % save the surrogate
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

%% 11. Surrogates and parameter sweeps for the 4th modal frequency (torsion) as a function of EI, GJ, Sxx, Szz (G2 test)

% This section generates many plots, which are saved rather than displayed on the screen
set(0, 'DefaultFigureVisible', 'off'); 

% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for EI) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for GJ) lower and upper uncertainty bound
InputOpts.Marginals(3).Type = 'Uniform';
InputOpts.Marginals(3).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for Sxx) lower and upper uncertainty bound
InputOpts.Marginals(4).Type = 'Uniform';
InputOpts.Marginals(4).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for Szz) lower and upper uncertainty bound
% The uncertain variables are inputs for physical maps that output QIs
myInput = uq_createInput(InputOpts); 

plotsfolderName = 'frq_4_g2_uq'; 
mkdir(plotsfolderName)

% Description of the physical model for UQLab
ModelOpts.mFile = 'model_G2_frq_4';
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

flag_parfor = false;            % can we run the physical model in parallel to build the training set? (True/False)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
N_train_increment = 8;          % we will increment the training set size until we reach convergence
N_train_max = 50;               % training budget (i.e., maximum number of training points allowed)
% run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
% recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
flag_test_for_mean_and_sigma = false;
flag_test_set = true;           % will a test set be generated for further surrogate validation?
load('groundTests\testData\SJD_groundTestData.mat', 'exprData');
experimental_data_set = exprData{1, 3}.frequencies(4, :);   

% Plots generator for parameter sweeps for the uncertain variables
inputs_name = ["EI scaling factor", "GJ scaling factor", "Sxx scaling factor", "Szz scaling factor"];  % list of the names of the uncertain variables
outputs_name = ["4th freq(Hz) at 0deg", "4th freq(Hz) at 10deg", "4th freq(Hz) at 20deg", "4th freq(Hz) at 30deg", "4th freq(Hz) at 60deg", "4th freq(Hz) at 90deg"]; % list of the names of the QIs
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate', MetaOpts.MetaType);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
plotsfolderName = 'frq_4_g2_uq';   
mkdir(plotsfolderName, 'plots_uq');
tic;
surrogates =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, plotsfolderName, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
totalTime = toc;
fprintf('Total surrogate building time: %.4f seconds\n', totalTime);
elementToSave = surrogates;
save(fullfile(plotsfolderName, 'surrogates_frq_4_g2.mat'), 'elementToSave'); % save the surrogate
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

%% 12. Surrogates and parameter sweeps for the 1st modal frequency (OOP bending) as a function of EI (G2 test)

% This section generates many plots, which are saved rather than displayed on the screen
set(0, 'DefaultFigureVisible', 'off'); 

% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for EI) lower and upper uncertainty bound
% The uncertain variables are inputs for physical maps that output QIs
myInput = uq_createInput(InputOpts); 

plotsfolderName = 'frq_1_g2_dim_red_uq'; 
mkdir(plotsfolderName)

% Description of the physical model for UQLab
ModelOpts.mFile = 'model_G2_frq_1_dim_red';
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

flag_parfor = false;            % can we run the physical model in parallel to build the training set? (True/False)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
N_train_increment = 8;          % we will increment the training set size until we reach convergence
N_train_max = 50;               % training budget (i.e., maximum number of training points allowed)
% run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
% recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
flag_test_for_mean_and_sigma = false;
flag_test_set = true;           % will a test set be generated for further surrogate validation?
load('groundTests\testData\SJD_groundTestData.mat', 'exprData');
experimental_data_set = exprData{1, 3}.frequencies(1, :);   

% Plots generator for parameter sweeps for the uncertain variables
inputs_name = ["EI scaling factor"];  % list of the names of the uncertain variables
outputs_name = ["1st freq(Hz) at 0deg", "1st freq(Hz) at 10deg", "1st freq(Hz) at 20deg", "1st freq(Hz) at 30deg", "1st freq(Hz) at 60deg", "1st freq(Hz) at 90deg"]; % list of the names of the QIs
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate', MetaOpts.MetaType);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
plotsfolderName = 'frq_1_g2_dim_red_uq';   
mkdir(plotsfolderName, 'plots_uq');
tic;
surrogates =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, plotsfolderName, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
totalTime = toc;
fprintf('Total surrogate building time: %.4f seconds\n', totalTime);
elementToSave = surrogates;
save(fullfile(plotsfolderName, 'surrogates_frq_1_g2_dim_red.mat'), 'elementToSave'); % save the surrogate
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

%% 13. Surrogates and parameter sweeps for the 2nd modal frequency (IP bending) as a function of EI, GJ (G2 test)

% This section generates many plots, which are saved rather than displayed on the screen
set(0, 'DefaultFigureVisible', 'off'); 

% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for EI) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for GJ) lower and upper uncertainty bound
% The uncertain variables are inputs for physical maps that output QIs
myInput = uq_createInput(InputOpts); 

plotsfolderName = 'frq_2_g2_dim_red_uq'; 
mkdir(plotsfolderName)

% Description of the physical model for UQLab
ModelOpts.mFile = 'model_G2_frq_2_dim_red';
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

flag_parfor = false;            % can we run the physical model in parallel to build the training set? (True/False)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
N_train_increment = 8;          % we will increment the training set size until we reach convergence
N_train_max = 50;               % training budget (i.e., maximum number of training points allowed)
% run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
% recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
flag_test_for_mean_and_sigma = false;
flag_test_set = true;           % will a test set be generated for further surrogate validation?
load('groundTests\testData\SJD_groundTestData.mat', 'exprData');
experimental_data_set = exprData{1, 3}.frequencies(2, :);   

% Plots generator for parameter sweeps for the uncertain variables
inputs_name = ["EI scaling factor", "GJ scaling factor"];  % list of the names of the uncertain variables
outputs_name = ["2nd freq(Hz) at 0deg", "2nd freq(Hz) at 10deg", "2nd freq(Hz) at 20deg", "2nd freq(Hz) at 30deg", "2nd freq(Hz) at 60deg", "2nd freq(Hz) at 90deg"]; % list of the names of the QIs
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate', MetaOpts.MetaType);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
plotsfolderName = 'frq_2_g2_dim_red_uq';   
mkdir(plotsfolderName, 'plots_uq');
tic;
surrogates =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, plotsfolderName, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
totalTime = toc;
fprintf('Total surrogate building time: %.4f seconds\n', totalTime);
elementToSave = surrogates;
save(fullfile(plotsfolderName, 'surrogates_frq_2_g2_dim_red.mat'), 'elementToSave'); % save the surrogate
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

%% 14. Surrogates and parameter sweeps for the 3rd modal frequency (2nd OOP bending) as a function of EI (G2 test)

% This section generates many plots, which are saved rather than displayed on the screen
set(0, 'DefaultFigureVisible', 'off'); 

% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for EI) lower and upper uncertainty bound
% The uncertain variables are inputs for physical maps that output QIs
myInput = uq_createInput(InputOpts); 

plotsfolderName = 'frq_3_g2_dim_red_uq'; 
mkdir(plotsfolderName)

% Description of the physical model for UQLab
ModelOpts.mFile = 'model_G2_frq_3_dim_red';
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

flag_parfor = false;            % can we run the physical model in parallel to build the training set? (True/False)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
N_train_increment = 8;          % we will increment the training set size until we reach convergence
N_train_max = 50;               % training budget (i.e., maximum number of training points allowed)
% run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
% recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
flag_test_for_mean_and_sigma = false;
flag_test_set = true;           % will a test set be generated for further surrogate validation?
load('groundTests\testData\SJD_groundTestData.mat', 'exprData');
experimental_data_set = exprData{1, 3}.frequencies(3, :);   

% Plots generator for parameter sweeps for the uncertain variables
inputs_name = ["EI scaling factor"];  % list of the names of the uncertain variables
outputs_name = ["3rd freq(Hz) at 0deg", "3rd freq(Hz) at 10deg", "3rd freq(Hz) at 20deg", "3rd freq(Hz) at 30deg", "3rd freq(Hz) at 60deg", "3rd freq(Hz) at 90deg"]; % list of the names of the QIs
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate', MetaOpts.MetaType);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
plotsfolderName = 'frq_3_g2_dim_red_uq';  
mkdir(plotsfolderName, 'plots_uq');
tic;
surrogates =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, plotsfolderName, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
totalTime = toc;
fprintf('Total surrogate building time: %.4f seconds\n', totalTime);
elementToSave = surrogates;
save(fullfile(plotsfolderName, 'surrogates_frq_3_g2_dim_red.mat'), 'elementToSave'); % save the surrogate
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

%% 15. Surrogates and parameter sweeps for the 4th modal frequency (torsion) as a function of GJ, Sxx, Szz (G2 test)

% This section generates many plots, which are saved rather than displayed on the screen
set(0, 'DefaultFigureVisible', 'off'); 

% Description of the uncertain variables for UQLab
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for GJ) lower and upper uncertainty bound
InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0.8, 1.2]; % (multiplicative scaling factor for Sxx) lower and upper uncertainty bound
% The uncertain variables are inputs for physical maps that output QIs
myInput = uq_createInput(InputOpts); 

plotsfolderName = 'frq_4_g2_dim_red_uq'; 
mkdir(plotsfolderName)

% Description of the physical model for UQLab
ModelOpts.mFile = 'model_G2_frq_4_dim_red';
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

flag_parfor = false;            % can we run the physical model in parallel to build the training set? (True/False)
seed = 100;                     % seed for reproducibility due to randomness in sampling the training set
N_train_increment = 8;          % we will increment the training set size until we reach convergence
N_train_max = 50;               % training budget (i.e., maximum number of training points allowed)
% run a test to check if surrogates are actually faster than classical MC for mean and sigma estimation  
% recommended only for cheap models (to find the true mean and sigma, we need a large MC with the physical model) 
flag_test_for_mean_and_sigma = false;
flag_test_set = true;           % will a test set be generated for further surrogate validation?
load('groundTests\testData\SJD_groundTestData.mat', 'exprData');
experimental_data_set = exprData{1, 3}.frequencies(4, :);   

% Plots generator for parameter sweeps for the uncertain variables
inputs_name = ["GJ scaling factor", "Sxx scaling factor"];  % list of the names of the uncertain variables
outputs_name = ["4th freq(Hz) at 0deg", "4th freq(Hz) at 10deg", "4th freq(Hz) at 20deg", "4th freq(Hz) at 30deg", "4th freq(Hz) at 60deg", "4th freq(Hz) at 90deg"]; % list of the names of the QIs
N_outputs = length(outputs_name);             % number of quantities of interest (QIs)
descriptive_title_for_plots = sprintf('%s surrogate', MetaOpts.MetaType);
N_eval = 100;                                                        % number of discretisation points for each uncertain variable (for plots)
plotsfolderName = 'frq_4_g2_dim_red_uq';   
mkdir(plotsfolderName, 'plots_uq');
tic;
surrogates =  surrogates_uq(MetaOpts, N_outputs, N_train_increment, N_train_max, flag_parfor, seed, plotsfolderName, flag_test_for_mean_and_sigma, flag_test_set); % Generates training points and builds the surrogates 
totalTime = toc;
fprintf('Total surrogate building time: %.4f seconds\n', totalTime);
elementToSave = surrogates;
save(fullfile(plotsfolderName, 'surrogates_frq_4_g2_dim_red.mat'), 'elementToSave'); % save the surrogate
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