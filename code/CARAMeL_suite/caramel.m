function [data_struct, ML_model] = ...
    caramel(phenotype_structure, interaction_structure, mode, varargin)

% DESCRIPTION: 
% Main script for using the CARAMeL algorithm to train a machine 
% learning (ML) model to predict drug interaction outcomes in condition-
% specific environments. Includes both training and prediction modes for 
% model development. 
% 
% STEPS: 
% 1. Parse through inputs
% 2. Ascertain input compatibility
% 3. Filter and process data
% 4. Define joint profiles (ML input)
% 5. Train or predict for the given interaction dataset
% 
% Author:   Carolina H. Chung
% Contact:  chechung@umich.edu

% I/O/U
%{
REQUIRED INPUTS: 
    1. phenotype_structure:     structure containing 4 fields
       a. data:         numeric matrix for drug and media flux profiles 
                        rows = features, columns = conditions
       b. features:             cell array of feature names 
                                (e.g. reaction names)
       c. conditions:           cell array of condition names 
                                (e.g. drug or media names)
       d. H:                    double array of drug and media entropies

    2. interaction_structure: 	structure containing at least 2 fields
                                'names' is required; choose between 
                                'scores' and 'class'
       a. names:        cell array of interacting drug and/or media names
                        each row corresponds to a single interaction
          -> note:      if list contains combinations of different lengths,
                        leave empty cells for combinations with smaller 
                        number of drugs
       b. scores:       double array of interaction scores that correspond
                        to entries in 'names'
       c. class:        cell or categorical array of interaction classes
                        that correspond to entries in 'names'
          -> accepted 
             values:    'Synergy', 'Additivity', and 'Antagonism'
       d. time:         numeric array of time values for each drug/media
                        treatment in a combination (relevant for 
                        sequential interactions)

    3. mode:                    string specifying mode for algorithm 
       -> choose between:       {'train', 'predict'}

OPTIONAL INPUTS: 
    1.  Key:                    cell array containing key information 
                                that matches conditions (from
                                phenotype_structure) to name values for  
                                interactions (from interaction_structure)
        -> note:                make sure 1st column corresponds to entries
                                in interaction names in structure variable
        -> default:             Key = []

    2.  MLtype:                 string specifying which ML method to use
        -> default:             MLtype = 'regRF'
        -> choose amongst:      {'regRF', 'classRF', 'fitrensemble', 
                                'fitcensemble', 'fitrsvm', 'fitcsvm'}

    3.  MLmodel:                variable to use as the ML model when
                                mode = 'predict' 
        -> default:             MLmodel = []

    4.  Verbose:                boolean specifying whether to print 
                                messages while algorithm is run
        -> default:             Verbose = false

    5.  regRFtrainOptions:      structure variable containing 
                                additional parameters for training the 
                                model based on 'regRF'
        -> default:             regRFtrainOptions = struct('importance',1)

    6. classRFtrainOptions:     structure variable containing 
                                additional parameters for training the 
                                model based on 'classRF'
        -> default:             classRFtrainOptions = struct('importance',1)

    7. classRFpredictOptions:   structure variable containing 
                                additional parameters for generating 
                                predictions based on 'classRF'
        -> default:             classRFpredictOptions = []

OTHER FUNCTION PARAMETERS: 
    1. fitrensemble:    any additional parameters available with the
                        built-in MATLAB function 'fitrensemble'
                        (when MLtype = 'fitrensemble')

    2. fitcensemble:    any additional parameters available with the
                        built-in MATLAB function 'fitcensemble'
                        (when MLtype = 'fitcensemble')

    3. fitrsvm:         any additional parameters available with the
                        built-in MATLAB function 'fitrsvm'
                        (when MLtype = 'fitrsvm')

    4. fitcsvm:         any additional parameters available with the
                        built-in MATLAB function 'fitcsvm'
                        (when MLtype = 'fitcsvm')

    5. predict:         any additional parameters available with the
                        built-in MATLAB function 'predict'
                        (when MLtype = 'fitrensemble', 'fitcensemble', 
                        'fitrsvm', or 'fitcsvm' and mode = 'predict')

    6. corr:            any additional parameters available with the
                        built-in MATLAB function 'corr'

OUTPUTS: 
    1. data_struct:         structure variable containing information
                            processed from function inputs; includes 
                            inputs used for ML model development

    2. ML_model:            variable for the ML model that is developed
                            (mode = 'train') or model that was used to 
                            generate predictions (mode = 'predict')

EXAMPLE USAGE: 
    1. Train a model given PS as the phenotype data structure and IS as the 
       structure variable with the interaction data: 
       >> [~, model] = caramel(PS, IS, 'train'); 

    2. Determine predictions for new interactions (provided through a 
       structure variable called TS) given a trained model (model): 
       >> predictions = caramel(PS, TS, 'predict', 'MLmodel', model); 
%}

%% PARSE THROUGH INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create an inputParser object
    p = inputParser; 
    if rem(length(varargin), 2) ~= 0
        error('Imbalanced number of name-value pair arguments provided.')
    end
    
    % Define accepted ML functions and class types
    regML = {'regRF','fitrensemble','fitrsvm'};
    classML = {'classRF','fitcensemble','fitcsvm'}; 
    classes = {'Synergy','Additivity','Antagonism'};
    
    % Define required input validation functions
    isValidPhenotype = @(x) ...
        isstruct(x) && numel(fieldnames(x)) == 4 && ...
        all(isfield(x, {'data','features','conditions','H'})) && ...
        size(x.data, 1) == numel(x.features) && ...
        size(x.data, 2) == numel(x.conditions) && ...
        size(x.data, 2) == numel(x.H);
    isValidInteraction = @(x) ...
        isstruct(x) && isfield(x, 'names') && ...
        (isfield(x, 'scores') && size(x.names, 1) == numel(x.scores)) || ...
        (isfield(x, 'class') && size(x.names, 1) == numel(x.class));
    isValidMode = @(x) any(validatestring(x, {'train', 'predict'}));
    
    % Define required input variables
    addRequired(p, 'phenotype_structure', isValidPhenotype)
    addRequired(p, 'interaction_structure', isValidInteraction)
    addRequired(p, 'mode', isValidMode)
    
    % Define optional input validation functions 
    isValidKey = @(x) iscell(x) && size(x, 2) == 2 && ...
        any(ismember(phenotype_structure.conditions, x));
    isValidMLtype = @(x) any(validatestring(x, horzcat(regML, classML)));
    isValidMLmodel = @(x) ~isempty(x); 
    isValidBoolean = @(x) islogical(x) || ...
        (isnumeric(x) && (x == 0 || x == 1));
    
    % Define optional input parameters
    addParameter(p, 'Key', [], isValidKey)
    addParameter(p, 'MLtype', 'regRF', isValidMLtype)
    addParameter(p, 'MLmodel', [], isValidMLmodel)
    addParameter(p, 'Verbose', false, isValidBoolean)
    addParameter(p, 'regRFtrainOptions', struct('importance', 1), @isstruct)
    addParameter(p, 'classRFtrainOptions', struct('importance', 1), @isstruct)
    addParameter(p, 'classRFpredictOptions', [], @isstruct)
    
    % Parse through inputs
    p.KeepUnmatched = true; 
    parse(p, phenotype_structure, interaction_structure, mode, varargin{:})
    if ~isempty(fieldnames(p.Unmatched))
        extraParams = p.Unmatched; 
        if logical(p.Results.Verbose)
            disp('Extra name-value pair arguments provided:')
            disp(p.Unmatched)
        end
    else
        extraParams = {};
    end
    
    % Define inputs by short names
    ps                  = phenotype_structure;
    is                  = interaction_structure;
    key                 = p.Results.Key;
    MLtype              = p.Results.MLtype;
    MLmodel             = p.Results.MLmodel;
    verbose             = logical(p.Results.Verbose);
    regRFtrainOpt       = p.Results.regRFtrainOptions;
    classRFtrainOpt     = p.Results.classRFtrainOptions;
    classRFpredictOpt   = p.Results.classRFpredictOptions;
    
    % Define list of optional input parameters from external functions
    fitrensembleParams = {'Method','NumLearningCycles','Learners',...
        'NPrint','NumBins','CategoricalPredictors','PredictorNames',...
        'ResponseName','ResponseTransform','CrossVal','CVPartition',...
        'Holdout','KFold','Leaveout','Weights','FResample','Replace',...
        'Resample','LearnRate','OptimizeHyperparameters',...
        'HyperparameterOptimizationOptions'};
    fitcensembleParams = {'Method','NumLearningCycles','Learners',...
        'NPrint','NumBins','CategoricalPredictors','PredictorNames',...
        'ResponseName','CrossVal','CVPartition','Holdout','KFold',...
        'Leaveout','ClassNames','Cost','Prior','ScoreTransform',...
        'Weights','FResample','Replace','Resample','LearnRate',...
        'RatioToSmallest','MarginPrecision','RobustErrorGoal',...
        'RobustMarginSigma','RobustMaxMargin','NPredToSample',...
        'OptimizeHyperparameters','HyperparameterOptimizationOptions'};
    fitrsvmParams = {};
    fitcsvmParams = {}; 
    predictParams = {'Learners','UseObsForLearner'};
    paramList = {fitrensembleParams, fitcensembleParams, ...
        fitrsvmParams, fitcsvmParams, predictParams}; 
    f = {'fitrensemble','fitcensemble','fitrsvm','fitcsvm', 'predict'};
    
    % Find any extra parameters and allocate to appropriate function
    if ~isempty(extraParams)
        extraNames = fieldnames(extraParams); 
        for k = 1:numel(paramList)
            [s, p_idx] = intersect(lower(extraNames), lower(paramList{k})); 
            if numel(p_idx) > 0
                eval([f{k} 'Varargin = cell(1, 2*numel(p_idx));'])
                for j = 1:numel(p_idx)
                    [i1, i2] = deal(2*(j-1)+1, 2*j); 
                    eval([f{k} 'Varargin{i1} = extraNames{p_idx(j)};'])
                    eval([f{k} 'Varargin{i2} = extraParams.' ...
                        extraNames{p_idx(j)}, ';'])
                end
                if verbose
                    fprintf('Extra input(s) passed to %s function. \n', ...
                        upper(f{k}))
                end
            else
                eval([f{k} 'Varargin = {};'])
            end
        end
    else
        for k = 1:numel(paramList)
            eval([f{k} 'Varargin = {};'])
        end
    end
    
    % Display input settings 
    if verbose
        if strcmpi(mode, 'train')
            disp('Initializing CARAMeL model training...')
        elseif strcmpi(mode, 'predict')
            disp('Predicting interaction outcomes...')
        end
        disp(strcat("Method chosen for RF: ", MLtype))
        disp('Running CARAMeL with the following user-defined inputs:')
        disp(p.Results)
        disp('Using default values for the following inputs:')
        disp(p.UsingDefaults')
    end

%% ASCERTAIN INPUT COMPATIBILITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Make sure required fields (based on MLtype) are present
    if ismember(MLtype, regML) && ~isfield(is, 'scores')
        error('Provide interaction scores for regression-based model.')
    end
    if ismember(MLtype, classML) && ~isfield(is, 'class')
        error('Provide interaction classes for classifier-based model.')
    end
    
    % Check if threshold information was provided
    if ismember(MLtype, regML) && ~isfield(is, 'threshold')
        if verbose
            disp(strcat("Provide threshold information to ", ...
                "classify interactions based on their score."))
        end
    end
    
%% FILTER AND PROCESS DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check interaction names against phenotype conditions
    drug_conditions = unique(is.names); 
    drug_conditions(cellfun(@isempty, drug_conditions)) = [];
    if ~all(ismember(drug_conditions, ps.conditions))
        if verbose
            disp('Not all interaction names match phenotype data.')
            fprintf('Attempting to match interaction and phenotype data...')
        end
        if isempty(key)
            error('Provide key data to match interactions to phenotypes.')
        elseif size(key, 2) ~= 2
            error('Key data does not include all data needed.')
        elseif ~any(ismember(drug_conditions, key(:)))
            error('Interaction names do not match provided key.')
        else
            [query_list, match_list] = deal(key(:, 1), key(:, 2)); 
        end
        % 	match interaction names to phenotype conditions
        [drugs, ~, idx] = intersect(drug_conditions, query_list, 'stable');
        combos = is.names;
        for i = 1:numel(drugs)
            combos(strcmpi(drugs{i}, is.names)) = match_list(idx(i));
        end
        if verbose
            fprintf('completed successfully. \n')
        end
    elseif all(ismember(drug_conditions, ps.conditions))
        if verbose
            disp('All interaction names accounted for in phenotype data.')
        end
        combos = is.names;
    end
    
    % Filter out incomplete interaction entries
    match_idx = ismember(combos, ps.conditions); 
    ix = sum(~cellfun(@isempty, is.names), 2) == sum(match_idx, 2);
    if sum(ix) < size(is.names, 1)
        warning('Not able to %s for all interactions.', mode)
    end
    interactionNames = combos(ix, :); 
    
    % Define actual interaction outcome data
    if ismember(MLtype, regML) 
        interactionScores = is.scores(ix);
        if isfield(is, 'class')
            interactionClass = is.class(ix); 
        elseif isfield(is, 'threshold')
            interactionClass = ...
                caramel_classify(interactionScores, is.threshold); 
        else
            interactionClass = [];
        end
    else
        interactionClass = is.class(ix); 
        if isfield(is, 'scores')
            interactionScores = is.scores(ix); 
        else
            interactionScores = [];
        end
    end
    
    % Account for time (if provided)
    if isfield(is, 'time')
        interactionTime = is.time(ix, :); 
    else
        interactionTime = []; 
    end
    
    % Extract phenotype information relevant to given interactions
    phenotypeConditions = unique(interactionNames(:), 'stable');
    %   account for empty cell
    phenotypeConditions(cellfun(@isempty, phenotypeConditions)) = [];
    [~, ~, idx] = intersect(phenotypeConditions, ps.conditions, 'stable'); 
    phenotypeArray = ps.data(:, idx);
    %   account for entropy 
    H = ps.H(idx)'; 
    %   define phenotype data variable
    t1 = table(vertcat(ps.features, 'Entropy'), 'VariableNames', {'Label'}); 
    col = unique(is.names(ix, :), 'stable'); 
    col(cellfun(@isempty, col)) = []; 
    t2 = array2table([phenotypeArray; H'], 'VariableNames', col'); 
    phenotypeData = horzcat(t1, t2); 

%% DEFINE JOINT PROFILES (ML INPUT) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ML input variable
    X = nan(2*size(phenotypeArray, 1) + 3, size(interactionNames, 1));
    %   for each interaction
    for i = 1:size(interactionNames, 1)
        % find conditions of interest
        [~, ~, idx] = intersect(interactionNames(i, :), ...
            phenotypeConditions, 'stable'); 
        % define sigma features
        sigma = sum(phenotypeArray(:, idx), 2) .* (2/numel(idx));
        % define delta + time features
        if isempty(interactionTime) || sum(interactionTime(i, :)) == 0 
            delta = sum(phenotypeArray(:, idx), 2) == 1;
            time = 0; 
        else
            p_data = phenotypeArray(:, idx); 
            t_idx = find(~isnan(interactionTime(i, :))); 
            time = sum(interactionTime(i, t_idx(1:end-1))); 
            delta = diff(interactionTime(i, t_idx) .* p_data(:, t_idx), ...
                numel(t_idx) - 1, 2) ./ sum(interactionTime(i, t_idx));
%             delta = diff(p_data(:, t_idx), sum(t_idx) - 1, 2) ./ time;
        end
        % define entropy features
        X(1:end - 2, i) = [sigma; delta; time];
        X(end - 1, i) = mean(H(idx)); 
        X(end, i) = sum(H(idx)); 
    end
    %   define feature labels
    labels = vertcat(strcat('sigma-', ps.features), strcat('delta-', ...
        ps.features), {'time','entropy-mean','entropy-sum'}'); 
    %   make sure labels are unique
    labels = matlab.lang.makeUniqueStrings(labels);
    %   define joint profile variable
    t1 = table(labels, 'VariableNames', {'Feature'}); 
    joinLabels = join(is.names(ix, :), '_'); 
    if numel(unique(joinLabels)) < numel(joinLabels)
        joinLabels = matlab.lang.makeUniqueStrings(joinLabels);
    end
    t2 = array2table(X, 'VariableNames', joinLabels); 
    jointProfiles = horzcat(t1, t2); 
    %   transpose input (for function compatibility)
    X = X';  

%% TRAIN OR PREDICT FOR THE GIVEN INTERACTION DATASET %%%%%%%%%%%%%%%%%%%%%

    % Define additional inputs (if MLmode = regRF or classRF)
    if ismember(MLtype, {'regRF','classRF'}) && strcmpi(mode, 'train') 
        ntree = 500; mtry = max(floor(size(X, 2)/3), 1);
        if verbose
            fprintf('\t Setting to defaults %d trees and mtry=%d \n', ...
                ntree, mtry);
        end
    end
    
    % Train or predict
    switch mode
        case 'train'
            switch MLtype
                case 'regRF'
                    ML_model = regRF_train(X, interactionScores, ntree, ...
                        mtry, regRFtrainOpt);
                    ML_model.features = labels; 
                case 'classRF'
                    y = zeros(size(interactionClass));
                    for i = 1:numel(classes)
                        y(ismember(interactionClass, classes{i})) = i;
                    end
                    ML_model = classRF_train(X, y, ntree, mtry, ...
                        classRFtrainOpt);
                    ML_model.features = labels; 
                case 'fitrensemble' 
                    ML_model = fitrensemble(X, interactionScores, ...
                        'PredictorNames', labels, fitrensembleVarargin{:});
                case 'fitcensemble'
                    ML_model = fitcensemble(X, interactionClass, ...
                        'PredictorNames', labels, 'ClassNames', classes, ...
                        fitcensembleVarargin{:});
                case 'fitrsvm'
                    ML_model = fitrsvm(X, interactionScores, ...
                        'PredictorNames', labels, fitrsvmVarargin{:}); 
                case 'fitcsvm'
                    ML_model = fitrsvm(X, interactionClass, ...
                        'PredictorNames', labels, fitcsvmVarargin{:}); 
            end
        case 'predict'
            if isempty(MLmodel)
                error('Provide trained ML model to generate predictions.')
            end
            switch MLtype
                case 'regRF'
                    predScores = regRF_predict(X, MLmodel);
                case 'classRF'
                    [yhat, votes] = classRF_predict(X, MLmodel, ...
                        classRFpredictOpt);
                    predClass = cell(size(yhat));
                    for i = 1:numel(classes)
                        predClass(yhat == i) = classes(i);
                    end
                    predScores = votes ./ ntree; 
                case 'fitrensemble'
                    predScores = predict(MLmodel, X, predictVarargin{:});
                case 'fitcensemble'
                    [predClass, predScores] = predict(MLmodel, X, ...
                        predictVarargin{:});
                case 'fitrsvm'
                    predScores = predict(MLmodel, X, predictVarargin{:});
                case 'fitcsvm'
                    predClass = predict(MLmodel, X, predictVarargin{:});
                    predScores = [];
            end
    end
    
    % Define ML method used for results
    if ismember(MLtype, regML)
        MLmethod = 'regression'; 
    elseif ismember(MLtype, classML)
        MLmethod = 'classification'; 
    end
    
    % Classify predictions (if ML method is regression-based)
    if strcmp(MLmethod, 'regression') && strcmpi(mode, 'predict')
        if isfield(is, 'threshold')
            predClass = caramel_classify(predScores, is.threshold); 
        else
            predClass = [];
        end
    end
    
    % Finalize ML portion based on function mode
    if strcmpi(mode, 'train')
        data_struct = struct(...
            'interactionNames', {is.names(ix, :)}, ...
            'interactionScores', interactionScores, ...
            'interactionClass', {interactionClass}, ...
            'phenotypeData', phenotypeData, ...
            'jointProfiles', jointProfiles, ...
            'MLmethod', MLmethod); 
        if verbose
            disp('Training complete. Exiting function now.')
        end
    elseif strcmpi(mode, 'predict')
        ML_model = MLmodel;
        data_struct = struct(...
            'interactionNames', {is.names(ix, :)}, ...
            'interactionScores', interactionScores, ...
            'interactionClass', {interactionClass}, ...
            'phenotypeData', phenotypeData, ...
            'jointProfiles', jointProfiles, ...
            'predScores', predScores, ...
            'predClass', {predClass}, ...
            'MLmethod', MLmethod); 
        if verbose
            disp('Predictions determined. Exiting function now.')
        end
    end

end