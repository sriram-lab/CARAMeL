function [cv_struct, best_model, idx] = ...
    caramel_crossValidation(phenotype_structure, interaction_structure, ...
    varargin)

% DESCRIPTION:
% This script conducts a k-fold cross-validation for training and
% testing CARAMeL models. 
% 
% STEPS:
% 1. Parse through inputs
% 2. Ascertain input compatibility
% 3. Conduct k-fold cross-validation
% 
% Author:   Carolina H. Chung
% Contact:  chechung@umich.edu
    
% I/O/U
%{
REQURED INPUTS: 
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

OPTIONAL INPUTS: 
    1. Kfold:                   positive integer specifying the k value 
                                for k-fold cross-validation
       -> default:              Kfold = 10

    2. Indices:                 numerical array of k-fold indices (length
                                must match number of drug combinations in
                                'interaction_struct')
       -> default:              Indices = crossvalind('Kfold', y, k)

    3. ApplyKfold:              Boolean array specifying which data points
                                to include for k-fold partitioning
       -> default:              ApplyKfold = true(size(y))

    4. DisplayOutput:           boolean specifying whether to show 
                                results and plots for model performance
       -> default:              DisplayOutput = false)

OTHER FUNCTION PARAMETERS: 
    1. caramel:          any additional parameters available with the
                        'caramel' function
    
OUTPUTS: 
    1. cv_struct:       structure containing the following fields:
        a. data:        array containing data structures from caramel
        b. model:       array containing models built for each k-fold
        c. results:     array containing results structures from caramel

    2. best_model:      best performing CARAMeL model (based on correlation)

    3. idx:             numerical array of k-fold indices

EXAMPLE USAGE:  
    1. Conduct 10-fold cross-validation given phenotype structure (PS) and 
       interaction structure (IS) and return best performing model: 
       >> [CV_results, best_model] = caramel_crossValidation(PS, IS); 

    2. Conduct 5-fold cross-validation for first 50 interaction entries: 
       >> applyk = false(size(IS.scores)); 
       >> applyk(1:50) = true(50, 1); 
       >> CV_results = caramel_crossValidation(PS, IS, 'Kfold', 5, ...
                                              'ApplyKfold', applyk); 
    3. Store k-fold indices from (1) and re-run cross-validation: 
       >> [CV_results, idx] = caramel_crossValidation(PS, IS); 
       >> CV_results2 = caramel_crossValidation(PS, IS, 'Indices, idx);
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
    
    % Define required input variables
    addRequired(p, 'phenotype_structure', isValidPhenotype)
    addRequired(p, 'interaction_structure', isValidInteraction)
    
    % Define optional input validation functions 
    isValidKfold = @(x) isnumeric(x) && ...
        numel(x) == 1 && rem(x, 1) == 0 && x > 0; 
    isValidBoolean = @(x) islogical(x) || ...
        (isnumeric(x) && (x == 0 || x == 1));
    isValidTimer = @(x) any(validatestring(x, {'on', 'off'}));
    
    % Define optional input parameters
    addParameter(p, 'Kfold', 10, isValidKfold)
    addParameter(p, 'Indices', [], isValidBoolean)
    addParameter(p, 'ApplyKfold', [], isValidBoolean)
    addParameter(p, 'Timer', 'on', isValidTimer)
    addParameter(p, 'Verbose', true, isValidBoolean)
    addParameter(p, 'DisplayOutput', false, isValidBoolean)
    
    % Parse through inputs
    p.KeepUnmatched = true; 
    parse(p, phenotype_structure, interaction_structure, varargin{:})
    if ~isempty(fieldnames(p.Unmatched))
        extraParams = p.Unmatched; 
        if p.Results.Verbose
            disp('Extra name-value pair arguments provided:')
            disp(p.Unmatched)
        end
    else
        extraParams = [];
    end
    
    % Define list of optional input parameters from external functions
    caramelParams = {'Key', 'MLtype', 'MLmodel', 'Verbose', ...
        'regRFtrainOptions', 'classRFtrainOptions', 'classRFpredicOptions'};
    
    % Find any extra parameters and allocate to CARDIG function
    if ~isempty(extraParams)
        extraNames = fieldnames(extraParams); 
        [~, p_idx] = intersect(lower(extraNames), lower(caramelParams)); 
        if numel(p_idx) > 0
            caramelVarargin = cell(1, 2*numel(p_idx)); 
            for j = 1:numel(p_idx)
                [i1, i2] = deal(2*j-1, 2*j); 
                caramelVarargin{i1} = extraNames{p_idx(j)};
                eval(sprintf('caramelVarargin{i2} = extraParams.%s;', ...
                   extraNames{p_idx(j)}))
            end
            if strcmpi(caramelVarargin, 'MLtype')
                idx = find(strcmpi(caramelVarargin, 'MLtype')); 
                MLtype = caramelVarargin{idx + 1}; 
            else
                MLtype = 'regRF'; 
            end
        else
            caramelVarargin = {}; MLtype = 'regRF'; 
        end
    else
        caramelVarargin = {}; MLtype = 'regRF'; 
    end
    
    % Define inputs by name
    ps      = p.Results.phenotype_structure; 
    is      = p.Results.interaction_structure; 
    k       = p.Results.Kfold; 
    idx     = p.Results.Indices; 
    applyk  = p.Results.ApplyKfold; 
    timer   = p.Results.Timer; 
    verbose = logical(p.Results.Verbose); 
    display = logical(p.Results.DisplayOutput); 
    
    % Start timer (if prompted)
    if strcmpi(timer, 'on')
        tStart = tic; 
    end
    
    % Display message
    if verbose
        if k < size(is.names, 1)
            fprintf('Initializing %d-fold cross-validation... \n', k)
        else
            disp('Initializing leave-one-out cross-validation...')
        end
    end

%% ASCERTAIN INPUT COMPATIBILITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Stop function if k > # of observations
    if k > size(is.names, 1)
        error('k is greater than the number of observations.')
    end
    
    % Define/check applyk parameter
    if isempty(applyk)
        applyk = true(size(is.scores)); 
    elseif size(applyk) ~= size(is.names, 1)
        error('Provide valid input for ApplyKfold parameter.')
    end
    [a, b] = deal(is.names(~applyk, :), is.names(applyk, :)); 
    [x, y] = deal(is.scores(~applyk), is.scores(applyk)); 
    
    % Make sure required fields (based on MLtype) are present
    if ismember(MLtype, regML) && ~isfield(is, 'scores')
        error('Provide interaction scores for regression-based model.')
    end
    if ismember(MLtype, classML) && ~isfield(is, 'class')
        error('Provide interaction classes for classifier-based model.')
    end
    
    % Specify classification thresholds (if needed)
    if ismember(MLtype, regML) && ~isfield(is, 'threshold')
        is.threshold = quantile(is.scores, [0.25 0.75]); 
    end
        
%% CONDUCT K-FOLD CROSS-VALIDATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Determine indices from k (if needed)
    if isempty(idx)
        idx = crossvalind('Kfold', y, k); 
    end
    
    % k-fold cross-validation
    predScores = nan(size(y)); predClass = cell(size(y)); 
    if ismember(MLtype, regML)  
        r = -1; 
    else
        accuracy = 0; 
    end
    if verbose
        if k < size(is.names, 1)
            progressbar(sprintf('Conducting %d-fold cross-validation...', k))
        else
            progressbar('Conducting leave-one-out cross-validation...')
        end
    end
    for i = 1:k
        % Define training structure
        train_struct = struct(...
            'names', {vertcat(a, b(idx ~= i, :))}, ...
            'scores', [x; y(idx ~= i)]); 
        % Define testing structure
        test_struct = struct(...
            'names', {b(idx == i, :)}, ...
            'scores', y(idx == i), 'threshold', is.threshold); 
        % Train caramel Model
        [~, model] = ...
            caramel(ps, train_struct, 'train', caramelVarargin{:}); 
        % Predict interaction scores using GEM-ML model
        data = caramel(ps, test_struct, 'predict', 'MLmodel', model, ...
            caramelVarargin{:});
        predScores(idx == i) = data.predScores; 
        if ~isempty(data.predClass)
            predClass(idx == i) = data.predClass; 
        end
        % Store importance data (if better than before)
        results = caramel_assessPerf(data); 
        if ismember(MLtype, regML) && results.corrStat > r
            best_model = model; 
        elseif ismember(MLtype, classML) && results.accuracy > accuracy
            best_model = model; 
        end
        % Update progress bar (if prompted)
        if verbose
            progressbar(i/k)
        end
    end
    
    % Define data structure to assess performance
    data_struct = caramel_processInteractions(is, 'Indices', applyk); 
    data_struct.predScores = predScores; 
    data_struct.predClass = predClass; 
    data_struct.MLmethod = data.MLmethod; 

    % Assess performance
    cv_struct = caramel_assessPerf(data_struct, ...
        'ShowPlots', display, 'DisplayResults', display); 
    
    % End message
    if verbose 
        fprintf('%d-fold cross-validation completed. \n', k)
    end
    
    % End timer (if needed)
    if strcmpi(timer, 'on')
        tEnd = toc(tStart);
        [m, s] = deal(floor(tEnd/60), round(rem(tEnd, 60)));
        fprintf('Elapsed time is %d minute(s) and %d second(s) \n', m, s)
    end
    
end