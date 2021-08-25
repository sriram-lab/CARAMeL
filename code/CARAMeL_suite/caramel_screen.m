function t = caramel_screen(ps, drug_list, model, varargin)

% DESCRIPTION: 
% This function generates outcome predictions for a given list of drugs and
% growth conditions. 
% 
% STEPS: 
% 1. Parse through inputs
% 2. Ascertain input compatibility
% 3. Define combinations + generate predictions
% 4. Construct output variable
%
% Author:   Carolina H. Chung
% Contact:  chechung@umich.edu
  
% I/O/U
%{
REQUIRED INPUTS: 
    1. ps:              phenotype structure
    2. drug_list:       cell array of drugs to screen through
    3. model:           trained CARAMeL model
    
OPTIONAL INPUTS: 
    1. Conditions:      cell array of growth conditions to account for
       -> type:         cell array
       -> default:      [];

    2. Order:           intraction order to assess
       -> type:         integer or array of integers
       -> default:      Order = 2 (pairwise interactions)

    3. Sequential:      Boolean to predict sequential interactions
       -> type:         logical
       -> default:      Sequential = true

    4. Duration:        treatment duration to use for sequential entries
       -> note:         number of elements must match Order value
       -> type:         array of numeric entries
       -> default:      Duration = [14 1]
    
OUTPUTS: 
    1. predictions:     structure variable containing model predictions and
                        other related information

EXAMPLE USAGE: 
    1. Generate predictions using default parameters given 'drug_list':
       >> p = caramel_screen(ps, drug_list, model); 

    2. Generate three-way drug predictions given 'conditions' and 14-14-1 
       treatment durations: 
       >> p = caramel_screen(ps, drug_list, model, ...
                            'Condition', conditions, ...
                            'Order', 3, 'Duration', [14 14 1]); 
%}

%% PARSE THROUGH INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create an inputParser object
    p = inputParser;
    if rem(length(varargin), 2) ~= 0    % for additional parameters
        error('Imbalanced number of name-value pair arguments provided.')
    end
    
    % Define custom input validation function(s)
    isValidPhenotype = @(x) ...
        isstruct(x) && numel(fieldnames(x)) == 4 && ...
        all(isfield(x, {'data','features','conditions','H'})) && ...
        size(x.data, 1) == numel(x.features) && ...
        size(x.data, 2) == numel(x.conditions) && ...
        size(x.data, 2) == numel(x.H);
    isValidList = @(x) ~isempty(x) && (iscell(x) || isstring(x));
    isValidMLmodel = @(x) ~isempty(x); 
    isValidInteger = @(x) isnumeric(x) && rem(x, 1) == 0;
    isValidBoolean = @(x) islogical(x) || ...
        (isnumeric(x) && (x == 1 || x == 0)); 
    isValidDuration = @(x) isnumeric(x) && numel(x) > 1; 
    
    % Define required input variables
    addRequired(p, 'ps', isValidPhenotype)
    addRequired(p, 'drug_list', isValidList)
    addRequired(p, 'model', isValidMLmodel)
    
    % Define required input variables
    addParameter(p, 'Conditions', [], isValidList)
    addParameter(p, 'Order', 2, isValidInteger)
    addParameter(p, 'Sequential', true, isValidBoolean)
    addParameter(p, 'Duration', [14 7], isValidDuration)
    addParameter(p, 'Verbose', true, isValidBoolean)
    
    % Parse through inputs
    p.KeepUnmatched = true; 
    parse(p, ps, drug_list, model, varargin{:})
    if ~isempty(fieldnames(p.Unmatched))
        extraParams = p.Unmatched; 
        if logical(p.Results.Verbose)
            disp('Extra name-value pair arguments provided:')
            disp(p.Unmatched)
        end
    else
        extraParams = {};
    end
    
    % Define inputs by name
    ps      = p.Results.ps;
    d       = p.Results.drug_list; 
    model   = p.Results.model;
    c       = p.Results.Conditions; 
    order   = p.Results.Order; 
    seq     = logical(p.Results.Sequential); 
    t       = p.Results.Duration; 
    verbose = p.Results.Verbose; 
    
    % Define list of optional input parameters from external functions
    caramelParams = {'Key', 'MLtype', 'MLmodel', 'Verbose', ...
        'regRFtrainOptions', 'classRFtrainOptions', 'classRFpredictOptions'};
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
    paramList = {caramelParams, fitrensembleParams, fitcensembleParams, ...
        fitrsvmParams, fitcsvmParams, predictParams}; 
    f = {'caramel','fitrensemble','fitcensemble','fitrsvm','fitcsvm', 'predict'};
    
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
    mainVarargin = cat(1, caramelVarargin, fitrensembleVarargin, ...
        fitcensembleVarargin, fitrsvmVarargin, fitcsvmVarargin, ...
        predictVarargin);
    
    % Display input settings 
    if verbose
        disp('Screening interaction outcomes...')
    end

%% ASCERTAIN INPUT COMPATIBILITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Ascertain compatibility between Order and Duration
    if seq
        if ~exist('order', 'var') || ~exist('t', 'var')
            error('Sequential set to true but Order or Duration not provided.')
        end
        if size(t, 2) ~= order
            error('Duration array size does not match Order specified.')
        end
    end

%% DEFINE COMBINATIONS + GENERATE PREDICTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Define input 
    d_combos = nchoosek(d, order); 
    if isempty(c)
        combos = d_combos; 
    else
        combos = horzcat(repmat(d_combos, numel(c), 1), ...
            repelem(c, size(d_combos, 1), 1)); 
    end
    X = struct('names', {combos}, 'scores', zeros(size(combos, 1), 1));     
    
    % Determine simultaneous predictions
    if verbose
        disp('Predicting simultaneous interactions...')
    end
    r = caramel(ps, X, 'predict', 'MLmodel', model, mainVarargin{:}); 
    Y = r.predScores; 
    
    % Determine sequential predictions (if prompted)
    if seq
        % define time
        if isempty(c)
            X.time = repmat(t, numel(X.scores), 1); 
        else
            X.time = repmat([t NaN], numel(X.scores), 1); 
        end
        % define sequences
        p = permn(1:order, order); 
        idx = all(diff(sort(p, 2), 1, 2), 2); 
        p(~idx, :) = [];
        s = join(string(p), '')'; 
        % generate predictions
        Y = [Y nan(numel(Y), numel(s))];
        for i = 1:size(p, 1)
            Xseq = X; 
            if isempty(c)
                Xseq.names = X.names(:, p(i, :)); 
            else
                Xseq.names = X.names(:, [p(i, :) end]);
            end
            if verbose
                fprintf('Predicting sequential interactions (%d/%d)... \n', ...
                    i, size(p, 1))
            end
            r = caramel(ps, Xseq, 'predict', 'MLmodel', model, mainVarargin{:});
            Y(:, i + 1) = r.predScores; 
        end
    end
    
%% DEFINE OUTPUT

    % Define combination table
    t1 = array2table(combos); 
    col = strcat('D', string(1:order)); 
    if isempty(c)
        t1.Properties.VariableNames = col;
    else
        t1.Properties.VariableNames = horzcat(col, 'Media'); 
    end
    
    % Define prediction table
    if seq
        t2 = array2table(Y);
        t2.Properties.VariableNames = horzcat('sim', strcat('seq_', s)); 
    else
        t2 = table(Y, 'VariableNames', {'sim'}); 
    end
    
    % Define output
    t = horzcat(t1, t2);

end