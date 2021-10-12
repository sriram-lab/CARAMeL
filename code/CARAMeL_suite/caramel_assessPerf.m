function results_struct = caramel_assessPerf(data_struct, varargin)

% DESCRIPTION:
% This script assessed CARAMeL model performance. 
% 
% STEPS:
% 1. Parse through inputs
% 2. Ascertain input compatibility
% 3. Assess correlation (if MLmethod = regression)
% 4. Assess confusion matrix (if MLmethod = classification or if possible)
% 5. Assess ROC curves (if MLmethod = classification or if possible)
% 6. Display output (if prompted)
% 
% Author:   Carolina H. Chung
% Contact:  chechung@umich.edu
    
% I/O/U
%{
REQURED INPUTS: 
    1. data_struct:         structure variable returned after predicting
                            interaction outcomes using a trained CARAMeL
                            model

OPTIONAL INPUTS: 
    1. Threshold:           2-element numeric vector specifying thresholds
                            to classify interactions based on their scores
       -> default:          Threshold = []

    2. InteractionType:     string specifying which class to determine
                            further classification-based performance stats
       -> choose from:      {'Synergy', 'Additivity', 'Antagonism'}
       -> default:          InteractionType = []

    3. Verbose:             boolean specifying whether to print messages
                            while algorithm is run
       -> default:          Verbose = false

    4. ShowPlots:           boolean specifying whether to show plots that
                            visualize model performance
       -> default:          ShowPlots = false

    5. DisplayResults:      boolean specifying whether to display final
                            results of CARAMeL model predictions
       -> default:          DisplayResults = false

OTHER FUNCTION PARAMETERS: 
    1. corr:                any additional parameters available with the
                            built-in MATLAB function 'corr'
    
OUTPUTS: 
    1. results_struct:      structure variable containing final results of 
                            CARAMeL model performance

EXAMPLE USAGE:  
    1. Return performance results for a CARAMeL model: 
       >> results = caramel_assessPerf(data_struct); 

    2. Same as (1), but also determine synergy-specific stats: 
       >> results = caramel_assessPerf(data_struct, ...
                                      'InteractionType', 'Synergy'); 
%}

%% PARSE THROUGH INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create an inputParser object
    p = inputParser; 
    if rem(length(varargin), 2) ~= 0    % for additional parameters
        error('Imbalanced number of name-value pair arguments provided.')
    end
    
    % Define interaction classes
    classes = {'Synergy','Additivity','Antagonism'};
    
    % Define custom input validation function(s)
    isValidData = @(x) isstruct(x) && all(isfield(x, ...
        {'interactionNames','interactionScores','predScores', 'MLmethod'}));
    isValidThreshold = @(x) isnumeric(x) && numel(x) == 2; 
    isValidInterType = @(x) any(validatestring(x, classes));
    isValidBoolean = @(x) islogical(x) || ...
        (isnumeric(x) && (x == 1 || x == 0)); 
    
    % Define required input variables
    addRequired(p, 'data_struct', isValidData)
    
    % Define optional input parameters
    addParameter(p, 'Threshold', [], isValidThreshold)
    addParameter(p, 'InteractionType', [], isValidInterType)
    addParameter(p, 'Verbose', false, isValidBoolean)
    addParameter(p, 'ShowPlots', false, isValidBoolean)
    addParameter(p, 'DisplayResults', false, isValidBoolean)
    
    % Parse through inputs
    p.KeepUnmatched = true; 
    parse(p, data_struct, varargin{:})
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
    d               = p.Results.data_struct; 
    thresh          = p.Results.Threshold; 
    interType       = p.Results.InteractionType;
    verbose         = logical(p.Results.Verbose);
    showPlots       = logical(p.Results.ShowPlots); 
    dispResults     = logical(p.Results.DisplayResults);
    
    % Define list of optional input parameters from CORR function
    corrParams = {'Type','Rows','Tail'}; 
    
    % Allocate extra parameters to CORR function (if applicable)
    if ~isempty(extraParams)
        extraNames = fieldnames(extraParams); 
        [s, p_idx] = intersect(lower(extraNames), lower(corrParams)); 
        corrVarargin = cell(1, 2*numel(p_idx));
        for j = 1:numel(p_idx)
            [i1, i2] = deal(2*(j-1)+1, 2*j); 
            corrVarargin{i1} = extraNames{p_idx(j)};
            eval(['corrVarargin{i2} = extraParams.' corrVarargin{i1} ';'])
        end
        if ~ismember('type', s)
            corrVarargin{end+1} = 'Type'; 
            corrVarargin{end+1} = 'Spearman';
        end
        if verbose
            disp('Extra input(s) passed to CORR function.')
        end
    else
        corrVarargin{1} = 'Type'; 
        corrVarargin{2} = 'Spearman';
    end
    
%% ASCERTAIN INPUT COMPATIBILITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check data structure fields 
    if strcmpi(d.MLmethod, 'regression')
        if ~all(isfield(d, {'interactionScores','predScores'}))
            error('Actual or predicted interaction scores missing.')
        elseif numel(d.interactionScores) ~= numel(d.predScores)
            error('Actual and predicted score vector sizes do not match.')
        end
    else
        if ~all(isfield(d, {'interactionClass','predClass'}))
            error('Actual or predicted interaction classes missing.')
        elseif numel(d.interactionClass) ~= numel(d.predClass)
            error('Actual and predicted class vector sizes do not match.')
        end
    end

%% ASSESS CORRELATION (if MLmethod = regression) %%%%%%%%%%%%%%%%%%%%%%%%%%

    if strcmpi(d.MLmethod, 'regression')
        [corrStat, corrPvalue] = ...
            corr(d.interactionScores, d.predScores, corrVarargin{:});
    else
        corrStat = NaN; corrPvalue = NaN;
        if verbose
            disp('Correlation irrelevant (MLmethod = classification).')
        end
    end

%% ASSESS CONFUSION MATRIX (if MLmethod = classification or if possible) %%

    % Classify interactions (if possible)
    if strcmpi(d.MLmethod, 'regression') && isempty(d.interactionClass)
        if ~isempty(thresh)
            d.interactionClass = ...
                caramel_classify(d.interactionScores, thresh); 
            d.predClass = caramel_classify(d.predScores, thresh); 
        else
            if verbose
                fprintf(strcat("Specified ML method is regression and ", ...
                "threshold information was not provided. \n Exiting ", ...
                "without determining confusion matrix nor related stats.\n"))
            end
            results_struct = struct(...
            'interactionNames', {d.interactionNames}, ...
            'interactionScores', d.interactionScores, ...
            'predScores', d.predScores, ...
            'MLmethod', d.MLmethod, ...
            'corrStat', corrStat, ...
            'corrPvalue', corrPvalue);
        end
        if showPlots
            caramel_plot(results_struct, 'scatter'); 
        end
        if dispResults
            disp(results_struct)
        end
        return
    end

    % Determine confusion matrix
    try
        [confMat, order] = confusionmat(d.interactionClass, ...
            d.predClass, 'Order', classes);
        [~, confMatStats] = confusion.getValues(confMat); 
    catch
%         [confMat, order] = confusionmat(d.interactionClass, d.predClass); 
        confMatStats = [];
    end
%     [~, confMatStats] = confusion.getValues(confMat); 
    
    % Assess class-specific stats (if prompted)
    if ~isempty(interType) && ~isempty(confMatStats)
        index = find(strcmp(order, interType));
        try
            Accuracy = sum(confMatStats.AccuracyInTotal);
            Error = sum(confMatStats.ErrorInTotal);
            Sensitivity = confMatStats.Sensitivity(index);
            Specificity = confMatStats.Specificity(index);
            Precision = confMatStats.Precision(index);
            FalsePositiveRate = confMatStats.FalsePositiveRate(index);
            F1_score = confMatStats.F1_score(index);
            MatthewsCorrelationCoefficient = ...
                confMatStats.MatthewsCorrelationCoefficient(index);
        catch
            Accuracy = sum(confMatStats.AccuracyInTotal);
            Error = sum(confMatStats.ErrorInTotal);
            Sensitivity = NaN;
            Specificity = NaN;
            Precision = NaN;
            FalsePositiveRate = NaN;
            F1_score = NaN;
            MatthewsCorrelationCoefficient = NaN;
            if verbose
                warning('Could not determine stats for %s.', interType)
            end
        end
    end
    
%% ASSESS ROC CURVES (if MLmethod = classification or if possible) %%%%%%%%

    if strcmpi(d.MLmethod, 'regression')
        try
            [synergyAUC, antagonismAUC] = ...
                caramel_auroc(d.interactionClass, d.MLmethod, d.predScores); 
        catch
            [synergyAUC, antagonismAUC] = deal(NaN, NaN); 
            warning('Could not calculate AUROC values.')
        end
    else
        try
            [synergyAUC, antagonismAUC] = ...
                caramel_auroc(d.interactionClass, d.MLmethod, d.predClass); 
        catch
            [synergyAUC, antagonismAUC] = deal(NaN, NaN); 
            warning('Could not calculate AUROC values.')
        end
    end

%% DEFINE OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Account for missing data
    if ~isfield(d, 'interactionScores')
        d.interactionScores = [];
    end
    if ~isfield(d, 'predScores')
        d.predScores = [];
    end
    
    % Define results structure
    if ~isempty(interType)
        results_struct = struct(...
            'interactionNames', {d.interactionNames}, ...
            'interactionScores', d.interactionScores, ...
            'interactionClass', {d.interactionClass}, ...
            'predScores', d.predScores, ...
            'predClass', {d.predClass}, ...
            'MLmethod', d.MLmethod, ...
            'corrStat', corrStat, ...
            'corrPvalue', corrPvalue, ...
            'confMatStats', confMatStats, ...
            'accuracy', Accuracy, ...
            'error', Error, ...
            'sensitivity', Sensitivity, ...
            'specificity', Specificity, ...
            'precision', Precision, ...
            'falsePositiveRate', FalsePositiveRate, ...
            'f1_score', F1_score, ...
            'MCC', MatthewsCorrelationCoefficient, ...
            'synergyAUC', synergyAUC, ...
            'antagonismAUC', antagonismAUC);
    else
        results_struct = struct(...
            'interactionNames', {d.interactionNames}, ...
            'interactionScores', d.interactionScores, ...
            'interactionClass', {d.interactionClass}, ...
            'predScores', d.predScores, ...
            'predClass', {d.predClass}, ...
            'MLmethod', d.MLmethod, ...
            'corrStat', corrStat, ...
            'corrPvalue', corrPvalue, ...
            'confMatStats', confMatStats, ...
            'synergyAUC', synergyAUC, ...
            'antagonismAUC', antagonismAUC);
    end
    
%% DISPLAY OUTPUT (IF PROMPTED) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Plots
    if showPlots
        if strcmpi(d.MLmethod, 'regression')
            if all(ismember(classes, d.interactionClass))
                type = {'scatter','confusion','box','roc'};
            else
                type = {'scatter','confusion'};
            end
        else
            if all(ismember(classes, d.interactionClass))
                type = {'confusion','roc'};
            else
                type = {'confusion'};
            end
        end
        caramel_plot(results_struct, type);
    end
    
    % Results
    if dispResults
        disp(results_struct)
    end

end