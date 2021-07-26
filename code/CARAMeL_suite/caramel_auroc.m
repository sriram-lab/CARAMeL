function [synAUC, antAUC] = caramel_auroc(interClass, MLmethod, pred, ...
    varargin)

% DESCRIPTION: 
% This function returns the area under the receiver operating curve (AUROC)
% for classifying interactions as synergistic or antagonistic.
% 
% STEPS: 
% 1. Parse through inputs
% 2. Ascertain input compatibility
% 3. Determine AUROC values
%
% Author:   Carolina H. Chung
% Contact:  chechung@umich.edu
  
% I/O/U
%{
REQUIRED INPUTS: 
    1. interClass:  classes for drug interactions

    2. MLmethod:    string specifying what ML type was used
       -> accepted: {'regression','classification'}

    3. pred:        CARAMeL model predictions
       -> note:     provide 'predScores' if MLtype = 'regression'
                    provide 'predClass' if MLtype = 'classification'

OPTIONAL INPUTS: 
    1. Reverse:     Boolean value whether to flip directionality used to
                    classify synergy and antagonism based on scores
       -> note:     relevant when MLtype = 'regression'
       -> default:  Reverse = false
    
OUTPUTS: 
    1. synAUC:      area under the curve (AUC) for synergy ROC

    2. antAUC:      area under the curve (AUC) for antagonism ROC

EXAMPLE USAGE: 
    1. Determine synergy and antagonism AUROC for regression model: 
       >> [synAUC, antAUC] = caramel_auroc(interClass, 'regression', ...
                                          predScores);

    2. Determine synergy and antagonism AUROC for classification model: 
       >> [synAUC, antAUC] = caramel_auroc(interClass, 'classification', ...
                                          predClass);
%}

%% PARSE THROUGH INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create an inputParser object
    p = inputParser; 
    if rem(length(varargin), 2) ~= 0    % for additional parameters
        error('Imbalanced number of name-value pair arguments provided.')
    end
    
    % Define custom input validation function(s)
    isValidClass = @(x) iscell(x) || iscategorical(x) || ...
        isnumeric(x) || islogical(x);
    isValidMLmethod = @(x) ismember(x, {'regression', 'classification'});
    isValidPred = @(x) isValidClass(x) || isnumeric(x); 
    isValidBoolean = @(x) islogical(x) || ...
        (isnumeric(x) && (x == 0 || x == 1)); 
    
    % Define required input variables
    addRequired(p, 'interClass', isValidClass)
    addRequired(p, 'MLmethod', isValidMLmethod)
    addRequired(p, 'pred', isValidPred)
    
    % Define optional input parameters
    addParameter(p, 'Reverse', false, isValidBoolean)
    
    % Parse through inputs
    parse(p, interClass, MLmethod, pred, varargin{:})
    
    % Define inputs by name
    interClass  = p.Results.interClass; 
    MLmethod    = p.Results.MLmethod;
    pred        = p.Results.pred; 
    Reverse     = p.Results.Reverse; 

%% ASCERTAIN INPUT COMPATIBILITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check interClass and predScores sizes
    if numel(interClass) ~= numel(pred)
        error('Class and prediction vector sizes do not match.')
    end
    
    % Check interClass contents
    if ~any(ismember({'Synergy','Antagonism'}, interClass))
        error('Synergy and Antagonism classes not present in class vector.')
    end

%% DETERMINE AUROC VALUES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Define scores based on model predictions provided
    if strcmpi(MLmethod, 'classification')
        if all(ismember({'Synergy','Antagonism'}, pred))
            [gIdx, gNames] = grp2idx(interClass); 
            scores = nan(size(pred)); 
            for i = 1:max(gIdx)
                idx = strcmpi(pred, gNames{i}); 
                scores(idx) = i; 
            end
            synIdx = find(strcmp(gNames, 'Synergy')); 
            antIdx = find(strcmp(gNames, 'Antagonism')); 
            addIdx = find(~ismember(gNames, {'Synergy','Antagonism'})); 
        elseif any(~ismember({'Synergy','Antagonism'}, pred))
            warning('Too difficult, quitting now')
            return
        end
    else 
        scores = pred; 
    end

    % Define scores for perfcurve function
    if strcmpi(MLmethod, 'regression') && Reverse
        syn_scores = rescale(scores); 
        ant_scores = rescale(-scores); 
    elseif strcmpi(MLmethod, 'regression')
        ant_scores = rescale(scores); 
        syn_scores = rescale(-scores); 
    elseif strcmpi(MLmethod, 'classification')
        syn_scores = scores(:, synIdx) - ...
            max(scores(:, addIdx), scores(:, antIdx)); 
        ant_scores = scores(:, antIdx) - ...
            max(scores(:, addIdx), scores(:, synIdx)); 
    end
    
    % Determine AUC values
    if ismember('Synergy', interClass)
        [~, ~, ~, synAUC] = perfcurve(interClass, syn_scores, 'Synergy');
    else
        warning('Count not determine AUROC for Synergy.')
        synAUC = NaN; 
    end
    if ismember('Antagonism', interClass)
        [~, ~, ~, antAUC] = perfcurve(interClass, ant_scores, 'Antagonism');
    else
        warning('Count not determine AUROC for Antagonism.')
        antAUC = NaN; 
    end

end