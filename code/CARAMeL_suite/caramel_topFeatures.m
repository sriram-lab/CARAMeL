function t_importance = caramel_topFeatures(model, varargin)

% DESCRIPTION: 
% This function extracts the top features for a CARAMeL model.
% 
% STEPS: 
% 1. Parse through inputs
% 2. Ascertain input compatibility
% 3. Re-order and extract features based on importance
%
% Author:   Carolina H. Chung
% Contact:  chechung@umich.edu
  
% I/O/U
%{
REQUIRED INPUTS: 
    1. model:           CARAMeL model
    
OPTIONAL INPUTS: 
    1. N:               number of top features to return 
       -> default:      N = []

    2. Variance:        amount of variance explained by returned features
       -> note:         provide a double within the (0, 1) range
       -> default:      Variance = 0.95

    3. Measure:         type of measure used to determine importance score
       -> note:         only relevant when MLmethod = 'regRF' or 'classRF'
       -> valid values: {'Accuracy','MSE'}
       -> default:      Measure = 'MSE'
    
OUTPUTS: 
    1. t_importance:    table variable storing the top features and
                        corresponding importance scores
       -> Feature:      top feature names
       -> Score:        top importance scores
       -> Explained:    cumulative varaiance explained
       -> Rank:         rank of top features
          -> note:      only returned when Variance is used for extraction

EXAMPLE USAGE: 
    1. Extract features using defaults given 'model':  
       >> t_importance = caramel_topFeatures(model); 

    2. Extract top 50 features: 
       >> t_importance = caramel_topFeatures(model, 'N', 50); 

    3. Extract top features explaining 75% of the variance: 
       >> t_importance = caramel_topFeatures(model, 'Variance', 0.75); 

    4. Extract features based on Accuracy measure: 
       >> t_importance = caramel_topFeatures(model, 'Measure', 'Accuracy'); 
%}

%% PARSE THROUGH INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create an inputParser object
    p = inputParser; 
    if rem(length(varargin), 2) ~= 0    % for additional parameters
        error('Imbalanced number of name-value pair arguments provided.')
    end
    
    % Define custom input validation function(s)
    isValidModel = @(x) contains(class(x), {'Ensemble','SVM'}) || ... 
        (isstruct(x) && all(isfield(x, {'importance', 'features'})));  
    isValidInteger = @(x) isnumeric(x) && rem(x, 1) == 0; 
    isValidVariance = @(x) isnumeric(x) && x > 0 && x < 1;
    isValidMeasure = @(x) any(validatestring(x, {'MSE', 'Accuracy'}));
    
    % Define required input variables
    addRequired(p, 'model', isValidModel)
    
    % Define optional input parameters
    addParameter(p, 'N', [], isValidInteger)
    addParameter(p, 'Variance', 0.95, isValidVariance)
    addParameter(p, 'Measure', 'MSE', isValidMeasure)
    
    % Parse through inputs
    parse(p, model, varargin{:})
    
    % Define inputs by name
    model       = p.Results.model; 
    N           = p.Results.N; 
    variance    = p.Results.Variance; 
    measure     = p.Results.Measure; 

%% ASCERTAIN INPUT COMPATIBILITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check model fields (for regRF or classRF models)
    if all(isfield(model, {'importance','features'}))
        if size(model.importance, 1) ~= numel(model.features)
            error('Model importance and feature sizes do not match.')
        end
    end

%% RE-ORDER AND EXTRACT FEATURES BASED ON IMPORTANCE %%%%%%%%%%%%%%%%%%%%%%

    % Create table and re-order
    if all(isfield(model, {'importance','features'}))
        Feature = model.features; 
        if strcmpi(measure, 'Accuracy')
            Score = model.importance(:, 1); 
        elseif strcmpi(measure, 'MSE')
            Score = model.importance(:, 2); 
        end 
    else
        Feature = model.PredictorNames(:); 
        Score = predictorImportance(model)';
    end
    t = table(Feature, Score); 
    t = sortrows(t, 'Score', 'descend'); 
    
    % Extract top features
    if isempty(N)
        t.Explained = cumsum(t.Score) / sum(t.Score); 
        idx = find(t.Explained > variance, 1); 
    else
        idx = N; 
    end
    t_importance = t(1:idx, :); 
    t_importance.Rank = transpose(1:idx); 

end