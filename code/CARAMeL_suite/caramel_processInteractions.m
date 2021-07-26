function data_struct = caramel_processInteractions(interaction_struct, ...
    varargin)

% DESCRIPTION:
% This script processes interaction data into a standardized format. 
% 
% STEPS:
% 1. Parse through inputs
% 2. Process interaction data
% 
% Author:   Carolina H. Chung
% Contact:  chechung@umich.edu
    
% I/O/U
%{
REQURED INPUTS: 
    1. interaction_structure: 	structure containing at least 2 fields
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
    1. Indices:                 numerical array of indices for which 
                                interactions to determine output
       -> default:              Indices = true(size(is.names, 1), 1)
    
OUTPUTS: 
    1. data_struct:         structure variable containing information
                            processed from function inputs; includes 
                            inputs used for ML model development

EXAMPLE USAGE:  
    1. Conduct 10-fold cross-validation given phenotype structure (PS) and 
       interaction structure (IS) and return best performing model: 
       >> [CV_results, best_model] = caramel_crossValidation(PS, IS); 
%}

%% PARSE THROUGH INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create an inputParser object
    p = inputParser; 
    if rem(length(varargin), 2) ~= 0
        error('Imbalanced number of name-value pair arguments provided.')
    end
    
    % Define required input validation functions
    isValidInteraction = @(x) ...
        isstruct(x) && isfield(x, 'names') && ...
        (isfield(x, 'scores') && size(x.names, 1) == numel(x.scores)) || ...
        (isfield(x, 'class') && size(x.names, 1) == numel(x.class));
    
    % Define required input variables
    addRequired(p, 'interaction_struct', isValidInteraction)
    
    % Define optional input validation functions 
    isValidIndices = @(x) (isnumeric(x) || islogical(x)) && ...
        numel(x) == size(interaction_struct.names, 1); 
    
    % Define optional input parameters
    addParameter(p, 'Indices', [], isValidIndices)
    
    % Parse through inputs
    parse(p, interaction_struct, varargin{:})
    
    % Define inputs by name
    is      = p.Results.interaction_struct; 
    idx     = p.Results.Indices; 
        
%% PROCESS INTERACTION DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Define idx if it's empty
    if isempty(idx)
        idx = true(size(is.names, 1), 1); 
    end
    j = max(sum(~cellfun(@isempty, is.names(idx, :)), 2)); 
    
    % Define interaction outcome data
    if isfield(is, 'scores')
        interactionScores = is.scores(idx);
    end
    if isfield(is, 'class')
        interactionClass = is.class(idx); 
    elseif all(isfield(is, {'scores','threshold'}))
        interactionClass = caramel_classify(is.scores(idx), is.threshold); 
    else
        interactionClass = [];
    end
    
    % Define output
    data_struct = struct(...
        'interactionNames', {is.names(idx, 1:j)}, ...
        'interactionScores', interactionScores, ...
        'interactionClass', {interactionClass});     

end