function [phenotypeData, jointProfiles, interactionData] = ...
    caramel_featurize(phenotype_structure, interaction_structure, varargin)

% DESCRIPTION: 
% Script for determining features for a CARAMeL model. Comprised of sigma, 
% delta, entropy, and time scores.
% 
% STEPS: 
% 1. Parse through inputs
% 2. Ascertain input compatibility
% 3. Filter and process data
% 4. Define joint profiles (CARAMeL features)
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

OPTIONAL INPUTS: 
    1.  Key:                    cell array containing key information 
                                that matches conditions (from
                                phenotype_structure) to name values for  
                                interactions (from interaction_structure)
        -> note:                make sure 1st column corresponds to entries
                                in interaction names in structure variable
        -> default:             Key = []

    2.  Verbose:                boolean specifying whether to print 
                                messages while algorithm is run
        -> default:             Verbose = false

OUTPUTS: 
    1. phenotypeData:       table variable containing the phenotype data 
                            used to determine joint profiles (CARAMeL
                            features)

    2. jointProfiles:       table variable containing CARAMeL feature names
                            in the first column and the feature information
                            for all queried drug combinations

EXAMPLE USAGE: 
    1. Determine CARAMeL features and corresponding phenotype datagiven PS 
       as the phenotype data structure and IS as the drug interaction data 
       structure: 
       >> [phenotypeData, jointProfiles] = caramel_featurize(PS, IS); 

    2. Do the same as example 1, but provide a Key to match between 
       phenotype conditions and drug names and only return features: 
       >> [~, jointProfiles] = caramel_featurize(PS, IS, 'Key', Key); 
%}

%% PARSE THROUGH INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create an inputParser object
    p = inputParser; 
    if rem(length(varargin), 2) ~= 0
        error('Imbalanced number of name-value pair arguments provided.')
    end
    
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
    isValidKey = @(x) iscell(x) && size(x, 2) == 2 && ...
        any(ismember(phenotype_structure.conditions, x));
    isValidBoolean = @(x) islogical(x) || ...
        (isnumeric(x) && (x == 0 || x == 1));
    
    % Define optional input parameters
    addParameter(p, 'Key', [], isValidKey)
    addParameter(p, 'Verbose', false, isValidBoolean)
    
    % Parse through inputs
    parse(p, phenotype_structure, interaction_structure, varargin{:})
    
    % Define inputs by short names
    ps                  = phenotype_structure;
    is                  = interaction_structure;
    key                 = p.Results.Key;
    verbose             = logical(p.Results.Verbose);

%% ASCERTAIN INPUT COMPATIBILITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Make sure at least one drug interaction is accounted by PS
    if isempty(key) && ~any(ismember(unique(is.names), ps.conditions)) 
        error('No intersection between interaction and phenotype data.')
    elseif ~isempty(key) && ~any(ismember(ps.conditions, key)) && ...
            ~any(ismember(unique(is.names), key))
        error('No intersection between interaction / phenotype / key data.')
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
        warning('Not able to use all interactions.')
    end
    interactionNames = combos(ix, :); 

    % Filter other interaction data fields 
    if isfield(is, 'scores')
        interactionScores = is.scores(ix, :); 
    else
        interactionScores = []; 
    end
    if isfield(is, 'class')
        interactionClass = is.class(ix, :); 
    else
        interactionClass = []; 
    end
    if isfield(is, 'time')
        interactionTime = is.time(ix, :); 
    else
        interactionTime = []; 
    end

    % Re-define interaction data
    interactionData = struct(); 
    interactionData.names = interactionNames; 
    interactionData.scores = interactionScores; 
    interactionData.class = interactionClass; 
    interactionData.time = interactionTime; 
    
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

end