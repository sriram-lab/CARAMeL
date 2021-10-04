function gem_table = caramel_features2gem(t_importance, gem, data, varargin)

% DESCRIPTION: 
% This function maps top CARAMeL model features to GEM reaction data.
% 
% STEPS: 
% 1. Parse through inputs
% 2. Ascertain input compatibility
% 3. Map top features to GEM data
% 4. Display message regarding other feature types (if present)
%
% Author:   Carolina H. Chung
% Contact:  chechung@umich.edu
  
% I/O/U
%{
REQUIRED INPUTS: 
    1. t_importance:    top features from CARAMeL model
       -> note:         output from 'extract_importance' function

    2. gem:             GEM model to map to
    
OPTIONAL INPUTS: 
    1. Normalize:       Boolean whether to normalize total importance
                        scores in output
       -> default:      Normalize = true
    
OUTPUTS: 
    1. gem_table:       table variable storing GEM information
                        corresponding to top CARAMeL model features
       -> rxn:          top feature names
       -> rxnName:      top importance scores
       -> subSystem:    cumulative varaiance explained
       -> Score:        total importance score

EXAMPLE USAGE: 
    1. Map top features in t_importance to GEM (w/ normalized scores):  
       >> gem_table = caramel_features2gem(t_importance, gem); 

    2. Map top features in t_importance to GEM (w/t normalizing scores):  
       >> gem_table = caramel_features2gem(t_importance, gem, ...
                                          'Normalize', false); 
%}

%% PARSE THROUGH INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create an inputParser object
    p = inputParser; 
    if rem(length(varargin), 2) ~= 0    % for additional parameters
        error('Imbalanced number of name-value pair arguments provided.')
    end
    
    % Define custom input validation function(s)
    isValidImp = @(x) istable(x) && all(ismember({'Feature','Score'}, ...
        x.Properties.VariableNames)); 
    isValidGEM = @(x) isstruct(x) && isfield(x, 'rxns'); 
    isValidData = @(x) isstruct(x) && all(isfield(x, ...
        {'interactionScores', 'jointProfiles'}));
    
    % Define required input variables
    addRequired(p, 't_importance', isValidImp)
    addRequired(p, 'gem', isValidGEM)
    addRequired(p, 'data', isValidData)
    
    % Parse through inputs
    parse(p, t_importance, gem, data)
    
    % Define inputs by name
    t       = p.Results.t_importance; 
    gem     = p.Results.gem; 
    data    = p.Results.data; 

%% ASCERTAIN INPUT COMPATIBILITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check whether feature names match to GEM reactions
    sd_features = t.Feature(contains(t.Feature, {'sigma','delta'})); 
    rxns = cellfun(@(x) x(11:end), sd_features, 'UniformOutput', false); 
    if ~any(ismember(rxns, gem.rxns))
        error('None of the top features map to GEM reactions.')
    end
    
    % Check that interaction number matches number of joint profiles
    features = data.jointProfiles.Feature; 
    p = table2array(data.jointProfiles(:, 2:end)); 
    if size(p, 2) ~= numel(data.interactionScores)
        error('Interaction number and joint profile number do not match.')
    end
    if ~any(ismember(features, sd_features))
        error('None of top features match joint profile features.')
    end
    
    % Check if any interactions are sequential
    if any(p(strcmpi(features, 'time'), :) > 0)
        ix = p(strcmpi(features, 'time'), :) > 0; 
        scores = struct(); 
        scores.sim = data.interactionScores(~ix); 
        scores.seq = data.interactionScores(ix); 
    else
        scores = data.interactionScores; 
    end

%% MAP TOP FEATURES TO GEM DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Map features to GEM reactions
    [rxn, ~, idx] = intersect(rxns, gem.rxns, 'stable'); 
    
    % Extract other reaction information (if available)
    try
        rxnName = gem.rxnNames(idx); 
    catch
        warning('Could not map features to reaction names.')
    end
    try
        subSystem = gem.subSystems(idx); 
    catch
        warning('Could not map features to reaction subsystems.')
    end
    
    % Determine total importance scores
    Pvalue = nan(size(rxn)); Direction = cell(size(rxn)); 
    for i = 1:numel(rxn)
        ix1 = find(endsWith(features, rxn{i})); 
        if isstruct(scores)
            [psim, pseq] = deal(1, 1); 
        else
            psim = 1; 
        end
        for j = 1:numel(ix1)
            ix2 = logical(p(ix1(j), :)); 
            % accout for both
            if isstruct(scores) 
                [ix3, ix4] = deal(ix2(~ix), ix2(ix)); 
                [x1, y1] = deal(scores.sim(ix3), scores.sim(~ix3)); 
                [x2, y2] = deal(scores.seq(ix4), scores.seq(~ix4)); 
                [~, p1, ~, stats1] = ttest2(x1, y1, 'VarType', 'unequal'); 
                [~, p2, ~, stats2] = ttest2(x2, y2, 'VarType', 'unequal');
                if p1 < psim
                    psim = p1; 
                end
                if p2 < pseq
                    pseq = p2; 
                end
            % simultaneous only    
            else 
                [x, y] = deal(scores(ix2), scores(~ix2)); 
                [~, p, ~, stats] = ttest2(x, y, 'VarType', 'unequal'); 
                if p < psim
                    psim = p; 
                end
            end
        end
        if isstruct(scores)
            % both significant
            if (psim < 0.05) && (pseq < 0.05)
                Pvalue(i) = max([psim pseq]); 
                if (stats1.tstat > 0) && (stats2.tstat > 0)
                    Direction(i) = {'A/R'}; 
                elseif (stats1.tstat < 0) && (stats2.tstat < 0)
                    Direction(i) = {'Sy/Se'}; 
                elseif (stats1.tstat < 0) && (stats2.tstat > 0)
                    Direction(i) = {'Sy/R'}; 
                else
                    Direction(i) = {'A/Se'}; 
                end
            % simultaneous only
            elseif (psim < 0.05)
                Pvalue(i) = psim; 
                if stats1.tstat > 0
                    Direction(i) = {'A'}; 
                else
                    Direction(i) = {'Sy'}; 
                end
            % sequential only
            elseif (pseq < 0.05)
                Pvalue(i) = pseq; 
                if stats2.tstat > 0
                    Direction(i) = {'R'}; 
                else
                    Direction(i) = {'Se'}; 
                end
            % neither
            else
                Pvalue(i) = min([psim pseq]); 
                Direction(i) = {'NS'}; 
            end
        else
            Pvalue(i) = psim; 
            if Pvalue(i) < 0.05
                if stats.tstat < 0
                    Direction(i) = {'Sy'}; 
                elseif stats.tstat > 0
                    Direction(i) = {'A'}; 
                end
            else
                Direction(i) = {'NS'}; 
            end
        end
%         idx = logical(p(endsWith(features, rxn{i}), :));
%         if size(idx, 1) > 1
%             idx = logical(sum(idx)); 
%         end
%         x = data.interactionScores(idx); 
%         y = data.interactionScores(~idx); 
%         [~, Pvalue(i), ~, stats] = ttest2(x, y, 'VarType', 'unequal'); 
%         if Pvalue(i) < 0.05
%             if stats.tstat < 0
%                 Direction(i) = {'Synergy'}; 
%             elseif stats.tstat > 0
%                 Direction(i) = {'Antagonism'}; 
%             end
%         else
%             Direction(i) = {'NS'}; 
%         end
    end
    
    % Define output
    if exist('rxnName', 'var') && exist('subSystem', 'var')
        gem_table = table(rxn, rxnName, subSystem, Pvalue, Direction); 
    elseif exist('rxnName', 'var') && ~exist('subSystem', 'var')
        gem_table = table(rxn, rxnName, Pvalue, Direction); 
    elseif ~exist('rxnName', 'var') && exist('subSystem', 'var')
        gem_table = table(rxn, subSystem, Pvalue, Direction); 
    else
        gem_table = table(rxn, Pvalue, Direction); 
    end
    gem_table = sortrows(gem_table, 'Pvalue', 'ascend'); 
    gem_table.Rank = transpose(1:numel(Pvalue));

%% DISPLAY MESSAGE REGARDING OTHER FEATURE TYPES (IF PRESENT) %%%%%%%%%%%%%

    % Entropy
    if any(contains(t.Feature, 'entropy'))
        disp('Entropy detected as a top feature:')
        disp(t(contains(t.Feature, 'entropy'), :))
    end
    
    % Time
    if any(strcmp(t.Feature, 'time'))
        disp('Time detected as a top feature:')
        disp(t(strcmp(t.Feature, 'time'), :))
    end

end