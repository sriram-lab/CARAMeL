function stat_table = caramel_anova(d_struct, g_table, varargin)

% DESCRIPTION: 
% This function conducts a one-way ANOVA test for all top reactions
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
    isValidDstr = @(x) isstruct(x) && all(isfield(x, {'phenotypeData', ...
        'interactionNames','interactionScores','interactionClass'})); 
    isValidGtab = @(x) istable(x) && ismember({'rxn'}, ...
        x.Properties.VariableNames); 
    isValidGEM = @(x) isstruct(x) && isfield(x, 'rxns'); 
    isValidType = @(x) any(validatestring(x, {'sim', 'seq'}));
    
    % Define required input variables
    addRequired(p, 'd_struct', isValidDstr)
    addRequired(p, 'g_table', isValidGtab)
    
    % Define required input variables
    addParameter(p, 'GEM', [], isValidGEM)
    addParameter(p, 'Type', 'sim', isValidType)
    
    % Parse through inputs
    parse(p, d_struct, g_table, varargin{:})
    
    % Define inputs by name
    sD      = p.Results.d_struct; 
    tG      = p.Results.g_table; 
    gem     = p.Results.GEM; 
    type    = p.Results.Type; 

%% ASCERTAIN INPUT COMPATIBILITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Define phenotype data parts
    tP = sD.phenotypeData; 
    [pLab, pCond] = deal(tP.Label, tP.Properties.VariableNames(2:end)); 

    % Check that all conditions are present in phenotype data
    conditions = unique(sD.interactionNames, 'stable'); 
    conditions(cellfun(@isempty, conditions)) = [];
    if ~all(ismember(conditions, pCond))
        error('Not all conditions represented in phenotype data.')
    end
    
    % Check that all reactions are present in phenotype data
    if sum(contains(pLab, tG.rxn)) < numel(tG.rxn)
        error('Not all GEM reactions represented in phenotype data.')
    end
    
    % Define phenotype data 
    pData = table2array(tP(contains(pLab, tG.rxn), 2:end)); 
    pLab = pLab(contains(pLab, tG.rxn)); 
    pRxn = cellfun(@(x) x(5:end), pLab, 'UniformOutput', false); 

%% MAP TOP FEATURES TO GEM DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Iterate through each interaction
    [names, scores] = deal(sD.interactionNames, sD.interactionScores); 
    pstats = [];
%     progressbar('Conducting one-way ANOVA...')
    for i = 1:numel(tG.rxn)
        idx = find(strcmp(pRxn, tG.rxn(i)));
        for j = 1:numel(idx)
            d = zeros(size(names)); 
            for k = 1:size(d, 1)
                [~, ~, cIdx] = intersect(names(k, :), pCond, 'stable');
                d(k, 1:numel(cIdx)) = pData(idx(j), cIdx);
            end
            if strcmpi(type, 'sim')
                g = sum(d, 2); 
            elseif strcmpi(type, 'seq')
                g = zeros(size(d, 1), 1);  
                g(d(:, 1) == 1 & d(:, 2) == 0) = 1; 
                g(d(:, 1) == 0 & d(:, 2) == 1) = 2; 
                g(sum(d, 2) == 2) = 3; 
            end
            sDir = {pLab{idx(j)}(1:3)}; 
            if numel(unique(g)) == 1
                continue
            end
            [~, ~, stats] = anova1(scores, g, 'off'); figure
%             subplot(4, 3, i)
            results = multcompare(stats, 'Display', 'on'); 
            if any(results(:, end) < 0.05) && isempty(pstats)
                rIdx = results(:, end) < 0.05; 
                rxn = repmat(tG.rxn(i), sum(rIdx), 1);
                if ~isempty(gem) && ...
                        all(isfield(gem, {'rxnNames','subSystems'}))
                    gIdx = strcmp(gem.rxns, tG.rxn(i)); 
                    rxnName = repmat(gem.rxnNames(gIdx), sum(rIdx), 1); 
                    subSystem = repmat(gem.subSystems(gIdx), sum(rIdx), 1); 
                end
                dir = repmat(sDir, sum(rIdx), 1); 
                n = repmat(numel(rIdx), sum(rIdx), 1); 
                gN = repmat(numel(unique(g)), sum(rIdx), 1); 
                pstats = results(rIdx, :); 
                if strcmpi(type, 'seq') || any(g == 0)
                    pstats(:, 1:2) = pstats(:, 1:2) - 1; 
                end
            elseif any(results(:, end) < 0.05)
                rIdx = results(:, end) < 0.05; 
                rxn = vertcat(rxn, repmat(tG.rxn(i), sum(rIdx), 1));
                if ~isempty(gem) && ...
                        all(isfield(gem, {'rxnNames','subSystems'}))
                    gIdx = strcmp(gem.rxns, tG.rxn(i)); 
                    rxnName = vertcat(rxnName, ...
                        repmat(gem.rxnNames(gIdx), sum(rIdx), 1)); 
                    subSystem = vertcat(subSystem, ...
                        repmat(gem.subSystems(gIdx), sum(rIdx), 1)); 
                end
                dir = vertcat(dir, repmat(sDir, sum(rIdx), 1)); 
                n = [n; repmat(numel(rIdx), sum(rIdx), 1)]; 
                gN = [gN; repmat(numel(unique(g)), sum(rIdx), 1)]; 
                p = results(rIdx, :); 
                if strcmpi(type, 'seq') || any(g == 0)
                    p(:, 1:2) = p(:, 1:2) - 1; 
                end
                pstats = [pstats; p]; 
            end
        end
%         progressbar(i / numel(tG.rxn))
    end
    if strcmpi(type, 'sim')
        g = repmat({'Synergy'}, size(pstats, 1), 1); 
        g(pstats(:, 4) > 0) = {'Antagonism'}; 
    else
        g = repmat({'Sensitivity'}, size(pstats, 1), 1); 
        g(pstats(:, 4) > 0) = {'Resistance'}; 
    end
    
    % Define output
    if ~isempty(gem) && all(isfield(gem, {'rxnNames', 'subSystems'}))
        t1 = table(rxn, rxnName, subSystem, dir, n, g, gN); 
    else
        t1 = table(rxn, dir, n, g, gN); 
    end
    t2 = array2table(pstats, 'VariableNames', {'g1','g2','l','x','u','p'}); 
    stat_table = horzcat(t1, t2); 
    stat_table = sortrows(stat_table, {'rxn', 'p'}, {'ascend', 'ascend'}); 

end