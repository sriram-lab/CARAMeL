function top_subsystems = caramel_rankSubsystems(gem_table, gem, varargin)

% DESCRIPTION: 
% This function returns a list of sub-systems represented in GEM reactions
% associated with top CARAMeL model features.
% 
% STEPS: 
% 1. Parse through inputs
% 2. Ascertain input compatibility
% 3. Rank sub-systems
%
% Author:   Carolina H. Chung
% Contact:  chechung@umich.edu
  
% I/O/U
%{
REQUIRED INPUTS: 
    1. gem_table:       GEM reactions associated with top CARAMeL features
       -> note:         output from 'caramel_features2gem' function

    2. gem:             GEM model 
    
OUTPUT: 
    1. gem_table:       table variable of sub-systems represented by 
                        GEM reactions associated with top CARAMeL features
       -> subSystem:    cumulative varaiance explained
       -> Count:        number of reactions belonging to sub-system
       -> Ratio:        ratio of Count to total number of sub-system 
                        reactions in GEM
       -> Pvalue:       p-value of hypergeometric test for sub-system

EXAMPLE USAGE: 
    1. Return list of sub-systems represented by top GEM reactions:  
       >> top_subsystems = caramel_rankSubsystems(gem_table, gem); 
%}

%% PARSE THROUGH INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create an inputParser object
    p = inputParser; 
    
    % Define custom input validation function(s)
    isValidGEMtable = @(x) istable(x) && ismember({'subSystem'}, ...
        x.Properties.VariableNames); 
    isValidGEM = @(x) isstruct(x) && isfield(x, 'subSystems');  
    
    % Define required input variables
    addRequired(p, 'gem_table', isValidGEMtable)
    addRequired(p, 'gem', isValidGEM)
    
    % Define optional input validation functions
    g = {'Sy', 'Se', 'Sy/Se', 'Sy/R', 'A', 'R', 'A/R', 'A/Se', 'all'}; 
    isValidGroup = @(x) any(validatestring(x, g)); 
    isValidBoolean = @(x) islogical(x) || ...
        (isnumeric(x) && (x == 0 || x == 1));
    
    % Define optional input parameters
    addParameter(p, 'Group', 'all', isValidGroup)
    addParameter(p, 'AdjustPval', true, isValidBoolean)
    
    % Parse through inputs
    parse(p, gem_table, gem, varargin{:})
    
    % Define inputs by name
    t       = p.Results.gem_table; 
    gem     = p.Results.gem; 
    group   = p.Results.Group; 
    padj    = p.Results.AdjustPval; 
    
%% ASCERTAIN INPUT COMPATIBILITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check that GEM table sub-systems match GEM sub-systems
    if ~any(ismember(t.subSystem, gem.subSystems))
        error('Sub-systems in GEM table and GEM do not match.')
    end
    
    % Filter out non-significant portion
    t = t(t.Pvalue < 0.05, :); 
    
    % Define data to work with based on group
    if ~strcmpi(group, 'all')
        t = t(strcmpi(t.Direction, group), :); 
    end
    
%% RANK SUB-SYSTEMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Return unique number of sub-systems
    subSystem = unique(t.subSystem); 
    
    % Determine Count and Ratio
    Count = nan(size(subSystem)); Ratio = Count; Pvalue = Count; 
    T = numel(gem.subSystems); 
    for i = 1:numel(subSystem)
        Count(i) = sum(strcmp(t.subSystem, subSystem{i})); 
        N = sum(strcmp(gem.subSystems, subSystem{i})); 
        Ratio(i) = round(Count(i) / N, 2); 
        Pvalue(i) = 1 - hygecdf(Count(i), T, N, size(t, 1)); 
    end

    % Adjust p-values (if prompted)
    if padj
        Pvalue = mafdr(Pvalue, 'BHFDR', true); 
    end
    
    % Define output
    top_subsystems = table(subSystem, Count, Ratio, Pvalue); 
    top_subsystems = sortrows(top_subsystems, ...
        {'Pvalue','Ratio'}, {'ascend','descend'});
    
end