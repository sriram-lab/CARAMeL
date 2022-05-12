function boxPlot = sig_boxplot(x, g, varargin)

% function boxPlot = sig_boxplot(data, groups, group_order, plot_points, point_label)

% Function description: 
% This function generates boxplots that include results for significance
% comparison.
%
% Function steps: 
% 1. Parse through inputs
% 2. Ascertain input compatibility
% 3. Generate grouped box plots

% I/O
%{
REQUIRED INPUTS: 
    1. x: vector of numerical data to plot

    2. g: vector of group labels for data, specified as: 
          a. Categorical, logical, or numeric vector
             e.g. [1 1 1 2 1 3 3 2 2]
          b. Character, string, or cell array of character vectors
             e.g. {'A' 'A' 'A' 'B' 'A' 'C' 'C' 'B' 'B'}

    -> Note: x, y, and g must all have the same length

OPTIONAL INPUTS: 
    1. GroupOrder: cell of unique group labels set in desired order
       a. Default: GroupOrder = unique(g)

    2. PlotPoints: Boolean to specify if points should be plotted
       a. Default: PlotPoints = true

    3. Metadata:   include other information for plot tooltip; specify as:
       a. Structure variable, where field names correspond to the category
          name and the field value contains category information for all 
          data points (i.e. field value length = number of data points)
          e.g. structure('Group Type', group_type)
       b. Default: Metadata = []

OUTPUTS: 
    1. boxPlot: object handle for output boxplots
%}

%% PARSE THROUGH INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create an inputParser object
    p = inputParser; 
    if rem(length(varargin),2) ~= 0
        error('Imbalanced number of name-value pair arguments provided.')
    end
    
    % Define required input validation function
    isValidGroup = @(x) iscategorical(x) || isnumeric(x) || iscell(x)||...
        ischar(x) || islogical(x);
    
    % Define required input variables
    addRequired(p, 'x', @isnumeric)
    addRequired(p, 'g', isValidGroup) 
    
    % Define optional input validation function
    isValidColor = @(x) ischar(x) || isnumeric(x) || iscell(x); 
    isValidBool = @(x) islogical(x) || isnumeric(x); 
    isValidLabel = @(x) ischar(x) || isstring(x) || iscell(x);
    
    % Define optional input parameters
    addParameter(p, 'GroupOrder', unique(g, 'stable'), isValidGroup)
    addParameter(p, 'PlotPoints', false, isValidBool)
    addParameter(p, 'MarkerSize', 36, @isnumeric)
    addParameter(p, 'MarkerFaceColor', 'k', isValidColor)
    addParameter(p, 'MarkerEdgeColor', 'k', isValidColor)
    addParameter(p, 'MarkerFaceAlpha', 0.2, @isnumeric)
    addParameter(p, 'XLabel', 'X', isValidLabel)
    addParameter(p, 'YLabel', 'Y', isValidLabel)
    addParameter(p, 'Metadata', [], @isstruct)
    
    % Parse through inputs
    parse(p, x, g, varargin{:})
    
    % Define inputs by name
    x        = p.Results.x; 
    g        = p.Results.g; 
    gorder   = p.Results.GroupOrder; 
    ppoints  = logical(p.Results.PlotPoints); 
    mkSz    = p.Results.MarkerSize; 
    mkFC    = p.Results.MarkerFaceColor; 
    mkEC    = p.Results.MarkerEdgeColor; 
    mkFA    = p.Results.MarkerFaceAlpha; 
    xlab     = p.Results.XLabel; 
    ylab     = p.Results.YLabel; 
    md       = p.Results.Metadata;

%% ASCERTAIN INPUT COMPATIBILITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Don't plot anything if x is empty
    if isempty(x) 
        warning('An empty input vector was provided. Exiting function.')
        if nargout > 0
            boxPlot = []; 
        end
        return
    end
    
    % Check x, g, and gorder
    if numel(x) ~= numel(g)
        error('Size of group variable does not match input data.');
    end
    
    if numel(mkFC) == 1
        mkFC = repmat(mkFC, numel(g), 1); 
        mkFA = repmat(mkFA, numel(g), 1); 
    end

    if isnumeric(gorder)
        gorder = cellstr(num2str(gorder)); 
    end
    
%     % Define mkFC and mkFA (if not defined)
%     if ~isempty(scatter2Varargin) && isfield(scatter2Varargin, 'Metadata')
%         idx = find(strcmpi('Metadata', scatter2Varargin)) + 1; 
%         s = scatter2Varargin{idx}; s.Group = g; 
%     else
%         scatter2Varargin = cat(2, scatter2Varargin, {'Metadata', ...
%             struct('Group', {g})}); 
%     end

%     % Check metadata
%     if ~isempty(md)
%         mdFields = fieldnames(md); 
%         for i = 1:length(mdFields)
%             if eval(sprintf('length(md.%1$s) ~= numel(x)', mdFields{i}))
%                 error([mdFields{i} ' field of metadata does not match ' ...
%                     'input data.'])
%             end
%         end
%     end
    
%     % Convert boolean parameters to logical variables
%     if isnumeric(ppoints)
%         ppoints = logical(ppoints); 
%     end

%% GENERATE GROUPED BOX PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Calculate pairwise p-values
    if iscell(g) || ischar(g) || iscategorical(g)
        pair_groups = nchoosek(cellstr(unique(g)), 2);
    elseif isnumeric(g) || islogical(g)
        pair_groups = nchoosek(unique(g), 2); 
    end
    p_vals = nan(size(pair_groups, 1), 1);
    for i = 1:size(pair_groups, 1)
        try
            vals1 = x(g == pair_groups(i, 1));
            vals2 = x(g == pair_groups(i, 2));
        catch
            vals1 = x(strcmp(g, pair_groups{i,1}));
            vals2 = x(strcmp(g, pair_groups{i,2}));
        end
        [~, p_vals(i)] = ttest2(vals1, vals2, 'VarType', 'unequal');
        if p_vals(i) >= 0.05
            p_vals(i) = NaN; 
        end
    end
    idx = ~isnan(p_vals);
    
    % Generate box plots
    boxPlot = boxplot(x, g, 'GroupOrder', gorder);
    set(boxPlot, 'LineWidth' , 2)
    
    % Add scattered points (if prompted)
    if ppoints
        hold on
        xCenter = 1:length(gorder);
        spread = 0.5; % 0 = no spread; 0.5 = random spread within bounds
        for i = 1:length(gorder)
            try
                ix = g == gorder{i};
            catch
                ix = strcmp(g, gorder{i});
            end
            points = x(ix);
            xcoord = rand(size(points)) * spread - (spread/2) + xCenter(i);
            c = char(unique(mkFC(ix), 'stable'))'; 
            if ~isempty(md)
                mdn = fieldnames(md); md_ix = struct();
                for j = 1:numel(mdn)
                    eval(sprintf('md_ix.%s = md.%s(ix);', mdn{j}, mdn{j}))
                end
            end
            s = gscatter2(xcoord, points, mkFC(ix), 'MarkerSize', mkSz, ...
                'MarkerFaceColor', c, 'MarkerEdgeColor', mkEC, ...
                'MarkerFaceAlpha', unique(mkFA(ix)), 'Legend', false, ...
                'XLabel', xlab, 'YLabel', ylab, 'Metadata', md_ix); 
%             s = scatter(rand(size(points)) * spread - (spread/2) + ...
%                 xCenter(i), points, 'ko', 'MarkerFaceColor', 'k', ...
%                 'MarkerFaceAlpha', 0.2);
%             s.DataTipTemplate.DataTipRows(1).Label = xlab; 
%             s.DataTipTemplate.DataTipRows(2).Label = ylab; 
%             row = dataTipTextRow('Group', g(ix)); 
%             s.DataTipTemplate.DataTipRows(end + 1) = row; 
%             if ~isempty(md)
%                 mdFields = fieldnames(md); 
%                 for j = 1:length(mdFields)
%                     row = dataTipTextRow(mdFields{j}, ...
%                         eval(sprintf('md.%s(ix)', mdFields{j})));
%                     s.DataTipTemplate.DataTipRows(end + 1) = row;
%                 end
%             end
        end
        set(gca, 'LineWidth' ,1)
    end
    
    % Add significance annotation
    pairs = cell(size(pair_groups, 1), 1);
    for i = 1:size(pair_groups, 1)
        pairs{i} = pair_groups(i,:);
    end
    sigstar(pairs(idx), p_vals(idx)); 
    axis 'auto y'; hold off; xlabel(xlab); ylabel(ylab); 
    
end