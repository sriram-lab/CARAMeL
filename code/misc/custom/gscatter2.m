function h = gscatter2(x, y, g, varargin)

% DESCRIPTION: 
% This function creates a grouped 2D scatter plot given X, Y, and grouping
% data. Optional inputs are available for plot customization as well 
% (refer to I/O section below). 
%
% STEPS:
% 1. Parse through inputs
% 2. Ascertain input compatibility
% 3. Generate grouped scatter plot

% I/O
%{
REQUIRED INPUTS: 
    1. x: x values, specified as a vector

    2. y: y values, specified as a vector

    3. g: grouping variable, specified as: 
       a. Categorical, logical, or numeric vector
          e.g. [1 1 1 2 1 3 3 2 2]
       b. Character, string, or cell array of character vectors
          e.g. {'A' 'A' 'A' 'B' 'A' 'C' 'C' 'B' 'B'}

    -> Note: x, y, and g must all have the same length

OPTIONAL INPUTS: 

    1. Colormap: colormap to use for picking random group colors
       a. Refer to 'colormap' MATLAB documentation for options
       b. Default: Colormap = 'hsv'

    2. MarkerSymbol: marker symbol, specified as: 
       a. Char array of symbols
          - Specify one character to use the same symbol for all markers
            e.g. 'o'
          - Specify a char array of same length as unique values in g to 
            use different edge colors for each group 
            e.g. 'o+*'
       b. Refer to 'MarkerType' input description in 'scatter' function for
          all available options
       c. Default: MarkerSymbol = 'o'

    3. MarkerSize: marker size, specified as: 
       a. Num array of size values
          - Specify one size value to use the same size for all markers
            e.g. 12
          - Specify a num array of same length as unique values in g to 
            use different marker sizes for each group 
            e.g. [10 11 12]
       b. Default: MarkerSize = 36

    4. MarkerFaceColor: marker face color, specified as: 
       a. Char array of color names (same length as unique values in g)
          e.g. 'rbg'
       b. Three column matrix of RGB triplets (number of rows = number of
          unique values in g)
          e.g. [1 0 0; 0 1 0; 0 0 1]
       c. Default: MarkerFaceColor = []

    5. MarkerEdgeColor: marker edge color, specified as: 
       a. Char array of color names
          - Specify one character to use the same color for all markers
            e.g. 'k'
          - Specify a char array of same length as unique values in g to 
            use different edge colors for each group 
            e.g. 'rbg'
       b. Three column matrix of RGB triplets 
          - Specify single vector to use the same color for all markers
            e.g. [1 0 0]; 
          - Specity an m x 3 matrix, where m = number of unique values in g
            e.g. [1 0 0; 0 1 0; 0 0 1]
       c. Default: MarkerEdgeColor = 'none'

    6. MarkerFaceAlpha: affects color transparency; specify as: 
       a. Num array of alpha values (between 0 and 1)
          - Specify one alpha value to use the same alpha for all markers
            e.g. 0.5
          - Specify a num array of same length as unique values in g to 
            use different alpha values for each group 
            e.g. [0.2 0.5 1]
       b. Default: MarkerFaceAlpha = 1

    7. LineWidth: line width for marker edges; specify as: 
       a. Num array of width values 
          - Specify one width value to use the same edge width for all 
            markers
            e.g. 2
          - Specify a num array of same length as unique values in g to 
            use different edge width values for each group 
            e.g. [2 4 6]
       b. Default: LineWidth = 0.5
    8. Filled: option to fill the interior of the markers
       a. Use this option with markers that have a face (e.g. 'o')
       b. Default = Filled = 'filled'

    9. Legend: specify whether to include a plot legend or not
       a. Use a boolean (false/true) or numeric (0/1) value
       b. Default: Legend = 'true'

    10. LegendPosition: specify legend position
        a. Refer to 'legend' MATLAB documentation for options
        b. Default: LegendPosition = 'best'

    11. LegendTitle: specify a title for the legend
        a. Can be specified as a cell, char, or string
           e.g. 'Legend Title'
        b. Default: LegendTitle = []

    12. XLabel: specify a label for the x-axis; also used in plot tooltip
        a. Can be specified as a cell, char, or string
           e.g. 'X'
        b. Default: XLabel = 'X'

    13. YLabel: specify a label for the y-axis; also used in plot tooltip
        a. Can be specified as a cell, char, or string
           e.g. 'Y'
        b. Default: XLabel = 'Y'

    14. Metadata: include other information for plot tooltip; specify as:
        a. Structure variable, where field names correspond to the category
           name and the field value contains category information for all 
           data points (i.e. field value length = number of data points)
           e.g. structure('Group Type', group_type)
        b. Default: Metadata = []
    
OUTPUTS: 
    1. h: graphics handles, returned as an array of Scatter objects. Each 
       Scatter object corresponds to a group in g
%}

%% PARSE THROUGH INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create an inputParser object
    p = inputParser; 
    if rem(length(varargin),2) ~= 0
        error('Imbalanced number of name-value pair arguments provided.')
    end
    
    % Define required input validation function
    isValidGroup = @(x) iscategorical(x) || isnumeric(x) || iscell(x)||...
        ischar(x) || isstring(x) || islogical(x);
    
    % Define required input variables
    addRequired(p, 'x', @isnumeric)
    addRequired(p, 'y', @isnumeric)
    addRequired(p, 'g', isValidGroup) 
    
    % Define optional input validation function
    isValidColor = @(x) ischar(x) || isnumeric(x) || iscell(x); 
    isValidLabel = @(x) ischar(x) || isstring(x) || iscell(x);
    isValidBoolean = @(x) islogical(x) || (isnumeric(x) && ...
        (x == 0 || x == 1)); 
    
    % Define optional input parameters
    addParameter(p, 'Colormap', 'hsv', isValidLabel)
    addParameter(p, 'MarkerSymbol', 'o', @ischar)
    addParameter(p, 'MarkerSize', 36, @isnumeric)
    addParameter(p, 'MarkerFaceColor', [], isValidColor)
    addParameter(p, 'MarkerEdgeColor', 'none', isValidColor)
    addParameter(p, 'MarkerFaceAlpha', 1, @isnumeric)
    addParameter(p, 'LineWidth', 0.5, @isnumeric)
    addParameter(p, 'Filled', 'filled', isValidLabel); 
    addParameter(p, 'ShowAxes', false, isValidBoolean)
    addParameter(p, 'Legend', true, @(x) islogical(x) || isnumeric(x))
    addParameter(p, 'LegendPosition', 'best', isValidLabel)
    addParameter(p, 'LegendNumColumns', 1, @isnumeric)
    addParameter(p, 'LegendTitle', [], isValidLabel)
    addParameter(p, 'XLabel', 'X', isValidLabel)
    addParameter(p, 'YLabel', 'Y', isValidLabel)
    addParameter(p, 'LabelAxes', true, @(x) islogical(x) || isnumeric(x))
    addParameter(p, 'Metadata', [], @isstruct)
    
    % Parse through inputs
    parse(p, x, y, g, varargin{:})
    
    % Define inputs by name
    x = p.Results.x; 
    y = p.Results.y; 
    g = p.Results.g; 
    cm = p.Results.Colormap; 
    mkSym   = p.Results.MarkerSymbol; 
    mkSz    = p.Results.MarkerSize; 
    mkFC    = p.Results.MarkerFaceColor; 
    mkEC    = p.Results.MarkerEdgeColor; 
    mkFA    = p.Results.MarkerFaceAlpha; 
    lw      = p.Results.LineWidth; 
    fill    = p.Results.Filled; 
    axes    = p.Results.ShowAxes; 
    leg     = p.Results.Legend; 
    legPos  = p.Results.LegendPosition; 
    legTit  = p.Results.LegendTitle; 
    legCol  = p.Results.LegendNumColumns; 
    xlab    = p.Results.XLabel; 
    ylab    = p.Results.YLabel; 
    label   = p.Results.LabelAxes; 
    md      = p.Results.Metadata; 
    
%% ASCERTAIN INPUT COMPATIBILITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Don't plot anything if either x or y is empty
    if isempty(x) || isempty(y)
        warning('An empty input vector was provided. Exiting function.')
        if nargout > 0
            h = []; 
        end
        return
    end
    
    % Check x and y
    if numel(x) ~= numel(y)
        error('Sizes for x and y do not match.');
    end

    % Check g (group variable)
    if length(g) ~= numel(x)
        error('Size of group variable does not match x and y.');
    end
    
    % Check metadata
    if ~isempty(md)
        mdFields = fieldnames(md); 
        for i = 1:length(mdFields)
            if eval(sprintf('length(md.%1$s) ~= numel(x)', mdFields{i}))
                error([mdFields{i} ' field of metadata does not match ' ...
                    'x nor y data.'])
            end
        end
    end
    
    % Convert boolean parameters to logical variables
    if isnumeric(leg)
        leg = logical(leg); 
    end
    if isnumeric(label)
        label = logical(label); 
    end

%% GENERATE GROUPED SCATTER PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Convert group indexing
    [gIdx, gLabels] = grp2idx(g); 
    n = max(gIdx); 
    
    % Assign colors based on colormap (if custom colors not provided)
    if isempty(mkFC)
        mkFC = feval(cm, n);
    end
    
    % Repeat marker attributes according to max number of groups
    if numel(mkSym) == 1                            % marker symbol
        mkSym = repelem(mkSym, n, 1); 
    elseif numel(mkSym) ~= n
        error('Input for marker symbol is invalid.')
    end
    if numel(mkSz) == 1                             % marker size
        mkSz = repelem(mkSz, n, 1); 
    elseif numel(mkSz) ~= n
        error('Input for marker size is invalid.')
    end
    if (ischar(mkFC) && numel(mkFC) ~= n) || ...    % marker face color
            (isnumeric(mkFC) && size(mkFC,1) ~= n) 
        error('Input for marker face color is invalid.')
    elseif ischar(mkFC) && numel(mkFC) == n
        mkFC = mkFC'; 
    end
    if numel(mkEC) == 1                             % marker edge color
        mkEC = repelem(mkEC, n, 1); 
    elseif strcmpi(mkEC, 'none')
        mkEC = repmat(mkEC, n, 1); 
    elseif numel(mkEC) ~= n
        error('Input for marker edge color is invalid.')
    end
    if numel(mkFA) == 1                             % marker face alpha
        mkFA = repelem(mkFA, n, 1); 
    elseif numel(mkFA) ~= n
        error('Input for marker face alpha is invalid.')
    end
    if numel(lw) == 1                               % line width
        lw = repelem(lw, n, 1); 
    elseif numel(lw) ~= n
        error('Input for line width is invalid.')
    end
    
    % Generate scatter plot
    h = cell(n,1); 
    for i = 1:n
        if i ~= 0
            hold on
        end
        s = scatter(x(gIdx == i), y(gIdx == i), ...
            mkSz(i), mkFC(i,:),  mkSym(i), fill, ...
            'MarkerEdgeColor', mkEC(i,:), 'LineWidth', lw(i)); 
        s.MarkerFaceAlpha = mkFA(i); 
        s.DataTipTemplate.DataTipRows(1).Label = xlab; 
        s.DataTipTemplate.DataTipRows(2).Label = ylab; 
        row = dataTipTextRow('Group', repmat(gLabels(i),sum(gIdx == i),1)); 
        s.DataTipTemplate.DataTipRows(end+1) = row; 
        if ~isempty(md)
            for j = 1:length(mdFields)
                row = dataTipTextRow(mdFields{j}, ...
                    eval(sprintf('md.%1$s(gIdx == i)', mdFields{j})));
                s.DataTipTemplate.DataTipRows(end+1) = row;
            end
        end
        h{i} = s; 
    end
    hold off
    
    % Add axis lines (if prompted)
    if axes
        vline(0, 'k'); hline(0, 'k')
    end
    
    % Change axis labels (if specified)
    if label
        xlabel(xlab); ylabel(ylab); 
    end
        
    % Add legend (if prompted)
    if leg 
        lgd = legend(gLabels, 'Location', legPos, 'NumColumns', legCol); 
        if ~isempty(legTit)
            title(lgd, legTit)
        end
    end
    
end