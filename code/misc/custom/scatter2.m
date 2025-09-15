function s = scatter2(x, y, varargin)

% Function description: 
% This function creates a custom 2D scatter plot given X, and Y. 
% Optional inputs are available for plot customization as well 
% (refer to I/O section below). 

% I/O
%{
REQUIRED INPUTS: 
    1. x: x values, specified as a vector

    2. y: y values, specified as a vector

    -> Note: x and y must have the same length

OPTIONAL INPUTS: 
    1. MarkerSymbol: marker symbol, specified as: 
       a. Char array of symbols
          - Specify one character to use the same symbol for all markers
            e.g. 'o'
          - Specify a char array of same length as unique values in g to 
            use different edge colors for each group 
            e.g. 'o+*'
       b. Refer to 'MarkerType' input description in 'scatter' function for
          all available options
       c. Default: MarkerSymbol = 'o'

    2. MarkerSize: marker size, specified as: 
       a. Num array of size values
          - Specify one size value to use the same size for all markers
            e.g. 12
          - Specify a num array of same length as unique values in g to 
            use different marker sizes for each group 
            e.g. [10 11 12]
       b. Default: MarkerSize = 36

    3. MarkerFaceColor: marker face color, specified as: 
       a. Char array of color names (same length as unique values in g)
          e.g. 'rbg'
       b. Three column matrix of RGB triplets (number of rows = number of
          unique values in g)
          e.g. [1 0 0; 0 1 0; 0 0 1]
       c. Default: MarkerFaceColor = []

    4. MarkerEdgeColor: marker edge color, specified as: 
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

    5. MarkerFaceAlpha: affects color transparency; specify as: 
       a. Num array of alpha values (between 0 and 1)
          - Specify one alpha value to use the same alpha for all markers
            e.g. 0.5
          - Specify a num array of same length as unique values in g to 
            use different alpha values for each group 
            e.g. [0.2 0.5 1]
       b. Default: MarkerFaceAlpha = 1

    6. LineWidth: line width for marker edges; specify as: 
       a. Num array of width values 
          - Specify one width value to use the same edge width for all 
            markers
            e.g. 2
          - Specify a num array of same length as unique values in g to 
            use different edge width values for each group 
            e.g. [2 4 6]
       b. Default: LineWidth = 0.5

    7. Filled: option to fill the interior of the markers
       a. Use this option with markers that have a face (e.g. 'o')
       b. Default: Filled = 'filled'

    8. ShowAxes: boolean to show axis lines or not
       a. Default: ShowAxes = true

    9. XLabel: specify a label for the x-axis; also used in plot tooltip
       a. Can be specified as a cell, char, or string
          e.g. 'X'
       b. Default: XLabel = 'X'

    10. YLabel: specify a label for the y-axis; also used in plot tooltip
       a. Can be specified as a cell, char, or string
          e.g. 'Y'
       b. Default: XLabel = 'Y'

    11. Metadata: include other information for plot tooltip; specify as:
        a. Structure variable, where field names correspond to the category
           name and the field value contains category information for all 
           data points (i.e. field value length = number of data points)
           e.g. structure('Group Type', group_type)
        b. Default: Metadata = []
    
OUTPUTS: 
    1. s: graphics handle, returned as a Scatter object. 
%}

%% PARSE THROUGH INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create an inputParser object
    p = inputParser; 
    if rem(length(varargin),2) ~= 0
        error('Imbalanced number of name-value pair arguments provided.')
    end
    
    % Define required input variables
    addRequired(p, 'x', @isnumeric)
    addRequired(p, 'y', @isnumeric)
    
    % Define optional input validation function
    isValidColor = @(x) ischar(x) || isnumeric(x); 
    isValidAlpha = @(x) isnumeric(x) && x >= 0 && x <= 1; 
    isValidLabel = @(x) ischar(x) || isstring(x) || iscell(x);
    isValidBoolean = @(x) islogical(x) || (isnumeric(x) && ...
        (x == 0 || x == 1)); 
    
    % Define optional input parameters
    addParameter(p, 'MarkerSymbol', 'o', @ischar)
    addParameter(p, 'MarkerSize', 36, @isnumeric)
    addParameter(p, 'MarkerFaceColor', 'k', isValidColor)
    addParameter(p, 'MarkerEdgeColor', 'none', isValidColor)
    addParameter(p, 'MarkerFaceAlpha', 1, isValidAlpha)
    addParameter(p, 'LineWidth', 0.5, @isnumeric)
    addParameter(p, 'Filled', 'filled', isValidLabel); 
    addParameter(p, 'ShowAxes', false, isValidBoolean)
    addParameter(p, 'XLabel', 'X', isValidLabel)
    addParameter(p, 'YLabel', 'Y', isValidLabel)
    addParameter(p, 'Metadata', [], @isstruct)
    
    % Parse through inputs
    parse(p, x, y, varargin{:})
    
    % Define inputs by name
    x       = p.Results.x; 
    y       = p.Results.y; 
    mkSym   = p.Results.MarkerSymbol; 
    mkSz    = p.Results.MarkerSize; 
    mkFC    = p.Results.MarkerFaceColor; 
    mkEC    = p.Results.MarkerEdgeColor; 
    mkFA    = p.Results.MarkerFaceAlpha; 
    lw      = p.Results.LineWidth; 
    fill    = p.Results.Filled; 
    axes    = p.Results.ShowAxes; 
    xlab    = p.Results.XLabel; 
    ylab    = p.Results.YLabel; 
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

%% GENERATE GROUPED SCATTER PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Generate scatter plot
    s = scatter(x, y, mkSz, mkFC, mkSym, fill, ...
        'MarkerEdgeColor', mkEC, 'LineWidth', lw); 
    s.MarkerFaceAlpha = mkFA; 
    s.DataTipTemplate.DataTipRows(1).Label = xlab; 
    s.DataTipTemplate.DataTipRows(2).Label = ylab; 
    if ~isempty(md)
        for j = 1:length(mdFields)
            row = dataTipTextRow(mdFields{j}, ...
                eval(sprintf('md.%1$s', mdFields{j})));
            s.DataTipTemplate.DataTipRows(end+1) = row;
        end
    end
    
    % Add axis lines (if prompted)
    if axes
        vline(0, 'k'); hline(0, 'k')
    end
    
    % Add axis labels
    xlabel(xlab); ylabel(ylab)
    
end