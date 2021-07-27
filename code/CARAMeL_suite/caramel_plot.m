function plots = caramel_plot(result_structure, plot_type, varargin)

% DESCRIPTION: 
% This function generates plots that assess CARAMeL model performance. 
% 
% STEPS: 
% 1. Parse through inputs
% 2. Ascertain input compatibility
% 3. Generate histogram 
% 4. Generate scatter plot
% 5. Generate confusion chart
% 6. Generate box plots
% 7. Generate ROC curves
% 8. Define output
%
% Author:   Carolina H. Chung
% Contact:  chechung@umich.edu
  
% I/O/U
%{
REQUIRED INPUTS: 
    1. result_structure:    results from CARAMeL model
    2. plot_type:           array specifying which plot types to generate
       -> type:             cell array with at least one element
       -> accepted values:  {'hist','scatter','confusion','box','roc'}
    
OPTIONAL INPUTS: 
    1. ClassThresh:         threshold values to classify drug interactions
       a. 1st element represents threshold for synergistic interactions
       b. 2nd element represents threshold for antagonistic interactions
       -> type:             numerical vector with two elements
       -> default:          ClassThresh = []

    2. Figure:              boolean to create a new figure with each plot
       -> type:             logical
       -> default:          Figure = true

    3. Subplot:             boolean to generate figures as subplots
       -> type:             logical 
       -> default:          Subplot = false

    4. Layout:              vector specifying subplot layout to use
       -> type:             numeric vector with two elements
       -> default:          Layout = []

    5. CorrType:            method to calculate correlation b/t scores
       -> type:             string or char
       -> note:             refer to built-in `corr` function for options
       -> default:          CorrType = 'spearman'
    
OUTPUTS: 
    1. plots:   structure containing plots stored as variables for further
                desired modifications for data visualization
       a. histPlot:         variable storing histogram of predicted scores
       b. scatterPlot:      variable storing scatter plot of actual
                            drug-pair interaction scores vs. predicted 
                            drug-pair interaction scores 
       c. confusionChart:   variable storing confusion chart 
       d. boxPlot:          variable storing box plots of predicted score
                            distribution according to true drug-pair
                            classification
       e. ROCplot:          variable storing ROC curves that assess the 
                            model's performance in correctly classifying
                            synergistic and antagonistic interactions

EXAMPLE USAGE: 
    1. Display scatter, confusion, and ROC plots: 
       >> caramel_plot(result_structure, {'scatter', 'confusion', 'roc'})

    2. Display confusion and ROC plots in single figure with 1x2 layout: 
       >> caramel_plot(result_structure, {'confusion', 'roc'}, ...
                      'Subplot', true, 'Layout', [1, 2])
%}

%% PARSE THROUGH INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Create an inputParser object
    p = inputParser; 
    if rem(length(varargin), 2) ~= 0    % for additional parameters
        error('Imbalanced number of name-value pair arguments provided.')
    end
    
    % Define custom input validation function(s)
    isValidResult = @(x) isstruct(x) && isfield(x, 'interactionNames'); 
    isValidType = @(x) all(ismember(x, ...
        {'hist','scatter','confusion','box','roc'})); 
    isValidClassThresh = @(x) isnumeric(x) && numel(x) == 2; 
    isValidBoolean = @(x) islogical(x) || ...
        (isnumeric(x) && (x == 1 || x == 0)); 
    isValidLayout = @(x) isnumeric(x) && numel(x) == 2 && ~any(rem(x, 1));
    isValidPosition = @(x) isnumeric(x) && numel(x) == 2 && all(x <= 1); 
    isValidIndices = @(x) islogical(x) || isnumeric(x); 
    
    % Define required input variables
    addRequired(p, 'result_structure', isValidResult)
    addRequired(p, 'plot_type', isValidType)
    
    % Define optional input parameters
    addParameter(p, 'ClassThresh', [], isValidClassThresh)
    addParameter(p, 'Figure', true, isValidBoolean)
    addParameter(p, 'Subplot', false, isValidBoolean)
    addParameter(p, 'Layout', [], isValidLayout)
    addParameter(p, 'Text', true, isValidBoolean)
    addParameter(p, 'Position', [0.8 0.1], isValidPosition)
    addParameter(p, 'Indices', [], isValidIndices)
    
    % Parse through inputs
    parse(p, result_structure, plot_type, varargin{:})
    
    % Define inputs by name
    res     = p.Results.result_structure; 
    type    = lower(p.Results.plot_type); 
    thresh  = p.Results.ClassThresh; 
    fig     = logical(p.Results.Figure); 
    sub     = logical(p.Results.Subplot); 
    layout  = p.Results.Layout; 
    txt     = logical(p.Results.Text); 
    pos     = p.Results.Position; 
    idx     = p.Results.Indices;
    
%% ASCERTAIN INPUT COMPATIBILITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Make sure inputs for each plot type are provided
    %   histogram
    if ismember('hist', type)
        if ~isfield(res, 'predScores')
            error('Provide predicted scores to generate histogram.')
        end
    end
    if ismember('scatter', type)
        if ~isfield(res, 'interactionScores') || ~isfield(res, 'predScores')
            error(strcat("Provide both actual and predicted scores ", ...
            "to generate scatter plot.")) 
        end
    end
    %   confusion plots
    if ismember('confusion', type)
        if ~isfield(res, 'MLmethod')
            error('Provide ML method.')
        elseif strcmpi(res.MLmethod, 'classification')
            if ~all(isfield(res, {'interactionClass','predClass'}))
                error('Provide both actual and predicted classes.')
            end
        elseif strcmpi(res.MLmethod, 'regression')
            if ~any(isfield(res, {'interactionClass','predClass'}))
                warning('Attempting to classify scores...')
                if isempty(thresh)
                    error('Provide threshold data for classification.')
                else
                    if ~isfield(res, 'interactionClass') && ...
                            ~isfield(res, 'interactionScores')
                        error('Provide actual interaction scores.')
                    elseif ~isfield(res, 'interactionClass') && ...
                            isfield(res, 'interactionScores')
                        res.interactionClass = ...
                            caramel_classify(res.interactionScores, thresh);
                    end
                    if ~isfield(res, 'predClass') && ...
                            ~isfield(res, 'predScores')
                        error('Provide predicted interaction scores.')
                    elseif ~isfield(res, 'predClass') && ...
                            isfield(res, 'predScores')
                        res.predClass = ...
                            caramel_classify(res.predScores, thresh);
                    end
                end
            end
        end
    end
    %   box plots
    if ismember('box', type)
        if ~isfield(res, 'predScores')
            error('Provide predicted scores to generate box plots.')
        end
        if ~isfield(res, 'interactionClass') 
            if isempty(thresh) || ~isfield(res, 'interactionScores')
                error(strcat("Insufficient information provided to ", ...
                    "classify drug interactions."))
            else
                res.interactionClass = ...
                    caramel_classify(res.interactionScores, thresh); 
            end
        end
    end
    %   ROC plots
    if ismember('roc', type)
        if ~isfield(res, 'predScores')
            error('Provide predicted scores.')
        end
        if ~isfield(res, 'MLmethod')
            error('Provide ML method.')
        elseif strcmpi(res.MLmethod, 'classification')
            if ~isfield(res, 'interactionClass')
                error('Provide both actual interaction classes.')
            end
        elseif strcmpi(res.MLmethod, 'regression')
            if ~isfield(res, 'interactionClass')
                warning('Attempting to classify actual scores...')
                if isempty(thresh)
                    error('Provide threshold data for classification.')
                elseif ~isfield(res, 'interactionScores')
                    error('Provide actual interaction scores.')
                else
                    res.interactionClass = ...
                        caramel_classify(res.interactionScores, thresh);
                end
            end
        end
    end
    
    % Define indexing variable (if sub = true)
    if sub
        fig = false; figure; k = 0; 
        % Make sure subplot layout is valid
        if prod(layout) < numel(type)
            error('Invalid subplot layout provided.')
        end
    end
    
    % Define data
    if all(idx == 0 | idx == 1)
        idx = logical(idx); 
    end
    if isempty(idx)
        idx = true(size(res.interactionScores)); 
    end
    names = res.interactionNames(idx, :); 
    scores = res.interactionScores(idx); 
    class = res.interactionClass(idx); 
    predScores = res.predScores(idx); 
    predClass = res.predClass(idx); 
    [r, p] = corr(scores, predScores, 'type', 'spearman'); 
    
%% GENERATE HISTOGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Histogram of predicted scores + distribution fit based on kernel
    if ismember('hist', type)
        if sub
            k = k + 1; subplot(layout(1), layout(2), k)
        elseif fig
            figure
        end
        % Generate histogram
        histPlot = histfit(predScores, 10, 'kernel'); 
        xlabel('Prediction'); ylabel('Count')
        title('Distribution of predicted interaction scores')
    else
        histPlot = [];
    end
    
%% GENERATE SCATTER PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Scatter plot of actual vs. predicted scores
    if ismember('scatter', type)
        if sub
            k = k + 1; subplot(layout(1), layout(2), k)
        elseif fig
            figure
        end
        % Generate scatter plot
        m = struct('Interaction', {join(names, ' + ')}); 
        scatterPlot = scatter2(scores, predScores, 'Metadata', m, ...
            'ShowAxes', false, 'MarkerFaceAlpha', 0.5, ...
            'XLabel', 'Actual', 'YLabel', 'Predicted'); 
        title('Correlation Assessment')
        % Add fit line
        fit_line = lsline; fit_line.Color = 'r'; fit_line.LineWidth = 2; 
        % Add correlation text (if prompted)
        if txt
            text(pos(1), pos(2), sprintf('R = %4.2f \n p = %4.2e', r, p), ...
                'Units', 'normalized')
        end
    else
        scatterPlot = [];
    end
    
%% GENERATE CONFUSION CHART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Confusion chart of interaction classes
    if ismember('confusion', type)
        if sub
            k = k + 1; subplot(layout(1), layout(2), k)
        elseif fig
            figure
        end
        [confMat, order] = confusionmat(class, predClass);
        % Define class order
        if all(ismember({'Synergy', 'Additivity', 'Antagonism'}, order))
            order = ordinal(order); 
            order = reorderlevels(order, ...
                {'Synergy', 'Additivity', 'Antagonism'});
        elseif all(ismember({'Synergy', 'Additivity'}, order))
            order = ordinal(order); 
            order = reorderlevels(order, {'Synergy','Additivity'});
        elseif all(ismember({'Additivity', 'Antagonism'}, order))
            order = ordinal(order); 
            order = reorderlevels(order, {'Additivity', 'Antagonism'});
        end
        % Generate confusion chart
        confusionChart = confusionchart(confMat, order); 
        confusionChart.Title = 'Drug Pair Classification';
    else
        confusionChart = [];
    end
    
%% GENERATE BOX PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Box plots of predicted score distribution for each class
    if ismember('box', type)
        if sub
            k = k + 1; subplot(layout(1), layout(2), k)
        elseif fig
            figure
        end
        % Define class order
        order = unique(class); 
        if all(ismember({'Synergy', 'Additivity', 'Antagonism'}, order))
            order = {'Synergy', 'Additivity', 'Antagonism'};
        elseif all(ismember({'Synergy', 'Additivity'}, order))
            order = {'Synergy','Additivity'};
        elseif all(ismember({'Additivity', 'Antagonism'}, order))
            order = {'Additivity', 'Antagonism'};
        end
        % Generate box plots
        m = struct('Interaction', {join(names, '+')}); 
        boxPlot = sig_boxplot(predScores, class, ...
            'GroupOrder', order, 'Metadata', m);
        xlabel('True Class'); ylabel('Predicted Scores')
        title('Predicted Score Distribution')
    else
        boxPlot = [];
    end
    
%% GENERATE ROC CURVES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ROC curves to determine how well model classified interactions
    if ismember('roc', type)
        if sub
            k = k + 1; subplot(layout(1), layout(2), k)
        elseif fig
            figure
        end
        % Define inputs for ROC curves
        if strcmpi(res.MLmethod, 'regression')
            syn_scores = rescale(-predScores); 
            ant_scores = rescale(predScores); 
        elseif strcmpi(res.MLmethod, 'classification')
            scores = res.predScores; 
            syn_scores = scores(:,1) - max(scores(:,2), scores(:,3)); 
            ant_scores = scores(:,3) - max(scores(:,1), scores(:,2)); 
        end
        [x1, y1, ~, auc1] = perfcurve(class, syn_scores, 'Synergy');
        [x2, y2, ~, auc2] = perfcurve(class, ant_scores, 'Antagonism');
        % Generate ROC curves
        ROCplot = plot(x1, y1, x2, y2, 'LineWidth', 2);
        ref = refline(1, 0); ref.Color = 'k'; ref.LineStyle = '--'; 
        legend(sprintf('S, AUROC = %4.2f', auc1), ...
            sprintf('A, AUROC = %4.2f', auc2), ...
            'Location', 'southeast')
        xlabel('1 - Specificity'); ylabel('Sensitivity')
        title('Model Performance Assessment')
    else
        ROCplot = [];
    end
    
%% DEFINE OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clear plots
    plots = struct(...
        'histPlot', histPlot, ...
        'scatterPlot', scatterPlot, ...
        'confusionChart', confusionChart, ...
        'boxPlot', boxPlot, ...
        'ROCplot', ROCplot); 
    
end