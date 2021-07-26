function varargout = extract_metadata(file, query, match, extract, varargin)

    % DESCRIPTION 
    % This function extracts metadata for a genome-scale metabolic model
    % (GEM) given an annotation table. 
    % 
    % STEPS 
    % 1. Parse through inputs
    % 2. Ascertain input compatibility
    % 3. Extract metadata of interest
    % 
    % Author:   Carolina H. Chung
    % Created:  June 6, 2019
    % Updates:  September 20, 2020 (Carolina H. Chung)
    
    % I/O
    %{
    REQUIRED INPUTS: 
        1. file:        filename for reference data/table (.xls/.xlsx file)
           -> type:     string/char with .xls or .xlsx extension
    
        2. query:       list of unique queries to compare agains reference 
           -> type:     string/char cell array of query list
           -> dim:      m x 1 (m = number of query entries)
    
        3. match:       column name in reference to match queries with 
           -> type:     string/char
    
        4. extract:     column name(s) to extract data from based on query
           -> type:     single string/char or string/char array
           -> dim:      1 x n (n = number of columns to extract data from)
    
    OPTIONAL INPUTS:
        1. READTABLE:   any additional inputs for 'readtable' function can
                        be included as Name-Value pairs when reading 'file'
    
    OUTPUTS:
        1. varargout:   variable-length list of outputs 
           -> type:     cell or numeric vector(s)
           -> dim:      m x 1 x n (m = query length, 
                                   n = # of columns to extract from)
           -> note:     returns 'query' entries for those that do not have
                        a match in the reference table
    
    EXAMPLE USAGE: 
        1. Extract from one column labeled 'EntrezID': 
           >> output = extract_metadata(reference_table.xlsx, ...
                                        gene_list, 'gene', 'EntrezID')
    
        2. Extract from two columns labeled 'EntrezID' and 'KeggID': 
           >> output = extract_metadata(reference_table.xlsx, ...
                                        gene_list, 'gene', ...
                                        {'EntrezID', 'KeggID'})
    
        3. Specify additional parameter for 'readtable': 
           >> output = extract_metadata(reference_table.xlsx, ...
                                        gene_list, 'gene', ...
                                        {'EntrezID', 'KeggID'}, ...
                                        'Sheet', 'gene_data')
    %}
    
%% PARSE THROUGH INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Create an inputParser object
    p = inputParser; 
    if rem(length(varargin), 2) ~= 0    % for additional parameters
        error('Imbalanced number of name-value pair arguments provided.')
    end
    
    % Define custom input validation function(s)
    isValidName = @(x) isstring(x) || ischar(x); 
    isValidCell = @(x) iscell(x) && length(x) >= 1 && ...
        (isstring(x{1}) || ischar(x{1})); 
    isValidList = @(x) isValidName(x) || isValidCell(x); 
    
    % Define required input variables
    addRequired(p, 'file', isValidName)
    addRequired(p, 'query', isValidCell)
    addRequired(p, 'match', isValidName)
    addRequired(p, 'extract', isValidList)
    
    % Parse through inputs
    p.KeepUnmatched = true; 
    parse(p, file, query, match, extract, varargin{:})
    
    % Define inputs by name
    file =      p.Results.file; 
    query =     p.Results.query; 
    match =     p.Results.match; 
    extract =   p.Results.extract;  
    
    % Define list of optional input parameters from external functions
    readtableParams = {'FileType','ReadVariableNames','ReadRowNames', ...
        'TreatAsEmpty','TextType','DatetimeType',...
        'PreserveVariableNames','Delimiter','HeaderLines','Format',...
        'EmptyValue','MultipleDelimsAsOne','CollectOutput',...
        'CommentStyle','ExpChars','EndOfLine','DateLocale','Encoding',...
        'DurationType','HexType','BinaryType','Sheet','Range','UseExcel'}; 
    
    % Find any optional parameters and allocate to appropriate function
    if ~isempty(fieldnames(p.Unmatched))
        disp('Extra name-value pair arguments provided:')
        disp(p.Unmatched)
        extraParams = p.Unmatched; 
        extraNames = fieldnames(extraParams); 
        % READTABLE function
        [~, p_idx] = intersect(extraNames, readtableParams);
        if numel(p_idx) > 0
            readtableVarargin = cell(1, 2*numel(p_idx));
            readtableVarargin(1:2:end-1) = extraNames(p_idx);
            i = 0; 
            for j = 2:2:length(readtableVarargin)
                i = i + 1; 
                readtableVarargin(j) = ...
                    {eval(sprintf('extraParams.%1$s', extraNames{i}))};
            end
            disp(['Extra name-value pair argument(s) passed to ', ...
                'READTABLE function.'])
        else
            readtableVarargin = {};
        end
    end
    
%% ASCERTAIN INPUT COMPATIBILITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Define reference table
    reference = readtable(file, 'PreserveVariableNames', true, ...
        readtableVarargin{:}); 
    column_names = reference.Properties.VariableNames; 

    % Make sure that reference table contains match and extract inputs
    if ~ismember(match, column_names)
        error("Input for 'match' not found in reference table.")
    end
    if isstring(extract) || ischar(extract)
        if ~ismember(extract, column_names)
            error("Not all entries in 'extract' found in reference table.")
        end
    else
        if ~all(ismember(extract, column_names))
            error("Not all entries in 'extract' found in reference table.")
        end
    end
    
    % Make sure all entries in 'query' are unique
    if numel(unique(query)) ~= length(query)
        error("Entries in 'query' not all unique.")
    end
 
%% EXTRACT METADATA OF INTEREST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Define output object
    varargout = cell(nargout, 1); 

    % Define match and extract column indices
    match_idx = strcmpi(fieldnames(reference), match); 
    [~, ~, ext_idx] = intersect(extract, fieldnames(reference), 'stable'); 
    
    % Find query entries present in reference
    [~, query_idx, ref_idx] = ...
        intersect(query, table2array(reference(:, match_idx)), 'stable'); 
    miss_idx = find(~ismember(1:1:numel(query), query_idx)); 
    
    % Extract metadata
    for i = 1:nargout
        if ~isempty(miss_idx)
            varargout{i} = cell(numel(query), 1); 
            varargout{i}(query_idx) = ...
                table2array(reference(ref_idx, ext_idx(i))); 
            varargout{i}(miss_idx) = query(miss_idx); 
        else
            varargout{i} = ...
                table2array(reference(ref_idx, ext_idx(i))); 
        end
    end

end
