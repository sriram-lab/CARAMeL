function [phenotype_data, phenotype_labels, conditions] = ...
    process_chemgen_v2(fname, z)

    % DESCRIPTION 
    % This function processes chemogenomic data to form a binary matrix
    % indicating conditions where KO strains were sensitive or resistant
    % to treatment. 
    % 
    % STEPS 
    % 1. Input processing
    % 2. Load data
    % 3. Convert gene IDs to standard IDs
    % 4. Transform data
    % 5. Define output phenotype variables
    % 
    % Author:   Sriram Chandrasekaran 
    % Created:  2018-10-23
    % Updated:  2021-05-21 (Carolina H. Chung)

    % I/O
    %{
    OPTIONAL INPUTS: 
        1. fname:   filename for chemogenomic data
        2. z:       threshold value to define significant effect on fitness
                    (default: z = 2)
    
    OUTPUTS:
        1. phenotype_data:      binary matrix containing 1 where E. coli
                                was sensitive or resistant to a treatment
        2. phenotype_labels:    labels (i.e. genes) for phenotype data
        3. conditions:          list of conditions (i.e. treatments)
    %}
    
%% INPUT PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(fname)
        fname = 'ecoli_phenotype_data_cell.xlsx'; % (Nichols et al.)
    end
    if ~exist('z','var') || isempty(z)
        z = 2;
    end
    
%% LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [phenotype_num, txt] = xlsread(fname);   % numerical data
    probelist = txt(2:end, 1);               % gene names (ECK numbers)
    conditions = txt(1, 2:end)';             % list of conditions

%% CONVERT GENE IDs TO STANDARD IDs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('ecoli_annotation_data','genenames_array','genenames_bnums')
    clear plist
    plist = cell(size(probelist));
    for i = 1:length(probelist)
        tem = regexp(probelist{i}, 'ECK[\d]*-', 'split');
        try
            plist(i) = tem(2);
        catch
            plist(i) = tem(1);
        end
    end
    plist = regexprep(plist, '''', '');
    [ix, pos] = ismember(upper(plist), upper(genenames_array));
    plist_bnums = plist; plist_bnums(ix) = genenames_bnums(pos(ix));
    plist_bnums = matlab.lang.makeUniqueStrings(plist_bnums);
    
%% TRANSFORM DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sensitive strains
    sensitive_phenotype_num = phenotype_num < -z; 
    idx = sum(sensitive_phenotype_num, 2) == 0;
    sensitive_data = sensitive_phenotype_num(~idx,:);
    sensitive_labels = plist_bnums(~idx);
    % Resistant strains
    resistant_phenotype_num = phenotype_num > z; 
    idx = sum(resistant_phenotype_num,2) == 0;
    resistant_data = resistant_phenotype_num(~idx,:);
    resistant_labels = plist_bnums(~idx);

%% DEFINE OUTPUT PHENOTYPE VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    phenotype_data = [sensitive_data; resistant_data];
    phenotype_labels = [sensitive_labels; resistant_labels];
    
end
