function [phenotype_data, phenotype_labels] = ...
    process_transcriptome_tb(expression_data, ...
    expression_data_rowlabels, z)

    % DESCRIPTION 
    % This function transforms TB transcriptomic data into a logical matrix
    % indicating conditions where genes were differentially expressed. 
    % 
    % STEPS 
    % 1. Define main variables
    % 2. Transform data
    % 3. Define output phenotype variables
    
    % DEFINE MAIN VARIABLES %
    phenotype_num = expression_data;
    plist_bnums = expression_data_rowlabels;
    if ~exist('z','var') || isempty(z)
        z = 2;
    end
    
    % TRANSFORM DATA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORIGINAL METHOD
%     % Down-regulated genes 
%     for i = 1:length(conditions)
%         te = plist_bnums(phenotype_num(:,i) < -z);
%         cell_z1_list_t(1:length(te),i) = te;
%         lte(i) = length(te);
%     end
%     cell_z1_list_t = regexprep(cell_z1_list_t,'''','');
%     cell_z1_list_t = regexprep(cell_z1_list_t,'[]','');
%     %phenotype_labels0 = unique(cell_z1_list_t(:));
%     phenotype_labels0 = unique(cell_z1_list_t (~cellfun(@isempty, cell_z1_list_t)));
%     clear nicholslistix_t
%     for i = 1:length(conditions)
%          cell_list_column = cell_z1_list_t(:,i);
%              nicholslistix_t(:,i) = ismember(phenotype_labels0,cell_list_column(~cellfun('isempty',cell_list_column)));
%     %    nicholslistix_t(:,i) = ismember(phenotype_labels0,cell_z1_list_t(:,i));
%     end
%     %phenotype_data = nicholslistix_t;
%     % Up-regulated genes
%     for i = 1:length(conditions)
%         te1 = plist_bnums(phenotype_num(:,i) > z);
%         cell_z1_list_t1(1:length(te1),i) = te1;
%         lte1(i) = length(te1);
%     end
%     cell_z1_list_t1 = regexprep(cell_z1_list_t1,'''','');
%     %phenotype_labels1 = unique(cell_z1_list_t1(:));
%     phenotype_labels1 = unique(cell_z1_list_t1(~cellfun(@isempty, cell_z1_list_t1)));
%     clear nicholslistix_t1
%     for i = 1:length(conditions)
%              cell_list_column = cell_z1_list_t1(:,i);
%              nicholslistix_t1(:,i) = ismember(phenotype_labels1,cell_list_column(~cellfun('isempty',cell_list_column)));
%     %    nicholslistix_t1(:,i) = ismember(phenotype_labels1,cell_z1_list_t1(:,i));
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW METHOD
    % Down-regulated genes
    down_phenotype_num = phenotype_num < -z; 
    idx = sum(down_phenotype_num,2) == 0;
    down_data = down_phenotype_num(~idx,:);
    down_labels = plist_bnums(~idx);
    % Up-regulated genes
    up_phenotype_num = phenotype_num > z; 
    idx = sum(up_phenotype_num,2) == 0;
    up_data = up_phenotype_num(~idx,:);
    up_labels = plist_bnums(~idx);
    
    % DEFINE OUTPUT PHENOTYPE VARIABLES % 
%     phenotype_data = [nicholslistix_t;nicholslistix_t1];
%     phenotype_labels = [phenotype_labels0;phenotype_labels1];
    phenotype_data = [down_data; up_data];
    phenotype_labels = [down_labels; up_labels];

end