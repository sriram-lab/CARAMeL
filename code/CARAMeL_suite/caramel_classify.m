function class = caramel_classify(scores, threshold, reverse)

% DESCRIPTION: 
% This function classifies combination therapies into Synergy, Additivity, 
% or Antagonism given interaction score and threshold information. 
%
% Author:   Carolina H. Chung
% Contact:  chechung@umich.edu
    
% I/O/U
%{
INPUTS: 
    1. scores:      numerical vector of interaction scores
    2. threshold:	2-element vector containing threshold values 
        
OUTPUT: 
    1. class:       ordinal array of classification labels 

EXAMPLE USAGE: 
    1. Classify interaction scores (scores) given a threshold (thresh): 
       >> class = caramel_classify(scores, thresh);
%}

    % Account for reverse (if not provided)
    if ~exist('reverse', 'var') || isempty(reverse)
        reverse = false; 
    end
    
    % Determine classes
    if reverse
        class = repmat({'Additivity'}, size(scores)); 
        class(scores < threshold(1)) = {'Antagonism'}; 
        class(scores > threshold(2)) = {'Synergy'}; 
    else
        class = repmat({'Additivity'}, size(scores)); 
        class(scores < threshold(1)) = {'Synergy'}; 
        class(scores > threshold(2)) = {'Antagonism'}; 
    end

end