%% Confusion matrix
%%
%%
% 
% * Developer Er.Abbas Manthiri S 
% * Date  25-12-2016
% * Mail Id: abbasmanthiribe@gmail.com
% * Reference
% * <http://www.dataschool.io/simple-guide-to-confusion-matrix-terminology/ Dataschool>
% * <https://en.wikipedia.org/wiki/Confusion_matrix Wikipedia>
% 
clc
clear all
close all
warning off all
rng('default')

%% Proof
disp('Running Proof....')
n=100;m=4;
actual=round(rand(1,n)*m);
[c_matrixp,Result]= confusion.getMatrix(actual,actual);

disp('Getting Values')
Accuracy=Result.Accuracy
Error=Result.Error
Sensitivity=Result.Sensitivity
Specificity=Result.Specificity
Precision=Result.Precision
FalsePositiveRate=Result.FalsePositiveRate
F1_score=Result.F1_score
MatthewsCorrelationCoefficient=Result.MatthewsCorrelationCoefficient
Kappa=Result.Kappa




%% Multiclass demo
disp('_____________Multiclass demo_______________')
disp('Runing Multiclass confusionmat')
n=100;m=2;
actual=round(rand(1,n)*m);
predict=round(rand(1,n)*m);
[c_matrix,Result,RefereceResult]= confusion.getMatrix(actual,predict);
%
% %DIsplay off
% % [c_matrix,Result,RefereceResult]= confusionmat(actual,predict,0)

%% Two Class  demo
disp('____________Two Class  demo________________')
disp('Running Simple Confusionmat...')
n=100;m=1;
actual=round(rand(1,n)*m);
predict=round(rand(1,n)*m);
% [c_matrix,Result]= confusionmat(actual,predict)
[c_matrix,Result]= confusion.getMatrix(actual,predict);

%% Get Calculation using confusion matrix
disp('____________Get Calculation using confusion matrix________________')
n=5;
c_matrix=randi([20,40],[n,n]);
disp('confusion matrix generated')
disp(c_matrix)
disp('Running Calcualtion...')
[Result,RefereceResult]=confusion.getValues(c_matrix);
disp(Result)
disp(RefereceResult)


