% Function for wrapping necessary variables for the logistic classifier.
% Takes the uber array and a couple other functions to output a struct with
% all desired predecision touch variables 

%INPUT: uber array
%OUTPUT: V array which is a struct for all predecision touch variables

function [V] = classifierWrapper(uberarray)

for rec = 1:length(uberarray)
    
    array=uberarray{rec};
    
    V(rec).varNames = {'theta','velocity','amplitude','setpoint','phase','deltaKappa'};
    
    [~ ,prelixGo, prelixNoGo, ~ ,~ ,~] = assist_predecisionVar(array);
    varx=[1:6];
    for f=1:length(varx)
        [~, ~,V(rec).var.hit{f}, V(rec).var.miss{f}, V(rec).var.FA{f}, V(rec).var.CR{f},~,~] = assist_vardistribution(array,varx(f),prelixGo,prelixNoGo,[-25:0],[5:25]);
        if f == 1
            V(rec).touchNum.hit=cellfun(@numel,V(rec).var.hit{1}); % count number of predecision touches using theta vals
            V(rec).touchNum.miss=cellfun(@numel,V(rec).var.miss{1});
            V(rec).touchNum.FA=cellfun(@numel,V(rec).var.FA{1});
            V(rec).touchNum.CR=cellfun(@numel,V(rec).var.CR{1});
        end
        V(rec).var.hit{f}   =      cell2mat(V(rec).var.hit{f}); %predecision variables for hit trials
        V(rec).var.miss{f}  =      cell2mat(V(rec).var.miss{f});
        V(rec).var.FA{f}    =      cell2mat(V(rec).var.FA{f});
        V(rec).var.CR{f}    =      cell2mat(V(rec).var.CR{f});
        
        
    end
end