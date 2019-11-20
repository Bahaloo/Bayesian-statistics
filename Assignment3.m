clc
clear all

prompt = 'enter the sequence of the results, '' '', ''+'', ''+-+'', ...:';
testResults = input(prompt);
if testResults == ' '; testResults = ''; end
%%
keySet =   {':(',    ':)',   '+|:(',   '+|:)',   '-|:(',    '-|:)'};
valueSet = [0.001    0.999    0.99      0.05     1-0.99     1-0.05];
P = containers.Map(keySet,valueSet);

P('+') = P('+|:)')*P(':)')+P('+|:(')*P(':(');   % = 0.0509
P('-') = P('-|:)')*P(':)')+P('-|:(')*P(':(');   % = 1-P('+')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(testResults)
    if testResults(i) == '+'; seq(i) = 1;
    else seq(i) = 0;
    end
end
testResults = str2num(testResults);
Evidence = (P('+'))^sum(seq)*(P('-'))^(length(seq)-sum(seq));
P(':(|+') = P('+|:(')*P(':(')/P('+');  %  = 0.0194
% suppose we get a seq like 11011 or 110 or 011 ...
for i = 1:length(testResults) 
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

