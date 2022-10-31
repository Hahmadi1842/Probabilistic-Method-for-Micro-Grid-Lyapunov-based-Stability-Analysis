%% Perform mode selection on LTI system
System = G; % Define System to reduce
UpperCutoffFrequency = 131.7541098582893;
LowerCutoffFrequency = 7.832287517424242;
 
% Create option set for freqsep command
Options = freqsepOptions();
 
% Select modes between lower and upper cutoff frequencies
ReducedSystem = freqsep(System,UpperCutoffFrequency,Options);
[~,ReducedSystem] = freqsep(ReducedSystem,LowerCutoffFrequency,Options);
 
% Create comparison plot
bode(System,ReducedSystem);
