% Global Sensitivity Analysis Via Optimal Transport
% MatLab files for reproducing figures

mypause=@()pause(1); % @()drawnow % @()snapnow
% seed the random number generator with a fixed value for reproducability
%% Figure 1,2 relocated to supplementary material
%% Figure 4
ishigamifunction = @(x)sin(x(:,1)).*(1.0+0.1*x(:,3).^4)+7.0*(sin(x(:,2))).^2;
rng(hex2dec('ABADFEED'));  
Univariate_TestWass22_Ishigami_20230426

mypause
%% Figure 5
%a
rng(hex2dec('ABADFEED'));   
Univariate_TestWass22_50corr999short20230425_2
mypause
%b (needs 32Gb of free memory)
rng(hex2dec('ABADFEED'));  
Univariate_TestWass22_50corr9999short20230425_2
mypause
%% Figure 6
rng(hex2dec('ABADFEED'));  
MultGaussComparisonNormalizedPartitionsBuresEntropic
mypause
%% Figure 7
rng(hex2dec('ABADFEED'));  
MultGaussComparisonNormalizedSmallSimplex
mypause