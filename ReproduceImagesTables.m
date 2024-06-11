% Global Sensitivity Analysis Via Optimal Transport
% MatLab files for reproducing figures

%mypause=@()pause(); % @()drawnow % @()snapnow
% seed the random number generator with a fixed value for reproducability
%% Figure 1
mainPartitions
% only keep fig2
h=figure(1);delete(h);
h=figure(3);delete(h);
h=figure(4);delete(h);
disp('press any key to continue')
pause;
%% Figure 2
ishigamifunction = @(x)sin(x(:,1)).*(1.0+0.1*x(:,3).^4)+7.0*(sin(x(:,2))).^2;
rng(hex2dec('ABADFEED'));  
Univariate_TestWass22_Ishigami_20230426
disp('press any key to continue')
pause
%% Figure 3
%a
rng(hex2dec('ABADFEED'));   
Univariate_TestWass22_50corr999short20230425_2
disp('press any key to continue')
pause
%b (needs 32Gb of free memory)
rng(hex2dec('ABADFEED'));  
Univariate_TestWass22_50corr9999short20230425_2
disp('press any key to continue')
pause
%% Figure 4
rng(hex2dec('ABADFEED'));  
MultGaussComparisonNormalizedPartitionsBuresEntropic
disp('press any key to continue')
pause
%% Figure 5
rng(hex2dec('ABADFEED'));  
MultGaussComparisonNormalizedSmallSimplex
% Table 3
disp('average computational times per solve')
T./[ 2 5 7 8 10]

disp('press any key to continue')
pause
%% ATO simulator
% Figure 6a
ATOmainPartitions % reads data from 'ATOdataInvyProfitTime2.mat'
% only keep fig2
h=figure(1);delete(h);
h=figure(3);delete(h);
h=figure(4);delete(h);
h=figure(5);delete(h);
disp('press any key to continue')
pause;
% Figure 6b
Univariate_ATO_Profit_20230429 % runs wassersi and bwsi on datafile
% Figure 7b: from ATO8132_Sink_Swap_Bures.xlsx
ATO8132Simplex
%Figure 7 was produced by postprocessing results from the subroutine ATO8132Simplex.m in Excel. In turn the ATO8132Simplex.m loads inputs and outputs data from ATO_OT_8192.mat. The Excel file is ATO8132_Sink_Swap_Bures.xlsx. All files are uploaded in the Github.

