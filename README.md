# OTsensitivity
Probabilistic Sensitivity Analysis with Optimal Transport 

This is a collection of Matlab subroutines accompanying the work: "Probabilistic Sensitivity via Optimal Transport", by E. Borgonovo, A. Figalli, E. Plischke and G. Savarè, 2024.

1.	Overview

The folder contains 36 files. Of these, 33 files are MATLAB scripts (.m), two are .mat files and one is an Excel file used to produce Figure 7 postprocessing the output of one of the MATLAB scripts. 

2.	Codes

The main subroutines are:

bigwassersteintest.m. 

The boiler plate is [results]=bigwassersteintest(x,y,M,selection), where x is the input/feature data, y is the output/target data, M is the cardinality of the partition and “selection” is a list of names among a set of OT solvers (or a bitfield). 
The available choices are:
'simplex': revised transport simplex algorithm by Luenberger Ye 2016
'sinkhorn': Sinkhorn-Knopp algorithm proposed by Cuturi 2013
'swap': implementation of Puccetti 2017
'hungarian': Andrej Lopatin’s version of the Kuhn Munkres algorithm
'habr': approximate solution of Jaroslav Habr (1961) 
'lawler': Lawler’s tree-based approach for the Hungarian problem, Date Nagi 2016
'bures': Bues Wasserstein (Frechét distance) approximation
'gradient': Accelerated gradient descent for OT, An Lai Gu 2022
'divergence': symmetric Sinkhorn algorithm using Sinkhorn divergences, Jean Feydy 2020
'sinkscale': downscaled version of the Sinkhorn algorithm
'mack': Bradford / Mack min cost flow algorithm
'revsimplex': Revised simplex of Luenberger Ye 2016 ported from the transport.R package
'otpython': earth mover’s distance, LEMON solver via python OT library
'sinkmem': small-memory footprint Sinkhorn 

The output of this subroutine comprises:
A list containing for each specified solver ‘type’ the separations Wtype (use mean for OT sensitivity) and the timing Ttype

***********************
Please make sure that the following list of auxiliary solvers is available on github, if we want to offer full support:
transsimp2.m
sinkfastCost.m
wasserstein_swap.m
hungarian.m
habr.m
ungarisch.m
bwsi.m
wassgrad2.m
sinkhorn4.m
lapmack.m
revsimplex.m
sinkfast.m
bicount.m
*********************
“bigwassersteintestEpsilon.m” and “bigwassersteintestEpsilon2.m” are the same as bigwassersteintest.m. They just differ because alternative values for the penalty in the Sinkhorn function have been used, in order to test the sensitivity of results to this parameter, as per the last paragraph of Section 5.

bwsi.m: This subroutine estimates the Wasserstein-Bures metric and the entropic Wasserstein Bures metric-
The script also contains two functions, B=mysqrtm(A) and b=tracesqrtm(A), which are two original functions for computing the square root of a matrix and its trace.

ReproduceImagesTables.m: Contains the code to reproduce the Figures and Tables in the paper. It requires
mainPartitions.m
Univariate_TestWass22_Ishigami_20230426.m
Univariate_TestWass22_50corr999short20230425_2.m
Univariate_TestWass22_50corr9999short20230425_2.m
MultGaussComparisonNormalizedPartitionsBuresEntropic.m
MultGaussComparisonNormalizedSmallSimplex.m
ATOmainPartitions.m
ATOdataInvyProfitTime2.mat
Univariate_ATO_Profit_20230429.m
ATO8132Simplex.m
wassersi.m
bwsi.m

All these solvers are original implementations by the authors.
In addition, the subroutine sobolpoints.m is a quasirandom Monte Carlo generator. The subroutine deltamim.m is a fast implementation to compute the so-called \delta-importance measure of Borgonovo (2007).

3.	Reproduction of Figures
In order to run the subroutines, please follow these steps:
1)	Download all files in the “MATLAB” folder of your pc (or create a subfolder in the MATLAB folder). Be sure to save them in the same folder.
2)	To reproduce all files simultaneously, in the Matlab Prompt type “ReproduceImagesTables”.
You can also choose to reproduce specific parts.

Note: Figures 1 to 6 in the paper were generated directly by MATLAB. Therefore running the corresponding section in ReproduceImagesTables yields directly these Figures.

Figure 7 was produced by postprocessing results from the subroutine ATO8132Simplex.m in Excel. In turn the ATO8132Simplex.m loads inputs and outputs data from ATO_OT_8192.mat. The Excel file is ATO8132_Sink_Swap_Bures.xlsx. All files are uploaded in the Github.

1.	Computational requirements

The calculations were performed on a PC with the following features: 64GB RAM. With such memory, calculations did not last over 1 hour for any of the examples (the majority run in less than 10 minutes), However, we cannot exclude that, with lower RAM capacity available, the times might be longer.

