Constant extension (naked DNA)
1. convertConstantExtension_nakedDNA.m
Script to convert raw data from constant-extension experiments with naked DNA; plots clamped tether extension-turn relation without topo II before clamp & extension, total turns applied, and topo II relaxation rate over time during clamp

2. f_3piece.m
Function defining 3-piece fit to the naked DNA extension-turn relation

Repeated winding (naked DNA)
1. convertRepeatedWinding_nakedDNA.m
Script to convert raw data from repeated winding experiments with naked DNA; plots extension-turn relations without topo II and extensions over time during repeated winding experiments with topo II

2. f_3piece.m
Function defining 3-piece fit to the naked DNA extension-turn relation

Pre-buckled relaxation
1. calculatePrebuckledTurns.m
Script to convert raw data from pre-buckled relaxation experiments; plots extension-turn relations without topo II and extensions over time during pre-buckled relaxation experiments with topo II; also calculates & plots turn state reached by topo II activity vs. duration of topo II activity in the pre-buckled regime 

2. alignOneTrace.m
Aligns relaxation events from one trace so that the buckling transition is reached at t = 0 & averages over time

3. alignAllTraces.m
Aligns relaxation events from all traces so that the buckling transition is reached at t = 0 & averages over time; used for direct conversion of pre-buckled relaxation data

4. f_3piece.m
Function defining 3-piece fit to the naked DNA extension-turn relation

Constant extension (chromatin)
1. convertConstantExtension_chromatin.m
Script to convert raw data from constant-extension experiments with chromatin; plots clamped tether extension-turn relation without topo II before clamp & extension, total turns applied, and topo II relaxation rate over time during clamp

2. f_5piece.m
Function defining 5-piece fit to the chromatin extension-turn relation

3. fit5piece.m
Script to fit chromatin extension-turn data to the 5-piece fit function; returns (turn, extension) coordinates of 4 points that uniquely define the fit

4. HCpara.m
Script to obtain extension-turn relation parameters (max extension, center, buckling positions, and slopes of the buckled regimes) from the 5-piece fit

5. getNucQuality.m
Script to determine quality of a chromatin fiber from its extension-turn relation

Repeated winding (chromatin)
1. convertRepeatedWinding_chromatin.m
Script to convert raw data from repeated winding experiments with chromatin; plots extension-turn relations without topo II and extensions over time during repeated winding experiments with topo II

2. f_5piece.m
Function defining 5-piece fit to the chromatin extension-turn relation

3. fit5piece.m
Script to fit chromatin extension-turn data to the 5-piece fit function; returns (turn, extension) coordinates of 4 points that uniquely define the fit

4. HCpara.m
Script to obtain extension-turn relation parameters (max extension, center, buckling positions, and slopes of the buckled regimes) from the 5-piece fit

5. getNucQuality.m
Script to determine quality of a chromatin fiber from its extension-turn relation

Source for other MATLAB codes to run our scripts
CircleFitByPratt.m: https://www.mathworks.com/matlabcentral/fileexchange/22643-circle-fit-pratt-method
