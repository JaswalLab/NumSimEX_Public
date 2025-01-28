INPUT: Fittings
Fittings require 3 files to be run:
1) An uncertainty file (for which this file is a template)
2) A simulation file (see the other template)
3) A file of experimental HXMS data.  This data may be formatted in either of two ways:
	i) A narrow dataset (Three columns; "Time", "m/z', and "Intensity")
	ii) A series of pairs of rows; the number of pairs is the number of timepoints collected, the first row of each pair consists of m/z values; the second, of intensity values
The input file must be formatted as follows. Please note that you should not include the parentheses in the actual file, just the uncertainties of the parameters described in the parentheses, separated by spaces. 
Line numbers are given as a function of the number of states simulated.
Here is the general input text file format to run a fitting:
---------------------------------
1|	(Ignored) Fitting notes
2|	0   k12 k13... k1n
3|	k21 0	k23... k2n
	k31 k32 0  ... k3n
	.
	.
	.
n+1|	kn1 kn2	kn3... 0
n+2|	(# fully protected exchanging sites in state 1 - native) ('' state 2) ... ('' state n - 1 - last intermediate before fully unfolded)
n+3| 	(kbreathe)
n+4| 	(Optimization type - "Coordinate", "Nelder_Mead", or "ABCSMC")
n+5|	If "Coordinate" was chosen:
				(max # iterations) 
				(starting # choices per parameter to be fit) 
				(percentage changes for coordinate search) 
				(number of threads to run concurrently)
		If "Nelder_Mead" was chosen:
				(fit tolerance) 
				(max # iterations) 
				(starting # choices per parameter to be fit) 
				(percentage changes for initial simplex) 
				(number of threads to run concurrently)
		If "ABCSMC" was chosen:
				(prior type for rate constants - "unif" or "logunif")
				(maximum number of intermediate distributions)  
				(desired number of points in final distribution)
				(proportion of variance to retain)
				(initial epsilon cutoff)
				(final epsilon cutoff)
				(number of threads to run concurrently)
---------------------------------

The simulation input file gives the starting values for parameters. This can include ALL values, even ones that will not be estimated, as whether or not a parameter is estimated is set by the uncertainties file (could be called the fitting input file). 
How the uncertainties file works:
In the uncertainties file, to keep a parameter fixed, set its uncertainty value to 0. 
For coordinate descent, for each estimated parameter value, multiple starting options are explored. The uncertainty value for each parameter helps determine how many and which starting values are explored. 
Under the coordinate descent options in the uncertainties file, 2 relevant values are set for this. One sets the number of starting choices per parameter to be estimated, which is constant for all parameters to be estimated. 

(The other sets the percentage change allowed during the coordinate search.) 
Consider the following values for the number of protected sites in the native state:
14 as the starting value in the siminput file
3 as the uncertainty value in the uncertainties file
4 as the number of starting choices for each parameter in the coordinate descent portion of the uncertainties file. 
The algorithm will then construct an interval from 11 to 17 (coming from 14 +/- 3) and choose 4 equidistant points within the interval as the starting points for coordinate descent (in this case, the values 11, 13, 15, 17 would be the starting points). (From those starting points, when the algorithm needs to move, it will move 10%).  
For our purposes, we often set the number of starting parameter values to be considered as 3. This means our uncertainties value gives us the starting points easily as starting value +/- uncertainty value. Three starting points in the example above would mean threads starting at 11, 14, and 17. 


