INPUT: Simulations
The input file for simulations must be formatted as follows. Please note that you should not include the parentheses in the actual file, just the numbers described in the parentheses separated by spaces. 
Line numbers are given as a function of the number of states simulated. For some lines, such as n+7, descriptions of entries are given on multiple lines - please enter all selections on the same line.
Here is the general input text file format to run a simulation:
---------------------------------
1|	(Ignored) Simulation notes
2| 	(number of states) (number of proteins)
3|	0   k12 k13... k1n
4|	k21 0	k23... k2n
	k31 k32 0  ... k3n
	.
	.
	.
n+2|	kn1 kn2	kn3... 0
n+3|	(kint1) (kint2) ... (kintn) 
n+4|  	(# fully protected exchanging sites in state 1 - native) ('' state 2) ... ('' state n - fully unfolded) 
n+5|	(timepoint 1) (timepoint 2) ... (timepoint n)
n+6|	(seed) (number of exchangers) (KbreathingRate) (mass) (charge state) (width of Gaussians) (resolution)
n+7|	(Method of choosing which kints correspond to exchanging sites - "Default", "Smallest", "Largest", "Random_Sample", "Located_Middle", "Located_Ends", or "Manual") 
		(Method of choosing which sites are protected in which states - "Default", or "Manual")
		(initial isotope saturation of the protein, "D" or "H") 
		(normalization type - "Height" or "Area") 
		((optional) Goodness of Fit measure - defaults to "Standard", "Inverse_Time" can be chosen)
Use the following lines only if "Manual" is chosen in line n+4, entering 1 if site is exposed; 0, if site is buried and cannot exchange; or -1, if the site is subject to breathing  
n+8|	(State 1, Site 1) (State 1, Site 2) (State 1, Site 3) ... (State 1, Site (#ofexchangers))
n+9|	(State 2, Site 1) (State 2, Site 2) (State 2, Site 3) ... (State 2, Site (#ofexchangers))
		.
		.
		.
2n+7| 	(State n, Site 1) (State n, Site 2) (State n, Site 3) ... (State n, Site (#ofexchangers))
---------------------------------
If "Manual" is chosen as the kint selection method in n+7, please enter only the kints for those sites that are exchanging in line n+3.
In line n+4, the native state (state 1) should have the greatest number of protected sites.
All lines past n+7 (or past 2n+7, if "Manual" was chosen in line n+4) are ignored and can be used for notes.