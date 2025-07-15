# Files source 
`GitHub\NumericalSimulationsMain\Wagawal_NumSim\ExampleNumSim\December2024_Resubmission\Experiments`

# How NumSimEx was run for doi:10.1002/pro.70045. PMID: 39865386; PMCID: PMC11761709.  (TwoStateNoBr as example) 

## Step 1: Creation of input files
- Simulation Input file (e.g. `SimInput2NoBr.txt`)
- Fitting Uncertainties file (e.g. `FittingUC2NoBr.txt`)
- Experimental HXMS data file (e.g. `B2M_data.txt`)
- Output file location (e.g. `2NoBr_Fit_Output.txt`)

## Step 2: Running a fit with NumSimEX
Fittings can be run from the command line with the following code: <br> 
`java -jar JaswalNS.jar SimInput.txt FitInput.txt ExpData.txt FitOutput.txt` <br>
or, in the case of TwoStateNoBr, <br>
`java -jar JaswalNS.jar SimInput2NoBr.txt FittingUC2NoBr.txt B2M_data.txt 2NoBr_Fit_Output.txt` <br>
Note: the argument `v` can be added to the end of each command to enable verbose printing, such as `detailed-fitting-results-2022-05-04.txt`. For example: <br>
`java -jar JaswalNS.jar SimInput2NoBr.txt FittingUC2NoBr.txt B2M_data.txt 2NoBr_Fit_Output.txt v`

## Step 3: Run a simulation with the parameters output from the fitting 

### Make a copy of the sim input file
E.g., `SimInput2NoBrCheck.txt` from `SimInput2NoBr.txt"

### Update sim input parameters
Copy optimized parameter values from the fit output file (e.g. `2NoBr_Fit_Output.txt`) into the correct places in the new sim input file (e.g., `SimInput2NoBrCheck.txt`). 

### Run the jar file again with simulation settings
Simulations can be run from the command line with the following code: <br>
`java -jar JaswalNS.jar UpdatedSimInput.txt SimOutput.txt ExpData.txt` <br>
or, in the case of TwoStateNoBr, <br>
`java -jar JaswalNS.jar SimInput2NoBrCheck.txt SimOutput2NoBr.txt ExpData.txt` <br>
