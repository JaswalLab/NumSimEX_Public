package numericalsimulations;

import java.io.*;
import java.util.*;
import java.nio.file.*;
import java.util.concurrent.atomic.AtomicInteger;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;

public class Fitting {

//Parameters for debugging/progress tracking
private boolean debug = false; //Whether to print out debugging statements
private final boolean progresstracker = true; //Whether to print out statements that show program progress
private boolean verbose = true; //Whether to print additional output files
private long starttime; //Time of creation of this object
private long start; //Gets modified for conceptual step time measurements

//Parameters for multithreading
private int numthreads; //Total number of threads the fitting will run - corresponds to the number of coordinate descent starts
private int maxconcurrentthreads = 1; //Maximum number of threads to be run concurrently
private Thread[] threads; //Stores the active daughter threads that run FittingRunnables

//Parameters for file/directory locations
private String originalsimlocation; //Location of sim input file
private String uncertaintylocation; //Location of uncertainties file
private String expdatalocation; //Location of experimental data - either format (narrow or alternating rows)
private String threaddirectory; //Location where the daughter threads will output their results
private String finaloutputlocation; //Location where the fitting summary file will be written

//Parameters for simulation
private Random randomgenerator; //Random object used to generate Unif(0,1), allows for a global seed
private int numberofstates; //Number of protein folding states (ie native + intermediates + fully unfolded)
private int numberofproteins; //Number of proteins in the simulation (see simproteins)
private int[] numberprotected; //Number of exchanging sites protected in each state
private int numberofexchangers; //Number of sites observed to exchange over the timescale of the experiment
private double[][] ratematrix; //Rate matrix of exchange between folding states
private double kbreathe; //Universal breathing rate constant for amide sites on the protein
private double[] kints; //List of kints (intrinsic rate constatnts of hydrogen exchange at each of the amide sites)
private double[][] protectionmatrix; //Assigns state specific protection values on site-by-site basis
private String protectionmethod; //Selected method to choose how to protect sites across sites
private String kintmethod; //Selected method to choose which kints to use in the simulation
private ParameterSet initialps; //Template for the starting parameter sets that get sent to FittingRunnable threads
private double[] timepoints; //List of time points at which to simulate
private double timepointdiscrepancy = 0.0005; //max allowed difference between two time points to be considered the same
private int ntps; //Number of time points
private double initialmass; //Uncharged, undeuterated mass of the protein
private int chargestate; //Charge state at which the data was collected
private double gaussianwidth; //Full Width at Half Max of the reference peak
private int resolution; //Number of m/z points at which simulations are evaluated
private double directionofexchange = -1; //1 if H to D, -1 otherwise (adding or subtracting mass)
private double M0; //Used by calculateSpectra function to provide the equivalent of a starting mass
private final double mdeuterium = 1.00627; //Mass, in daltons, of a deuterium
private boolean normalizetoheight; //Normalizes output intensity to either area or height
private double[][] expdata; //Narrow format of experimental data (Time, m/z, Intensity)
private double[] mzs; //If simulation run to compare to experimental data, the experimental data's m/zs; otherwise, estimated
private String gofmethod = "Standard"; //Way to calculate goodness of fit

//Parameters for fitting
private double[][] deltaratematrix; //Input uncertainties in values of the folding state rate matrix
private double[] deltaprotected; //Input uncertainties in the number of protected sites in each state
private double deltakbreathe; //Uncertainty in the universal breathing rate constant, kbreathe
private String optimizationmethod; //Algorithm to be used to fit simulations to the experimental data
private boolean iscoordinate; //If the algorithm to be used is Coordinate Descent
private boolean isneldermead; //If the algorithm to be used is Nelder-Mead
private boolean isabcsmc; //If the algorithm to be used is ABC SMC
private int startingchoicesperparameter; //Number of starting choices per parameter to be given
private int maxiterations; //Maximum number of iterations either algorithm can take
private double percentchange; //Percent change in each of the parameters as part of coordinate descent
private int ndist; //Maximum number of distributions to create as part of ABCSMC
private int currentdist; //Current distribution in ABCSMC
private int nfinalpoints; //Number of points to retain from each sample from each distribution
private double epsilon0; //The first epsilon cutoff to use
private double epsilonT; //The last epsilon cutoff to use
private double epsilont; //The epsilon cutoff of the current distribution
private double epsilonprev; //The epsilon cutoff the previous distribution
private double propvariance; //The proportion of variance to retain for each parameter when perturbing
private boolean[] isfit; //Which of the possible parameters that can be fit, is to be fit, as specified by an uncertainty of 0
private int numparametersfit; //Count of the paramters that are fit, as indicated in 'isfit'
private ParameterSet[] startingsets; //List of all the starting parameter sets, each is the starting point for coordinate descent
private double[][] allsimresults; //Table of parameters for all simulations run, along with thread number and gof, and weight
private double[][] prevsimresults; //Table of parameters for the last set of simulations run, along with thread number, gof, and weight
private double[][] prevsimcovariates; //Table of fit parameters in the previous population
private double[][] priordistvariances; //Covariance matrix of fit parameters in the previous population
private MultivariateNormalDistribution perturbationkernellikelihood; //Used to find likelihood (height) of the perturbation kernel
private int gofcolumn; //The column index of gof in allsimresults and prevsimresults
private int weightcolumn; //The column index of weight in allsimresults and prevsimresults
private double[] finalparameters; //The best parameters identified by fitting
private double fittolerance; //Used for the Nelder-Mead ccmponent of fitting
private double[][] P0params; //the distributional parameters of the P0 distribution
private String priortype = "unif"; //The distribution of the prior (support has same endpoints regardless)

    public Fitting (String siminput, String fittinginput, String experimentaldata, String finaloutput, boolean shoulddebug, 
            boolean shouldverbose){ //Constructor using file paths of input/output files
        //Sim Input: location of simulation input, 
        //Fitting Input: location of fiting uncertainties/algorithmic parameters, 
        //Experimental Data: location of exp data input, File Directory: location of thread results, 
        //Final Output: location of final output, 
        //shoulddebug: whether to debug, verbose: whether to create additional output files
        starttime = start = System.nanoTime();
        debug = shoulddebug;
        verbose = shouldverbose;
        setImportLocations(siminput, fittinginput, experimentaldata, finaloutput);
        if(debug || progresstracker){
            System.out.println("The fitting program has completed step 1: Make a Fitting object (" 
                    + (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
        }
        
        start = System.nanoTime();
        importFittingFiles();
        if(debug || progresstracker){
            System.out.println("The fitting program has completed step 2: Import all data from the specified locations (" 
                    + String.format("%.3f", (((double )System.nanoTime() - start)/1000000000)) + " seconds)");
        }
        
        if(isabcsmc){
            ABCSMC();
        }
        else{
            start = System.nanoTime();
            initializeParameterSets(); 
            if(debug || progresstracker){
                System.out.println("The fitting program has completed step 3: Generate all initial sets of parameters for " 
                        + optimizationmethod + " (" + (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
            }

            start = System.nanoTime();
            initializeThreads();
            if(debug || progresstracker){
                System.out.println("The fitting program has completed step 4: Initialize and run the " + numthreads + 
                        " threads to parallelize the work ("  + 
                        String.format("%.1f", (((double )System.nanoTime() - start)/(1000000000 * 3600))) + " seconds)");
            }

            start = System.nanoTime();
            importThreadResults(0, numthreads, 0);
            if(debug || progresstracker){
                System.out.println("The fitting program has completed step 5: Read in the output of each of the " + numthreads + 
                        " threads ("  + (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
            }

            start = System.nanoTime();
            additionalAnalysis(); 
            if(debug || progresstracker){
                System.out.println("The fitting program has completed step 6: Conduct additional analysis on the results of fitting (" + 
                        (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
            }

            start = System.nanoTime();
            printResults(); 
            if(debug || progresstracker){
                System.out.println("The fitting program has completed step 7: Write output of fitting to '" + finaloutputlocation + 
                        "' (" + (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
            }
        }
    }

    //Stores locations of input files
    private void setImportLocations(String siminput, String fittinginput, String experimentaldata, String finaloutput){
       originalsimlocation = siminput;
       uncertaintylocation = fittinginput;
       expdatalocation = experimentaldata;
       if(Paths.get(finaloutput).getParent() != null){
           threaddirectory = Paths.get(finaloutput).getParent() + FileSystems.getDefault().getSeparator();
       }else{
           threaddirectory = "";
       }
       finaloutputlocation = finaloutput;
    }

    //Read in and set simulation parameters, experimental data, and uncertainties/algorithmic parameters
    private void importFittingFiles(){ 
        importSim();  
        importExpData();
        importUncertainties();
        analyzeUncertainties();        
    }

    private void calculateNumParameters(){ //Calculate total number of parameters fit
        for(int i = 0; i < isfit.length; i++){
            if(isfit[i]){
                numparametersfit++;
            }
        }
    }

    private void initializeParameterSets(){ //Create the initial parameter sets on which to run coordinate descent
        initializeSimPS();
        startingsets = new ParameterSet[(int)Math.pow(startingchoicesperparameter, numparametersfit)];
        int[][] combinatorialmatrix = initializeCombinatorialMatrix(startingchoicesperparameter, numparametersfit);
        //nrows is number of choices per parameter, ncols is number of parameters fit
        double[][] startingparameterchoices = initializeStartingChoices(); 
        for(int i = 0; i < startingsets.length; i++){
            ParameterSet currentps = initialps.deepCopy();
            currentps.setRandomGenerator(runif() * 1000000000);
            int parameterfit = 0;
            int position = 0;
            //Set of rate matrix parameters
            for(int k = 0; k < numberofstates; k++){
                for(int j = 0; j < numberofstates; j++){
                    if(k != j){
                        if(isfit[position]){
                            currentps.setRateMatrixPosition(k, j, 
                                    startingparameterchoices[combinatorialmatrix[i][parameterfit]][parameterfit]);
                            parameterfit++;

                        }
                        position++;
                    }
                }
            }

            //Set of protection information
            for(int j = 0; j < numberofstates-1; j++){
                if(isfit[position]){
                    currentps.setNumProtected(j, (int) startingparameterchoices[combinatorialmatrix[i][parameterfit]][parameterfit]);
                    parameterfit++;
                }
                position++;
            }

            //Kbreathe
            if(isfit[position]){
                currentps.setKbreathe(startingparameterchoices[combinatorialmatrix[i][parameterfit]][parameterfit]);
                parameterfit++;
            }
            position++;
            startingsets[i] = currentps;
        }
    }

    private void initializeSimPS(){ //Creates initial parameter set using various parameters
        //Contains the integer parameters: number of proteins, charge state, and resolution
        int[] inputintparams = new int[]{numberofproteins, chargestate, resolution};
        //Contains the double parameters: mass, gaussian width, initial isotype saturation, and kbreathe
        double[] inputdoubleparams = new double[]{initialmass, gaussianwidth, directionofexchange, kbreathe};
        //Contains whether or not to normalze to height
        boolean[] inputbooleanparams = new boolean[]{normalizetoheight};
        initialps = new ParameterSet(protectionmatrix, ratematrix, kints, inputintparams, inputdoubleparams, inputbooleanparams, 
                timepoints, randomgenerator, mzs, gofmethod);
    }

    //Creates an enumeration of every element of the direct product group Zn*Zn*...*Zn 
    //where n is zmod, and the number of groups is numzmod
    private int[][] initializeCombinatorialMatrix(int zmod, int numzmod){ 
        int[][] directproductgroup = new int[(int) Math.pow(zmod, numzmod)][numzmod];
        for(int i = 0; i < directproductgroup.length; i++){ //Each row is an element of the direct product group
            for(int j = 0; j < directproductgroup[0].length; j++){
                directproductgroup[i][j] = ((int) Math.floor((i)/Math.pow(zmod, (numzmod - j - 1))))%zmod;
            }
        }
        return directproductgroup;
    }

    private double[][] initializeStartingChoices(){ //Creates table with a set of possible starting values for each parameter fit
        double[][] startingchoices = new double[startingchoicesperparameter][numparametersfit];
        int position = 0;
        int currentfitparam = 0;
        //Set of rate matrix parameters
        for(int k = 0; k < numberofstates; k++){
            for(int j = 0; j < numberofstates; j++){
                if(k != j){
                    if(isfit[position]){
                        if(startingchoicesperparameter == 1){
                            startingchoices[0][currentfitparam] = ratematrix[k][j];
                        }
                        else{
                            for(int i = 0; i < startingchoicesperparameter; i++){
                                startingchoices[i][currentfitparam] = (ratematrix[k][j] - deltaratematrix[k][j]) + 
                                        ((2*deltaratematrix[k][j])/(startingchoicesperparameter-1))*i;
                            }
                        }
                        currentfitparam++;
                   }
                   position++;
                }
            }
        }

        //Set of protection information
        for(int j = 0; j < numberofstates-1; j++){
            if(isfit[position]){
                if(startingchoicesperparameter == 1){
                    startingchoices[0][currentfitparam] = numberprotected[j];
                }
                else{
                    for(int i = 0; i < startingchoicesperparameter; i++){
                        startingchoices[i][currentfitparam] = (numberprotected[j] - deltaprotected[j]) + 
                                (deltaprotected[j] + 
                                Math.min(deltaprotected[j], numberofexchangers-numberprotected[j]))/(startingchoicesperparameter-1)*i;
                    }
                }
                currentfitparam++;
           }
           position++;
        }

        //Kbreathe
        if(isfit[position]){
            if(startingchoicesperparameter == 1){
                startingchoices[0][currentfitparam] = kbreathe;
            }
            else{
                for(int i = 0; i < startingchoicesperparameter; i++){
                    startingchoices[i][currentfitparam] = (kbreathe - deltakbreathe) + 
                            (2*deltakbreathe)/(startingchoicesperparameter-1)*i;
                }
            }
        }
        position++;
        
        if(debug){
            System.out.println("The possible initial values for each parameter to be fit were the following: ");
            for(int i = 0 ; i < startingchoices.length; i++){
                for(int j = 0 ; j < startingchoices[0].length; j++){
                    System.out.print(startingchoices[i][j] + " ");
                }
                System.out.println();
            }
        }
        return startingchoices;
    }

    private void ABCSMC(){ //Runs the ABCSMC algorithm
        start = System.nanoTime();
        long temptime = System.nanoTime();
        numthreads = ndist*maxconcurrentthreads;
        initializeSimPS();
        initializeOtherABCInputs();
        while(currentdist < ndist && epsilonprev > epsilonT){
            samplePoints(currentdist);
            analyzePoints(currentdist);
            updateEpsilon();
            currentdist++;
        }
        if(debug || progresstracker){
            System.out.println("The fitting program has completed step 3: Evaluate all distributions (" + 
                    (int) (((double )System.nanoTime() - temptime)/1000000000) + " seconds)");
        }
        start = System.nanoTime();
        exportFullABC();
        if(debug || progresstracker){
            System.out.println("The fitting program has completed step 4: Write output of fitting to '" + finaloutputlocation + 
                    "' ("  + String.format("%.3f", (((double )System.nanoTime() - start)/1000000000)) + " seconds)");
        }
    }

    //Updates the epsilon cutoff for the next population to be the median epsilon of the current population
    private void updateEpsilon(){ 
        Arrays.sort(allsimresults, (a, b) -> Double.compare(b[gofcolumn], a[gofcolumn]));
        epsilonprev = epsilont;
        if(nfinalpoints % 2 == 0){
            epsilont = (Math.abs(1-allsimresults[nfinalpoints/2-1][gofcolumn]) + 
                    Math.abs(1-allsimresults[nfinalpoints/2][gofcolumn]))/2;
        }else{
            epsilont = Math.abs(1-allsimresults[nfinalpoints/2][gofcolumn]);
        }
        epsilont = Math.max(epsilont, epsilonT);
    }

    private void initializeOtherABCInputs(){
        epsilont = epsilon0;
        epsilonprev = epsilon0+1;
        allsimresults = new double[1][(int) Math.pow(numberofstates, 2) + 3];
        priordistvariances = new double[numparametersfit][numparametersfit];
        if(debug){
            System.out.println("The ends of the support of the prior distribution were: ");
            for(int i = 0; i < P0params.length; i++){
                for(int j = 0; j < P0params[0].length; j++){
                    System.out.print(P0params[i][j] + " ");
                }
                System.out.println();
            }
        }
    }
    
    private void samplePoints(int dist){
        start = System.nanoTime();
        threads = new Thread[maxconcurrentthreads];
        AtomicInteger currentpop = new AtomicInteger(0);
        for(int i = 0; i < maxconcurrentthreads; i++){
            //Create new thread on a Runnable object
            Thread newfittingthread;
            newfittingthread = new Thread(new ABCSMCRunnable(currentdist, initialps.deepCopy(), expdata, isfit, P0params, priortype, 
                    allsimresults, priordistvariances, propvariance, currentpop, epsilont, nfinalpoints, dist*maxconcurrentthreads+i, 
                    threaddirectory, new Random(randomgenerator.nextLong()), debug));
            //Store the thread to access it
            threads[i] = newfittingthread;
            //Start thread
            threads[i].start(); 
        }
        for(int i = 0; i < maxconcurrentthreads; i++){ //Wait for each of the threads to finish
            try{
                threads[i].join(); 
            }
            catch(Exception e){
                System.out.println("Thread " + i + " is still running");
            }
        }
        if(debug || progresstracker){
            System.out.println("The fitting program has completed distance evaluation for population " + (dist+1) + " of " + ndist);
        }
    }
    
    private void analyzePoints(int population){
        if(currentdist != 0){
            prevsimresults = new double[allsimresults.length][allsimresults[0].length];
            for(int i = 0; i < allsimresults.length; i++){
                for(int j = 0; j < allsimresults[0].length; j++){
                    prevsimresults[i][j] = allsimresults[i][j];
                }
            }
        }
        importThreadResults(population*maxconcurrentthreads, population*maxconcurrentthreads + maxconcurrentthreads, 1);
        Arrays.sort(allsimresults, (a, b) -> Double.compare(b[weightcolumn], a[weightcolumn]));
        makeSampleWeights();
        exportPopulationResults(population);
        calculateVariances();
    }

    private void calculateVariances(){ //calculates the sample variance of each of the parameters from the distribution
        //private double[][] prevsimcovariates; //Table of fit parameters in the previous population
        //private double[][] priordistvariances; //Covariance matrix of fit parameters in the previous population
        prevsimcovariates = new double[nfinalpoints][numparametersfit];
        int col = 0;
        for(int i = 0; i < isfit.length; i++){
            if(isfit[i]){
                for(int j = 0; j < nfinalpoints; j++){
                    prevsimcovariates[j][col] = allsimresults[j][i+1];
                }
                col++;
            }
        }
        priordistvariances = new Covariance(prevsimcovariates, true).getCovarianceMatrix().getData();
        if(debug){
            System.out.println("The covariates (values of fit parameters) from the previous population were: ");
            for(int i = 0; i < prevsimcovariates.length; i++){
                for(int j = 0; j < prevsimcovariates[0].length; j++){
                    System.out.print(prevsimcovariates[i][j] + " ");
                }
                System.out.println();
            }
            System.out.println("The covariances of the fit parameters in the previous populations were: ");
            for(int i = 0; i < priordistvariances.length; i++){
                for(int j = 0; j < priordistvariances[0].length; j++){
                    System.out.print(priordistvariances[i][j] + " ");
                }
                System.out.println();
            }
        }
    }
    
    private void exportPopulationResults(int population){
        //Allows for output of results of every thread
        String fulloutputlocation = threaddirectory + "Populationresults" + population + ".txt";
        try{
            PrintWriter writer = new PrintWriter(new FileWriter(fulloutputlocation), true); //Makes txt file, adds variable headings
            writer.print("Population ");
            writer.print("Thread ");
            //Set of rate matrix parameters
            for(int k = 0; k < numberofstates; k++){
                for(int j = 0; j < numberofstates; j++){
                    if(k != j){
                        writer.print("k" + (k+1) + "" + (j+1) + " ");
                    }
                }
            }
            for(int i = 0; i < numberofstates-1; i++){
                writer.print("NumberProtectedinState" + (i+1) + " ");
            }
            writer.print("kbreathe ");
            writer.print("GoF ");
            writer.print("Weight ");
            writer.println("Epsilon");
            makeSampleWeights();
            for(int i = 0; i < allsimresults.length; i++){
                writer.print(population + " ");
                for(int j = 0; j < allsimresults[0].length; j++){
                    writer.print(allsimresults[i][j] + " ");
                }
                writer.println(epsilont);
            }
            writer.close();
        }
        catch(Exception e){
            System.out.println("There was an error writing out the full results");
        }
    }
    
    private void makeSampleWeights(){ //give the retained points appropriate weights
        for(int i = 0; i < allsimresults.length; i++){
            if(i >= nfinalpoints){
                allsimresults[i][weightcolumn] = 0;
            }
        }
        double normfactor = 0;
        if(currentdist == 0){
            for(int i = 0; i < nfinalpoints; i++){
                normfactor = normfactor + allsimresults[i][weightcolumn];
            }
        }else{
            double[][] paramvariances = new double[priordistvariances.length][priordistvariances[0].length];
            for(int i = 0; i < paramvariances.length; i++){
                for(int j = 0; j < paramvariances[0].length; j++){
                    paramvariances[i][j] = Math.sqrt(propvariance)*priordistvariances[i][j];
                }
            }
            perturbationkernellikelihood = new MultivariateNormalDistribution(new double[paramvariances.length], paramvariances);
            for(int i = 0; i < nfinalpoints; i++){
                allsimresults[i][weightcolumn] = heightP0(i)/heightDenominator(i);
                normfactor = normfactor + allsimresults[i][weightcolumn];
            }
        }
        for(int i = 0; i < nfinalpoints; i++){
            allsimresults[i][weightcolumn] = allsimresults[i][weightcolumn]/normfactor;
        }
    }

    private double heightP0(int row){ //Height of the prior distribution
        //double[][] P0params the distributional parameters of the P0 distribution
        //row 0 the mins of the uniform distribution
        //row 1 the maxes of the uniform distribution
        //ncolumns total number of parameters
        //isfit whether each parameter is fit
        double height = 1;
        int position = 0;
        double[] parameters = allsimresults[row];
        //Set of rate matrix parameters
        for(int k = 0; k < numberofstates; k++){
            for(int j = 0; j < numberofstates; j++){
                if(k != j){
                    if(isfit[position]){
                        if(priortype.equalsIgnoreCase("unif")){
                            height = height*dunif(P0params[0][position], P0params[1][position], parameters[position]);
                        }else{
                            height = height*dlogunif(P0params[0][position], P0params[1][position], parameters[position]);
                        }
                    }
                    position++;
                }
            }
        }

        //Set of protection information
        for(int j = 0; j < numberofstates-1; j++){
            if(isfit[position]){
                height = height*dunif(P0params[0][position], P0params[1][position], parameters[position]);
            }
            position++;
        }

        //Kbreathe
        if(isfit[position]){
            if(priortype.equalsIgnoreCase("unif")){
                height = height*dunif(P0params[0][position], P0params[1][position], parameters[position]);
            }else{
                height = height*dlogunif(P0params[0][position], P0params[1][position], parameters[position]);
            }
        }
        return(height);
    }

    private double heightDenominator(int row){ //See Toni 2009 for explanation of this denominator
        double sum = 0;
        for(int i = 0; i < nfinalpoints; i++){
            double[] perturbeddistances = new double[numparametersfit];
            double kernelheight = 0;
            int position = 0;
            int parameterfit = 0;
            //Set of rate matrix parameters
            for(int k = 0; k < numberofstates; k++){
                for(int j = 0; j < numberofstates; j++){
                    if(k != j){
                        if(isfit[position]){
                            perturbeddistances[parameterfit] = allsimresults[row][position+1]-prevsimresults[i][position+1];
                            parameterfit++;
                        }
                        position++;
                    }
                }
            }

            //Set of protection information
            for(int j = 0; j < numberofstates-1; j++){
                if(isfit[position]){
                    perturbeddistances[parameterfit] = allsimresults[row][position+1]-prevsimresults[i][position+1];
                    parameterfit++;
                }
                position++;
            }

            //Kbreathe
            if(isfit[position]){
                perturbeddistances[parameterfit] = allsimresults[row][position+1]-prevsimresults[i][position+1];
            }
            kernelheight = perturbationkernellikelihood.density(perturbeddistances);
            double weightedkernelheight = prevsimresults[i][weightcolumn]*kernelheight;
            sum = sum + weightedkernelheight;
        }
        return(sum);
    }

    private double dunif(double min, double max, double x){
        double density = 1.0/Math.abs(max - min);
        return(density);
    }
    
    private double dlogunif(double tempmin, double tempmax, double x){
        double logmin = Math.log(Math.min(tempmin, tempmax));
        double logmax = Math.log(Math.max(tempmin, tempmax));
        double density = 1.0/(x*(logmax-logmin));
        return(density);
    }
       
    private void exportFullABC(){
        ArrayList<String> abcresults = new ArrayList<>();
        for(int i = 0; i < currentdist; i++){
            Scanner sc = new Scanner("");
            String filename = threaddirectory + "Populationresults" + i + ".txt";
            File inputfile = new File(filename);
            try{
                sc = new Scanner(inputfile);
            }
            catch(Exception e){
                System.out.println("Couldn't make scanner object");
            }
            if(i != 0){
                sc.nextLine(); //ignore variable header
            }
            while(sc.hasNextLine()){
                abcresults.add(sc.nextLine());
            }
            sc.close();
        }
        try{
            PrintWriter writer = new PrintWriter(new FileWriter(finaloutputlocation), true); //Makes txt file, adds variable headings
            for(int i = 0; i < abcresults.size(); i++){
                writer.println(abcresults.get(i));
            }
            writer.close();
            numthreads = currentdist*maxconcurrentthreads;
            for(int i = 0; i < numthreads; i++){
                String filename = threaddirectory + "Threadresults" + i + ".txt";
                File inputfile = new File(filename);
                if(!inputfile.delete()){
                    System.out.println("Failed to delete '" + filename + "'");
                }
            }
            for(int i = 0; i < currentdist; i++){
                String filename = threaddirectory + "Populationresults" + i + ".txt";
                File inputfile = new File(filename);
                if(!inputfile.delete()){
                    System.out.println("Failed to delete '" + filename + "'");
                }
            }
        }
        catch(Exception e){
            System.out.println("There was an error writing out the full results");
        }
    }
    
    private void initializeThreads(){ //Creates a thread for each starting parameter set to run coordinate descent
        numthreads = startingsets.length;
        if(debug){
            System.out.println("The thread initialization method for " + numthreads + " threads was reached");
        }    
        threads = new Thread[maxconcurrentthreads];
        for(int i = 0; i < numthreads/maxconcurrentthreads; i++){
            for(int j = 0; j < maxconcurrentthreads; j++){
                //Create new thread on a Runnable object
                Thread newfittingthread;
                if(iscoordinate){
                    newfittingthread = new Thread(new CoordinateDescentRunnable(startingsets[i*maxconcurrentthreads+j], expdata, isfit, 
                            i*maxconcurrentthreads+j, maxiterations, percentchange, threaddirectory));
                }
                else{//change to 'else if(isneldermead)', if you want to add new fitting algorithm
                    newfittingthread = new Thread(new NelderMeadRunnable(startingsets[i*maxconcurrentthreads+j], expdata, isfit, 
                            numparametersfit, percentchange, maxiterations, fittolerance, 
                            i*maxconcurrentthreads+j, threaddirectory, debug));
                }
                //Store the thread to access it
                threads[j] = newfittingthread;
                //Start thread
                threads[j].start(); 
            }
            for(int j = 0; j < maxconcurrentthreads; j++){ //Wait for each of the threads to finish
                try{
                    threads[j].join(); 
                }
                catch(Exception e){
                    System.out.println("Thread " + i + " is still running");
                }
            }
            System.out.println("----------------------------------------------------------------------");
            System.out.println("Estimated time remaining : " + 
                    String.format("%.1f", ((((double)(System.nanoTime() - starttime))*
                            (numthreads - (i+1)*maxconcurrentthreads))/((i+1)*maxconcurrentthreads*(double)1000000000*60))) + 
                    " minutes");
            System.out.println("----------------------------------------------------------------------");
        }
        //Repeat process above for the last set of threads
        for(int i = 0; i < numthreads-(numthreads/maxconcurrentthreads)*maxconcurrentthreads; i++){
            Thread newfittingthread;
            if(iscoordinate){
                newfittingthread = new Thread(
                        new CoordinateDescentRunnable(startingsets[(numthreads/maxconcurrentthreads)*maxconcurrentthreads+i], 
                                expdata, isfit, i+(numthreads/maxconcurrentthreads)*maxconcurrentthreads, 
                                maxiterations, percentchange, threaddirectory));
            }
            else{//change to 'else if(isneldermead)', if you want to add new fitting algorithm
                newfittingthread = new Thread(
                        new NelderMeadRunnable(startingsets[(numthreads/maxconcurrentthreads)*maxconcurrentthreads+i], 
                                expdata, isfit, numparametersfit, percentchange, maxiterations, fittolerance, 
                                i+(numthreads/maxconcurrentthreads)*maxconcurrentthreads, threaddirectory, debug));
            }
            threads[i] = newfittingthread;
            threads[i].start();
        }
        for(int i = 0; i < numthreads-(numthreads/maxconcurrentthreads)*maxconcurrentthreads; i++){
            try{
                threads[i].join();
            }
            catch(Exception e){
                System.out.println("Thread " + i + " is still running");
            }
        }
        if(debug){
            System.out.println("All " + numthreads + " threads have been run");
        }  
    }

    private void importThreadResults(int startindex, int endindex, int includeweight){ //Imports each of files created by the threads
        ArrayList<String> simresults = new ArrayList<>();
        for(int i = startindex; i < endindex; i++){
            try{
                Scanner sc = new Scanner("");
                String filename = threaddirectory + "Threadresults" + i + ".txt";
                File inputfile = new File(filename);
                sc = new Scanner(inputfile);
                //System.out.println(filename);
                sc.nextLine(); //ignore variable header
                while(sc.hasNextLine()){
                    simresults.add(sc.nextLine());
                }
                sc.close();
                }
            catch(Exception e){
                System.out.println("Couldn't make scanner object");
            }
        }
        allsimresults = new double[simresults.size()][(int) Math.pow(numberofstates, 2) + 2 + includeweight];
        Scanner sc = new Scanner("");
        for(int i = 0; i < simresults.size(); i++){
            sc = new Scanner(simresults.get(i));
            int col = 0;
            while(sc.hasNextDouble()){
                allsimresults[i][col] = sc.nextDouble();
                col++;
            }
        }
        sc.close();
        if(debug){
            System.out.println("The top of the simulation results set merged from the threads is: ");
            for(int i = 0; i < 5; i++){
                for(int j = 0; j < allsimresults[0].length; j++){
                    System.out.print(allsimresults[i][j] + " ");
                }
                System.out.println();
            }
        }
    }

    private void additionalAnalysis(){ //As of now, just finds best performing set of parameters from all threads run
        finalparameters = new double[allsimresults[0].length];
        int locationofbestresult = 0;
        gofcolumn = (int) Math.pow(numberofstates, 2) + 1;
        double bestgofsofar = allsimresults[0][gofcolumn];
        for(int i = 0; i < allsimresults.length; i++){
            if(allsimresults[i][gofcolumn] > bestgofsofar){
                locationofbestresult = i;
                bestgofsofar = allsimresults[i][gofcolumn];
            }
        }
        for(int i = 0; i < allsimresults[0].length; i++){
            finalparameters[i] = allsimresults[locationofbestresult][i];
        }
        if(debug){
            System.out.println("The best parameters were: ");
            for(int i = 0; i < finalparameters.length; i++){
                System.out.print(finalparameters[i] + " ");
            }
            System.out.println();
        }
        if(verbose){
            runValidationSim();
        }
    }
    
    private void runValidationSim(){
        ParameterSet finalps = initialps.deepCopy();
        int position = 1;
        //Set of rate matrix parameters
        for(int k = 0; k < numberofstates; k++){
            for(int j = 0; j < numberofstates; j++){
                if(k != j){
                    finalps.setRateMatrixPosition(k, j, finalparameters[position]);
                    position++;
                }
            }
        }
        //Set of protection information
        for(int j = 0; j < numberofstates-1; j++){
            finalps.setNumProtected(j, (int) finalparameters[position]);
            position++;
        }
        //Kbreathe
        finalps.setKbreathe(finalparameters[position]);
        
        Sim validationSim = new Sim(finalps);
        double[][] simdata = validationSim.getResults();
        PrintWriter validationwriter;
        try{
            validationwriter = new PrintWriter(new FileWriter(threaddirectory + "best-sim-results-" + 
                    java.time.LocalDate.now() + ".txt"), true);
            validationwriter.println("Time mz Intensity Source");
            for(int i = 0; i < expdata.length; i++){
                validationwriter.println(expdata[i][0] + " " + expdata[i][1] + " " + expdata[i][2] + " Experiment");
            }
            for(int i = 0; i < simdata.length; i++){
                validationwriter.println(simdata[i][0] + " " + simdata[i][1] + " " + simdata[i][2] + " Simulation");
            }
            validationwriter.close();
        }catch(Exception e){
            System.out.println("There was an error saving a simulation with the best parameters found");
        }
        if(verbose){
            printGraphCode();
        }
    }
    
    private void printGraphCode(){ //Make an R file to graph the results
        String rlocation = threaddirectory + "graphingsimcode.R";
        String pnglocation = "graphsim.png";
        PrintWriter graphingcodewriter;
        try{
            graphingcodewriter = new PrintWriter(new FileWriter(rlocation), true);
            graphingcodewriter.println("# Code to graph " + "best-sim-results-" + java.time.LocalDate.now() + ".txt" + " \n");
            graphingcodewriter.println("library(tidyverse) \n"
                                    + "library(grDevices) \n"
                                    + "coldata <- read.table('" + "best-sim-results-" + java.time.LocalDate.now() + 
                                        ".txt" + "', header = TRUE)%>% \n"
                                    + "  mutate(Time = as.factor(signif(Time, digits = 5)), \n"
                                    + "  Source = as.factor(Source)) \n"
                                    + "simplot <- ggplot(coldata, aes(x = mz, y = Intensity, color = Source)) + \n"
                                    + "  geom_line(size = 3) + \n"
                                    + "  facet_wrap(~ Time) + \n"
                                    + "  theme_bw() + \n"
                                    + "  theme(axis.title.y = element_blank(), \n"
                                    + "        axis.ticks.y = element_blank(), \n"
                                    + "        axis.text.y = element_blank(), \n"
                                    + "        axis.text.x = element_text(angle = 90), \n"
                                    + "        panel.grid = element_blank(), \n"
                                    + "        legend.position = 'top', \n"
                                    + "        text = element_text(size = 60) \n"
                                    + "  ) + \n"
                                    + "  xlab('m/z') \n"
                                    + "png('" + pnglocation + "', width = 3840, height = 2160) \n"
                                    + "simplot \n"
                                    + "dev.off() ");
            graphingcodewriter.close();
        }catch(Exception e){
            System.out.println("There was an error saving a simulation with the best parameters found");
        }
    }

    private void printResults(){ //Write out summary file about the fitting conducted and best parameter set
        PrintWriter finalwriter;
        try{
            finalwriter = new PrintWriter(new FileWriter(finaloutputlocation), true);
            finalwriter.println("This fitting analysis finished " + new java.util.Date() + " and took " + 
                    String.format("%.1f", (((double )System.nanoTime() - starttime)/1000000000)) + " seconds");
            finalwriter.println("The fitting was performed on the following experimental data: " + expdatalocation + 
                    " using " + numthreads + " threads");
            finalwriter.println("Within the uncertainties given, the best parameters " + 
                    "describing the protein's HXMS behavior were as follows: ");
            finalwriter.println();         
    
            int position = 1;
            //Set of rate matrix parameters
            for(int k = 0; k < numberofstates; k++){
                for(int j = 0; j < numberofstates; j++){
                    if(k != j){
                        finalwriter.println("k" + (k+1) + "" + (j+1) + " = " + finalparameters[position]);
                        position++;
                    }
                }
            }

            //Set of protection information
            for(int i = 0; i < numberofstates-1; i++){
                finalwriter.println("The number of protected sites in state " + (i+1) + " = " + finalparameters[position]);
                position++;
            }

            //Kbreathe
            finalwriter.println("kbreathe = " + finalparameters[position]);
            position++;
            finalwriter.println();
            finalwriter.println("The goodness of fit measure for this fitting was " + finalparameters[position]);
            position++;
            finalwriter.println("The fitting method used was " + optimizationmethod);
            finalwriter.close();
        }
        catch(Exception e){
            System.out.println("There was an error writing out the fitting results");
        }
        if(verbose){
            //Allows for output of results of every thread
            String fulloutputlocation = threaddirectory + "detailed-fitting-results-" + java.time.LocalDate.now() + ".txt";
            try{
                PrintWriter writer = new PrintWriter(new FileWriter(fulloutputlocation), true); //Makes txt file, adds variable headings
                writer.print("Thread ");
                //Set of rate matrix parameters
                for(int k = 0; k < numberofstates; k++){
                    for(int j = 0; j < numberofstates; j++){
                        if(k != j){
                            writer.print("k" + (k+1) + "" + (j+1) + " ");
                        }
                    }
                }
                for(int i = 0; i < numberofstates-1; i++){
                    writer.print("NumberProtectedinState" + (i+1) + " ");
                }
                writer.print("kbreathe ");
                writer.println("GoF");
                for(int i = 0; i < allsimresults.length; i++){
                    for(int j = 0; j < allsimresults[0].length-1; j++){
                        writer.print(allsimresults[i][j] + " ");
                    }
                    writer.print(allsimresults[i][allsimresults[0].length-1]);
                    writer.println();
                }
                writer.close();
            }
            catch(Exception e){
                System.out.println("There was an error writing out the full results");
            }
        }
        for(int i = 0; i < numthreads; i++){
            String filename = threaddirectory + "Threadresults" + i + ".txt";
            File inputfile = new File(filename);
            if(!inputfile.delete()){
                System.out.println("Failed to delete '" + filename + "'");
            }
        }
    }
    
    private void setSeed(double seed){ //Sets seed for the simulation
        randomgenerator = new Random();
        randomgenerator.setSeed((long) seed);
    }

    private void fixTimepoints(){  //Arrange the timepoints in order, gives error if time points occur at unallowable times
        Arrays.sort(timepoints);
        if(timepoints[0] <= 0){
            System.out.println("At least one of the time points entered is too small (<= 0).  Please fix this and try again");
        }
    }
    
    //Select subset of kints to use from the provided list, using desired selection criteria
    private void chooseKints(String kintmethod){ 
        if(kints.length == numberofexchangers){                
        }
        //Manual selection of exchanging sites, requires the correct number of kints to be input in the input file
        else if(kintmethod.equalsIgnoreCase("Manual")){ 
            System.out.println("An incorrect number of kints were inputted; please try again");
        }
        else if (kintmethod.equalsIgnoreCase("Default")){ //Takes a random sample of the smallest 50% of the kints
            if(kints.length > 2*numberofexchangers){
                int halfway = kints.length/2;
                Arrays.sort(kints);
                ArrayList<Integer> templocations = sample(numberofexchangers, halfway);
                ArrayList<Double> tempkints = new ArrayList<>(numberofexchangers);
                for(int i = 0; i < numberofexchangers; i++){
                    tempkints.add(kints[templocations.get(i)]);
                }
                kints = new double[numberofexchangers];
                for(int i = 0; i < numberofexchangers; i++){
                    kints[i] = tempkints.get(i);
                }  
            }else{
                ArrayList<Integer> templocations = sample(numberofexchangers, kints.length);
                ArrayList<Double> tempkints = new ArrayList<>(numberofexchangers);
                for(int i = 0; i < numberofexchangers; i++){
                    tempkints.add(kints[templocations.get(i)]);
                }
                kints = new double[numberofexchangers];
                for(int i = 0; i < numberofexchangers; i++){
                    kints[i] = tempkints.get(i);
                } 
            }
        }
        else if (kintmethod.equalsIgnoreCase("Smallest")){ //Uses the smallest magnitude kints inputted
            Arrays.sort(kints);
            ArrayList<Double> tempkints = new ArrayList<>(numberofexchangers);
            for(int i = 0; i < numberofexchangers; i++){
                tempkints.add(kints[i]);
            }
            kints = new double[numberofexchangers];
            for(int i = 0; i < numberofexchangers; i++){
                kints[i] = tempkints.get(i);
            }
        }
        else if (kintmethod.equalsIgnoreCase("Largest")){ //Uses the largest magnitude kints inputted
            Arrays.sort(kints);
            ArrayList<Double> tempkints = new ArrayList<>(numberofexchangers);
            for(int i = kints.length-1; i > kints.length-numberofexchangers-1; i--){
                tempkints.add(kints[i]);
            }
            kints = new double[numberofexchangers];
            for(int i = 0; i < numberofexchangers; i++){
                kints[i] = tempkints.get(i);
            }
        }
        else if (kintmethod.equalsIgnoreCase("Random_Sample")){ //Takes a random sample of all inputted kints
            ArrayList<Integer> templocations = sample(numberofexchangers, kints.length);
            ArrayList<Double> tempkints = new ArrayList<>(numberofexchangers);
            for(int i = 0; i < numberofexchangers; i++){
                tempkints.add(kints[templocations.get(i)]);
            }
            kints = new double[numberofexchangers];
            for(int i = 0; i < numberofexchangers; i++){
                kints[i] = tempkints.get(i);
            } 
        }
        //Takes the kints located in the middle of the protein (or at least the middle of the kints as inputted)
        else if (kintmethod.equalsIgnoreCase("Located_Middle")){ 
            int begin = (int) ((kints.length - numberofexchangers)/2);
            ArrayList<Double> tempkints = new ArrayList<>(numberofexchangers);
            for(int i = begin; i < numberofexchangers+begin; i++){
                tempkints.add(kints[i]);
            }
            kints = new double[numberofexchangers];
            for(int i = 0; i < numberofexchangers; i++){
                kints[i] = tempkints.get(i);
            }
        }
        //Takes the kints located at the ends of the proteins (or at least the ends of the kints as inputted)
        else if (kintmethod.equalsIgnoreCase("Located_Ends")){  
            int end = (numberofexchangers/2);
            int begin = kints.length-numberofexchangers/2;
            if(numberofexchangers % 2 != 0){
                begin--;
            }
            ArrayList<Double> tempkints = new ArrayList<>(numberofexchangers);
            for(int i = 0; i < end; i++){
                tempkints.add(kints[i]);
            }
            for(int j = begin; j < kints.length; j++){
                tempkints.add(kints[j]);
            }
            kints = new double[numberofexchangers];
            for(int k = 0; k < numberofexchangers; k++){
                kints[k] = tempkints.get(k);
            }
        }
        else{
            System.out.println("No kint selection method was chosen; please try again");
        }
        if(debug){
            System.out.println("The kints chosen are the following:");
            for(int i = 0; i < kints.length; i++){
                System.out.print(kints[i] + " ");
            }
            if(kints.length != numberofexchangers){
                System.out.println("The wrong number of kints were selected");
            }
        }
    }
    
    //Makes the protection matrix object, which gives information on whether a site is 
    //completely protected, breathing, or completely exposed
    private void initializeProtectionMatrix(String protectionmethod, Scanner inputfile){ 
        protectionmatrix = new double[numberofstates][numberofexchangers];
        if(protectionmethod.equalsIgnoreCase("Manual")){ //Manual specification of the protection scheme
            for(int i = 0; i < numberofstates; i++){
                for(int j = 0; j < numberofexchangers; j++){
                    protectionmatrix[i][j] = inputfile.nextInt();
                }
            }
            if(debug){
                System.out.println("The protection matrix was set manually");
            }
            inputfile.close();
        }
        //By default, native interior is protected, intermediate exterior is protected, and unfolded is random
        else if(protectionmethod.equalsIgnoreCase("Default")){  
            //Selecting r sites from n total to be protected
            int r = numberprotected[0];
            int n = numberofexchangers;            
            /* interior
            int begin = ((n - r)/2);
            int end = ((n + r)/2);
            
            //exterior
            int begin = (r/2);
            int end = (n-r/2);
            */
            int begin = (r/2);
            int end = begin + n - r;
            for(int i = 0; i < n; i++){ //The exterior sites are protected by default in the native state
                if(i >= begin && i < end){
                   protectionmatrix[0][i] = -1; 
                }
                else{
                    protectionmatrix[0][i] = 0;
                } 
            }
            for(int i = 1; i < numberofstates-1; i++){   //In the intermediate states:
                r = numberprotected[i];
                n = numberofexchangers;  
                if(i % 2 != 0){ //If we have one intermediate state, this is the default; interior is protected
                    begin = ((n - r)/2);
                    end = ((n + r)/2);
                    for(int j = 0; j < n; j++){ //The exterior sites are protected by default in the native state
                        if(j >= begin && j < end){
                           protectionmatrix[i][j] = 0; 
                        }
                        else{
                            protectionmatrix[i][j] = 1;
                        } 
                    }
                }
                else{ 
                    //Exterior is protected; this allows intermediate states to alternate whether protected sites are interior/exterior
                    begin = (r/2);
                    end = (n-r/2);
                    for(int j = 0; j < n; j++){ //The exterior sites are protected by default in the native state
                        if(j >= begin && j < end){
                           protectionmatrix[i][j] = 1; 
                        }
                        else{
                            protectionmatrix[i][j] = 0;
                        } 
                    }
                }
            }
            r = numberprotected[numberofstates-1];
            n = numberofexchangers; 
            ArrayList<Integer> templocations = sample(r, n);  
            for(int j = 0; j < n; j++){ //A random sample of all sites are set as fully protected; the rest exchange without breathing
                if(templocations.contains(j)){
                    protectionmatrix[numberofstates-1][j] = 0; 
                }
                else{
                    protectionmatrix[numberofstates-1][j] = 1;
                }
            }
            inputfile.close();
        }
        else{
            System.out.println("The protection matrix could not be specified; please try again");
            inputfile.close();
        }
        for(int i = 0; i < numberofstates; i++){  //Change the -1s to kbreathe for later use
            for(int j = 0; j < numberofexchangers; j++){
                if(protectionmatrix[i][j] == -1){
                    protectionmatrix[i][j] = kbreathe;
                }
            }
        }
        if(debug){
            System.out.println("The protection matrix (adjusted by kbreathe) was set to be the following: ");
            for(int i = 0; i < numberofstates; i++){
                System.out.print("For state " + i + ": ");
                for(int j = 0; j < numberofexchangers; j++){
                    System.out.print(protectionmatrix[i][j] + " ");
                }
                System.out.println();
            }
            for(int i = 0; i < numberofstates; i++){
                int count = 0;
                for(int j = 0; j < numberofexchangers; j++){
                    if(protectionmatrix[i][j] == 0){
                        count++;
                    }
                }
                System.out.println("The number of effectively protected sites in state " + i + " is " + count);
            }
        }
        inputfile.close();
    }
    
    private ArrayList<Integer> sample(int r, ArrayList<Integer> from){  //Generates a random sample from an integer arraylist
        ArrayList<Integer> toreturn = new ArrayList<>();
        ArrayList<Integer> temp = sample(r, from.size());
        for(int i = 0; i < r; i++){
            toreturn.add(from.get(temp.get(i)));
        }
        return toreturn;
    }
    
    private ArrayList<Integer> sample(int r, int n){  //generates a random sample of size r from a pool of size n without replacement
        ArrayList<Integer> toreturn = new ArrayList<>();
        while(toreturn.size() < r){
            int draw = (int) (runif()*n);
            boolean alreadydrawn = false;
            for(int i = 0; i < toreturn.size(); i++){
                if(toreturn.get(i) == draw){
                    alreadydrawn = true;
                }
            }
            if(!alreadydrawn){
                toreturn.add(draw);
            }
        }
        Collections.sort(toreturn);
        return toreturn;
    }

    public double runif(){ //Generates a random draw from a Uniform(0, 1) random variable, complying with global seed
        return randomgenerator.nextDouble();
    }
    
    //Generates a random draw from a Uniform(min, max) random variable, complying with global seed
    public double runif(double min, double max){ 
        double draw = Math.abs(max - min) * runif() + Math.min(min, max);
        return draw;
    }

    private void fixProtection(){ //Checks if the protection values input are allowable
        for(int i = 1; i < numberprotected.length-1; i++){
            if(numberprotected[i] < 0){
                System.out.println("A negative number of protected sites was entered for at least one state.  " + 
                        "Please fix this and try again");
            }
            else if(numberprotected[i] > numberofexchangers){
                System.out.println("A number of protected sites greater than the number of exchanging sites was entered " + 
                        "for at least one state.  Please fix this and try again");
            }
        }
    }
    
    private void importSim(){ //Imports values from a simulation input file
        try{
            File inputfile = new File(originalsimlocation);
            Scanner sc = new Scanner("");
            try{
                sc = new Scanner(inputfile);
            }
            catch(Exception e){
                System.out.println("A scanner could not be attached to the simulation input file");
                sc.close();
            }
                        
            //Line 1: Ignore as of now
            sc.nextLine();
            if(debug){
                System.out.println("Line 1 of the input file was ignored correctly");
            }

            //Line 2: #states and proteins
            numberofstates = sc.nextInt();
            gofcolumn = (int) Math.pow(numberofstates, 2) + 1;
            weightcolumn = gofcolumn+1;
            numberofproteins = sc.nextInt();
            if(debug){
                System.out.println("Read in line 2 as : ");
                System.out.println("Number of states is : " + numberofstates);
                System.out.println("Number of proteins is : " + numberofproteins);
            }

            //Lines 3 through n+2: Rate matrix
            ratematrix = new double[numberofstates][numberofstates];
            for(int i = 0; i < numberofstates; i++){
                for(int j = 0; j < numberofstates; j++){
                    ratematrix[i][j] = sc.nextDouble();
                }
                ratematrix[i][i] = 0;
            }
            if(debug){
                System.out.println("Read in lines 3 through n+2 as : ");
                System.out.println("The protein folding transition matrix : ");
                for(int i = 0; i < numberofstates; i++){
                    for(int j = 0; j < numberofstates; j++){
                        System.out.print(ratematrix[i][j]+ " ");
                    }
                    System.out.println();
                }                
            }

            //Line n + 3: kint list
            ArrayList<Double> tempkints = new ArrayList<>();
            sc.nextLine();

            Scanner kintsc = new Scanner(sc.nextLine());
            
            while(kintsc.hasNextDouble()){
                tempkints.add(kintsc.nextDouble());
            }
            kintsc.close();
            kints = new double[tempkints.size()];
            for(int i = 0; i < kints.length; i++){
                kints[i] = tempkints.get(i);
            }
            
            if(debug){
                System.out.println("Read in line n+3 as : ");
                System.out.println("The kints entered are: ");
                for(int i = 0; i < kints.length; i++){
                    System.out.print(kints[i] + " ");
                }
                System.out.println();
            }

            //Line n + 4: # protected sites in each state
            numberprotected = new int[numberofstates];
            for(int i = 0; i < numberofstates; i++){
                numberprotected[i] = sc.nextInt();
            }
            if(debug){
                System.out.println("Read in line n+4 as : ");
                for(int i = 0; i < numberofstates; i++){
                    System.out.println("The number of protected sites in state " + i + " is : " + numberprotected[i]);
                }
            }

            //Line n + 5: time points of simulation
            ArrayList<Double> temptps = new ArrayList<>();    
            sc.nextLine();
            Scanner tpsc = new Scanner(sc.nextLine());
            
            while(tpsc.hasNextDouble()){
                temptps.add(tpsc.nextDouble());
            }
            
            tpsc.close();
            
            timepoints = new double[temptps.size()];
            for(int i = 0; i < timepoints.length; i++){
                timepoints[i] = temptps.get(i);
            }
            fixTimepoints();
            ntps = timepoints.length;
            if(debug){
                System.out.println("Read in line n+5 as : ");
                for(int i = 0; i < ntps; i++){
                    System.out.print("Time point " + i + " is ");
                    System.out.println(timepoints[i]);
                }
                System.out.println("The total number of time points is " + ntps);
            }

            //Line n + 6: Seed, nex, kbreathe, mass, charge state, width of gaussian, resolution of output
            setSeed(sc.nextDouble());
            numberofexchangers = sc.nextInt();
            if(numberofexchangers > kints.length){
                System.out.println("There are too many exchanging sites for the number of kints you have entered; please try again");
                sc.close();
                return;
            }

            kbreathe = sc.nextDouble();
            initialmass = sc.nextDouble();
            chargestate = sc.nextInt();
            gaussianwidth = sc.nextDouble(); //Akin to the concept of full width at half maximum for experimental HXMS gaussians
            resolution = sc.nextInt();

            if(debug){
                System.out.println("Read in line n+6 as : ");
                System.out.println("The number of exchangers is : " + numberofexchangers);
                System.out.println("The breathing rate constant is : " + kbreathe);
                System.out.println("The initial mass of the protein is : " + initialmass);
                System.out.println("The charge state of the protein is : " + chargestate);
                System.out.println("The gaussian width was set to be : " + gaussianwidth);
                System.out.println("The resolution was set to be : " + resolution);
            }

            //Line n + 7: Method selection
            kintmethod = sc.next();
            chooseKints(kintmethod);

            protectionmethod = sc.next();   

            if(sc.next().charAt(0) == 'H'){
                directionofexchange = 1;
            }
            
            M0 = initialmass + chargestate;
            if(directionofexchange < 0){
                M0 = M0 + numberofexchangers*mdeuterium;  //Change shift location of start if going D to H
            }

            String normalizationmethod = sc.next();
            if(normalizationmethod.equalsIgnoreCase("Height")){
                normalizetoheight = true;
            }else if(normalizationmethod.equalsIgnoreCase("Area")){
                normalizetoheight = false;
            }else{
                if(debug){
                    System.out.println("No normalization method was selected, defaulted to normalization by area");
                }
            }

            if(!sc.hasNextDouble() && sc.hasNext()){
                gofmethod = sc.next();
                if(debug && gofmethod.equalsIgnoreCase("Standard")){
                    System.out.println("The goodness of fit measure was selected as an unweighted time integral");
                }
                else if(debug && gofmethod.equalsIgnoreCase("Inverse_Time")){
                    System.out.println("The goodness of fit measure was selected as a time integral weighted by 1/time");
                }
            }

            if(debug){
                System.out.println("The kint selection method was chosen to be : " + kintmethod);
                System.out.println("The protection method of sites at different folding states was chosen to be : " + protectionmethod);
                if(directionofexchange == 1){
                    System.out.println("The HXMS timecourse is going from H to D");
                }else{
                    System.out.println("The HXMS timecourse is going from D to H");
                }
                System.out.print("The normalization will be conducted to make ");
                if(normalizetoheight){
                    System.out.print("maximum height");
                }else{
                    System.out.print("total area under the curve");
                }
                System.out.println(" for each timepoint equal 1");
                System.out.println("The goodness of fit method chosen was " + gofmethod);
            }

            //Lines n + 8 through 2n + 7 (Optional manual specification of protection matrix)
            initializeProtectionMatrix(protectionmethod, sc);
            sc.close();                      
        }
        catch(Exception except){
            System.out.println("Unable to open file '" + originalsimlocation + "'"); 
        }
        fixProtection();
    }

    //Imports experimental HXMS data, converts to narrow format, and normalizes to total area at each time point = 0
    private void importExpData(){ 
        File inputfile = new File(expdatalocation);
        Scanner sc = new Scanner("");
        try{
            sc = new Scanner(inputfile);
        }
        catch(Exception e){
            System.out.println("Couldn't make scanner object");
        }
        if(sc.hasNextDouble()){ //If the experimental data is alternating rows of m/z and intensity values
            double[][] tempexpdata = new double[ntps*2][resolution];
            for(int i = 0; i < ntps*2; i++){ //double it so that each tp has 2 rows, 1 for m/z, 1 for intensity
                for(int j = 0; j < resolution; j++){
                    tempexpdata[i][j] = sc.nextDouble();
                }
            }
            mzs = tempexpdata[0];
            expdata = new double[resolution * ntps][3];
            for(int i = 0; i < ntps; i++){
                for(int j = 0; j < resolution; j++){
                    expdata[i*resolution+j][0] = timepoints[i]; //time
                    expdata[i*resolution+j][1] = tempexpdata[2*i][j]; //mz
                    expdata[i*resolution+j][2] = tempexpdata[2*i+1][j]; //intensity
                }
            }
        }
        else{ //If the experimental data is in the (preferred) tall format
            expdata = new double[resolution * ntps][3];
            sc.nextLine();
            int nlines = 0;
            while(sc.hasNextLine()){
                nlines++;
                sc.nextLine();
            }
            if(debug){
                System.out.println("The experimental data has " + nlines + " data points");
            }
            if(nlines%resolution != 0){
                System.out.println("The resolution of the experimental data does not match the sim input file");
            }
            try{
                sc = new Scanner(inputfile);
            }
            catch(Exception e){
                System.out.println("Couldn't make scanner object");
            }
            sc.nextLine();
            double[][] tempexpdata = new double[nlines][3];
            for(int i = 0; i < nlines; i++){
                for(int j = 0; j < 3; j++){
                    tempexpdata[i][j] = sc.nextDouble();
                }
            }
            int tempposition = 0;
            for(int i = 0; i < ntps; i++){
                if(debug){
                    System.out.println("The data for time point " + i + " (" + timepoints[i] + ") is being imported, starting at line "+
                            (tempposition + 2) + " of the data file");
                }
                while(Math.abs(tempexpdata[tempposition][0] - timepoints[i]) > timepointdiscrepancy){
                    if(debug){
                        System.out.println("Experimental time " + tempexpdata[tempposition][0] + " was not equal to " + 
                                " sim input time " + timepoints[i]);
                    }
                    tempposition = tempposition + resolution;
                }
                for(int j = 0; j < resolution; j++){
                    for(int k = 0; k < 3; k++){
                        expdata[i*resolution + j][k] = tempexpdata[j+tempposition][k];
                    }
                }
                tempposition = tempposition + resolution;
            }
            mzs = new double[resolution];
            for(int i = 0; i < resolution; i++){
                mzs[i] = expdata[i][1];
            }
        }
        sc.close();

        double tparea = 0;
        double tempdeltamz = mzs[1] - mzs[0];
        for(int i = 0; i < ntps; i++){
            tparea = 0;
            for(int j = 0; j < resolution-1; j++){
                tparea = tparea + tempdeltamz * (expdata[i*resolution + j][2] + expdata[i*resolution + j + 1][2]);
            }
            tparea = tparea/2;
            for(int j = 0; j < resolution; j++){
                expdata[i*resolution + j][2] = expdata[i*resolution + j][2]/tparea;
            }
        }
        if(debug){
            System.out.println("The top of the experimental data imported was: ");
            System.out.println("Time | m/z | Intensity");
            for(int i = 0; i < 5; i++){
                for(int j = 0; j < 3; j++){
                    System.out.print(expdata[i][j] + " ");
                }
                System.out.println();
            }
            System.out.println("The bottom of the experimental data imported was: ");
            System.out.println("Time | m/z | Intensity");
            for(int i = ntps*resolution - 5; i < ntps*resolution; i++){
                for(int j = 0; j < 3; j++){
                    System.out.print(expdata[i][j] + " ");
                }
                System.out.println();
            }
        }
    }

    private void importUncertainties(){ //Imports values from file containing uncertainties and algorithmic parameters
        File inputfile = new File(uncertaintylocation);
        Scanner sc = new Scanner("");
        try{
            sc = new Scanner(inputfile);
        }
        catch(Exception e){
            System.out.println("Couldn't make scanner object");
        }
        //Line 1: notes, ignored
        sc.nextLine();
        if(debug){
            System.out.println("Correctly ignored line 1 of the uncertainties file");
        }
        
        //Lines 2 through n+1: uncertainties in rates between states
        deltaratematrix = new double[numberofstates][numberofstates];
        for(int i = 0; i < numberofstates; i++){
            for(int j = 0; j < numberofstates; j++){
                deltaratematrix[i][j] = sc.nextDouble();
            }
            deltaratematrix[i][i] = 0;
        }
        if(debug){
            System.out.println("Read in lines 2 through n+1 of the uncertainties file as: ");
            System.out.println("The uncertainties in the rate matrix are : ");
            for(int i = 0; i < numberofstates; i++){
                for(int j = 0; j < numberofstates; j++){
                    System.out.print(deltaratematrix[i][j]+ " ");
                }
                System.out.println();
            }
                            
        }
        
        //Line n+2: uncertainties of the number of fully protected sites in each state
        deltaprotected = new double[numberofstates-1];
        for(int i = 0; i < numberofstates-1; i++){
            deltaprotected[i] = sc.nextDouble();
        }
        if(debug){
            System.out.println("Read in line n+2 of the uncertainties file as: ");
            for(int i = 0; i < numberofstates-1; i++){
                System.out.println("The uncertainty in the number of sites protected in state " + i + " is " + deltaprotected[i]);
            }
        }
        
        //Line n+3: uncertainty in kbreathe
        deltakbreathe = sc.nextDouble();
        if(debug){
            System.out.println("Read in line n+3 of the uncertainties file as: ");
            System.out.println("The uncertainty in kbreathe is : " + deltakbreathe);
        }
        
        //Line n+4: optimization method
        optimizationmethod = sc.next();
        
        if(debug){
            System.out.println("Read in line n+4 as the following: ");
            System.out.println("The fitting method chosen is " + optimizationmethod);
        }
        
        //Line n+5: optimization hyperparameters
        if(optimizationmethod.equalsIgnoreCase("Coordinate") || optimizationmethod.equalsIgnoreCase("C")){
            iscoordinate = true;
            maxiterations = sc.nextInt();
            startingchoicesperparameter = sc.nextInt();
            percentchange = sc.nextDouble()/100;
            maxconcurrentthreads = sc.nextInt();
            if(debug){
                System.out.println("Read in line n+5 as the following: ");
                System.out.println("The maximum number of coordinate descent iterations is " + maxiterations);
                System.out.println("The starting number of choices per parameter fit is " + startingchoicesperparameter);
                System.out.println("The proportion change used for coordinate descent is " + percentchange);
            }
        }
        else if(optimizationmethod.equalsIgnoreCase("Nelder_Mead") || optimizationmethod.equalsIgnoreCase("NM")){
            isneldermead = true;
            fittolerance = sc.nextDouble();
            maxiterations = sc.nextInt();
            startingchoicesperparameter = sc.nextInt();
            percentchange = sc.nextDouble()/100;
            maxconcurrentthreads = sc.nextInt();
            if(debug){
                System.out.println("Read in line n+5 as the following: ");
                System.out.println("The tolerance to be used is " + fittolerance);
                System.out.println("The maximum number of Nelder-Mead iterations is " + maxiterations);
                System.out.println("The starting number of choices per parameter fit is " + startingchoicesperparameter);
                System.out.println("The proportion change used for the vertices of the initial simplex is " + percentchange);
            }
        }
        else if(optimizationmethod.equalsIgnoreCase("ABCSMC") || optimizationmethod.equalsIgnoreCase("ABC")){
            isabcsmc = true;
            priortype = sc.next();
            ndist = sc.nextInt() + 1; //add the prior and posterior distributions
            nfinalpoints = sc.nextInt();
            propvariance = sc.nextDouble();
            epsilon0 = sc.nextDouble();
            epsilonT = sc.nextDouble();
            maxconcurrentthreads = sc.nextInt();
            if(debug){
                System.out.println("The prior distribution was chosen to be " + priortype + " for all rate constants");
                System.out.println("Including the prior and posterior distributions, up to " + ndist + " distributions will be used");
                System.out.println("At each step, " + nfinalpoints + " points will be retained");
                System.out.println("The epsilon cutoffs at each step will start with " + epsilon0 + " and end with " + epsilonT);
            }
        }
        else{
            System.out.println("No optimization algorithm was selected. Please enter an optimization method and try again.");
        }
        if(debug){
            System.out.println("The program will use " + maxconcurrentthreads + " concurrent threads");
        }
        
        sc.close();
    }

    private void analyzeUncertainties(){ //Based on the uncertainties entered, sets parameters as fit or not fit
        isfit = new boolean[(int) Math.pow(numberofstates, 2)];
        P0params = new double[2][(int) Math.pow(numberofstates, 2)];
        int position = 0;
        //Set of rate matrix parameters
        for(int i = 0; i < numberofstates; i++){
            for(int j = 0; j < numberofstates; j++){
                if(i != j){
                    if(deltaratematrix[i][j] != 0){
                        isfit[position] = true;
                        deltaratematrix[i][j] = Math.min(Math.abs(deltaratematrix[i][j]), Math.abs(ratematrix[i][j]));
                        P0params[0][position] = Math.max(0, ratematrix[i][j] - deltaratematrix[i][j]);
                        P0params[1][position] = ratematrix[i][j] + deltaratematrix[i][j];
                    }
                    else{
                        isfit[position] = false;
                        P0params[0][position] = ratematrix[i][j];
                        P0params[1][position] = ratematrix[i][j];
                    }
                    position++;
                }
            }
        }

        //Set of protection information
        for(int i = 0; i < numberofstates-1; i++){
            if(deltaprotected[i] != 0){
                isfit[position] = true;
                deltaprotected[i] = Math.min(numberofexchangers-Math.abs(getNumProtected(i)), 
                        Math.min(Math.abs(deltaprotected[i]), Math.abs(getNumProtected(i))));
                P0params[0][position] = getNumProtected(i) - deltaprotected[i];
                P0params[1][position] = getNumProtected(i) + deltaprotected[i];
            }
            else{
                isfit[position] = false;
                P0params[0][position] = getNumProtected(i);
                P0params[1][position] = getNumProtected(i);
            }
            position++;
        }
        
        //Kbreathe
        if(deltakbreathe != 0){
            isfit[position] = true;
            deltakbreathe = Math.min(deltakbreathe, kbreathe);
            P0params[0][position] = Math.max(0, kbreathe - deltakbreathe);
            P0params[1][position] = kbreathe + deltakbreathe;
        }
        else{
            isfit[position] = false;
            P0params[0][position] = kbreathe;
            P0params[1][position] = kbreathe;
        }
        position++;
        calculateNumParameters();
    }
    
    private int getNumProtected(int state){
        int numprotected = 0;
        for(int i = 0; i < protectionmatrix[state].length; i++){
            if(protectionmatrix[state][i] == 0){
                numprotected++;
            }
        }
        return numprotected;
    }
}