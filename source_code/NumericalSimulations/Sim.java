package numericalsimulations;

import java.io.*;
import java.util.*;
import flanagan.math.Matrix;
import java.nio.file.FileSystems;
import java.nio.file.Paths;

public class Sim{

//Parameters for debugging/timing 
private long start = System.nanoTime(); //Measures time of creation of Sim Object
private boolean debug = false; //Show selected important internal values and calculations
private boolean verbose = false; //Create R file to graph simulation results
private boolean progresstracker = true; //Shows the progress of just the 8 main steps of the simulation

//Parameters for simulations component
private Random randomgenerator; //Random object used to generate Unif(0, 1), allows for a global seed
private int numberofstates; //Number of protein folding states (ie native + intermediates + fully unfolded)
private int numberofproteins; //Number of proteins in the simulation (see simproteins)
private int numberofexchangers; //Number of sites observed to exchange over the timescale of the experiment
private int[] numberprotected; //Number of sites completely protected at each state
private double[][] ratematrix; //Rate matrix of exchange between folding states
private double[][] protectionmatrix; //Assigns state specific protection values on site-by-site basis
//Condenses information from kints and protection to a matrix of effective kex at each site at each folding state
private double[][] finalmatrix; 
private double[] kints; //List of kints (intrinsic rate constatnts of hydrogen exchange at each of the amide sites)
private double[] cdf; //Cumulative Distribution Function of equilibrium probability distribution, used to sample from equilibrium states
private double[] timepoints; //List of time points at which to simulate
private double timepointdiscrepancy = 0.0005; //max allowed difference between two time points to be considered the same
private int ntps; //Number of time points
private double endtime; //Last timepoint desired from simulation
private double kbreathe; //Breathing rate of the protein
private ArrayList<Protein> simproteins; //List of protein objects to be used in the simulation
private String kintmethod; //Selected method to choose which kints to use in the simulation
private String protectionmethod; //Selected method to choose how to protect sites across sites

//Parameters for convolution to HXMS timecourse
private double[][] exchangetimecourse; //Contains probability distributions for number of amide sites exchanged conditioning on time
private double initialmass; //Standard mass of a protein when all amide sites are occupied with protons
private int chargestate; //Charge state from mass spectrometer
private char isosat; //Direction of exchange as a character
private int directionofexchange = -1; //1 if H to D, -1 otherwise (adding or subtracting mass)
private final double mdeuterium = 1.00627; //Mass difference, in daltons, of a deuterium compared to a proton
private int resolution; //Number of m/z values at which the simulations are evaluated
private double gaussianwidth; //Full Width at Half Max of the reference peak
private double gaussiansd; //Used in default m/z interval calculation
private double gaussianfactor; //Used in convolution of gaussians with simulation output
private double startvalue; //Smallest m/z value calculated in output
private double endvalue; //Largest m/z value calculated in output
private double deltamz; //Constant difference between sequential m/z values
private double[] mzs; //If simulation run to compare to experimental data, the experimental data's m/zs; otherwise, estimated
private double M0; //Used by calculateSpectra function to provide the equivalent of a starting mass
private boolean normalizetoheight; //Normalizes output intensity to either area or height
private double[][] simHXMSds; //Dataset that gets output in narrow format (Time, m/z, Intensity)

//Parameters for comparison with other simulations or experimental data
private double gof = 0; //Goodness of fit measure
private String gofmethod = "Standard"; //Way to calculate goodness of fit
private double[][] expdata; //Narrow format of experimental data (Time, m/z, Intensity)
    
    //Running simulation directly from an input file
    public Sim(String inputlocation, String outputlocation, boolean shoulddebug, boolean shouldverbose){ 
        debug = shoulddebug;
        verbose = shouldverbose;
        setParameters(inputlocation); //Extract values from input file & set seed
        generateMZs(true);
        if(debug || progresstracker){
            System.out.println("The simulation program has completed step 1: Set the parameters from the input file                  " +
                    "               (" + (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
        }
        
        start = System.nanoTime();
        makeProteins(); //Make the proteins in the simulations
        if(debug || progresstracker){
            System.out.println("The simulation program has completed step 2: Make the proteins in the simulation                     " +
                    "               (" + (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
        }
        
        start = System.nanoTime();
        generateStateTransitions(); //Generate history of transitions between states for the proteins
        if(debug || progresstracker){
            System.out.println("The simulation program has completed step 3: Find the global folding events for each protein         " +
                    "               (" + (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
        }
        
        start = System.nanoTime();
        findExchangeTimes(); //Finds the timepoints when each of the exchange sites on each protein irrevocably exchange their hydrogens
        if(debug || progresstracker){
            System.out.println("The simulation program has completed step 4: Find when each exchange site on the protein has exchanged"+
                    "              (" + (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
        }
        
        start = System.nanoTime();
        mzTimepoints(); //From the exchange timepoints, make discrete m/z distribution
        if(debug || progresstracker){
            System.out.println("The simulation program has completed step 5: Calculate distribution of # protons exchanged for " +
                    "specified timepoints (" + (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
        }
        
        start = System.nanoTime();
        mzSmooth(); //Apply smoothing function to convolute to data you might get out of a mass spec
        if(debug || progresstracker){
            System.out.println("The simulation program has completed step 6: Make a continuous distribution from the " +
                    "discrete distribution          (" + (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
        }
        
        start = System.nanoTime();        
        simHXMSds = normalizeTimecourse(simHXMSds); //Normalize m/z distribution to area or height = 1
        if(debug || progresstracker){
            System.out.println("The simulation program has completed step 7: Normalize intensities                                   " +
                    "               (" + (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
        }
        
        start = System.nanoTime();
        output(outputlocation); //Export out m/z distribution at the desired timepoints
        if(verbose){
            printGraphCode(outputlocation);
        }
        if(debug || progresstracker){
            System.out.println("The simulation program has completed step 8: Write out the simulation results in a narrow dataset    " +
                    "               (" + (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
        }
        
        
    }
    
    public Sim(String inputlocation, String expdatalocation, long seed, int nproteins){ //Special constructor to test timing
        setParameters(inputlocation); //Extract values from input file & set seed
        if(!"".equals(expdatalocation)){
            importExpData(expdatalocation);
            generateMZs(false);
        }else{
            generateMZs(true);
        }
        setSeed(seed);
        numberofproteins = nproteins;
        makeProteins(); //Make the proteins in the simulations
        generateStateTransitions(); //Generate history of transitions between states for the proteins
        findExchangeTimes(); //Finds the timepoints when each of the exchange sites on each protein irrevocably exchange their hydrogens
        mzTimepoints(); //From the exchange timepoints, make discrete m/z distribution
        mzSmooth(); //Apply smoothing function to convolute to data you might get out of a mass spec      
        simHXMSds = normalizeTimecourse(simHXMSds); //Normalize m/z distribution to area or height = 1
        if(!"".equals(expdatalocation)){
            GOF(expdata);
        }
    }
    
    public Sim(ParameterSet setofparameters){ //Running simulation from a ParameterSet object (used for fittings)
        setParameters(setofparameters);
        makeProteins();
        generateStateTransitions();
        findExchangeTimes();
        mzTimepoints();
        mzSmooth();
        simHXMSds = normalizeTimecourse(simHXMSds);
    }
    
    //Running simulation with additional experimental data, uses experimental m/z and reports GoF automatically
    public Sim(String inputlocation, String outputlocation, String expdatalocation, boolean shoulddebug, boolean shouldverbose){ 
        debug = shoulddebug;
        verbose = shouldverbose;
        setParameters(inputlocation); //Extract values from input file & set seed
        if(debug || progresstracker){
            System.out.println("The simulation program has completed step 1:  Set the parameters from the input file                 " +
                    "                (" + (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
        }
        
        start = System.nanoTime();
        importExpData(expdatalocation);
        if(debug || progresstracker){
            System.out.println("The simulation program has completed step 2:  Read in experimental data to compare to the simulation " +
                    "                (" + (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
        }
        
        start = System.nanoTime();
        generateMZs(false);
        makeProteins(); //Make the proteins in the simulations
        if(debug || progresstracker){
            System.out.println("The simulation program has completed step 3:  Make the proteins in the simulation                    " +
                    "                (" + (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
        }
        
        start = System.nanoTime();
        generateStateTransitions(); //Generate history of transitions between states for the proteins
        if(debug || progresstracker){
            System.out.println("The simulation program has completed step 4:  Find the global folding events for each protein        " +
                    "                (" + (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
        }
        
        start = System.nanoTime();
        findExchangeTimes(); //Finds the timepoints when each of the exchange sites on each protein irrevocably exchange their hydrogens
        if(debug || progresstracker){
            System.out.println("The simulation program has completed step 5:  Find when each exchange site on the protein " +
                    "has exchanged              (" + (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
        }
        
        start = System.nanoTime();
        mzTimepoints(); //From the exchange timepoints, make discrete m/z distribution
        if(debug || progresstracker){
            System.out.println("The simulation program has completed step 6:  Calculate distribution of # protons exchanged for " +
                    "specified timepoints (" + (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
        }
        
        start = System.nanoTime();
        mzSmooth(); //Apply smoothing function to convolute to data you might get out of a mass spec
        if(debug || progresstracker){
            System.out.println("The simulation program has completed step 7:  Make a continuous distribution from the " +
                    "discrete distribution          (" + (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
        }  
        
        start = System.nanoTime();
        simHXMSds = normalizeTimecourse(simHXMSds); //Normalize m/z distribution to area or height = 1
        if(debug || progresstracker){
            System.out.println("The simulation program has completed step 8:  Normalize intensities                                  " +
                    "                (" + (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
        }
        
        start = System.nanoTime();
        output(outputlocation); //Export out m/z distribution at the desired timepoints
        if(verbose){
            printGraphCode(outputlocation, expdatalocation);
        }
        if(debug || progresstracker){
            System.out.println("The simulation program has completed step 9:  Write out the simulation results in a narrow dataset   " +
                    "                (" + (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
        }

        start = System.nanoTime();
        if(debug || progresstracker){
            System.out.println("The simulation program has completed step 10: Calculate goodness of fit value of " + 
                    String.format("%.3f", GOF(expdata)) + " with the experimental data    (" + 
                    (int) (((double )System.nanoTime() - start)/1000000) + " ms)");
        }
    }
    
    private void setParameters(String inputlocation){ //Read in simulation parameters from input file
        try{
            File inputfile = new File(inputlocation);
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
            setEquilibriumDist();  //Using the rate matrix, find folding state probability distribution

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
            endtime = timepoints[ntps-1];
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
            gaussiansd = gaussianwidth/chargestate; //Standard deviation of the gaussians used to convolute the simulation results
            double gaussianvariance = Math.pow(gaussiansd, 2);
            gaussianfactor = (-1 / (2 * gaussianvariance));
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

            isosat = sc.next().charAt(0);
            if(isosat == 'H'){
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
                if(isosat == 'H'){
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
            System.out.println("Unable to open file '" + inputlocation + "'"); 
        }
        fixProtection();
        initializeFinalMatrix();
    }
    
    private void setParameters(ParameterSet setofparameters){ //Set simulation parameters from ParameterSet object
        debug = false;
        progresstracker = false;
        numberofproteins = setofparameters.getNumberProteins();
        ratematrix = setofparameters.getRateMatrix();
        numberofstates = ratematrix.length;
        setEquilibriumDist();
        timepoints = setofparameters.getTimePoints();
        ntps = timepoints.length;
        endtime = timepoints[ntps-1];
        finalmatrix = setofparameters.getFinalMatrix();
        randomgenerator = setofparameters.getRandomGenerator();      
        numberofexchangers = finalmatrix[0].length;
        initialmass = setofparameters.getMass();
        chargestate = setofparameters.getChargeState();
        gaussianwidth = setofparameters.getGaussianWidth();
        gaussiansd = gaussianwidth/chargestate; //Standard deviation of the gaussians used to convolute the simulation results
        double gaussianvariance = Math.pow(gaussiansd, 2);
        gaussianfactor = (-1 / (2 * gaussianvariance));
        resolution = setofparameters.getResolution();
        gofmethod = setofparameters.getGoFMethod();
        M0 = initialmass + chargestate;
        directionofexchange = (int) setofparameters.getDirectionofExchange();
        if(directionofexchange < 0){
            M0 = M0 + numberofexchangers*mdeuterium;  //Change shift location of start if going D to H
        }       
        mzs = setofparameters.getMZ();
        deltamz = mzs[1] - mzs[0];
        normalizetoheight = false;
    }
    
    private void generateMZs(boolean mztype){ //generate m/z values from guesses
        mzs = new double[resolution];
        if(mztype){ //if the m/zs are to be automatically generated
            //(Fully hydrogen mass + charge state protons + ~qnorm(0.005)*gaussian sd) divided by charge = m/z
            startvalue = (initialmass + chargestate -16 * gaussiansd)/chargestate; 
            //Similar thing, but adding the mass of the additional neutrons before dividing by m/z
            endvalue = (initialmass + chargestate + numberofexchangers * mdeuterium + 16 * gaussiansd)/chargestate; 
            deltamz = (endvalue - startvalue)/resolution;
            for(int i = 0; i < resolution; i++){
                mzs[i] = startvalue + i*deltamz;
            }
        }
        else{ //if the m/zs are to be used from the experimental data
            for(int i = 0; i < resolution; i++){
                mzs[i] = expdata[i][1];
            }
            deltamz = mzs[1] - mzs[0];
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
        //Takes the kints located at the ends of the proteins (or at least the ends of the kints as input)
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
    
    private void initializeFinalMatrix(){ 
        //Using the kint and protection information, make a new matrix storing all this information in a compact format
        finalmatrix = new double[numberofstates][numberofexchangers];
        for(int i = 0; i < numberofstates; i++){
            for(int j = 0; j < numberofexchangers; j++){
                finalmatrix[i][j] = protectionmatrix[i][j]*kints[j];
            }
        }
        if(debug){
            System.out.println("The final information matrix was set to be the following: ");
            for(int i = 0; i < numberofstates; i++){
                System.out.print("For state " + i + ": ");
                for(int j = 0; j < numberofexchangers; j++){
                    System.out.print(finalmatrix[i][j] + " ");
                }
                System.out.println();
            }
        }
    }
    
    private void makeProteins(){ //Makes all the Protein objects in the simulation
        simproteins = new ArrayList<>(numberofproteins);
        for(int i = 0; i < numberofproteins; i++){
            simproteins.add(new Protein(finalmatrix, randomgenerator));
            //sampleFromEquilibrium() samples a protein's initial state from equilibrium state distribution
            //- not correct for intrinsically disordered proteins like alpha lytic protease
            //Set protein at time = 0 to a randomly assigned state from the equilibrium distribution
            simproteins.get(i).setStateChange(0, sampleFromEquilibrium()); 
        }
    }
    
    private void generateStateTransitions(){ //Find the folding state transitions for each protein in the simulation
        for(int i = 0; i < numberofproteins; i++){ //Going through each protein
            double time = 0; //starting at time 0
            while(time < endtime){ //During the timecourse, while the protein still has sites to exchange.  
                time = time + generateNewState(time, i);
            }
        }
    }
    
    private void findExchangeTimes(){ //Finds the HX exchange times for all sites on every protein in the simulation
        for(int i = 0; i < numberofproteins; i++){
            simproteins.get(i).findExchangers(timepoints);
        }
        if(debug){
            System.out.println("This is a sample protein's state history : ");
            int numtransitions = simproteins.get(0).getStateChangeTimes().size();
            System.out.println("Time | New State | # sites exchanged");
            for(int i = 0; i < numtransitions; i++){
                System.out.println(simproteins.get(0).getStateChangeTimes().get(i) + " " + 
                        simproteins.get(0).getStateHistory().get(i) + " " + 
                        simproteins.get(0).getNumberSitesExchanged(simproteins.get(0).getStateChangeTimes().get(i)));
            }
        }
    }
    
    //Finds the marginal distribution of the number of sites exchanged conditioned at each of the timepoints
    private void mzTimepoints(){ 
        exchangetimecourse = new double[ntps][numberofexchangers+1];
        for(int i = 0; i < ntps; i++){
            for(int j = 0; j < numberofproteins; j++){
                //If the protein has k number of sites exchanged at time i, 
                //increase the number of proteins in that exchanged number "bin" by 1
                exchangetimecourse[i][simproteins.get(j).getNumberSitesExchanged(timepoints[i])]++; 
            }
        }
        for(int i = 0; i < ntps; i++){
            for(int j = 0; j <= numberofexchangers; j++){
                //Normalize distribution of # of sites exchanged to be a valid probability distribution
                exchangetimecourse[i][j] = exchangetimecourse[i][j]/numberofproteins;  
            }
        }
        if(debug){
            System.out.println("The conditional probability mass function for the number of sites " + 
                    "exchanged conditioning on time is : ");
            System.out.println("Time | # sites exchanged | p(x)");
            for(int i = 0; i < ntps; i++){
                for(int j = 0; j <= numberofexchangers; j++){
                    System.out.println(timepoints[i] + " " + j + " " + exchangetimecourse[i][j]);
                }
            }
        }
    }
    
    private void mzSmooth(){ //Simulate the spectra for each m/z value and timepoint
        //Narrow data format for ease of analysis.  The column variables are Time, M/Z, and Intensity, respectively
        simHXMSds = new double[resolution*ntps][3]; 
        //Make Time column
        for(int i = 0; i < ntps; i++){
            for(int j = i*resolution; j < (i+1)*resolution; j++){
                simHXMSds[j][0] = timepoints[i];
            }
        }

        //Make m/z value column
        for(int i = 0; i < ntps; i++){
            for(int j = 0; j < resolution; j++){
                simHXMSds[i*resolution+j][1] = mzs[j];
            }
        }
                
        //Make intensity column   
        if(debug){
            System.out.println("The gaussian factor used to smooth the data is: " + gaussianfactor);
        }
        for(int i = 0; i < ntps; i++){
            double [] tpproportions = new double[numberofexchangers+1];
            for(int j = 0; j < numberofexchangers + 1; j++){
                tpproportions[j] = exchangetimecourse[i][j];
            }
            for(int k = i*resolution; k < (i+1)*resolution; k++){
                simHXMSds[k][2] = computeIntensity(tpproportions, simHXMSds[k][1]);
            }
        }
    }
    
    //Given an array indicating proportion of proteins with i number of sites exchanged and an m/z value, return the expected intensity
    private double computeIntensity(double[] distribution, double mz){ 
        double sum = 0;
        for(int i = 0; i <= numberofexchangers; i++){
            sum = sum + distribution[i] * evaluateNormal(mz, (M0 + directionofexchange * i * mdeuterium)/chargestate);
            //sum = sum + distribution[i] * Math.exp((-1) * 0.01*gaussianwidth *
            //Math.pow((mz - ((M0 + directionofexchange * i * mdeuterium)/chargestate)), 2)); 
            //John Strahan's function for convolution with gaussians
        }
        double finalsum = Math.sqrt(gaussianwidth/Math.PI)*sum; //Normalization
        return finalsum;
    }
    
    //Calculates the (unnormalized) height of the normal probability density function of Normal(mu, sqrt(2*gaussianfactor))
    private double evaluateNormal(double x, double mu){ 
        //Note that the variance of the normal distribution is gaussianvariance, calculated earlier
        return Math.exp(gaussianfactor * Math.pow(x - mu, 2));
    }
    
    //Normalizes intensities to either height = 1 or area = 1 at each timepoint
    private double[][] normalizeTimecourse(double[][] timecourse){  
        double[][] normalized = new double[timecourse.length][timecourse[0].length];
        
        if(normalizetoheight){
            for(int i = 0; i < ntps; i++){
                double maxintensity = 0;
                for(int j = i*resolution; j < (i+1)*resolution; j++){ //Find max intensity for a particular timepoint
                    if(timecourse[j][2] > maxintensity){
                        maxintensity = timecourse[j][2];
                    }             
                }
                for(int j = i*resolution; j < (i+1)*resolution; j++){ //Normalize intensities to have height 1
                    normalized[j][2] = timecourse[j][2]/maxintensity;              
                }
            }
        }else{
            double sum = 0;
            for(int i = 0; i < ntps; i++){
                sum = 0;
                for(int j = i*resolution; j < (i+1)*resolution-1; j++){ //Find total intensity for a particular timepoint
                    sum = sum + deltamz*(timecourse[j][2]+timecourse[j+1][2]);            
                }
                sum = sum/2;                             
                for(int j = i*resolution; j < (i+1)*resolution; j++){ //Normalize intensities to have area 1
                    normalized[j][2] = timecourse[j][2]/sum;              
                }
            }
        }
        //Manual deep copy of time and m/z columns
        for(int i = 0; i < ntps; i++){
            for(int j = i*resolution; j < (i+1)*resolution; j++){      
                normalized[j][0] = timecourse[j][0];   
                normalized[j][1] = timecourse[j][1];   
            }
        }
        /*
        System.out.println();
        for (int x = 0; x < normalized.length; x++){
            for (int y = 0; y < normalized[0].length; y++){
                System.out.print(normalized[x][y] + " ");
            }
            System.out.println();
        }
        */
        return normalized;
    }
    
    //Write out dataset (simHXMSds matrix) containing time, m/z, and normalized intensity column titles
    private void output(String outputlocation){  
        if(debug){
            System.out.println("The top of the data frame to be output is: ");
            System.out.println("Time | m/z | Intensity");
            for(int i = 0; i < 5; i++){
                for(int j = 0; j < 3; j++){
                    System.out.print(simHXMSds[i][j]+ " ");
                }
                System.out.println();
            }
        }
        try{
            PrintWriter writer = new PrintWriter(new FileWriter(outputlocation), true);
            writer.println("Time mz Intensity");
            for(int i = 0; i < resolution*ntps; i++){
                for(int j = 0; j < 2; j++){
                    writer.print(simHXMSds[i][j]+" ");
                }
                writer.println(simHXMSds[i][2]);
            }
            writer.close();
        }
        catch(Exception e) {
            System.out.println("Error writing out data");
        }
    }
    
    private double generateNewState(double currenttime, int currentprotein){ //Calculates waiting time and new state of given protein
        //See "Introduction to Stochastic Processes with R", by Robert P. Dobrow, page 270 for the best explanation
        //Essentially, there are two "clocks" for leaving the current state for non-end states 
        //(only one for native and fully unfolded states if all intermediates are on-pathway)
        //Whichever clock "goes off" first dictates the new folding state
        //The times of the alarm clocks are random draws from independent exponential distributions specified by the 2 
        //"off" rate constants (1 in the case of native or fully unfolded states)
        //Conceptually, the higher the rate constant, the faster the protein's folding state is likely to change, 
        //and the more likely it changes to that state 
        int currentstate = simproteins.get(currentprotein).getState();
        int newstate = -1; //Set value of newstate to an impossible value
        double tau = -1; //Make the time an impossible value
        try{
            for(int i = 0; i < numberofstates; i++){
                if(i != currentstate){ //We're transitioning to an accessible state other than the one we are currently in
                    if(currentstate < 0){
                        System.out.println("Aha got you! The cdf was ");
                        for(int j = 0; j < cdf.length; j++){
                            System.out.print(cdf[j] + " ");
                        }
                        System.out.println();
                    }
                    double temptau = rexp(ratematrix[currentstate][i]); //Generate random value from exponential distribution
                    if(temptau < tau || tau < 0){ //We want the state that we most quickly transition to, and the time it takes to do so
                        tau = temptau;
                        newstate = i;
                    }
                }
            }        
            //Record the details of transition between folding states
            simproteins.get(currentprotein).setStateChange(currenttime+tau, newstate); 
            return tau; 
        }
        catch(Exception e){
            e.printStackTrace();
            System.out.println("There was an error progressing forward to the next folding state using this rate matrix:");
            for(int i = 0; i < numberofstates; i++){
                for(int j = 0; j < numberofstates; j++){
                    System.out.print(ratematrix[i][j] + " ");
                }
                System.out.println();
            }
            System.out.println("The proposed new state was state " + newstate + " with switching time " + tau);
            System.exit(0);
            return(-1);
        }
    }
     
    public void setEquilibriumDist(){ //Use matrix of rate constants to make a discrete probability distribution from which to draw
        try{
            double[][] tempratematrix = new double[numberofstates+1][numberofstates];
            for(int i = 0; i < numberofstates; i++){
                double sum = 0;
                for(int j = 0; j < numberofstates; j++){
                    tempratematrix[i][j] = ratematrix[j][i]; //transpose the input rate matrix
                    sum = sum + ratematrix[i][j];
                }
                tempratematrix[numberofstates][i]=1; //Set the bottom row of the matrix to all have value 1
                tempratematrix[i][i] = -sum;
            }
            Matrix m = new Matrix(tempratematrix);
            double[] pmf = new double[numberofstates+1];
            pmf[numberofstates]=1;
            //This is equivalent to finding the eigenvector with eigenvalue 0.  
            //Due to mathemagic reasons, this eigenvector gives the normalized equilibrium population proportions
            pmf = m.solveLinearSet(pmf); 
            cdf = new double[numberofstates+1]; //Make cdf array
            for(int i = 1; i <= numberofstates; i++){
                double cutoff = 0;
                for(int j = 0; j < i; j++){
                    cutoff = cutoff+pmf[j];
                }
                cdf[i] = cutoff;
            }
            if(debug){
                System.out.println("The equilibrium proportions of the folding states was calculated to be: ");
                for(int i = 0; i < pmf.length; i++){
                    System.out.print(pmf[i] + " ");
                }
                System.out.println();
                System.out.println("The cumulative distribution function for the equilibrium proportions is: ");
                for(int i = 0; i < cdf.length; i++){
                    System.out.print(cdf[i] + " ");
                }
                System.out.println();
            }
        }
        catch(Exception e){
            cdf = new double[numberofstates+1];
            for(int i = 0; i < cdf.length; i++){
                cdf[i] = 1.0;
            }
            if(debug){
                System.out.println(e);
                e.printStackTrace();
                System.out.println("The equilibrium probability distribution could not be calculated (singular generator matrix)");
                System.out.println("Instead, the proteins were all initialized as native state with the following cdf:");
            }
            for(int i = 0; i < cdf.length; i++){
                System.out.print(cdf[i] + " ");
            }
            System.out.println();
        }
    }
    
    public int sampleFromEquilibrium(){  //Generate a random state given the equilibrium probability distribution
        try{
            double draw = runif();
            for(int i = 1; i <= numberofstates; i++){
                if(draw < cdf[i]){
                    return i-1;
                }
            }
        }
        catch(Exception e){
            //System.out.println("The inverse cdf method was conducted incorrectly");
        }
        return -1;
    }
    
    public double[][] getResults(){ //Extract the ds with time, m/z, and intensity, without column names
        return simHXMSds;
    }
    
    //Imports experimental HXMS data, converts to narrow format, and normalizes to total area at each time point = 0
    private void importExpData(String expdatalocation){ 
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
                    System.out.println("The data for time point " + i + " (" + timepoints[i] + 
                            ") is being imported, starting at line " + (tempposition + 2) + " of the data file");
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
    
    //Calculates and returns the GOF measure of this simulation with another simulation
    public double getGOF(double[][] experimentaldata){  
        //System.out.println("IN GETGOF METHOD");
        if(gof != 0){
            return gof;
        }
        else if(experimentaldata.length != resolution*ntps){
            System.out.println("The simulated data cannot be compared to the experimental data");
            return 1;
        }else{
            return GOF(experimentaldata);
        }
    }
    
    //If the simulation was created with an experimental data input
    public double getGOF(){ 
        return(gof);
    }

    private double GOF(double[][] experimentaldata){
        //Need to normalize simulation to area = 1 for goodness of fit calculation
        normalizetoheight = false;
        simHXMSds = normalizeTimecourse(simHXMSds);
        /*
        //Deviance method of calculating GOF
        double deviances[] = new double[ntps];
        double areatoadd = 0;
        for(int i = 0; i < ntps; i++){
            areatoadd = 0;
            for(int j = 0; j < resolution-1; j++){
                areatoadd = deltamz * experimentaldata[i*resolution + j][2] * 
        Math.log(experimentaldata[i*resolution + j][2]/simHXMSds[i*resolution + j][2]); //Integration wrt dmz
                if(!Double.isNaN(areatoadd)){
                    deviances[i] = deviances[i] + areatoadd;
                }
            }
            deviances[i] = deviances[i]/2;
            if(debug){
                double simsum = 0;
                double expsum = 0;
                for(int j = i*resolution; j < (i+1)*resolution-1; j++){ //Find total intensity for a particular timepoint
                    simsum = simsum + deltamz*(simHXMSds[j][2]+simHXMSds[j+1][2]);            
                    expsum = expsum + deltamz*(experimentaldata[j][2]+experimentaldata[j+1][2]);     
                }
                simsum = simsum/2;
                expsum = expsum/2;
                System.out.println("The area under the simulation curve at time point " + timepoints[i] + " was found to be " + simsum);
                System.out.println("The area under the experimental curve at time point " + timepoints[i] + 
        " was found to be " + expsum);
                System.out.println("The total deviance at time " + timepoints[i] + " was found to be " + deviances[i]);
            }
        }
        double deviance = 0;
        for(int i = 0; i < ntps-1; i++){
            deviance = deviance + (timepoints[i+1] - timepoints[i])*(deviances[i]+deviances[i+1]); //dt
        }
        deviance = (-1)*deviance/2;
        gof = deviance;
        return gof;
        */
        
        double avg = 0;
        double[] squerror = new double[ntps];
        double[] avgerror = new double[ntps];
        for(int i = 0; i < ntps; i++){ //For each time point
            avg = 2/(mzs[0] + mzs[resolution-1]); //find average height
            //System.out.println("finding average height: " + avg);
            for(int k = 0; k < resolution-1; k++){
                deltamz = experimentaldata[i*resolution+k+1][1] - experimentaldata[i*resolution+k][1];
                squerror[i] = squerror[i] + deltamz*(Math.pow((experimentaldata[i*resolution+k][2] - simHXMSds[i*resolution+k][2]), 2) +
                        Math.pow((experimentaldata[i*resolution+k+1][2] - simHXMSds[i*resolution+k+1][2]), 2));
                avgerror[i] = avgerror[i] + deltamz*(Math.pow((experimentaldata[i*resolution+k][2] - avg), 2) + 
                        Math.pow((experimentaldata[i*resolution+k+1][2] - avg), 2));
            }
        }
        double totalerror = 0;
        double totalavgerror = 0;
        if(gofmethod.equalsIgnoreCase("Standard")){
            for(int i = 0; i < ntps-1; i++){
                double a = timepoints[i];
                double b = timepoints[i+1];
                totalerror = totalerror + (b - a)*(squerror[i] + squerror[i+1]);
                totalavgerror = totalavgerror + (b - a)*(avgerror[i] + avgerror[i+1]);
            }
        }else{
            for(int i = 0; i < ntps-1; i++){
                double a = timepoints[i];
                double b = timepoints[i+1];
                totalerror = totalerror + Math.log(b/a)*(squerror[i] + squerror[i+1]);
                totalavgerror = totalavgerror + Math.log(b/a)*(avgerror[i] + avgerror[i+1]);
            }
        }
        if(debug){
            System.out.println("The total error is : " + totalerror);
            System.out.println("The total avg error is : " + totalavgerror);
        }
        gof = (1-totalerror/totalavgerror);
        return gof;
    }
        
    private ArrayList<Integer> sample(int r, int n){  //generates a random sample of size r from a pool of size n without replacement
        ArrayList<Integer> toreturn = new ArrayList<>(r);
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
    
    private void setSeed(double seed){ //Sets seed for the simulation
        randomgenerator = new Random();
        randomgenerator.setSeed((long) seed);
    }
    
    public double runif(){ //Generates random draw from Uniform(0, 1) random variable, using a global seed
        return randomgenerator.nextDouble();
    }
    
    public double rexp(double rate) {  //Generate a sample of size 1 from an exponential distribution with the specified rate parameter
        if(rate != 0){
            return -(Math.log(runif()) / rate);
        }else{
            return(endtime+1);
        }
    }
    
    private void fixTimepoints(){  //Arrange the timepoints in order, give error if 
        Arrays.sort(timepoints);
        if(timepoints[0] <= 0){
            System.out.println("At least one of the time points entered is too small (<= 0).  Please fix this and try again");
        }
    }
    
    private void fixProtection(){ //Checks if the protection values input are allowable
        for(int i = 1; i < numberprotected.length-1; i++){
            if(numberprotected[i] < 0){
                System.out.println("A negative number of protected sites was entered for at least one state. " + 
                        "Please fix this and try again");
            }
            else if(numberprotected[i] > numberofexchangers){
                System.out.println("A number of protected sites greater than the number of exchanging sites was entered " + 
                        "for at least one state.  Please fix this and try again");
            }
        }
    }
    
    private void printGraphCode(String outputlocation){ //Make an R file to graph the results
        String rlocation;
        if(Paths.get(outputlocation).getParent() != null){
            rlocation = Paths.get(outputlocation).getParent() + FileSystems.getDefault().getSeparator() + "graphingsimcode.R";
        }
        else{
            rlocation = "graphingsimcode.R";
        }
        String pnglocation = "graphsim.png";
        PrintWriter graphingcodewriter;
        try{
            graphingcodewriter = new PrintWriter(new FileWriter(rlocation), true);
            graphingcodewriter.println("# Code to graph " + Paths.get(outputlocation).getFileName() + " \n");
            graphingcodewriter.println("library(tidyverse) \n"
                                    + "library(grDevices) \n"
                                    + "simdata <- read.table('" + Paths.get(outputlocation).getFileName() + "', header = TRUE)%>% \n"
                                    + "  mutate(Time = as.factor(Time))\n"
                                    + "simplot <- ggplot(simdata, aes(x = mz, y = Intensity, color = Time)) + \n"
                                    + "  geom_line(size = 3) + \n"
                                    + "  facet_wrap(~ Time) + \n"
                                    + "  theme_bw() + \n"
                                    + "  theme(axis.title.y = element_blank(), \n"
                                    + "        axis.ticks.y = element_blank(), \n"
                                    + "        axis.text.y = element_blank(), \n"
                                    + "        axis.text.x = element_text(angle = 90), \n"
                                    + "        panel.grid = element_blank(), \n"
                                    + "        legend.position = 'none', \n"
                                    + "        text = element_text(size = 60) \n"
                                    + "  ) + \n"
                                    + "  xlab('m/z') \n"
                                    + "png('" + pnglocation + "', width = 3840, height = 2160) \n"
                                    + "simplot \n"
                                    + "dev.off() ");
            graphingcodewriter.close();
        }catch(Exception e){
            System.out.println("There was an error making an R file to graph this simulation");
        }
    }
    
    private void printGraphCode(String outputlocation, String expdatalocation){ //Make an R file to graph the results
        String rlocation;
        if(Paths.get(outputlocation).getParent() != null){
            rlocation = Paths.get(outputlocation).getParent() + FileSystems.getDefault().getSeparator() + "graphingsimcode.R";
        }
        else{
            rlocation = "graphingsimcode.R";
        }
        String pnglocation = "graphsim.png";
        PrintWriter graphingcodewriter;
        try{
            graphingcodewriter = new PrintWriter(new FileWriter(rlocation), true);
            graphingcodewriter.println("# Code to graph " + Paths.get(outputlocation).getFileName() + " \n");
            graphingcodewriter.println("library(tidyverse) \n"
                                    + "library(grDevices) \n"
                                    + "simdata <- read.table('" + Paths.get(outputlocation).getFileName() + "', header = TRUE)%>% \n"
                                    + "  mutate(Source = 'Simulation') \n"
                                    + "expdata <- read.table('" + Paths.get(expdatalocation).getFileName() + "', header = TRUE)%>% \n"
                                    + "  mutate(Source = 'Experiment') \n"
                                    + "coldata <- rbind(simdata, expdata)%>% \n"
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
            System.out.println("There was an error making an R file to graph this simulation");
        }
    }
}