package numericalsimulations;

import flanagan.math.*;
import java.io.FileWriter;
import java.io.PrintWriter;

public class NelderMeadRunnable implements Runnable {
    
    //Parameters for debugging/progress
    boolean debug = false; //Whether to print statements for debugging purposes
    boolean progresstracker = true; //Whether to print statements to give progress
    long start; //When this fitting runnable was started
    
    //Parameters for fitting
    ParameterSet current; //The starting ParameterSet object coordinate descent modifies
    int numberofstates; //Number of protein folding states
    int threadnumber; //This runnable's number, for documenation purposes
    double[][] expdata; //Narrow format of experimental data (Time, m/z, Intensity)
    String directory; //Location of folder to which the results will be output
    double[] startingparameters; //Starting parameters fed to the maximization algorithm
    double[] deltaparameters; //Uncertainties in the starting parameters
    boolean[] isfit; //Which of the possible parameters that can be fit, is to be fit, as specified by an uncertainty of 0
    int numparametersfit; //The number of parameters to be fit
    double minrateconstant = Math.pow(10, -20); //Minimum value for rate constants
    int maxiterations; //The maximum number of coordinate descent iterations that can be taken
    double changepercent; //The percent change in each of the parameters fit that coordinate descent tests
    double epsilon; //Convergence tolerance (maximum sd of the values of the function at apices of simplex) for stop
    
    //Parameters for moving data
    PrintWriter writer; //Allows for writing out of nelder mead process 
    
    public NelderMeadRunnable(ParameterSet startingset, double[][] experimentaldata, boolean[] whichparamsfit, int nfit, 
            double percentchange, int maxnmiterations, double ftol, int threadn, String directorypath, boolean todebug){
        start = System.nanoTime(); //Documents the starting time of this runnable
        current = startingset;
        expdata = experimentaldata;
        isfit = whichparamsfit;
        numparametersfit = nfit;
        changepercent = percentchange;
        maxiterations = maxnmiterations;
        epsilon = ftol;
        threadnumber = threadn;
        directory = directorypath;
        debug = todebug;
    }
    
    @Override
    public void run(){ //When the threads in Fitting are "started", this method is run
        try{
            if(debug | progresstracker){
                System.out.println("Thread " + threadnumber + " has been started");
            }
            Maximization maximizer = new Maximization();
            NelderMeadFunction function = new NelderMeadFunction();
            initializeFunction(maximizer, function);
            maximizer.nelderMead(function, startingparameters, deltaparameters, epsilon, maxiterations);
            if(debug | progresstracker){
                System.out.println("Thread " + threadnumber + " has finished (" + 
                        String.format("%.1f", (((double )System.nanoTime() - start)/1000000000)) + " seconds)");
            }
            writer.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.out.println(e);
            System.out.println("There was a problem running thread " + threadnumber);
            
        }
    }
    
    //Sets initial values, initial steps, algorithm constraints, and parameter constraints
    private void initializeFunction(Maximization m, NelderMeadFunction f){ 
        f.setParameterSet(current);
        f.setExperimentalData(expdata);
        f.setThreadNumber(threadnumber);
        f.setFitParameters(isfit);
        String filelocation = directory + "Threadresults" + threadnumber + ".txt";
        try{
            writer = new PrintWriter(new FileWriter(filelocation), true);
            f.makeStorage(writer);
        }
        catch(Exception e){
            System.out.println("Error: unable to make file '" + filelocation + "' in directory '" + directory + "'");
        }
        numberofstates = current.getNStates();
        startingparameters = new double[numparametersfit];
        deltaparameters = new double[numparametersfit];
        int position = 0;
        int positioninfunction = 0;
        //Set of rate matrix parameters
        for(int i = 0; i < numberofstates; i++){
            for(int j = 0; j < numberofstates; j++){
                if(i != j){
                    if(isfit[position]){
                        startingparameters[positioninfunction] = current.getRateMatrix(i, j);
                        m.addConstraint(positioninfunction, -1, minrateconstant); //can't have a rate constant below the set minimum
                        positioninfunction++;
                    }
                    position++;
                }
            }
        }

        //Set of protection information
        for(int i = 0; i < numberofstates-1; i++){
            if(isfit[position]){
                startingparameters[positioninfunction] = current.getNumProtected(i);
                m.addConstraint(positioninfunction, -1, 0); //can't have a negative number of protected sites
                m.addConstraint(positioninfunction, 1, current.getNumExchangers()); //can't have a negative number of protected sites
                positioninfunction++;
            }
            position++;
        }
        
        //Kbreathe
        
        if(isfit[position]){
            startingparameters[positioninfunction] = current.getKbreathe();
            m.addConstraint(positioninfunction, -1, 0); //can't have a negative kbreathe
            positioninfunction++;
        }
        
        for(int i = 0; i < startingparameters.length; i++){
            deltaparameters[i] = Math.max(changepercent*startingparameters[i], minrateconstant);
        }
        
        if(debug){
            System.out.println("The parameter arguments for thread " + threadnumber + " were: ");
            for(int i = 0; i < startingparameters.length; i++){
                System.out.print(startingparameters[i] + " ");
            }
            System.out.println();
            System.out.println("The step arguments for thread " + threadnumber + " were: ");
            for(int i = 0; i < deltaparameters.length; i++){
                System.out.print(deltaparameters[i] + " ");
            }
        }
    }
}
