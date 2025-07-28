package numericalsimulations;

import java.io.FileWriter;
import java.io.PrintWriter;

public class CoordinateDescentRunnable implements Runnable {

    //Parameters for debugging/progress
    boolean debug = false; //Whether to print statements for debugging purposes
    boolean progresstracker = true; //Whether to print statements to give progress
    long start; //When this fitting runnable was started

    //Parameters for fitting
    ParameterSet current; //The starting ParameterSet object coordinate descent modifies
    int numberofstates; //Number of protein folding states
    int threadnumber; //This runnable's number, for documenation purposes
    int maxiterations; //The maximum number of coordinate descent iterations that can be taken
    double changepercent; //The percent change in each of the parameters fit that coordinate descent tests
    String directory; //Location of folder to which the results will be output
    boolean[] isfit; //Which of the possible parameters that can be fit, is to be fit, as specified by an uncertainty of 0
    int numparametersfit; //Count of the paramters that are fit, as indicated in 'isfit'
    double[][] expdata; //Narrow format of experimental data (Time, m/z, Intensity)
    int iteration = 1; //Current iteration of coordinate descent
    int stabilitycounter = 0; //Allows coordinate descent to end earlier if no improvements are being made (local max of gof)

    //Parameters for moving data
    PrintWriter writer; //Allows for writing out of coordinate descent process    
    
    public CoordinateDescentRunnable(ParameterSet startingset, double[][] experimentaldata, boolean[] whichparamsfit, int threadn, 
            int maxcditerations, double percentchange, String filedirectory){ 
        start = System.nanoTime(); //Documents the starting time of this runnable
        current = startingset.deepCopy(); //Initialize starting parameters
        threadnumber = threadn; //Set thread number
        maxiterations = maxcditerations; //Set max number of Coordinate Descent iterations
        changepercent = percentchange; //The proportion by which each parameter is changed in coordinate descent
        directory = filedirectory; //Set output directory
        isfit = whichparamsfit; //Sets which parameters are fit in coordinate descent
        expdata = experimentaldata; //Save input experimental HXMS time course   
        if(debug){
            if(experimentaldata == null){
                System.out.println("The experimental data passed in to Thread " + threadnumber + " was null");
            }
            if(expdata == null){
                System.out.println("Thread " + threadnumber + "'s experimental data is null");
            }
            System.out.println("Thread " + threadnumber + " has been constructed");
        }
    }

    @Override
    public void run(){ //When the threads in Fitting are "started", this method is run
        try{
            if(debug | progresstracker){
                System.out.println("Thread " + threadnumber + " has been started");
            }
            numberofstates = current.getNStates(); 
            calculateNumParameters();  
            makeStorage();
            runCoordinateDescent();
            writer.close();
            if(debug | progresstracker){
                System.out.println("Thread " + threadnumber + " has finished (" + 
                        String.format("%.1f", (((double )System.nanoTime() - start)/1000000000)) + " seconds)");
            }
        }
        catch(Exception e){
            e.printStackTrace();
            System.out.println(e);
            System.out.println("There was a problem running thread " + threadnumber);
            
        }
    }
    
    //Writes a data file containing the information for the current threadnumber, the off-diagonal rate constants, 
    //the number of sites fully protected in not fully unfolded states, kbreathe, and goodness of fit
    private void makeStorage(){ 
        try{
            //Makes txt file, adds variable headings
            writer = new PrintWriter(new FileWriter(directory + "Threadresults" + threadnumber + ".txt"), true); 
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
        }
        catch(Exception e) {
            System.out.println("Error making storage file");
        }
    }
    
    //Called whenever a new simulation is run to store all values for that simulation, including goodness of fit
    private void writeSimToStorage(ParameterSet testedsim, double gof){ 
        writer.print(threadnumber + " ");
        double[][] tempratematrix = testedsim.getRateMatrix();
        //Set of rate matrix parameters
        for(int k = 0; k < numberofstates; k++){
            for(int j = 0; j < numberofstates; j++){
                if(k != j){
                    writer.print(tempratematrix[k][j] + " ");
                }
            }
        }
        for(int i = 0; i < numberofstates-1; i++){
            writer.print(testedsim.getNumProtected(i) + " ");
        }
        writer.print(testedsim.getKbreathe() + " ");
        writer.println(gof);
    }

    private void runCoordinateDescent(){ //Runs coordinate descent on a single set of starting parameters
        /*
        The basic idea behind coordinate descent is iterative improvement of the objective function (in our case, Goodness of Fit)
        The GoF of the starting parameter set is obtained
        For each parameter, a gof is obtained from runnning a simulation with the current parameter moved up a small increment or down 
        a small decrement (holding the rest of the parameters constant)
        If either of the increment or decrement increase gof, the parameter is set to the value that had better gof
        Otherwise, the parameter is kept at its original value and the stability counter increases
        This process is done for each parameter sequentially, until such point as the maximum number of iterations is reached or 
        changing the parameters hasn't improved the gof for two loops through all the parameters are fit
        Every simulation conducted has its parameters written out to the transient file for this particular Runnable
        */
        long starttime;
        ParameterSet psup;
        ParameterSet psdown;
        Sim simup;
        Sim simdown;
        Sim currentsim = new Sim(current);
        double gofcurrent = currentsim.getGOF(expdata);
        double gofup;
        double gofdown;
        if(debug){
            System.out.println("Thread " + threadnumber + " has started coordinate descent");
        }
        while(iteration <= maxiterations && stabilitycounter < 2*numparametersfit){ 
            starttime = System.nanoTime();
            int position = 0;
            // Evaluating rate matrix parameters
            for(int i = 0; i < numberofstates; i++){
                for(int j = 0; j < numberofstates; j++){
                    if(i != j){
                        if(isfit[position]){
                        //Evaluating what happens when the specified rate constant between folding states goes up
                        psup = current.deepCopy();
                        psup.incrementStateRateConstant(i, j, changepercent);
                        simup = new Sim(psup);
                        gofup = simup.getGOF(expdata);
                        writeSimToStorage(psup, gofup);
                        //Evaluating what happens when the specified rate constant between folding states goes down
                        psdown = current.deepCopy();
                        psdown.incrementStateRateConstant(i, j, (-1.0)*changepercent);
                        simdown = new Sim(psdown);
                        gofdown = simdown.getGOF(expdata);
                        writeSimToStorage(psdown, gofdown);

                        if(gofcurrent >= gofup && gofcurrent >= gofdown){
                            stabilitycounter++;
                        }
                        else{
                            stabilitycounter = 0;
                            if(gofup > gofdown){
                                current = psup;
                                //simcurrent = simup;
                                gofcurrent = gofup;
                            }
                            else{
                                current = psdown;
                                //simcurrent = simdown;
                                gofcurrent = gofdown;
                            }
                        }
                    }
                    position++;
                    }
                }
            }

            // Evaluating number of protected sites
            for(int i = 0; i < numberofstates-1; i++){ //Number of fully protected sites in non-fully unfolded states
                if(isfit[position]){
                    //Evaluating what happens when the breathing rate constant goes up
                    //- (int) Math.max(changepercent*psup.getNumProtected(i), 1)
                    psup = current.deepCopy();
                    psup.setNumProtected(i, psup.getNumProtected(i) + (int) Math.max(changepercent*psup.getNumProtected(i), 1));
                    simup = new Sim(psup);
                    gofup = simup.getGOF(expdata);
                    writeSimToStorage(psup, gofup);
                    //Evaluating what happens when the breathing rate constant goes down
                    psdown = current.deepCopy();
                    psup.setNumProtected(i, psup.getNumProtected(i) + (-1) * (int) Math.max(changepercent*psup.getNumProtected(i), 1));
                    simdown = new Sim(psdown);
                    gofdown = simdown.getGOF(expdata);
                    writeSimToStorage(psdown, gofdown);

                    if(gofcurrent >= gofup && gofcurrent >= gofdown){
                        stabilitycounter++;
                    }
                    else{
                        stabilitycounter = 0;
                        if(gofup > gofdown){
                            current = psup;
                            //simcurrent = simup;
                            gofcurrent = gofup;
                        }
                        else{
                            current = psdown;
                            //simcurrent = simdown;
                            gofcurrent = gofdown;
                        }
                    }
                }
                position++;
            }
            if(isfit[position]){
                //Evaluating what happens when the breathing rate constant goes up
                psup = current.deepCopy();
                psup.incrementKbreathe(changepercent);
                simup = new Sim(psup);
                gofup = simup.getGOF(expdata);
                writeSimToStorage(psup, gofup);
                //Evaluating what happens when the breathing rate constant goes down
                psdown = current.deepCopy();
                psdown.incrementKbreathe((-1.0)*changepercent);
                simdown = new Sim(psdown);
                gofdown = simdown.getGOF(expdata);
                writeSimToStorage(psdown, gofdown);

                if(gofcurrent >= gofup && gofcurrent >= gofdown){
                    stabilitycounter++;
                }
                else{
                    stabilitycounter = 0;
                    if(gofup > gofdown){
                        current = psup;
                        //simcurrent = simup;
                        gofcurrent = gofup;
                    }
                    else{
                        current = psdown;
                        //simcurrent = simdown;
                        gofcurrent = gofdown;
                    }
                }
                position++;
            }
            if(debug | progresstracker){
                System.out.println("Thread " + threadnumber + " has finished iteration " + iteration + " of coordinate descent (" + 
                        String.format("%.1f", (((double )System.nanoTime() - starttime)/1000000000)) + " seconds)");
            }       
            iteration++;
        }
    }
    
    private void calculateNumParameters(){
        for(int i = 0; i < isfit.length; i++){
            if(isfit[i]){
                numparametersfit++;
            }
        }
    }
    


}