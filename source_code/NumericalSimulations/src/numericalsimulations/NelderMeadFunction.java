package numericalsimulations;

import flanagan.math.*;
import java.io.PrintWriter;

class NelderMeadFunction implements MaximizationFunction{
    
    //Parameters for debugging/progress
    boolean debug = false; //Whether to print statements for debugging purposes
    boolean progresstracker = true; //Whether to print statements to give progress
    long starttime; //When this fitting runnable was started
    
    //Parameters for fitting
    ParameterSet current; //The starting ParameterSet object Nelder Mead optimization modifies 
    double[][] expdata; //Narrow format of experimental data (Time, m/z, Intensity)
    int threadnumber; //This runnable's number, for documenation purposes
    int numberofstates; //Number of protein folding states
    boolean[] isfit; //Which of the possible parameters that can be fit, is to be fit
    int iteration = 1; //Current iteration of nelder mead (technically just a function evaluation)
    double bestgof = 0; //Current best goodness of fit achieved
    
    //Parameters for moving data
    PrintWriter writer; //Allows for writing out of Nelder Mead process  
    
    @Override
    public double function(double[] params){
        modifyParameterSet(params);
        Sim tempsim = new Sim(current);
        double gof = tempsim.getGOF(expdata);
        writeSimToStorage(current, gof);
        iteration++;
        return(gof);
    }
    
    public void setParameterSet(ParameterSet starting){
        this.current = starting;
        numberofstates = current.getNStates();
    }
    
    public void setThreadNumber(int threadn){
        threadnumber = threadn;
    }
    
    public void setFitParameters(boolean[] whichparamsfit){
        isfit = whichparamsfit;
    }
    
    private void modifyParameterSet(double[] params){
        int position = 0;
        int positioninfunction = 0;
        //Set of rate matrix parameters
        for(int k = 0; k < numberofstates; k++){
            for(int j = 0; j < numberofstates; j++){
                if(k != j){
                    if(isfit[position]){
                        current.setRateMatrixPosition(k, j, params[positioninfunction]);
                        positioninfunction++;
                    }
                    position++;
                }
            }
        }

        //Set of protection information
        for(int j = 0; j < numberofstates-1; j++){
            if(isfit[position]){
                current.setNumProtected(j, (int) params[positioninfunction]);
                positioninfunction++;
            }
            position++;
        }

        //Kbreathe
        if(isfit[position]){
            current.setKbreathe(params[positioninfunction]);
        }
    }
    
    public void setExperimentalData(double[][] experimentaldata){
        expdata = experimentaldata;
    }
    
    //Writes a data file containing the information for the current threadnumber, the off-diagonal rate constants, 
    //the number of sites fully protected in not fully unfolded states, kbreathe, and goodness of fit
    public void makeStorage(PrintWriter writer){ 
        try{
            //Makes txt file, adds variable headings
            this.writer = writer;
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
        if(gof > bestgof){
            if(debug || progresstracker){
                System.out.println("Thread " + threadnumber + " has finished iteration " + iteration + " of Nelder-Mead (gof = " + 
                        String.format("%.3f", gof) + ")");
            }
            bestgof = gof;
        }
    }
}
