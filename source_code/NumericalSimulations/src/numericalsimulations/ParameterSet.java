package numericalsimulations;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

public class ParameterSet{

private boolean debug = false; //Debugging parameter    
private double[][] protectionmatrix; //Matrix containing protection, and breathing information for each state
private double[][] ratematrix; //Matrix of rate constants
private double[] kints; //List of selected kints
private int[] intparams; //Contains the integer parameters: number of proteins, charge state, and resolution
private double[] doubleparams; //Contains the double parameters: mass, gaussian width, initial isotype saturation, and kbreathe
private double[] tps; //List of experimental timepoints
private boolean[] booleanparams; //Contains whether or not to normalze to height
private Random generator; //Used to set global seed
private double[] mzs; //MZ values
private final double minrateconstant = Math.pow(10, -20); //Minimum rate constant between states
private String gofmethod; //Goodness of fit algorithm
    
    public ParameterSet(double[][] inputprotectionmatix, double[][] inputratematrix, double[] inputkints, int[] inputintparams, 
            double[] inputdoubleparams, boolean[] inputbooleanparams, double[] timepoints, 
            Random globalgenerator, double[] mzvals, String goftouse){
        protectionmatrix = inputprotectionmatix;
        ratematrix = inputratematrix;
        kints = inputkints;
        intparams = inputintparams;
        doubleparams = inputdoubleparams;
        tps = timepoints;
        booleanparams = inputbooleanparams;
        generator = globalgenerator;
        mzs = mzvals;
        gofmethod = goftouse;
    }
    
    public String getGoFMethod(){
        return gofmethod;
    }
    public void setGoFMethod(String goftouse){
        gofmethod = goftouse;
    }
    public int getNStates(){ //Returns number of states
        return ratematrix.length;
    }

    public double[][] getProtectionMatrix(){ //Returns protection matrix
        return protectionmatrix;
    }
    
    public int getNumExchangers(){ //Returns protection matrix
        return protectionmatrix[0].length;
    }

    public double[] getMZ(){ //Returns m/z values at which to calculate intensities
        return mzs;
    }

    public void setMZ(double[] mzvals){ //Sets list of m/z values to input
        mzs = mzvals;
    }
    
    public double[][] getRateMatrix(){ //Returns the rate matrix (or generator matrix, Q) for transition between folding states
        return ratematrix;
    }
    
    //Returns the rate matrix (or generator matrix, Q) for transition between folding states
    public double getRateMatrix(int row, int col){ 
        return ratematrix[row][col];
    }

    public void setRateMatrixPosition(int row, int col, double value){ //Sets desired position of rate matrix to given value
        ratematrix[row][col] = Math.max(minrateconstant, value);
    }
    
    //Returns matrix containing rate constants of HX (kHX) for all sites in all global folding states
    public double[][] getFinalMatrix(){ 
        double[][] finalmatrix = new double[protectionmatrix.length][protectionmatrix[0].length];
        for(int i = 0; i < finalmatrix.length; i++){  //Change the -1s to kbreathe for later use
            for(int j = 0; j < finalmatrix[0].length; j++){
                if(protectionmatrix[i][j] == -1){
                    finalmatrix[i][j] = doubleparams[3]*kints[j];
                }
                else if(protectionmatrix[i][j] == 1){
                    finalmatrix[i][j] = kints[j];
                }
            }
        }
        return finalmatrix;
    }
    
    public int getNumberProteins(){ //Returns number of proteins to simulate
        return intparams[0];
    }
    
    public int getChargeState(){ //Returns charge state of experimental data
        return intparams[1];
    }
    
    public int getResolution(){ //Returns number of m/z values at which to calculate intensities at each time point
        return intparams[2];
    }
    
    public double[] getTimePoints(){ //Returns list of time points at which to calculate HX distribution
        return tps;
    }
    
    public double getMass(){ //Returns standard mass of the protein (not including mass from charged ions)
        return doubleparams[0];
    }
    
    public double getGaussianWidth(){ //Returns gaussian width of reference peak
        return doubleparams[1];
    }
    
    public double getDirectionofExchange(){ //Returns direction of exchange; -1 for D->H, 1 for H->D
        return doubleparams[2];
    }
    
    //Returns breathing rate constant (assumed to be constant across sites, to be combined with k intrinsic)
    public double getKbreathe(){ 
        return doubleparams[3];
    }

    public void setKbreathe(double newkbreathe){ //Sets breathing rate constant
        doubleparams[3] = Math.max(0, newkbreathe);
    }
    
    public Random getRandomGenerator(){ //Returns Random object used to implement global seed
        return generator;
    }
    
    public void setRandomGenerator(double seed){ //Sets Random object used to implement global seed
        generator = new Random();
        generator.setSeed((long) seed);
    }
    
    public boolean getNormalizationType(){ //Returns normalization type simulation should use
        return booleanparams[0];
    }
    
    public double incrementKbreathe(double percentchange){ //Modifies kbreathe by specified percentage up or down
        double nextvalue = doubleparams[3]*(1 + percentchange);
        this.setKbreathe(nextvalue);
        return(this.getKbreathe());
    }

    public void setNumProtected(int state, int newnum){ //Sets number of sites protected in a specified state, to the number desired
        int newprotected;
        if(newnum < 0){
            newprotected = 0;
        }
        else if(newnum > protectionmatrix[0].length){
            newprotected = protectionmatrix[0].length;
        }
        else if(newnum == getNumProtected(state)){
            return;
        }
        else{
            newprotected = newnum;
        }
        //Selecting r sites from n total to be protected
        int r = newprotected;
        int n = protectionmatrix[0].length;
        int begin;
        int end;
        
        /* interior
        int begin = ((n - r)/2);
        int end = ((n + r)/2);

        //exterior
        int begin = (r/2);
        int end = (n-r/2);
        */

        if(state == 0){
            begin = (r/2);
            end = begin + n - r;
            for(int i = 0; i < n; i++){ //The exterior sites are protected by default in the native state
                if(i >= begin && i < end){
                    protectionmatrix[0][i] = -1; 
                }
                else{
                    protectionmatrix[0][i] = 0;
                } 
            }
        }else if (state < getNStates() - 1){
            if(state % 2 != 0){ //If we have one intermediate state, this is the default; interior is protected
                begin = ((n - r)/2);
                end = ((n + r)/2);
                for(int j = 0; j < n; j++){ //The exterior sites are protected by default in the native state
                    if(j >= begin && j < end){
                        protectionmatrix[state][j] = 0; 
                    }
                    else{
                        protectionmatrix[state][j] = 1;
                    } 
                }
            }
            else{ //Exterior is protected; this allows intermediate states to alternate whether protected sites are interior/exterior
                begin = (r/2);
                end = (n-r/2);
                for(int j = 0; j < n; j++){ //The exterior sites are protected by default in the native state
                    if(j >= begin && j < end){
                        protectionmatrix[state][j] = 1; 
                    }
                    else{
                        protectionmatrix[state][j] = 0;
                    } 
                }
            }
        }else{
            ArrayList<Integer> templocations = sample(r, n);  
            for(int j = 0; j < n; j++){ //A random sample of all sites are set as fully protected; the rest exchange without breathing
                if(templocations.contains(j)){
                    protectionmatrix[getNStates()-1][j] = 0; 
                }
                else{
                    protectionmatrix[getNStates()-1][j] = 1;
                }
            }
        }
        if(debug){
            System.out.println("We tried to set " + r + " sites to be protected in state " + state);
            for(int i = 0; i < getNStates(); i++){
                int count = 0;
                for(int j = 0; j < n; j++){
                    if(protectionmatrix[i][j] == 0){
                        count++;
                    }
                }
                System.out.println(count + " sites are protected in state " + i);
            }
        }
    }

    public int getNumProtected(int state){ //Returns the number of sites protected in the given state
        int numprotected = 0;
        for(int i = 0; i < protectionmatrix[state].length; i++){
            if(protectionmatrix[state][i] == 0){
                numprotected++;
            }
        }
        return numprotected;
    }
    
    //Modifies folding rate constants by specified percentage up or down
    public double incrementStateRateConstant(int row, int column, double percentchange){ 
        double nextvalue = ratematrix[row][column]*(1 + percentchange);
        this.setRateMatrixPosition(row, column, nextvalue);
        return(this.getRateMatrix(row, column));
    }
    
    public ParameterSet deepCopy(){ //Returns a completely new copy of the current parameter set (not just a pointer to the same object)
        double[][] newprotectionmatrix = new double[protectionmatrix.length][protectionmatrix[0].length];
        for(int i = 0; i < protectionmatrix.length; i++){
            for(int j = 0; j < protectionmatrix[0].length; j++){
                newprotectionmatrix[i][j] = protectionmatrix[i][j];
            }
        }
        double[][] newratematrix = new double[ratematrix.length][ratematrix[0].length];
        for(int i = 0; i < ratematrix.length; i++){
            for(int j = 0; j < ratematrix[0].length; j++){
                newratematrix[i][j] = ratematrix[i][j];
            }
        }
        double[] newkints = new double[kints.length];
        for(int i = 0; i < kints.length; i++){
            newkints[i] = kints[i];
        }
        int[] newintparams = new int[intparams.length];
        for(int i = 0; i < intparams.length; i++){
            newintparams[i] = intparams[i];
        }
        double[] newdoubleparams = new double[doubleparams.length];
        for(int i = 0; i < doubleparams.length; i++){
            newdoubleparams[i] = doubleparams[i];
        }
        boolean[] newbooleanparams = new boolean[booleanparams.length];
        for(int i = 0; i < booleanparams.length; i++){
            newbooleanparams[i] = booleanparams[i];
        }
        double[] newtps = new double[tps.length];
        for(int i = 0; i < tps.length; i++){
            newtps[i] = tps[i];
        }
        double[] newmzs = new double[mzs.length];
        for(int i = 0; i < mzs.length; i++){
            newmzs[i] = mzs[i];
        }
        //note that the random object is later overridden in Fitting's initializeParameterSets method to enable multithreaded seeding
        return new ParameterSet(newprotectionmatrix, newratematrix, newkints, newintparams, 
                newdoubleparams, newbooleanparams, newtps, generator, newmzs, gofmethod);
    }
    
    @Override
    public String toString(){ //Prints out the rate matrix, protection matrix, and kbreathe for the current parameter set
        System.out.println("The rate matrix of this parameter set object is: ");
        for(int i = 0; i < ratematrix.length; i++){
            for(int j = 0; j < ratematrix[0].length; j++){
                System.out.print(ratematrix[i][j] + " ");
            }
            System.out.println();
        }
        System.out.println("The protection matrix of this parameter set object is: ");
        for(int i = 0; i < protectionmatrix.length; i++){
            for(int j = 0; j < protectionmatrix[0].length; j++){
                System.out.print(protectionmatrix[i][j] + " ");
            }
            System.out.println();
        }
        System.out.println("The kbreathe of this parameter set object is: " + getKbreathe());
        return "";
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
    
    public double runif(){ //Generates random draw from Uniform(0, 1) random variable, using a global seed
        return generator.nextDouble();
    }
}