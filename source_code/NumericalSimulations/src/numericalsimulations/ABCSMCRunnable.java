package numericalsimulations;
import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;
import java.io.FileWriter;
import java.io.PrintWriter;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.random.*;

public class ABCSMCRunnable implements Runnable{
    
    //Parameters for debugging/progress
    private boolean debug = false; //Whether to print statements for debugging purposes
    private boolean progresstracker = true; //Whether to print statements to give progress
    private long start; //When this fitting runnable was started

    private int numberofstates; //Number of protein folding states
    private double[][] expdata; //Narrow format of experimental data (Time, m/z, Intensity)
    private boolean[] isfit; //Which of the possible parameters that can be fit, is to be fit, as specified by an uncertainty of 0
    private double[][] deltaparameters; //The min and max of each parameter under the prior
    private String priortype; //The distribution (either unif or logunif) of the prior
    private double[][] distparams; //The parameters and weights for all accepted particles from the last distribution
    private int weightcol; //The column index of the weights in the distparams matrix
    private double[][] paramvariances; //The covariance matrix of parameters in the last distribution
    private double propvariance; //The proportion of variance to retain for each parameter when perturbing
    private Random randomgenerator; //Random object used by perturbation kernel, allows for a global seed
    private MultivariateNormalDistribution perturbationkernel; //Multivariate normal random variable particles are perturbed with

    private int threadnumber; //Thread number of the current object
    private ParameterSet currentps; //ParameterSet that gets changed and then tested
    private int currentdist; //The current distribution of ABCSMC
    private double currentgof; //The gof of the current particle tested
    private double[] sampledparams; //The parameters of the current particle
    private int desiredpopsize; //The number of accepted particles needed to end the sampling
    private double epsilon; //The distance of the current particle
    private double epsilont; //The epsilon cutoff for this population
    private int checkinterval = 25; //After this number of simulations, check that other threads haven't already finished the task
    private AtomicInteger currentpopsize; //The current number of accepted particles
    private String directory; //Where to write out the data

    //Parameters for moving data
    PrintWriter writer; //Allows for writing out of ABC SMC sample evaluation process 

    public ABCSMCRunnable(int dist, ParameterSet startingset, double[][] experimentaldata, boolean[] whichparamsfit, 
            double[][] priorlimits, String distofprior, double[][] allsimresults, double[][] priordistvariances, 
            double propretainedvariance, AtomicInteger pop, double epsiloncutoff, int goaln, int threadn, String directorypath, 
            Random generator, boolean shoulddebug){ 
        start = System.nanoTime();
        currentdist = dist;
        currentps = startingset;
        currentps.setRandomGenerator(generator.nextDouble());
        expdata = experimentaldata;
        isfit = whichparamsfit;
        deltaparameters = priorlimits;
        priortype = distofprior;
        distparams = allsimresults;
        paramvariances = new double[priordistvariances.length][priordistvariances[0].length];        
        propvariance = propretainedvariance;
        for(int i = 0; i < priordistvariances.length; i++){
            for(int j = 0; j < priordistvariances[0].length; j++){
                paramvariances[i][j] = Math.sqrt(propvariance)*priordistvariances[i][j];
            }
        }
        threadnumber = threadn;
        currentpopsize = pop;
        epsilont = epsiloncutoff;
        desiredpopsize = goaln;
        threadnumber = threadn;
        directory = directorypath;
        randomgenerator = generator;
        weightcol = distparams[0].length-1;
        debug = shoulddebug;
    }

    @Override
    public void run(){ //When the threads in Fitting are "started", this method is run
        try{
            if(debug | progresstracker){
                System.out.println("Thread " + threadnumber + " has been started");
            }
            makeStorage();
            makePerturber();
            while(currentpopsize.get() < desiredpopsize){
                epsilon = epsilont + 1;
                int i = 1;
                while(epsilon > epsilont){
                    try{
                        sample();
                        perturb();
                        test();
                        i++;
                    }
                    catch(Exception e){
                        System.out.println("There was an issue sampling, perturbing, and testing in thread " + threadnumber);
                        e.printStackTrace();
                        System.out.println(e);
                        System.exit(0);
                    }
                    if(checkinterval % i == 0){
                        if(currentpopsize.get() >= desiredpopsize){
                            break;
                        }
                    }
                }
            }
            writer.close();
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
    
    public void makePerturber(){ //Create multivariate normal distribution to perturb with
        if(currentdist != 0){
            try{
                perturbationkernel = new MultivariateNormalDistribution(new JDKRandomGenerator(randomgenerator.nextInt()), 
                        new double[paramvariances.length], paramvariances);
            }catch(Exception e){
                //This exception occurs when the covariance matrix is singular, i.e., not invertible
                //This exception occurs because of rounding errors, as the "invertible" matrices have inverses when tested in R
                e.printStackTrace();
                System.out.println(e);
                String temp = "The following covariance matrix to be used for multivariate normal random sampling was singular: \n";
                for(int i = 0; i < paramvariances.length; i++){
                    for(int j = 0; j < paramvariances[0].length; j++){
                        temp = temp + paramvariances[i][j] + " ";
                    }
                    temp = temp + "\n";
                }
                System.out.println(temp);
                System.exit(0);
            }
        }
        if(debug){
            System.out.println("The prior limits perturbed particles must fall in are : ");
            for(int i = 0; i < deltaparameters.length; i++){
                for(int j = 0; j < deltaparameters[0].length; j++){
                    System.out.print(deltaparameters[i][j] + " ");
                }
                System.out.println();
            }
        }
    }

    public void sample(){ //sample a point to test
        if(currentdist == 0){
            sampledparams = distparams[0];
            int position = 0;
            //Set of rate matrix parameters
            for(int k = 0; k < numberofstates; k++){
                for(int j = 0; j < numberofstates; j++){
                    if(k != j){
                        if(isfit[position]){
                            if(priortype.equalsIgnoreCase("unif")){
                                sampledparams[position+1] = runif(deltaparameters[0][position], deltaparameters[1][position]);
                            }else{
                                sampledparams[position+1] = rlogunif(deltaparameters[0][position], deltaparameters[1][position]);
                            }
                            currentps.setRateMatrixPosition(k, j, sampledparams[position+1]);
                        }
                        position++;
                    }
                }
            }

            //Set of protection information
            for(int j = 0; j < numberofstates-1; j++){
                if(isfit[position]){
                    sampledparams[position+1] = (int) runif(deltaparameters[0][position], deltaparameters[1][position]+1.0);
                    currentps.setNumProtected(j, (int) sampledparams[position+1]);
                }
                position++;
            }

            //Kbreathe
            if(isfit[position]){
                if(priortype.equalsIgnoreCase("unif")){
                    sampledparams[position+1] = runif(deltaparameters[0][position], deltaparameters[1][position]);
                }else{
                    sampledparams[position+1] = rlogunif(deltaparameters[0][position], deltaparameters[1][position]);
                }
                currentps.setKbreathe(sampledparams[position+1]);
            }
        }else{
            double draw = runif();
            double cumweight = 0;
            for(int i = 0; i< distparams.length; i++){
                cumweight = cumweight + distparams[i][weightcol];
                if(draw <= cumweight){
                    sampledparams = distparams[i];
                    if(debug){
                        System.out.println("Thread " + threadnumber + " sampled the following particle: ");
                        for(int j = 0; j < sampledparams.length; j++){
                            System.out.print(sampledparams[j] + " ");
                        }
                        System.out.println();
                    }
                    return;
                }
            }
        }
    }

    public void perturb(){ //Pertub the particle using a gaussian random walk
        if(currentdist != 0){
            boolean notinpriorsupport = true;
            //update currentps
            //need to know what the sample from
            while(notinpriorsupport){
                double[] perturbdistances = perturbationkernel.sample();
                int position = 0;
                int parameterfit = 0;
                //Set of rate matrix parameters
                for(int k = 0; k < numberofstates; k++){
                    for(int j = 0; j < numberofstates; j++){
                        if(k != j){
                            if(isfit[position]){
                                currentps.setRateMatrixPosition(k, j, sampledparams[position+1] + perturbdistances[parameterfit]);
                                parameterfit++;
                            }
                            position++;
                        }
                    }
                }

                //Set of protection information
                for(int j = 0; j < numberofstates-1; j++){
                    if(isfit[position]){
                        currentps.setNumProtected(j, (int) (sampledparams[position+1] + perturbdistances[parameterfit]));
                        parameterfit++;
                    }
                    position++;
                }

                //Kbreathe
                if(isfit[position]){
                    currentps.setKbreathe(sampledparams[position+1] + perturbdistances[parameterfit]);
                }
                notinpriorsupport = !inPriorSupport();
                if(debug && notinpriorsupport){
                    System.out.println("Thread " + threadnumber + 
                            " threw away a possible particle for not being in the prior's support");
                }
            }
            if(debug){
                System.out.println("Thread " + threadnumber + " then perturbed the particle to: ");
                currentps.toString();
            }
        }
    }

    private boolean inPriorSupport(){
        int position = 0;
        //Set of rate matrix parameters
        for(int k = 0; k < numberofstates; k++){
            for(int j = 0; j < numberofstates; j++){
                if(k != j){
                    if(currentps.getRateMatrix(k, j) < deltaparameters[0][position] || 
                            currentps.getRateMatrix(k, j) > deltaparameters[1][position]){
                        sample();
                        return(false);
                    }
                    position++;
                }
            }
        }

        //Set of protection information
        for(int j = 0; j < numberofstates-1; j++){
            if(currentps.getNumProtected(j) < deltaparameters[0][position] || 
                    currentps.getNumProtected(j) > deltaparameters[1][position]){
                sample();
                return(false);
            }
            position++;
        }

        //Kbreathe
        if(currentps.getKbreathe() < deltaparameters[0][position] || currentps.getKbreathe() > deltaparameters[1][position]){
            sample();
            return(false);
        }
        return(true);
    }

    public void test(){
        currentgof = new Sim(currentps).getGOF(expdata);
        epsilon = Math.abs(1-currentgof);
        double weight = 0;
        if(epsilon < epsilont){
            int nacceptedparticles = currentpopsize.get();
            if(nacceptedparticles < desiredpopsize){
                weight = 1;
                System.out.println(currentpopsize.incrementAndGet() + " out of " + desiredpopsize + 
                        " particles have been accepted (epsilon cutoff = " + epsilont + ")");
            }
        }
        writeSimToStorage(currentps, weight);
    }


    public double runif(){
        return randomgenerator.nextDouble();
    }

    public double runif(double tempmin, double tempmax){
        double min = Math.min(tempmin, tempmax);
        double max = Math.max(tempmin, tempmax);
        double x = min + runif()*(max - min);
        return x;
    }
    
    public double rlogunif(double tempmin, double tempmax){
        double logmin = Math.log(Math.min(tempmin, tempmax));
        double logmax = Math.log(Math.max(tempmin, tempmax));
        double logx = logmin + runif()*(logmax - logmin);
        double x = Math.exp(logx);
        return x;
    }

    //Writes a data file containing the information for the current threadnumber, the off-diagonal rate constants, 
    //the number of sites fully protected in not fully unfolded states, kbreathe, goodness of fit, and weight
    private void makeStorage(){ 
        try{
            numberofstates = currentps.getNStates(); 
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
            writer.print("GoF ");
            writer.println("Weight");
        }
        catch(Exception e) {
            System.out.println("Error making storage file");
        }
    }

    //Called whenever a new simulation is run to store all values for that simulation, including goodness of fit
    private void writeSimToStorage(ParameterSet testedsim, double weight){ 
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
        writer.print(currentgof + " ");
        writer.println(weight); //the true weight is calculated in the main thread
    }
}


