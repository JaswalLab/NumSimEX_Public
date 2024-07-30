package numericalsimulations;

import java.io.*;
import java.util.*;

public class NumericalSimulations {
   
    public static void main(String[] args) {  //Main method, allows for the creation of Sim or Fitting objects
        System.out.println("Welcome to the Jaswal Lab Hydrogen Exchange Mass Spectrometry numerical simulations program");
        long start = System.nanoTime(); //start timer
        ArrayList<String> filelocations = new ArrayList<>();
        boolean debug = false;
        boolean verbose = false;
        boolean benchmark = false;
        for(int i = 0; i < args.length; i++){
            if(args[i].equalsIgnoreCase("debug") || args[i].equalsIgnoreCase("d")){
                debug = true;
            }
            else if(args[i].equalsIgnoreCase("verbose") || args[i].equalsIgnoreCase("v")){
                verbose = true;
            }
            else if(args[i].equalsIgnoreCase("benchmark") || args[i].equalsIgnoreCase("b")){
                benchmark = true;
            }
            else{
                filelocations.add(args[i]);
            } 
        }
        int numberofentries = filelocations.size();
        for(int i = 0; i < numberofentries; i++){
            System.out.println("Filepath entry " + (i + 1) + " to the program was '" + filelocations.get(i) + "'");
        }
        Sim run;
        switch (numberofentries) {
            case 2:
                if(!benchmark){
                    run = new Sim(filelocations.get(0), filelocations.get(1), debug, verbose);
                }else{
                    benchmark(filelocations.get(0), filelocations.get(1), " ");
                }
                System.out.println("This program took " + String.format("%.3f", (((double )System.nanoTime() - start)/1000000000)) + 
                        " seconds");
                break;
            case 3:
                if(!benchmark){
                    run = new Sim(filelocations.get(0), filelocations.get(1), filelocations.get(2), debug, verbose);
                }else{
                    benchmark(filelocations.get(0), filelocations.get(1), filelocations.get(2));
                }                
                System.out.println("This program took " + String.format("%.3f", (((double )System.nanoTime() - start)/1000000000)) + 
                        " seconds");
                break;
            case 4:
                Fitting fit = new Fitting(filelocations.get(0), filelocations.get(1), filelocations.get(2), 
                        filelocations.get(3), debug, verbose);
                System.out.println("This program took " + String.format("%.3f", (((double )System.nanoTime() - start)/1000000000)) + 
                        " seconds");
                break;
            default:
                boolean simpleinput = true; //whether to have prompts to read in information for simulations and fittings
                if(simpleinput){
                    commandLineRun();
                }
                else{
                    System.out.println("Command line running was not selected");
                }
                break;
        } 
    }
    
    public static void commandLineRun(){ //Allows user to choose between command-line entry of either fitting or creation
        Scanner kb = new Scanner(System.in);
        System.out.println("Would you like to conduct a simulation or a fitting? (enter 'sim' or 'fit')");
        String programtype = kb.nextLine();
        try{
            if(programtype.equalsIgnoreCase("sim") | programtype.equalsIgnoreCase("simulation")){
                commandLineSim();
            }
            else if(programtype.equalsIgnoreCase("fit") | programtype.equalsIgnoreCase("fitting")){
                commandLineFit();
            }
            else{
                System.out.println("Please enter either 'sim' or fit' (without the '')");
            }
        }catch(Exception e){
            System.out.println("Uhoh, there appears to have been an error. Try running this program in debug mode.");
        }
        kb.close();
    }
    
    public static void commandLineSim(){ //Allows command-line creation of a simulation object (with or without experimental data input)
        Scanner kb = new Scanner(System.in);
        String siminputprompt = "Please enter the location of the input file - i.e., 'C:\\\\Users\\\\...\\\\TestFolder\\\\SimInput.txt'"
                + " on Windows or '/Users/Username/.../TestFolder/SimInput.txt' on OSX or enter 'help' for the file template";
        System.out.println(siminputprompt);
        String inputlocation = kb.nextLine();
        while(!(new File(inputlocation).isFile()) || inputlocation.equalsIgnoreCase("help") || inputlocation.equalsIgnoreCase("h")){
            if(inputlocation.equalsIgnoreCase("help") || inputlocation.equalsIgnoreCase("h")){
                printSimTemplate();
                System.out.println(siminputprompt);
            }else{
                System.out.println("Please enter a valid file location");
            }
            inputlocation = kb.nextLine();
        }
        
        System.out.println("Please enter the location of the output file - i.e., 'C:\\\\Users\\\\...\\\\TestFolder\\\\SimOutput.txt'" +
                " on Windows or '/Users/Username/.../TestFolder/SimInput.txt' on OSX");
        String outputlocation = kb.nextLine();
        while(outputlocation.equalsIgnoreCase("")){
            System.out.println("Please enter a valid file location");
            outputlocation = kb.nextLine();
        }
        
        System.out.println("If you have experimental data to which you'd like to compare this simulation, " +
                "please enter its file location; otherwise, please press enter");
        System.out.println("(i.e., 'C:\\\\Users\\\\...\\\\TestFolder\\\\ExperimentalData.txt' on Windows or " + 
                "'/Users/Username/.../TestFolder/ExperimentalData.txt' on OSX)");
        String expdatalocation = kb.nextLine();
        
        System.out.println("If you would like to debug this simulation, please enter 'y'; otherwise, please press enter");
        boolean debug = !kb.nextLine().equalsIgnoreCase("");
        
        System.out.println("If you would like to produce an R file to graph this simulation, please enter 'y';" + 
                " otherwise, please press enter");
        boolean verbose = !kb.nextLine().equalsIgnoreCase("");
        
        kb.close();
        try{
            long start = System.nanoTime();
            if(expdatalocation.equalsIgnoreCase("")){
                Sim run = new Sim(inputlocation, outputlocation, debug, verbose);
            }
            else{
                Sim run = new Sim(inputlocation, outputlocation, expdatalocation, debug, verbose);
            }
            System.out.println("The program took " + String.format("%.3f", (((double )System.nanoTime() - start)/1000000000)) + 
                    " seconds");
        }
        catch(Exception e){
            System.out.println("There appears to be a problem");
        }
    }
    
    public static void commandLineFit(){ //Allows command-line creation of a fitting object
        Scanner kb = new Scanner(System.in);
        String siminputprompt = "Please enter the location of the sim input file - i.e., " + 
                "'C:\\\\Users\\\\...\\\\TestFolder\\\\SimInput.txt' on Windows or "+
                "'/Users/Username/.../TestFolder/SimInput.txt' on OSX or enter 'help' for the file template";
        System.out.println(siminputprompt);
        String siminputlocation = kb.nextLine();
        while(!(new File(siminputlocation).isFile()) || siminputlocation.equalsIgnoreCase("help") || 
                siminputlocation.equalsIgnoreCase("h")){
            if(siminputlocation.equalsIgnoreCase("help") || siminputlocation.equalsIgnoreCase("h")){
                printSimTemplate();
                System.out.println(siminputprompt);
            }else{
                System.out.println("Please enter a valid file location");
            }
            siminputlocation = kb.nextLine();
        }
        
        String uncertaintiesprompt = "Please enter the location of the uncertainties file - i.e., " + 
                "'C:\\\\Users\\\\...\\\\TestFolder\\\\Uncertainties.txt' on Windows or " + 
                "'/Users/Username/.../TestFolder/Uncertainties.txt' on OSX or enter 'help' for the file template";
        System.out.println(uncertaintiesprompt);
        String uncertaintieslocation = kb.nextLine();
        while(!(new File(uncertaintieslocation).isFile()) || uncertaintieslocation.equalsIgnoreCase("help") || 
                uncertaintieslocation.equalsIgnoreCase("h")){
            if(uncertaintieslocation.equalsIgnoreCase("help") || uncertaintieslocation.equalsIgnoreCase("h")){
                printFitTemplate();
                System.out.println(uncertaintiesprompt);
            }else{
                System.out.println("Please enter a valid file location");
            }
            uncertaintieslocation = kb.nextLine();
        }
        
        System.out.println("Please enter the location of the experimental data file - i.e., " + 
                "'C:\\\\Users\\\\...\\\\TestFolder\\\\ExperimentalData.txt' on Windows or " + 
                "'/Users/Username/.../TestFolder/ExperimentalData.txt' on OSX");
        String expdatalocation = kb.nextLine();
        while(!(new File(expdatalocation).isFile())){
            System.out.println("Please enter a valid file location");
            expdatalocation = kb.nextLine();
        }
        System.out.println("Please enter the location of the output file - i.e., " + 
                "'C:\\\\Users\\\\...\\\\TestFolder\\\\SimOutput.txt' on Windows or " + 
                "'/Users/Username/.../TestFolder/SimOutput.txt' on OSX");
        String outputlocation = kb.nextLine();
        while(outputlocation.equalsIgnoreCase("")){
            System.out.println("Please enter a valid file location");
            outputlocation = kb.nextLine();
        }
        
        System.out.println("If you would like to debug this fitting (recommended only while using 1 thread), " + 
                "please enter 'y'; otherwise, please press enter");
        boolean debug = !kb.nextLine().equalsIgnoreCase("");
        
        System.out.println("If you would like to produce additional files (for plotting), " + 
                "please enter 'y'; otherwise, please press enter");
        boolean verbose = !kb.nextLine().equalsIgnoreCase("");
        
        kb.close();
        try{
            long start = System.nanoTime();
            Fitting newfitting = new Fitting(siminputlocation, uncertaintieslocation, expdatalocation, outputlocation, debug, verbose);
            System.out.println("The program took " + String.format("%.3f", (((double )System.nanoTime() - start)/1000000000)) + 
                    " seconds");
        }
        catch(Exception e){
            System.out.println("There appears to be a problem");
        }
    }
    
    //Test the duration of a set of simulations
    public static void benchmark(String inputlocation, String timeoutputlocation, String expdatalocation){ 
        int[] nproteins = new int[]{1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 
            10000, 20000, 50000, 100000, 200000, 500000, 1000000};
        System.out.println("Simulation speed benchmarking will be conducted on the following numbers of proteins:");
        for(int i = 0; i < nproteins.length; i++){
            System.out.print(nproteins[i] + " ");
        }
        System.out.println();
        Scanner kb = new Scanner(System.in);
        System.out.println("Please enter a seed");
        long seed = kb.nextLong();
        Random randomgenerator = new Random();
        randomgenerator.setSeed(seed);
        
        System.out.println("Please enter a number of replicates");
        int replicates = kb.nextInt();
        
        long[] seeds = new long[nproteins.length*replicates];
        for(int i = 0; i < seeds.length; i++){
            seeds[i] = randomgenerator.nextLong();
        }
        PrintWriter writer;
        try{
            writer = new PrintWriter(new FileWriter(timeoutputlocation), true);
            writer.println("n_proteins seed time gof");
        
            for(int i = 0; i < nproteins.length; i++){
                System.out.println("Starting simulations for " + nproteins[i] + " proteins");
                for(int j = 0; j < replicates; j++){
                    long start = System.nanoTime();
                    long simseed = seeds[i*replicates + j];
                    Sim tempsim = new Sim(inputlocation, expdatalocation, simseed, nproteins[i]);
                    writer.println(nproteins[i] + " " + simseed + " " + 
                            String.format("%.3f", (((double )System.nanoTime() - start)/1000000000)) + " " + tempsim.getGOF());
                }
            }
            writer.close();
        }
        catch(Exception e){
            System.out.println("Error making records of time");
        }
    }
    
    public static void printSimTemplate(){
       System.out.println(
        "INPUT: Simulations \n" + 
        "The input file for simulations must be formatted as follows. Please note that you should not include the parentheses in the " + 
                "actual file, just the numbers described in the parentheses separated by spaces. \n" + 
        "Line numbers are given as a function of the number of states simulated. For some lines, such as n+7, descriptions of entries "+
                "are given on multiple lines - please enter all selections on the same line. \n" + 
        "Here is the general input text file format to run a simulation: \n" + 
        "--------------------------------- \n" + 
        "1|	(Ignored) Simulation notes \n" + 
        "2| 	(number of states) (number of proteins) \n" + 
        "3|	0   k12 k13... k1n \n" + 
        "4|	k21 0	k23... k2n \n" + 
        "	k31 k32 0  ... k3n \n" + 
        "	. \n" + 
        "	. \n" + 
        "	. \n" + 
        "n+2|	kn1 kn2	kn3... 0 \n" + 
        "n+3|	(kint1) (kint2) ... (kintn) \n" + 
        "n+4|  	(# fully protected exchanging sites in state 1 - native) ('' state 2) ... ('' state n - fully unfolded) \n" + 
        "n+5|	(timepoint 1) (timepoint 2) ... (timepoint n) \n" + 
        "n+6|	(seed) (number of exchangers) (KbreathingRate) (mass) (charge state) (width of Gaussians) (resolution) \n" + 
        "n+7|	(Method of choosing which kints correspond to exchanging sites - 'Default', 'Smallest', 'Largest', 'Random_Sample', " + 
                "'Located_Middle', 'Located_Ends', or 'Manual') \n" + 
        "		(Method of choosing which sites are protected in which states - 'Default', or 'Manual') \n" + 
        "		(initial isotope saturation of the protein, 'D' or 'H') \n" + 
        "		(normalization type - 'Height' or 'Area') \n" + 
        "		((optional) Goodness of Fit measure - defaults to 'Standard', 'Inverse_Time' can be chosen) \n" + 
        "Use the following lines only if 'Manual' is chosen in line n+4, entering 1 if site is exposed; 0, if site is buried and " + 
                "cannot exchange; or -1, if the site is subject to breathing \n" + 
        "n+8|	(State 1, Site 1) (State 1, Site 2) (State 1, Site 3) ... (State 1, Site (#ofexchangers)) \n" + 
        "n+9|	(State 2, Site 1) (State 2, Site 2) (State 2, Site 3) ... (State 2, Site (#ofexchangers)) \n" + 
        "		. \n" + 
        "		. \n" + 
        "		. \n" + 
        "2n+7| 	(State n, Site 1) (State n, Site 2) (State n, Site 3) ... (State n, Site (#ofexchangers)) \n" + 
        "--------------------------------- \n" + 
        "If 'Manual' is chosen as the kint selection method in n+7, please enter only the kints for those sites that are exchanging " + 
                "in line n+3. \n" + 
        "In line n+4, the native state (state 1) should have the greatest number of protected sites. \n" + 
        "All lines past n+7 (or past 2n+7, if 'Manual' was chosen in line n+4) are ignored and can be used for notes. \n"
       );
    }
    
    public static void printFitTemplate(){
        System.out.println(
        "INPUT: Fittings \n" + 
        "Fittings require 3 files to be run: \n" + 
        "1) An uncertainty file (for which this file is a template) \n" + 
        "2) A simulation file (see the other template) \n" + 
        "3) A file of experimental HXMS data.  This data may be formatted in either of two ways: \n" + 
        "	i) A narrow dataset (Three columns; 'Time', 'm/z', and 'Intensity') \n" + 
        "	ii) A series of pairs of rows; the number of pairs is the number of timepoints collected, the first row of each pair " +
                "consists of m/z values; the second, of intensity values \n" + 
        "The input file must be formatted as follows. Please note that you should not include the parentheses in the actual file, " + 
                "just the uncertainties of the parameters described in the parentheses, separated by spaces.  \n" + 
        "Line numbers are given as a function of the number of states simulated. \n" + 
        "Here is the general input text file format to run a fitting: \n" + 
        "--------------------------------- \n" + 
        "1|	(Ignored) Fitting notes \n" + 
        "2|	0   k12 k13... k1n \n" + 
        "3|	k21 0	k23... k2n \n" + 
        "	k31 k32 0  ... k3n \n" + 
        "	. \n" + 
        "	. \n" + 
        "	. \n" + 
        "n+1|	kn1 kn2	kn3... 0 \n" + 
        "n+2|	(# fully protected exchanging sites in state 1 - native) ('' state 2) ... ('' state n - 1 - last intermediate before " +
                "fully unfolded) \n" + 
        "n+3| 	(kbreathe) \n" + 
        "n+4| 	(Optimization type - 'Coordinate', 'Nelder_Mead', or 'ABCSMC') \n" + 
        "n+5|	If 'Coordinate' was chosen: \n" + 
        "				(max # iterations)  \n" + 
        "				(starting # choices per parameter to be fit)  \n" + 
        "				(percentage changes for coordinate search)  \n" + 
        "				(number of threads to run concurrently) \n" + 
        "		If 'Nelder_Mead' was chosen: \n" + 
        "				(fit tolerance)  \n" + 
        "				(max # iterations)  \n" + 
        "				(starting # choices per parameter to be fit)  \n" + 
        "				(percentage changes for initial simplex)  \n" + 
        "				(number of threads to run concurrently) \n" + 
        "		If 'ABCSMC' was chosen: \n" + 
        "				(prior type for rate constants - 'unif' or 'logunif') \n" + 
        "				(maximum number of intermediate distributions)   \n" + 
        "				(desired number of points in final distribution) \n" + 
        "				(proportion of variance to retain) \n" + 
        "				(initial epsilon cutoff) \n" + 
        "				(final epsilon cutoff) \n" + 
        "				(number of threads to run concurrently) \n" + 
        "--------------------------------- \n"
        );
    }
}
