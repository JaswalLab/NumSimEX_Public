package numericalsimulations;

import java.util.*;

class Protein{

    ArrayList<Site> sites;  //Arraylist of Site objects that can exchange
    ArrayList<Double> statechangetimes = new ArrayList<>(); //Stores a list of when the protein changes folding state
    ArrayList<Integer> statehistory = new ArrayList<>(); //Stores a list of the folding states a protein has been in
    private boolean fullyexchanged = false; //Whether all the sites on the protein have already exchanged
    private final double[][] finalmatrix; //Matrix containing all information about kex (nstate rows by nex columns)
    private final int nex; //Number of exchangers a protein has
    
    
    //Constructor requiring matrix with all information about kex for each site at every state
    public Protein (double[][] infotable, Random globalrand){ 
        finalmatrix = infotable;
        nex = finalmatrix[0].length;
        sites = new ArrayList<>(nex);
        initializeSites(globalrand);
    }
    
    //Adds a new Site object to the stored arraylist for each of the exchanging sites we have
    private void initializeSites(Random globalrand){ 
        for(int i = 0; i < nex; i++){
            Site newsite = new Site(getColumn(i, finalmatrix), globalrand); //Sim object gets passed through to allow for a global seed
            sites.add(newsite);
        }
    }
    
    private double[] getColumn(int columnnum, double[][] from){ //Extracts kth column from an i*j matrix
        int nrows = from.length;
        double[] toreturn = new double[nrows];
        for(int i = 0; i < nrows; i++){
            toreturn[i] = from[i][columnnum];
        }
        return toreturn;
    }
    
    //Documents a transition in states with the new state and the time of the event
    public void setStateChange (double time, int stateTo){  
        statechangetimes.add(time);
        statehistory.add(stateTo);
    }
    
    public Site getSite(int whichSite){  //Returns the specified nth Site object
        return sites.get(whichSite);
    }
    
    public int getState(){ //Used in Sim to get the most recent state when folding analysis is conducted
        return statehistory.get(statehistory.size()-1);
    }
    
    public void findExchangers(double[] timepoints){ //Finds exchange times for each site in this protein
        setCheckingPoints(timepoints);
        double time = 0;
        double tau;
        int laststate;
        for(int i = 1; i < statehistory.size(); i++){
            if(!isFullyExchanged()){
                laststate = statehistory.get(i-1);
                tau = statechangetimes.get(i) - statechangetimes.get(i-1);
                //Iterate across the sites to randomly test exchanging for those who haven't already
                for(int j = 0; j < sites.size(); j++){ 
                    //The last parameter gives the state the protein last switched into
                    sites.get(j).exchangeCheck(time, tau, laststate); 
                }
                time = time + tau;
            }
            else{
                return;
            }
        }
    }
    
    //Combines the two lists of timepoints (one input from user, one generated as the protein changes folding states) 
    //and the lists' associated folding states
    public void setCheckingPoints(double[] inputtps) { 
        for(int i = 0; i < inputtps.length; i++){
            int simtpsindex = 0;
            while(statechangetimes.get(simtpsindex) < inputtps[i]){
                simtpsindex++;
            }
            statehistory.add(simtpsindex, statehistory.get(simtpsindex - 1));
            statechangetimes.add(simtpsindex, inputtps[i]);
        }
    }
    
    public ArrayList<Integer> getStateHistory () { //Returns arraylist of states in which this protein has been
        return statehistory;
    }
    
    public ArrayList<Double> getStateChangeTimes() { //Returns arraylist of the times at which the protein changed states
        return statechangetimes;
    }
    
    public boolean isFullyExchanged(){  //Checks if the protein has been fully exchanged yet to increase efficiency
        if(fullyexchanged){
            return true;
        }
        for(int i = 0; i < sites.size(); i++){
            if(!sites.get(i).getExchanged()){
                return false;
            }
        }
        fullyexchanged = true;
        return fullyexchanged;
    }
    
    //Used to create conditional pmf of number of sites exchanged given a particular time
    public int getNumberSitesExchanged (double time){  
        int sum = 0;
        for(int i = 0; i < nex; i++){
            if(sites.get(i).getExchanged() && sites.get(i).getExchangedTime() < time){
                sum = sum + 1;
            }
        }
        return sum;
    }
}