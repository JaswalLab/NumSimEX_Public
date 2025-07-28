package numericalsimulations;

import java.util.Random;

class Site{

    private boolean exchanged = false; //Whether the site has exchanged
    private double exchangetime; //Approximate time at which the site has exchanged
    private final double[] kex; //Kex at each state, incorporating both kint and kbreathe as necessary
    private Random randomgenerator; //Random object from Sim object in which this exchanger exists (Used to implement global seed)
    
    
    public Site(double[] kexs, Random globalrand){ //Constructor
        kex = kexs;
        randomgenerator = globalrand;
    }
    
    public boolean getExchanged(){  //Return whether or not the site has exchanged
        return exchanged;
    }
    
    public void setExchanged(){  //Set the site as exchanged
        exchanged = true;
    }
    
    public void setExchangedTime(double time){  //Set approximate exchange time
        exchangetime = time;
    }
    
    public double getExchangedTime(){  //Return approximate exchange time
        return exchangetime;
    }
    
    //Over a particular length of time tau, checks whether a amide site exchanges
    public void exchangeCheck(double time, double tau, int currentstate){ 
        if(exchanged){ //If the site has already exchanged, nothing is calculated
            return;
        }
        double draw = runif();
        //Second part of expression in 'if' statement gives the probability of exchanging
        if(draw < (1-Math.exp((-1) * kex[currentstate] * tau))){ 
            //Approximation of exchange time as occurring halfway between the last time it was unexchanged
            //and the first time it was exchanged
            setExchangedTime(time+tau/2); 
            setExchanged();
        }
    }
    
    public double runif(){
        return randomgenerator.nextDouble();
    }
}