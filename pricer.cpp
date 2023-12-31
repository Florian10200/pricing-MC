#ifndef pricer_cpp
#define pricer_cpp

#include"pricer.h"
#include<random>
#include<cmath>
#include<iostream>

using namespace std;


// General functions --------------------------------------------------------------------------------------------------

double normalCDF(double x) { //cdf of N(0,1)
   return 0.5 * erfc(-x * sqrt(0.5));
}

double generateN01() { //random number generator
    random_device rd; //initialisation with a seed
    mt19937 generator(rd());

    normal_distribution<double> distribution(0.0, 1.0); //normal distribution N(0,1)

    double nb = distribution(generator);

    return nb;
}

// Constructor of the class option -------------------------------------------------------------------------------------------

option::option(double S_t, double K, double vol, double mat, double tx) {
    this -> spot = S_t;
    this -> strike = K;
    this -> volatility = vol;
    this -> maturity = mat;
    this -> txinteret = tx;
};



// Functions to get the private values of the class option -------------------------------------------------------------

double option::get_strike() {
    return strike;
};

double option::get_spot() {
    return spot;
};

double option::get_volatility() {
    return volatility;
};

double option::get_maturity() {
    return maturity;
};

double option::get_tx() {
    return txinteret;
};



// European option pricing ------------------------------------------------------------------------------------------------

european_option::european_option(double S_t, double K, double vol, double mat, double tx):option(S_t, K, vol, mat, tx){}; //Constructor of the european option class

double european_option::pricing_european_call_BS() { //pricing european call with BS
    double T = get_maturity();
    double S0 = get_spot();
    double K = get_strike();
    double sigma = get_volatility();
    double r = get_tx();

    double d1 = (log(S0/K) + (r + pow(sigma, 2)/2.0)*T)/(sigma*sqrt(T)); //General formula by Black-Scholes model
    double d2 = d1 - sigma*sqrt(T);
    double Nd1 = normalCDF(d1);
    double Nd2 = normalCDF(d2);
    
    
    double prix_call_BS = S0 * Nd1 - K * exp(-r*T) * Nd2;

    return prix_call_BS;
};


double european_option::delta_european_call() {
    double T = get_maturity();
    double S0 = get_spot();
    double K = get_strike();
    double sigma = get_volatility();
    double r = get_tx();

    double d1 = (log(S0/K) + (r + pow(sigma, 2)/2.0)*T)/(sigma*sqrt(T));
    
    double Nd1 = normalCDF(d1);
  
    return Nd1;
};


double european_option::pricing_european_put_BS() { //pricing european put with BS
    double T = get_maturity();
    double S0 = get_spot();
    double K = get_strike();
    double sigma = get_volatility();
    double r = get_tx();

    double d1 = (log(S0/K) + (r + pow(sigma, 2)/2.0)*T)/(sigma*sqrt(T));
    double d2 = d1 - sigma*sqrt(T);
    double Nd1 = normalCDF(d1);
    double Nd2 = normalCDF(d2);
    
    
    double prix_put_BS = S0 * (Nd1 - 1) - K * exp(-r*T) * (Nd2 -1);
   
    return prix_put_BS;
};


double european_option::delta_european_put() {
    double T = get_maturity();
    double S0 = get_spot();
    double K = get_strike();
    double sigma = get_volatility();
    double r = get_tx();

    double d1 = (log(S0/K) + (r + pow(sigma, 2)/2.0)*T)/(sigma*sqrt(T));
    
    double Nd1 = normalCDF(d1);

    return Nd1 - 1;
    
};

double european_option::pricing_european_call_MC(int N) { //pricing european call with MC
    double T = get_maturity();
    double S0 = get_spot();
    double K = get_strike();
    double sigma = get_volatility();
    double r = get_tx();
    
    double ST; //price of the underlying asset at maturity
    double sumPayoff = 0.0; // We initialize a variable sumPayoff at 0 and it will contain the payoff of each path at maturity
    
    for (int i = 0; i < N; i++) { // For each path
        double z = generateN01();
        ST = S0*exp((r-pow(sigma, 2)/2.0)*T + sigma*sqrt(T)*z); // We estimate the spot price at maturity by supossing a lognormal distribution
        sumPayoff = sumPayoff + max(ST - K, 0.0); // The max refers to a payoff for a call option
    }
    
    double prix = sumPayoff/N*exp(-r*T); // We take the average and actualize the value with the continous discount factor

    return prix;
}


double european_option::pricing_european_put_MC(int N) { //pricing european put with MC : the same as for a call, just the payoff changes
    double T = get_maturity();
    double S0 = get_spot();
    double K = get_strike();
    double sigma = get_volatility();
    double r = get_tx();

    double ST; //price of the underlying asset at maturity
    double sumPayoff = 0.0;
    
    for (int i = 0; i < N; i++) {
        double z = generateN01();
        ST = S0*exp((r-pow(sigma, 2)/2.0)*T + sigma*sqrt(T)*z);
        sumPayoff = sumPayoff + max(K - ST, 0.0); // Payoff for a put option
    }
    
    double prix = sumPayoff/N*exp(-r*T);

    return prix;
}


// Lookback option pricing with Monte-Carlo -------------------------------------------------------------------------------

lookback_option::lookback_option(double S_t, double K, double vol, double mat, double tx):option(S_t, K, vol, mat, tx){}; //Constructor of the lookback option class

double lookback_option::lookback_pricingCallMC(int N) { //pricing lookback call with MC
    double T = get_maturity();
    double S0 = get_spot();
    double sigma = get_volatility();
    double r = get_tx();
    
    double dt = 1.0 / 365.0; //one step corresponds to a day : we actualize the value each day
    double sumPayoff = 0.0;
    double minSt; // We define the value of the minimum spot price of all the paths
    double payoff;
    
    for (int i = 0; i < N; i++) { //For each path, we calculate the minimum price that is taken by the asset
    
        double S_jdt = S0; //price of the underlying asset at time dt*j
        minSt = S0; //we reinitialise the value of the min at t=0 for each path
        for (int j = 0; j < T*365; ++j) { //for each path
            double z = generateN01();
            S_jdt *= exp((r-pow(sigma, 2)/2.0)*(dt) + sigma*sqrt(dt)*z); //The new price is just the former one multiplied by an exponential
            
            if (S_jdt < minSt) { // Condition to get the minimum
                minSt = S_jdt;
            }
        }
        
        payoff = max(S_jdt - minSt, 0.0); //max(price at maturity - minimum price, 0)
        sumPayoff = sumPayoff + payoff;
    }
    
    double prix = sumPayoff/N*exp(-r*T); // Average of the payoff and actualization of value with the discount factor

    return prix;
}


double lookback_option::lookback_pricingPutMC(int N) { //pricing lookbackput option with MC : quite similar than a call
    double T = get_maturity();
    double S0 = get_spot();
    double sigma = get_volatility();
    double r = get_tx();
    
    double dt = 1.0 / 365.0; 
    double sumPayoff = 0.0;
    double maxSt;
    double payoff;
    
    for (int i = 0; i < N; i++) { // For each path, we calculate the max price taken by the asset
        
        double S_jdt = S0; 
        maxSt = S0; 
        for (int j = 0; j < T*365; ++j) {
            double z = generateN01();
            S_jdt *= exp((r-pow(sigma, 2)/2.0)*(dt) + sigma*sqrt(dt)*z);
            
            if (S_jdt > maxSt) { // Condition to get the max
                maxSt = S_jdt;
            }
        }
        
        payoff = max(maxSt - S_jdt, 0.0); //max(maximum price - prix at maturity, 0)
        sumPayoff = sumPayoff + payoff;
    }
    
    double prix = sumPayoff/N*exp(-r*T);
 
    return prix;
}


// Asian option pricing with Monte-Carlo -----------------------------------------------------------------------------------

asian_option::asian_option(double S_t, double K, double vol, double mat, double tx):option(S_t, K, vol, mat, tx){}; //Constructor of the asian option class

double asian_option::asian_pricingCallMC(int N) { //pricing asian call with MC
    double T = get_maturity();
    double K = get_strike();
    double S0 = get_spot();
    double sigma = get_volatility();
    double r = get_tx();
    
    double dt = 1.0 / 365.0; // The same than a lookback option
    double sumPayoff = 0.0;
    double payoff;
    
    for (int i = 0; i < N; i++) { //For each path, we calculate the average price of the underlying asdet
        
        double sum_price_traj = 0;
        double S_jdt = S0; 
        
        for (int j = 0; j < T*365; ++j) {
            double z = generateN01();
            S_jdt *= exp((r-pow(sigma, 2)/2.0)*(dt) + sigma*sqrt(dt)*z); //Same modelisation of the underlying asset than in the lookback option
            sum_price_traj += S_jdt;
        }
        double average_price_traj = sum_price_traj/(T*365); //We calculate the average
        payoff = max(average_price_traj - K, 0.0); //max(average price  - strike price, 0)
        sumPayoff = sumPayoff + payoff;
    }
    
    double prix = sumPayoff/N*exp(-r*T); // Average of the payoffs + actualization at the present date

    return prix;
}



double asian_option::asian_pricingPutMC(int N) { //pricing asian put par MC
    double T = get_maturity();
    double K = get_strike();
    double S0 = get_spot();
    double sigma = get_volatility();
    double r = get_tx();
    
    double dt = 1.0 / 365.0; 
    double sumPayoff = 0.0;
    double payoff;

    
    for (int i = 0; i < N; i++) { 
        
        double sum_price_traj = 0;   
        double S_jdt = S0; 
        
        for (int j = 0; j < T*365; ++j) {
            double z = generateN01();
            S_jdt *= exp((r-pow(sigma, 2)/2.0)*(dt) + sigma*sqrt(dt)*z);
            sum_price_traj += S_jdt;
        }
        double average_price_traj = sum_price_traj/(T*365);
        payoff = max(K - average_price_traj, 0.0); //max(exercice price - average price, 0)
        sumPayoff = sumPayoff + payoff;
    }
    
    double prix = sumPayoff/N*exp(-r*T); // average of the payoffs + Actualization at the present date

    return prix;
}


// For Barrier option ----------------------------------------------------------------------------------------------------------

//Constructor of the barrier option class with 2 additional arguments
barrier_option::barrier_option(double S_t, double K, double vol, double mat, double tx, double Q, double b):option(S_t, K, vol, mat, tx){
    this -> q=Q;
    this -> B=b;
};

double barrier_option::get_B(){
    return B;
}

double barrier_option::get_q(){
    return q;
}


double barrier_option::PriceCall(){
    double T = get_maturity();
    double S0 = get_spot();
    double K = get_strike();
    double sig = get_volatility();
    double r = get_tx();
    double q = get_q();
    double B = get_B();
    double d1 = (log(S0/K)+(r-q+pow(sig,2)/2)*T)/(sig*sqrt(T));
    double d2 = d1-sig*sqrt(T);
    double Nd1=normalCDF(d1);
    double Nd2=normalCDF(d2);
    double C = S0*exp(-q*T)*Nd1-K*Nd2*exp(-r*T);
    return C;
}


double barrier_option::PricePut(){
    double T = get_maturity();
    double S0 = get_spot();
    double K = get_strike();
    double sig = get_volatility();
    double r = get_tx();
    double q = get_q();
    double B = get_B();
    double d1 = (log(S0/K)+(r-q+pow(sig,2)/2)*T)/(sig*sqrt(T));
    double d2 = d1-sig*sqrt(T);
    double Nd1=normalCDF(-d1);
    double Nd2=normalCDF(-d2);
    double P = -(S0*exp(-q*T)*Nd1-K*Nd2*exp(-r*T));
    return P;
}

double barrier_option::Call_UpAndIn_UpAndOut(){
    double T = get_maturity();
    double S0 = get_spot();
    double K = get_strike();
    double sig = get_volatility();
    double r = get_tx();
    double q = get_q();
    double B = get_B();
    double lam = (r-q+pow(sig,2)/2)/pow(sig,2);
    double C = PriceCall();
    double x1 = (log(S0/B)+lam*pow(sig,2)*T)/(sig*sqrt(T));
    double y = log(pow(B,2)/(S0*K))/(sig*sqrt(T))+lam*sig*sqrt(T);
    double y1 = (B/log(S0)+lam*pow(sig,2)*T)/(sig*sqrt(T));
    double Ny=normalCDF(-y);
    double Nyy=normalCDF(-y+sig*sqrt(T));
    double Ny1=normalCDF(-y1);
    double Ny11=normalCDF(-y1+sig*sqrt(T));
    double Nx1 = normalCDF(x1);
    double Nx11=normalCDF(x1-sig*sqrt(T));
    if (B > K){
        double cui=S0*exp(-q*T)*Nx1-Nx11*K*exp(-r*T)-S0*exp(-q*T)*pow(B/S0,2*lam)*(Ny-Ny1)+(Nyy-Ny11)*K*pow(B/S0,2*lam-2)*exp(-r*T);
        cout<<"The value of the CUI is "<< cui <<endl;
        cout<<"The value of the CUO is "<< (C-cui)<<endl;
    }
    else {
        cout<<"The value of the CUI is "<< C <<endl;
        cout<<"The value of the CUO is "<< 0.0<<endl;
    }
    return 0.0;
}

double barrier_option::Call_DownAndIn_DownAndOut(){
    double T = get_maturity();
    double S0 = get_spot();
    double K = get_strike();
    double sig = get_volatility();
    double r = get_tx();
    double q = get_q();
    double B = get_B();
    double lam = (r-q+pow(sig,2)/2)/pow(sig,2);
    double y = log(pow(B,2)/(S0*K))/(sig*sqrt(T))+lam*sig*sqrt(T);
    double Ny=normalCDF(y);
    double Ny2=normalCDF(y-sig*sqrt(T));
    double C = PriceCall();
    double x1 = (log(S0/B)+lam*pow(sig,2)*T)/(sig*sqrt(T));
    double y1 = (B/log(S0)+lam*pow(sig,2)*T)/(sig*sqrt(T));
    double Nx1 = normalCDF(x1);
    double Nx11=normalCDF(x1-sig*sqrt(T));
    double Ny1 = normalCDF(y1);
    double Ny11=normalCDF(y1-sig*sqrt(T));
    if (B <= K){
        double cdi= S0*exp(-q*T)*pow(B/S0,2*lam)*Ny-Ny2*K*pow(B/S0,2*lam-2)*exp(-r*T);
        cout<<"The value of the CDI is "<< cdi <<endl;
        cout<<"The value of the CDO is "<<(C-cdi)<<endl;
    }
    else {
        double cdo=S0*exp(-q*T)*Nx1-Nx11*K*exp(-r*T)-S0*exp(-q*T)*pow(B/S0,2*lam)*Ny1+Ny11*K*pow(B/S0,2*lam-2)*exp(-r*T);
        cout<<"The value of the CDI is "<< (C-cdo) <<endl;
        cout<<"The value of the CDO is "<< cdo <<endl;
    }
    return 0.0;

}


double barrier_option::Put_UpAndIn_UpAndOut(){
    double T = get_maturity();
    double S0 = get_spot();
    double K = get_strike();
    double sig = get_volatility();
    double r = get_tx();
    double q = get_q();
    double B = get_B();
    double lam = (r-q+pow(sig,2)/2)/pow(sig,2);
    double y = log(pow(B,2)/(S0*K))/(sig*sqrt(T))+lam*sig*sqrt(T);
    double Ny = normalCDF(-y);
    double Nyy = normalCDF(-y+sig*sqrt(T));
    double x1 = (log(S0/B)+lam*pow(sig,2)*T)/(sig*sqrt(T));
    double y1 = (B/log(S0)+lam*pow(sig,2)*T)/(sig*sqrt(T));
    double Nx1 = normalCDF(-x1);
    double Ny1 = normalCDF(-y1);
    double Nx11=normalCDF(-x1+sig*sqrt(T));
    double Ny11=normalCDF(-y1+sig*sqrt(T));
    double P = PricePut();
    if (B >= K){
        double pui = S0*(-1)*exp(-q*T)*pow(B/S0,2*lam)*Ny+Nyy*K*pow(B/S0,2*lam-2)*exp(-r*T);
        cout<<"The value of the PUI is "<< pui <<endl;
        cout<<"The value of the PUO is "<< (P-pui)<<endl;
    }
    else{
        double puo = S0*(-1)*exp(-q*T)*Nx1+Nx11*K*exp(-r*T)+S0*exp(-q*T)*pow(B/S0,2*lam)*Ny1-Ny11*K*pow(B/S0,2*lam-2)*exp(-r*T);
        cout<<"The value of the PUI is "<< (P-puo) <<endl;
        cout<<"The value of the PUO is "<< puo <<endl;
    }
    return 0.0;
}

double barrier_option::Put_DownAndIn_DownAndOut(){
    double T = get_maturity();
    double S0 = get_spot();
    double K = get_strike();
    double sig = get_volatility();
    double r = get_tx();
    double q = get_q();
    double B = get_B();
    double lam = (r-q+pow(sig,2)/2)/pow(sig,2);
    double P = PricePut();
    double x1 = (log(S0/B)+lam*pow(sig,2)*T)/(sig*sqrt(T));
    double y1 = (B/log(S0)+lam*pow(sig,2)*T)/(sig*sqrt(T));
    double y = log(pow(B,2)/(S0*K))/(sig*sqrt(T))+lam*sig*sqrt(T);
    double Nx1=normalCDF(-x1);
    double Nx11 = normalCDF(-x1+sig*sqrt(T));
    double Ny = normalCDF(y);
    double Ny1 = normalCDF(y1);
    double Nyy = normalCDF(y-sig*sqrt(T));
    double Ny11 = normalCDF(y1-sig*sqrt(T));
    if (B <= K){
        double pdi = S0*(-1)*exp(-q*T)*Nx1+Nx11*K*exp(-r*T)+S0*exp(-q*T)*pow(B/S0,2*lam)*(Ny-Ny1)-(Nyy-Ny11)*K*pow(B/S0,2*lam-2)*exp(-r*T);
        cout<<"The value of the PDI is "<< pdi <<endl;
        cout<<"The value of the PDO is "<< (P-pdi) <<endl;
    }
    else{
        cout<<"The value of the PDI is "<< P <<endl;
        cout<<"The value of the PDO is "<< 0.0 <<endl;
    }
    return 0.0;
}




#endif // pricer_cpp
