#ifndef option_pricing_cpp
#define option_pricing_cpp

#include"option_pricing.h"
#include<random>
#include<cmath>
#include<iostream>

using namespace std;


// Fonctions générales --------------------------------------------------

double normalCDF(double x) { //cdf de la loi normale centrée réduite
   return 0.5 * erfc(-x * sqrt(0.5));
}

double generateN01() { //générateur de nombres aléatoires
    random_device rd; //initialisation avec une graine
    mt19937 generator(rd());

    normal_distribution<double> distribution(0.0, 1.0); //distribution normale avec une moyenne de 0 et un écart-type de 1

    double nb = distribution(generator);

    return nb;
}

// Constructeur de la classe option, european_option, lookback_option, asian option, barrier option

option::option(double S_t, double K, double vol, double mat, double tx) {
    this -> spot = S_t;
    this -> strike = K;
    this -> volatility = vol;
    this -> maturity = mat;
    this -> txinteret = tx;
};

european_option::european_option(double S_t, double K, double vol, double mat, double tx):option(S_t, K, vol, mat, tx){};

lookback_option::lookback_option(double S_t, double K, double vol, double mat, double tx):option(S_t, K, vol, mat, tx){};

asian_option::asian_option(double S_t, double K, double vol, double mat, double tx):option(S_t, K, vol, mat, tx){};

barrier_option::barrier_option(double S_t, double K, double vol, double mat, double tx, double Q, double b):option(S_t, K, vol, mat, tx){
    this -> q=Q;
    this -> B=b;
};

// Fonction permettant de récupérer les valeurs des variables de la classe option

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



// Pricing d'options européenne

double european_option::pricing_european_call_BS() { //pricing european call par BS
    double T = get_maturity();
    double S0 = get_spot();
    double K = get_strike();
    double sigma = get_volatility();
    double r = get_tx();

    double d1 = (log(S0/K) + (r + pow(sigma, 2)/2.0)*T)/(sigma*sqrt(T));
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


double european_option::pricing_european_put_BS() { //pricing european put par BS
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

double european_option::pricing_european_call_MC(int N) { //pricing european call par MC
    double T = get_maturity();
    double S0 = get_spot();
    double K = get_strike();
    double sigma = get_volatility();
    double r = get_tx();
    
    double ST; //prix du sous jacent à maturité
    double sumPayoff = 0.0;
    
    for (int i = 0; i < N; i++) {
        double z = generateN01();
        ST = S0*exp((r-pow(sigma, 2)/2.0)*T + sigma*sqrt(T)*z);
        sumPayoff = sumPayoff + max(ST - K, 0.0);
    }
    
    double prix = sumPayoff/N*exp(-r*T);

    return prix;
}


double european_option::pricing_european_put_MC(int N) { //pricing european put par MC
    double T = get_maturity();
    double S0 = get_spot();
    double K = get_strike();
    double sigma = get_volatility();
    double r = get_tx();

    double ST; //prix du sous jacent à maturité
    double sumPayoff = 0.0;
    
    for (int i = 0; i < N; i++) {
        double z = generateN01();
        ST = S0*exp((r-pow(sigma, 2)/2.0)*T + sigma*sqrt(T)*z);
        sumPayoff = sumPayoff + max(K - ST, 0.0);
    }
    
    double prix = sumPayoff/N*exp(-r*T);

    return prix;
}


// Lookback option pricing avec Monte-Carlo


double lookback_option::lookback_pricingCallMC(int N) { //pricing lookback call par MC
    double T = get_maturity();
    double S0 = get_spot();
    double sigma = get_volatility();
    double r = get_tx();
    
    double dt = 1.0 / 365.0; //un pas correspond à un jour, on actualise le prix du sous jacent chaque jour
    double sumPayoff = 0.0;
    double minSt;
    double payoff;
    
    for (int i = 0; i < N; i++) { //pour chaque trajectoire on calcul le prix minimum atteint
    
        double S_jdt = S0; //prix du sous jacent à j*dt
        minSt = S0; //on réinitialise la valeur du prix min à t=0 pour chaque trajectoire
        for (int j = 0; j < T*365; ++j) {
            double z = generateN01();
            S_jdt *= exp((r-pow(sigma, 2)/2.0)*(dt) + sigma*sqrt(dt)*z);
            
            if (S_jdt < minSt) {
                minSt = S_jdt;
            }
        }
        
        payoff = max(S_jdt - minSt, 0.0); //max(prix à maturité - prix min, 0)
        sumPayoff = sumPayoff + payoff;
    }
    
    double prix = sumPayoff/N*exp(-r*T);

    return prix;
}


double lookback_option::lookback_pricingPutMC(int N) { //pricing lookbackput par MC
    double T = get_maturity();
    double S0 = get_spot();
    double sigma = get_volatility();
    double r = get_tx();
    
    double dt = 1.0 / 365.0; //un pas correspond à un jour, on actualise le prix du sous jacent chaque jour
    double sumPayoff = 0.0;
    double maxSt;
    double payoff;
    
    for (int i = 0; i < N; i++) { //pour chaque trajectoire on calcul le prix maximum atteint
        
        double S_jdt = S0; //prix du sous jacent à j*dt
        maxSt = S0; //on réinitialise la valeur du prix max à t=0 pour chaque trajectoire
        for (int j = 0; j < T*365; ++j) {
            double z = generateN01();
            S_jdt *= exp((r-pow(sigma, 2)/2.0)*(dt) + sigma*sqrt(dt)*z);
            
            if (S_jdt > maxSt) {
                maxSt = S_jdt;
            }
        }
        
        payoff = max(maxSt - S_jdt, 0.0); //max(prix max - prix à maturité, 0)
        sumPayoff = sumPayoff + payoff;
    }
    
    double prix = sumPayoff/N*exp(-r*T);
 
    return prix;
}


// Asian option pricing avec Monte-Carlo --------------------------------------------------


double asian_option::asian_pricingCallMC(int N) { //pricing lookback call par MC
    double T = get_maturity();
    double K = get_strike();
    double S0 = get_spot();
    double sigma = get_volatility();
    double r = get_tx();
    
    double dt = 1.0 / 365.0; //un pas correspond à un jour, on actualise le prix du sous jacent chaque jour
    double sumPayoff = 0.0;
    double payoff;
    
    for (int i = 0; i < N; i++) { //pour chaque trajectoire on calcul le prix moyen de l'actif sous-jacent
        
        double sum_price_traj = 0;
        double S_jdt = S0; //prix du sous jacent à j*dt initialisé à t=0
        
        for (int j = 0; j < T*365; ++j) {
            double z = generateN01();
            S_jdt *= exp((r-pow(sigma, 2)/2.0)*(dt) + sigma*sqrt(dt)*z);
            sum_price_traj += S_jdt;
        }
        double average_price_traj = sum_price_traj/(T*365);
        payoff = max(average_price_traj - K, 0.0); //max(prix moyen de la trajectoire - prix d'exercice, 0)
        sumPayoff = sumPayoff + payoff;
    }
    
    double prix = sumPayoff/N*exp(-r*T); // Moyenne des payoffs + actualisation à la date présente

    return prix;
}



double asian_option::asian_pricingPutMC(int N) { //pricing lookback put par MC
    double T = get_maturity();
    double K = get_strike();
    double S0 = get_spot();
    double sigma = get_volatility();
    double r = get_tx();
    
    double dt = 1.0 / 365.0; //un pas correspond à un jour, on actualise le prix du sous jacent chaque jour
    double sumPayoff = 0.0;
    double payoff;

    
    for (int i = 0; i < N; i++) { //pour chaque trajectoire on calcul le prix moyen de l'actif sous-jacent
        
        double sum_price_traj = 0;   
        double S_jdt = S0; //prix du sous jacent à j*dt initialisé à t=0
        
        for (int j = 0; j < T*365; ++j) {
            double z = generateN01();
            S_jdt *= exp((r-pow(sigma, 2)/2.0)*(dt) + sigma*sqrt(dt)*z);
            sum_price_traj += S_jdt;
        }
        double average_price_traj = sum_price_traj/(T*365);
        payoff = max(K - average_price_traj, 0.0); //max(prix moyen de la trajectoire - prix d'exercice, 0)
        sumPayoff = sumPayoff + payoff;
    }
    
    double prix = sumPayoff/N*exp(-r*T); // Moyenne des payoffs + actualisation à la date présente

    return prix;
}

// For Barrier option ------------------------------------------------

double barrier_option::get_B(){
    return B;
}

double barrier_option::get_q(){
    return q;
}




#endif // option_pricing_cpp
