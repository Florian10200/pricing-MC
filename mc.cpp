#include <iostream>
#include <math.h>
#include <random>
#include "mc.h"

using namespace std;



// FONCTIONS GENERALES --------------------------------------------------

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



//EUROPEAN OPTIONS --------------------------------------------------

europeanOption::europeanOption(double st, double sp, double vol, double mat, double tx) {
    this -> strike = st;
    this -> spot = sp;
    this -> volatility = vol;
    this -> maturity = mat;
    this -> txinteret = tx;}

double europeanOption::get_strike() {
    return strike;}

double europeanOption::get_spot() {
    return spot;}

double europeanOption::get_volatility() {
    return volatility;}

double europeanOption::get_maturity() {
    return maturity;}

double europeanOption::get_tx() {
    return txinteret;}


double europeanOption::pricing_european_call_BS() { //pricing european call par BS
    double T = maturity;
    double S0 = spot;
    double K = strike;
    double sigma = volatility;
    double r = txinteret;

    double d1 = (log(S0/K) + (r + pow(sigma, 2)/2.0)*T)/(sigma*sqrt(T));
    double d2 = d1 - sigma*sqrt(T);
    double Nd1 = normalCDF(d1);
    double Nd2 = normalCDF(d2);
    
    
    double prix_call_BS = S0 * Nd1 - K * exp(-r*T) * Nd2;

    cout << "prix du european call par BS est : ";
    
    return prix_call_BS;
}


double europeanOption::delta_european_call() {
    double T = maturity;
    double S0 = spot;
    double K = strike;
    double sigma = volatility;
    double r = txinteret;

    double d1 = (log(S0/K) + (r + pow(sigma, 2)/2.0)*T)/(sigma*sqrt(T));
    
    double Nd1 = normalCDF(d1);
    
    cout << "stratégie de réplication pour un call : delta = ";
    
    return Nd1;
}


double europeanOption::pricing_european_put_BS() { //pricing european put par BS
    double T = maturity;
    double S0 = spot;
    double K = strike;
    double sigma = volatility;
    double r = txinteret;

    double d1 = (log(S0/K) + (r + pow(sigma, 2)/2.0)*T)/(sigma*sqrt(T));
    double d2 = d1 - sigma*sqrt(T);
    double Nd1 = normalCDF(d1);
    double Nd2 = normalCDF(d2);
    
    
    double prix_put_BS = S0 * (Nd1 - 1) - K * exp(-r*T) * (Nd2 -1);

    cout << "prix du european put par BS est : ";
    
    return prix_put_BS;
}


double europeanOption::delta_european_put() {
    double T = maturity;
    double S0 = spot;
    double K = strike;
    double sigma = volatility;
    double r = txinteret;

    double d1 = (log(S0/K) + (r + pow(sigma, 2)/2.0)*T)/(sigma*sqrt(T));
    
    double Nd1 = normalCDF(d1);
    
    cout << "stratégie de réplication pour un put : delta = ";
    
    return Nd1 - 1;
}


double europeanOption::pricing_european_call_MC(int N) { //pricing european call par MC
    double T = maturity;
    double S0 = spot;
    double K = strike;
    double sigma = volatility;
    double r = txinteret;
    
    double ST; //prix du sous jacent à maturité
    double sumPayoff = 0.0;
    
    for (int i = 0; i < N; i++) {
        double z = generateN01();
        ST = S0*exp((r-pow(sigma, 2)/2.0)*T + sigma*sqrt(T)*z);
        sumPayoff = sumPayoff + max(ST - K, 0.0);
    }
    
    double prix = sumPayoff/N*exp(-r*T);
    
    cout << "prix du european call par MC est : ";
    
    return prix;
}


double europeanOption::pricing_european_put_MC(int N) { //pricing european put par MC
    double T = maturity;
    double S0 = spot;
    double K = strike;
    double sigma = volatility;
    double r = txinteret;

    double ST; //prix du sous jacent à maturité
    double sumPayoff = 0.0;
    
    for (int i = 0; i < N; i++) {
        double z = generateN01();
        ST = S0*exp((r-pow(sigma, 2)/2.0)*T + sigma*sqrt(T)*z);
        sumPayoff = sumPayoff + max(K - ST, 0.0);
    }
    
    double prix = sumPayoff/N*exp(-r*T);
    
    cout << "prix du european put par MC est : ";
    
    return prix;
}



//LOOKBACK (flottant) OPTIONS --------------------------------------------------

lookbackOption::lookbackOption(double sp, double vol, double mat, double tx) {
    this -> spot = sp;
    this -> volatility = vol;
    this -> maturity = mat;
    this -> txinteret = tx;}

double lookbackOption::get_spot() {
    return spot;}

double lookbackOption::get_volatility() {
    return volatility;}

double lookbackOption::get_maturity() {
    return maturity;}

double lookbackOption::get_tx() {
    return txinteret;}
    
    
double lookbackOption::pricing_lookback_call_MC(int N) { //pricing lookback call par MC
    double T = maturity;
    double S0 = spot;
    double sigma = volatility;
    double r = txinteret;
    
    double dt = 1.0 / 365.0; //un pas correspond à un jour, on actualise le prix du sous jacent chaque jour
    double S_jdt; //prix du sous jacent à j*dt
    double sumPayoff = 0.0;
    double minSt = 0.0;
    double payoff;
    
    for (int i = 0; i < N; i++) { //pour chaque trajectoire on calcul le prix minimum atteint
        
        double minSt = S0; //on réinitialise la valeur du prix min à t=0 pour chaque trajectoire
        for (int j = 0; j < T*365; ++j) {
            double z = generateN01();
            S_jdt = S0*exp((r-pow(sigma, 2)/2.0)*(j*dt) + sigma*sqrt(j*dt)*z);
            
            if (S_jdt < minSt) {
                minSt = S_jdt;
            }
        }
        
        payoff = max(S_jdt - minSt, 0.0); //max(prix à maturité - prix min, 0)
        sumPayoff = sumPayoff + payoff;
    }
    
    double prix = sumPayoff/N*exp(-r*T);
    
    cout << "prix du lookback call par MC est : ";
    
    return prix;
}


double lookbackOption::pricing_lookback_put_MC(int N) { //pricing lookbackput par MC
    double T = maturity;
    double S0 = spot;
    double sigma = volatility;
    double r = txinteret;
    
    double dt = 1.0 / 365.0; //un pas correspond à un jour, on actualise le prix du sous jacent chaque jour
    double S_jdt; //prix du sous jacent à j*dt
    double sumPayoff = 0.0;
    double maxSt = 0.0;
    double payoff;
    
    for (int i = 0; i < N; i++) { //pour chaque trajectoire on calcul le prix maximum atteint
        
        double minSt = S0; //on réinitialise la valeur du prix max à t=0 pour chaque trajectoire
        for (int j = 0; j < T*365; ++j) {
            double z = generateN01();
            S_jdt = S0*exp((r-pow(sigma, 2)/2.0)*(j*dt) + sigma*sqrt(j*dt)*z);
            
            if (S_jdt > maxSt) {
                maxSt = S_jdt;
            }
        }
        
        payoff = max(maxSt - S_jdt, 0.0); //max(prix max - prix à maturité, 0)
        sumPayoff = sumPayoff + payoff;
    }
    
    double prix = sumPayoff/N*exp(-r*T);
    
    cout << "prix du lookback put par MC est : ";
    
    return prix;
}
