#include <iostream>
#include <math.h>
#include <random>
#include "mc.h"

using namespace std;

double normalCDF(double x) { //cdf de la loi normale centrée réduite
   return 0.5 * erfc(-x * sqrt(0.5));
}

double generateN01() {
    //initialisation du générateur de nombres aléatoires avec une graine
    random_device rd;
    mt19937 generator(rd());

    //définition de la distribution normale avec une moyenne de 0 et un écart-type de 1
    normal_distribution<double> distribution(0.0, 1.0);

    double nb = distribution(generator);

    return nb;
}

Option::Option(double st, double sp, double vol, double mat, double tx) {
    this -> strike = st;
    this -> spot = sp;
    this -> volatility = vol;
    this -> maturity = mat;
    this -> txinteret = tx;}

double Option::get_strike() {
    return strike;}

double Option::get_spot() {
    return spot;}

double Option::get_volatility() {
    return volatility;}

double Option::get_maturity() {
    return maturity;}

double Option::get_tx() {
    return txinteret;}


double Option::pricing_european_call_BS() { //pricing european call par BS
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

    cout << "prix du european call par BS est : " << prix_call_BS << endl;
    
    cout << "delta_call : " << Nd1 << endl;
    
    return 0;
}

double Option::pricing_european_put_BS() { //pricing european put par BS
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

    cout << "prix du european put par BS est : " << prix_put_BS << endl;

    cout << "delta_put : " << Nd1 -1 << endl;
    
    
    return 0;
}

double Option::pricing_european_call_MC(int N) { //pricing european call par MC
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
    
    cout << "prix du european call par MC est : " << prix << endl;
    
    return 0;
}


double Option::pricing_european_put_MC(int N) { //pricing european put par MC
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
    
    cout << "prix du european put par MC est : " << prix << endl;
    
    return 0;
}


double Option::pricing_american_call_MC(int N) { //pricing american call par MC
    double T = maturity;
    double S0 = spot;
    double K = strike;
    double sigma = volatility;
    double r = txinteret;
    
    double dt = 1.0 / 365.0; //un pas correspond à un jour, on actualise le prix du sous jacent chaque jour
    double S_jdt; //prix du sous jacent à j*dt
    double sumPayoff = 0.0;
    double maxPayoff = 0.0;
    double payoff;
    
    for (int i = 0; i < N; ++i) { //pour chaque trajectoire on calcul le spot maximum
        
        double maxSpot = S0; //on réinitialise la valeur du maxSpot à 0 pour chaque trajectoire
        for (int j = 0; j < T*365; ++j) {
            double z = generateN01();
            S_jdt = S0*exp((r-pow(sigma, 2)/2.0)*(j*dt) + sigma*sqrt(j*dt)*z);
            
            if (S_jdt > maxSpot) {
                maxSpot = S_jdt;
            }
        }
        
        payoff = max(S_jdt - K, 0.0);
        sumPayoff = sumPayoff + payoff;
    }
    
    double prix = sumPayoff/N*exp(-r*T);
    
    cout << "prix du american call par MC est : " << prix << endl;
    
    return 0;
}


double Option::pricing_american_put_MC(int N) { //pricing american put par MC
    double T = maturity;
    double S0 = spot;
    double K = strike;
    double sigma = volatility;
    double r = txinteret;
    
    double dt = 1.0 / 365.0; //un pas correspond à un jour, on actualise le prix du sous jacent chaque jour
    double S_jdt; //prix du sous jacent à j*dt
    double sumPayoff = 0.0;
    double maxPayoff = 0.0;
    double payoff;
    
    for (int i = 0; i < N; ++i) { //pour chaque trajectoire on calcul le spot minimum
        
        double minSpot = S0; //on réinitialise la valeur du minSpot à 0 pour chaque trajectoire
        for (int j = 0; j < T*365; ++j) {
            double z = generateN01();
            S_jdt = S0*exp((r-pow(sigma, 2)/2.0)*(j*dt) + sigma*sqrt(j*dt)*z);
            
            if (S_jdt < minSpot) {
                minSpot = S_jdt;
            }
        }
        
        payoff = max(K - S_jdt, 0.0);
        sumPayoff = sumPayoff + payoff;
    }
    
    double prix = sumPayoff/N*exp(-r*T);
    
    cout << "prix du american put par MC est : " << prix << endl;    
    return 0;
}
