#ifndef option_pricing_h
#define option_pricing_h
#include<iostream>

using namespace std;

double normalCDF(double x);

double generateN01(); 


class option {
    private:
        double spot; //in $
        double strike; //in $
        double volatility; //in %
        double maturity; //in years
        double txinteret; //in %
    public:
        option(double S_t, double K, double vol, double maturity, double tx); //constructeur
        double get_strike(); //recupere la valeur du strike
        double get_spot(); //recupere la valeur du spot
        double get_volatility(); //recupere la valeur de la volatilite
        double get_maturity(); //recupere la valeur de la maturite
        double get_tx(); //recupere la valeur du tx d'interet
};        

// European option Pricing par les formules analytiques et par Monte

class european_option : public option {
    public:
        european_option(double S_t, double K, double vol, double maturity, double tx);
        
        double pricing_european_call_BS(); //pricing european call par BS
        double pricing_european_put_BS(); //pricing european put par BS
    
        double delta_european_call(); //calcul du delta pour un european call
        double delta_european_put(); //calcul du delta pour un european put
        
        double pricing_european_call_MC(int N); //pricing european call par MC, N le nombre de trajectoires
        double pricing_european_put_MC(int N); //pricing european put par MC, N le nombre de trajectoires
};

// For a lookback option avec Monte-Carlo


class lookback_option : public option {
    public:
        lookback_option(double S_t, double K, double vol, double maturity, double tx);
        
        double lookback_pricingCallMC(int N); //pricing lookback call par MC, N le nombre de trajectoires
        double lookback_pricingPutMC(int N); //pricing lookback put par MC, N le nombre de trajectoires
        
};


// For an asian option with Monte-Carlo

class asian_option : public option {
    public:
        asian_option(double S_t, double K, double vol, double maturity, double tx);
    
        double asian_pricingCallMC(int N); // Pricing d'une asian call option par MC, N lo nombre de trajectoire
        double asian_pricingPutMC(int N); // Pricing d'une asian put option par MC, N lo nombre de trajectoire
};


// For a barrier option_pricing_h

class barrier_option : public option {
    private:
        double q;
        double B;
    public:
        barrier_option(double S_t, double K, double vol, double maturity, double tx, double q, double B);
        double get_q();
        double get_B();
};











#endif // option_pricing_h
