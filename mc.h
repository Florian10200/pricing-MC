#include <iostream>
#include <math.h>
#include <random>

using namespace std;

double normalCDF(double x);

double generateN01(); 

class europeanOption { //european option, sans dividende
private:
    double strike; //in $
    double spot; //in $
    double volatility; //in %
    double maturity; //in years
    double txinteret; //in %
public:
    europeanOption (double st, double sp, double vol, double maturity, double tx); //constructeur
    double get_strike(); //recupere la valeur du strike
    double get_spot(); //recupere la valeur du spot
    double get_volatility(); //recupere la valeur de la volatilite
    double get_maturity(); //recupere la valeur de la maturite
    double get_tx(); //recupere la valeur du tx d'interet
    
    double pricing_european_call_BS(); //pricing european call par BS
    double pricing_european_put_BS(); //pricing european put par BS
    
    double delta_european_call(); //calcul du delta pour un european call
    double delta_european_put(); //calcul du delta pour un european put
    
    double pricing_european_call_MC(int N); //pricing european call par MC, N le nombre de trajectoires
    double pricing_european_put_MC(int N); //pricing european put par MC, N le nombre de trajectoires
};

class lookbackOption { //lookback flottant, sans dividende
private: //pas de strike
    double spot; //in $
    double volatility; //in %
    double maturity; //in years
    double txinteret; //in %
public:
    lookbackOption (double sp, double vol, double maturity, double tx); //constructeur d'un lookback flottant
    
    double get_strike(); //recupere la valeur du strike
    double get_spot(); //recupere la valeur du spot
    double get_volatility(); //recupere la valeur de la volatilite
    double get_maturity(); //recupere la valeur de la maturite
    double get_tx(); //recupere la valeur du tx d'interet
    
    double pricing_lookback_call_MC(int N); //pricing lookback call par MC, N le nombre de trajectoires
    double pricing_lookback_put_MC(int N); //pricing lookback put par MC, N le nombre de trajectoires
    
};
