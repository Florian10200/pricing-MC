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
    double getstrike(); //recupere la valeur du strike
    double getspot(); //recupere la valeur du spot
    double getvolatility(); //recupere la valeur de la volatilite
    double getmaturity(); //recupere la valeur de la maturite
    double gettx(); //recupere la valeur du tx d'interet
    
    double pricingCallBS(); //pricing european call par BS
    double pricingPutBS(); //pricing european put par BS
    
    double deltaCall(); //calcul du delta pour un european call
    double deltaPut(); //calcul du delta pour un european put
    
    double pricingCallMC(int N); //pricing european call par MC, N le nombre de trajectoires
    double pricingPutMC(int N); //pricing european put par MC, N le nombre de trajectoires
};

class lookbackOption { //lookback flottant, sans dividende
private: //pas de strike
    double spot; //in $
    double volatility; //in %
    double maturity; //in years
    double txinteret; //in %
public:
    lookbackOption (double sp, double vol, double maturity, double tx); //constructeur d'un lookback flottant
    
    double getstrike(); //recupere la valeur du strike
    double getspot(); //recupere la valeur du spot
    double getvolatility(); //recupere la valeur de la volatilite
    double getmaturity(); //recupere la valeur de la maturite
    double gettx(); //recupere la valeur du tx d'interet
    
    double pricingCallMC(int N); //pricing lookback call par MC, N le nombre de trajectoires
    double pricingPutMC(int N); //pricing lookback put par MC, N le nombre de trajectoires
    
};
