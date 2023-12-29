#pragma once
#include <iostream>
using namespace std;

double normalCDF(double x);

class barrierOption{
private:
    double spot;//initial price
    double strike; //strike price
    double maturity; //time until maturity
    double txinteret; //risk-free rate
    double volatility; //volatility
    double dividends; //Continuous return on dividends
    double barrier; //barrier level
public:
    //we code the getter
    barrierOption (double S0, double K, double T, double r, double sigma, double q, double B); //constructor
    double get_S0();
    double get_K();
    double get_T();
    double get_r();
    double get_sigma();
    double get_q();
    double get_B();

    double Call_UpAndIn_UpAndOut();
    double Call_DownAndIn_DownAndOut();
    double Put_UpAndIn_UpAndOut();
    double Put_DownAndIn_DownAndOut();
    double PriceCall();
    double PricePut();
};
