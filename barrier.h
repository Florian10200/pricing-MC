#pragma once
#include <iostream>
using namespace std;

double normalCDF(double x);

class barrierOption{
private:
    double S0;//initial price
    double K; //strike price
    double T; //time until maturity
    double r; //risk-free rate
    double sigma; //volatility
    double q; //Continuous return on dividends
    double B; //barrier level
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

    double Call_UpAndIn_UpAndOut();//fait
    double Call_DownAndIn_DownAndOut();//fait
    double Put_UpAndIn_UpAndOut();//fait
    double Put_DownAndIn_DownAndOut();
    double PriceCall();
    double PricePut();
};
