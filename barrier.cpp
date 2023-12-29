#include <iostream>
#include "barrier.h"
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <random>
using namespace std;

double normalCDF(double x) { //cdf de la loi normale centrée réduite
   return 0.5 * erfc(-x * sqrt(0.5));
}

double generateN01() { //gnérateur de nombres aléatoires
    random_device rd; //initialisation avec une graine
    mt19937 generator(rd());

    normal_distribution<double> distribution(0.0, 1.0); //distribution normale avec une moyenne de 0 et un écart-type de 1

    double nb = distribution(generator);

    return nb;
}

double barrierOption::simassetPrice(){
    double T = maturity;
    double S0 = spot;
    double K = strike;
    double sig = volatility;
    double r = txinteret;
    double q = dividends;
    double B = barrier;
    double dt = T/365;
    double z = sqrt(dt) * static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    double ST = S0*exp((r-pow(sig, 2)/2.0)*T + sig*z);
    return ST;
}

double barrierOption::MC2(string type, int N){
    double T = maturity;
    double S0 = spot;
    double K = strike;
    double sig = volatility;
    double r = txinteret;
    double q = dividends;
    double B = barrier;
    double total=0.0;
    for (int i = 0; i < N; ++i) {
        double ST = simassetPrice();
        double payoff = barrierPayoff(type, ST);
        total+=payoff;
    }
    double avg = total/static_cast<double>(N);
    cout<<"Le prix de l'option avec MC2 est "<< avg*exp(-r*T)<<endl;
    return 0.0;
}

double barrierOption::extremumTraj(string type) {
    double T = maturity;
    double S0 = spot;
    double K = strike;
    double sig = volatility;
    double r = txinteret;
    double q = dividends;
    double B = barrier;
    double maxP = B-1;
    double minP = B+1;
    double S_jdt;
    double dt = 1.0 / 365.0;
    /*double epsilon = rand() / (RAND_MAX + 1.0); // Nombre aléatoire entre 0 et 1
    double z = sqrt(-2.0 * log(epsilon)) * cos(2.0 * M_PI * epsilon); // Transformation de Box-Muller*/

    // Calcul du prix du sous-jacent selon le modèle de diffusion géométrique
    for (int j = 0; j < T*365; ++j) {
            double z = generateN01();
            S_jdt = S0*exp((r-pow(sig, 2)/2.0)*(j*dt) + sig*sqrt(j*dt)*z);
            if (S_jdt > maxP) {
                maxP = S_jdt;
            }
            else if (S_jdt < minP){
                minP = S_jdt;
            }
        }
    if (type=="CDI" or type=="CDO" or type=="PDI" or type=="PDO"){
        return minP;
    }
    else {
        return maxP;
    }
}

bool barrierOption::barrierTouched(string type, double z){
    double B = barrier;
    if (type=="CUI" or type=="PUI" or type=="PUO" or type=="CUO"){
        return z >= B;
    }
    else {
        return z <= B;
    }
}

double barrierOption::barrierPayoff(string type, double z){
    double T = maturity;
    double S0 = spot;
    double K = strike;
    double sig = volatility;
    double r = txinteret;
    double q = dividends;
    double B = barrier;
    double ST = S0*exp(r*T);
    if (type=="CUI" or type=="CDI"){
        if(barrierTouched(type, z)){
            return fmax(0, ST-K);
        }
        else {
            return 0.0;
        }
    }
    else if (type=="CUO" or type=="CDO"){
        if(barrierTouched(type, z)){
            return 0.0;
        }
        else {
            return fmax(0, ST-K);
        }
    }
    else if (type=="PUI" or type=="PDI"){
        if(barrierTouched(type, z)){
            return fmax(0, K-ST);
        }
        else {
            return 0.0;
        }
    }
    else if (type=="PUO" or type=="PDO"){
        if(barrierTouched(type, z)){
            return 0.0;
        }
        else {
            return fmax(0, K-ST);
        }
    }
    return 0.0;
}

double barrierOption::MCbarrier(string type, int N){
    double T = maturity;
    double S0 = spot;
    double K = strike;
    double sig = volatility;
    double r = txinteret;
    double q = dividends;
    double B = barrier;
    srand(static_cast<unsigned>(time(nullptr)));
    double SumPayoff = 0.0;
    for (int i = 0; i < N; ++i) {
        double simulatedPrice = extremumTraj(type);
        double payoff = barrierPayoff(type, simulatedPrice);
        SumPayoff += payoff;
    }
    double optionPrice = (SumPayoff / static_cast<double>(N)) * exp(-r * T);
    cout<<"The barrier option price using Monte-Carlo is "<< optionPrice<<endl;
    return 0.0;
}

barrierOption::barrierOption(double s, double k, double t, double rate, double sig, double Q, double b){
    this -> spot=s;
    this -> strike=k;
    this -> maturity=t;
    this -> txinteret=rate;
    this -> volatility=sig;
    this -> dividends=Q;
    this -> barrier=b;
}

double barrierOption::get_S0(){
    return spot;
}

double barrierOption::get_K(){
    return strike;
}

double barrierOption::get_T(){
    return maturity;
}

double barrierOption::get_B(){
    return barrier;
}

double barrierOption::get_r(){
    return txinteret;
}

double barrierOption::get_sigma(){
    return volatility;
}

double barrierOption::get_q(){
    return dividends;
}

double barrierOption::PriceCall(){
    double T = maturity;
    double S0 = spot;
    double K = strike;
    double sig = volatility;
    double r = txinteret;
    double q = dividends;
    double B = barrier;
    double d1 = (log(S0/K)+(r-q+pow(sig,2)/2)*T)/(sig*sqrt(T));
    double d2 = d1-sig*sqrt(T);
    double Nd1=normalCDF(d1);
    double Nd2=normalCDF(d2);
    double C = S0*exp(-q*T)*Nd1-K*Nd2*exp(-r*T);
    return C;
}

double barrierOption::PricePut(){
    double T = maturity;
    double S0 = spot;
    double K = strike;
    double sig = volatility;
    double r = txinteret;
    double q = dividends;
    double B = barrier;
    double d1 = (log(S0/K)+(r-q+pow(sig,2)/2)*T)/(sig*sqrt(T));
    double d2 = d1-sig*sqrt(T);
    double Nd1=normalCDF(-d1);
    double Nd2=normalCDF(-d2);
    double P = -(S0*exp(-q*T)*Nd1-K*Nd2*exp(-r*T));
    return P;
}

double barrierOption::Call_DownAndIn_DownAndOut(){
    double T = maturity;
    double S0 = spot;
    double K = strike;
    double sig = volatility;
    double r = txinteret;
    double q = dividends;
    double B = barrier;
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

double barrierOption::Call_UpAndIn_UpAndOut(){
    double T = maturity;
    double S0 = spot;
    double K = strike;
    double sig = volatility;
    double r = txinteret;
    double q = dividends;
    double B = barrier;
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

double barrierOption::Put_UpAndIn_UpAndOut(){
    double T = maturity;
    double S0 = spot;
    double K = strike;
    double sig = volatility;
    double r = txinteret;
    double q = dividends;
    double B = barrier;
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

double barrierOption::Put_DownAndIn_DownAndOut(){
    double T = maturity;
    double S0 = spot;
    double K = strike;
    double sig = volatility;
    double r = txinteret;
    double q = dividends;
    double B = barrier;
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
