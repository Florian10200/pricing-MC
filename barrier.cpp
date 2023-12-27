#include <iostream>
#include "barrier.h"
#include <cmath>
using namespace std;

double normalCDF(double x) { //cdf de la loi normale centrée réduite
   return 0.5 * erfc(-x * sqrt(0.5));
}

barrierOption::barrierOption(double s, double k, double t, double rate, double sig, double Q, double b){
    this -> S0=s;
    this -> K=k;
    this -> T=t;
    this -> r=rate;
    this -> sigma=sig;
    this -> q=Q;
    this -> B=b;
}

double barrierOption::get_S0(){
    return S0;
}

double barrierOption::get_K(){
    return K;
}

double barrierOption::get_T(){
    return T;
}

double barrierOption::get_B(){
    return B;
}

double barrierOption::get_r(){
    return r;
}

double barrierOption::get_sigma(){
    return sigma;
}

double barrierOption::get_q(){
    return q;
}



double barrierOption::PriceCall(){
    double d1 = (log(this->S0/this->K)+(this->r-this->q+pow(this->sigma,2)/2)*this->T)/(this->sigma*sqrt(this->T));
    double d2 = d1-this->sigma*sqrt(this->T);
    double Nd1=normalCDF(d1);
    double Nd2=normalCDF(d2);
    double C = this->S0*exp(-this->q*this->T)*Nd1-this->K*Nd2*exp(-this->r*this->T);
    return C;
}

double barrierOption::PricePut(){
    double d1 = (log(this->S0/this->K)+(this->r-this->q+pow(this->sigma,2)/2)*this->T)/(this->sigma*sqrt(this->T));
    double d2 = d1-this->sigma*sqrt(this->T);
    double Nd1=normalCDF(-d1);
    double Nd2=normalCDF(-d2);
    double P = -(this->S0*exp(-this->q*this->T)*Nd1-this->K*Nd2*exp(-this->r*this->T));
    return P;
}

double barrierOption::Call_DownAndIn_DownAndOut(){
    double lam = (this->r-this->q+pow(this->sigma,2)/2)/pow(this->sigma,2);
    double y = log(pow(this->B,2)/(this->S0*this->K))/(this->sigma*sqrt(this->T))+lam*this->sigma*sqrt(this->T);
    double Ny=normalCDF(y);
    double Ny2=normalCDF(y-this->sigma*sqrt(this->T));
    double C = PriceCall();
    double x1 = (log(this->S0/this->B)+lam*pow(this->sigma,2)*this->T)/(this->sigma*sqrt(this->T));
    double y1 = (this->B/log(this->S0)+lam*pow(this->sigma,2)*this->T)/(this->sigma*sqrt(this->T));
    double Nx1 = normalCDF(x1);
    double Nx11=normalCDF(x1-this->sigma*sqrt(this->T));
    double Ny1 = normalCDF(y1);
    double Ny11=normalCDF(y1-this->sigma*sqrt(this->T));
    if (this->B <= this->K){
        double cdi= this->S0*exp(-this->q*this->T)*pow(this->B/this->S0,2*lam)*Ny-Ny2*this->K*pow(this->B/this->S0,2*lam-2)*exp(-this->r*this->T);
        cout<<"The value of the CDI is "<< cdi <<endl;
        cout<<"The value of the CDO is "<<(C-cdi)<<endl;
    }
    else {
        double cdo=this->S0*exp(-this->q*this->T)*Nx1-Nx11*this->K*exp(-this->r*this->T)-this->S0*exp(-this->q*this->T)*pow(this->B/this->S0,2*lam)*Ny1+Ny11*this->K*pow(this->B/this->S0,2*lam-2)*exp(-this->r*this->T);
        cout<<"The value of the CDI is "<< (C-cdo) <<endl;
        cout<<"The value of the CDO is "<< cdo <<endl;
    }
    return 0.0;

}

double barrierOption::Call_UpAndIn_UpAndOut(){
    // variable d1,d2,d3,d4
    double lam = (this->r-this->q+pow(this->sigma,2)/2)/pow(this->sigma,2);
    double C = PriceCall();
    double x1 = (log(this->S0/this->B)+lam*pow(this->sigma,2)*this->T)/(this->sigma*sqrt(this->T));
    double y = log(pow(this->B,2)/(this->S0*this->K))/(this->sigma*sqrt(this->T))+lam*this->sigma*sqrt(this->T);
    double y1 = (this->B/log(this->S0)+lam*pow(this->sigma,2)*this->T)/(this->sigma*sqrt(this->T));
    double Ny=normalCDF(-y);
    double Nyy=normalCDF(-y+this->sigma*sqrt(this->T));
    double Ny1=normalCDF(-y1);
    double Ny11=normalCDF(-y1+this->sigma*sqrt(this->T));
    double Nx1 = normalCDF(x1);
    double Nx11=normalCDF(x1-this->sigma*sqrt(this->T));
    if (this->B > this->K){
        double cui=this->S0*exp(-this->q*this->T)*Nx1-Nx11*this->K*exp(-this->r*this->T)-this->S0*exp(-this->q*this->T)*pow(this->B/this->S0,2*lam)*(Ny-Ny1)+(Nyy-Ny11)*this->K*pow(this->B/this->S0,2*lam-2)*exp(-this->r*this->T);
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
    double lam = (this->r-this->q+pow(this->sigma,2)/2)/pow(this->sigma,2);
    double y = log(pow(this->B,2)/(this->S0*this->K))/(this->sigma*sqrt(this->T))+lam*this->sigma*sqrt(this->T);
    double Ny = normalCDF(-y);
    double Nyy = normalCDF(-y+this->sigma*sqrt(this->T));
    double x1 = (log(this->S0/this->B)+lam*pow(this->sigma,2)*this->T)/(this->sigma*sqrt(this->T));
    double y1 = (this->B/log(this->S0)+lam*pow(this->sigma,2)*this->T)/(this->sigma*sqrt(this->T));
    double Nx1 = normalCDF(-x1);
    double Ny1 = normalCDF(-y1);
    double Nx11=normalCDF(-x1+this->sigma*sqrt(this->T));
    double Ny11=normalCDF(-y1+this->sigma*sqrt(this->T));
    double P = PricePut();
    if (this->B >= this->K){
        double pui = this->S0*(-1)*exp(-this->q*this->T)*pow(this->B/this->S0,2*lam)*Ny+Nyy*this->K*pow(this->B/this->S0,2*lam-2)*exp(-this->r*this->T);
        cout<<"The value of the PUI is "<< pui <<endl;
        cout<<"The value of the PUO is "<< (P-pui)<<endl;
    }
    else{
        double puo = this->S0*(-1)*exp(-this->q*this->T)*Nx1+Nx11*this->K*exp(-this->r*this->T)+this->S0*exp(-this->q*this->T)*pow(this->B/this->S0,2*lam)*Ny1-Ny11*this->K*pow(this->B/this->S0,2*lam-2)*exp(-this->r*this->T);
        cout<<"The value of the PUI is "<< (P-puo) <<endl;
        cout<<"The value of the PUO is "<< puo <<endl;
    }
    return 0.0;
}

double barrierOption::Put_DownAndIn_DownAndOut(){
    double lam = (this->r-this->q+pow(this->sigma,2)/2)/pow(this->sigma,2);
    double P = PricePut();
    double x1 = (log(this->S0/this->B)+lam*pow(this->sigma,2)*this->T)/(this->sigma*sqrt(this->T));
    double y1 = (this->B/log(this->S0)+lam*pow(this->sigma,2)*this->T)/(this->sigma*sqrt(this->T));
    double y = log(pow(this->B,2)/(this->S0*this->K))/(this->sigma*sqrt(this->T))+lam*this->sigma*sqrt(this->T);
    double Nx1=normalCDF(-x1);
    double Nx11 = normalCDF(-x1+this->sigma*sqrt(this->T));
    double Ny = normalCDF(y);
    double Ny1 = normalCDF(y1);
    double Nyy = normalCDF(y-this->sigma*sqrt(this->T));
    double Ny11 = normalCDF(y1-this->sigma*sqrt(this->T));
    if (this->B <= this->K){
        double pdi = this->S0*(-1)*exp(-this->q*this->T)*Nx1+Nx11*this->K*exp(-this->r*this->T)+this->S0*exp(-this->q*this->T)*pow(this->B/this->S0,2*lam)*(Ny-Ny1)-(Nyy-Ny11)*this->K*pow(this->B/this->S0,2*lam-2)*exp(-this->r*this->T);
        cout<<"The value of the PDI is "<< pdi <<endl;
        cout<<"The value of the PDO is "<< (P-pdi) <<endl;
    }
    else{
        cout<<"The value of the PDI is "<< P <<endl;
        cout<<"The value of the PDO is "<< 0.0 <<endl;
    }
    return 0.0;
}
