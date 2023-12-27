#include <iostream>
#include "barrier.h"
#include <cmath>
using namespace std;

int main() {
    string type;
    double S0;
    double K;
    double T;
    double r;
    double sigma;
    double q;
    double B;
    cout<<"What type of barrier option do you want to price? (possible answer: CD, CU, PU, PD)";
    cin>> type;
    if (type=="CD"){
        cout<<"Enter the characteristics of the option(stock price, strike price, maturity, rate, volatility, yield, barrier):";
        cin>>S0;
        cin>>K;
        cin>>T;
        cin>>r;
        cin>>sigma;
        cin>>q;
        cin>>B;
        barrierOption option(S0, K, T, r, sigma, q, B);
        option.Call_DownAndIn_DownAndOut();
    }
    else if(type=="CU"){
        cout<<"Enter the characteristics of the option(stock price, strike price, maturity, rate, volatility, yield, barrier):";
        cin>>S0;
        cin>>K;
        cin>>T;
        cin>>r;
        cin>>sigma;
        cin>>q;
        cin>>B;
        barrierOption option(S0, K, T, r, sigma, q, B);
        option.Call_UpAndIn_UpAndOut();
    }
    else if(type=="PU"){
        cout<<"Enter the characteristics of the option(stock price, strike price, maturity, rate, volatility, yield, barrier):";
        cin>>S0;
        cin>>K;
        cin>>T;
        cin>>r;
        cin>>sigma;
        cin>>q;
        cin>>B;
        barrierOption option(S0, K, T, r, sigma, q, B);
        option.Put_UpAndIn_UpAndOut();
    }
    else if(type=="PD"){
        cout<<"Enter the characteristics of the option(stock price, strike price, maturity, rate, volatility, yield, barrier):";
        cin>>S0;
        cin>>K;
        cin>>T;
        cin>>r;
        cin>>sigma;
        cin>>q;
        cin>>B;
        barrierOption option(S0, K, T, r, sigma, q, B);
        option.Put_DownAndIn_DownAndOut();
    }
    else{
        cout<<"error";
    }
    /*barrierOption option(S0, K, T, r, sigma, q, B);

    option.Call_DownAndIn_DownAndOut();
    option.Call_UpAndIn_UpAndOut();
    option.Put_UpAndIn_UpAndOut();
    option.Put_DownAndIn_DownAndOut();*/
    return 0;
};
