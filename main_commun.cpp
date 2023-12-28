#include"option_pricing.h"
#include<iostream>

using namespace std;

int main(){
    
// For a european option -------------------------------------------------------

    european_option option1(50.0, 52.0, 0.2, 2.0, 0.05);
    
    double S0 = option1.get_spot();
    double K = option1.get_strike();
    double sigma = option1.get_volatility();
    double T = option1.get_maturity();
    double r = option1.get_tx();

    cout << "spot :" << S0 << endl;
    cout << "strike :" << K << endl;
    cout << "volatility :" << sigma << endl;
    cout << "maturity :" << T << endl;
    cout << "taux d'intéret :" << r << endl;

    cout << "pricing european call with BS : " << option1.pricing_european_call_BS() << endl;
    cout << "pricing european put with BS : " << option1.pricing_european_put_BS() << endl;
    
    cout << "delta for a european call : " << option1.delta_european_call() << endl;
    cout << "delat for a european put : " << option1.delta_european_put() << endl;
    
    cout << "pricing european call with Monte-Carlo : " << option1.pricing_european_call_MC(1000) << endl;
    cout << "pricing european put avec Monte-Carlo : " << option1.pricing_european_put_MC(1000) << endl;

// For a lookback option ---------------------------------------------------------------------

    lookback_option option2(50.0, 52.0, 0.2, 2.0, 0.05);
    
    cout<<"pricing lookback call option avec MC : "<<option2.lookback_pricingCallMC(100)<<endl;
    cout<<"pricing lookback put option avec MC : "<<option2.lookback_pricingPutMC(100)<<endl;


// For an asian option ----------------------------------------------------------------------

    asian_option option3(50.0, 52.0, 0.2, 2.0, 0.05);

    cout<<"pricing asian call option avec MC : "<<option3.asian_pricingCallMC(100)<<endl;
    cout<<"pricing asian put option avec MC : "<<option3.asian_pricingPutMC(100)<<endl;

// For a barrier option ------------------------------------------------------------------------

    barrier_option option4(50.0, 52.0, 0.2, 2.0, 0.05,33.3,44.4);
    



};
  
