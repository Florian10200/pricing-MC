#include"option_pricing.h"
#include<iostream>

using namespace std;

int main(){
    
// For a european option -----------------------------------------------------------------------------------------------------------------------------------------

    european_option option1(50.0, 52.0, 0.2, 0.5, 0.05);
    
    double S0 = option1.get_spot(); // We get the spot price at time 0
    double K = option1.get_strike(); // We get the strike price
    double sigma = option1.get_volatility(); // We get the volatility
    double T = option1.get_maturity(); // We get the maturity
    double r = option1.get_tx(); // We get the interest rate

    cout << "spot :" << S0 << endl; // We show the different values of the option
    cout << "strike :" << K << endl;
    cout << "volatility :" << sigma << endl;
    cout << "maturity :" << T << endl;
    cout << "interest rate :" << r << endl;

    cout << "pricing european call with BS : " << option1.pricing_european_call_BS() << endl; // We get the price for a call and a put with Black-Scholes formula
    cout << "pricing european put with BS : " << option1.pricing_european_put_BS() << endl;
    
    cout << "delta for a european call : " << option1.delta_european_call() << endl; // We calculate the delta
    cout << "delat for a european put : " << option1.delta_european_put() << endl;
    
    cout << "pricing european call with Monte-Carlo : " << option1.pricing_european_call_MC(1000) << endl; // We price them with Monte-Carlo as well
    cout << "pricing european put with Monte-Carlo : " << option1.pricing_european_put_MC(1000) << endl;

    
// For a lookback option ----------------------------------------------------------------------------------------------------------------------------------

    lookback_option option2(50.0, 52.0, 0.2, 0.5, 0.05); // The strike price will not be used for this pricing
    
    cout<<"pricing lookback call option with MC : "<<option2.lookback_pricingCallMC(1000)<<endl; // Pricing with Monte-Carlo for a call and a put
    cout<<"pricing lookback put option with MC : "<<option2.lookback_pricingPutMC(1000)<<endl;


// For an asian option --------------------------------------------------------------------------------------------------------------------------------

    asian_option option3(50.0, 52.0, 0.2, 0.5, 0.05);

    cout<<"pricing asian call option with MC : "<<option3.asian_pricingCallMC(1000)<<endl; // Pricing with Monte-Carlo for a call and a put
    cout<<"pricing asian put option with MC : "<<option3.asian_pricingPutMC(1000)<<endl;

    
// For a barrier option -------------------------------------------------------------------------------------------------------------------------------------

    barrier_option option4(50.0, 52.0, 0.2, 0.5, 0.05,0,10);
    
    cout<<option4.Call_UpAndIn_UpAndOut()<<endl;
    cout<<option4.Call_DownAndIn_DownAndOut()<<endl;
    cout<<option4.Put_UpAndIn_UpAndOut()<<endl;
    cout<<option4.Put_DownAndIn_DownAndOut()<<endl;
    



};
  
