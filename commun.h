#ifndef option_pricing_h
#define option_pricing_h
#include<iostream>

using namespace std;

double normalCDF(double x); // Cumulative distribution fonction of the N(0,1)
double generateN01(); // Generate a distribution of N(0,1)



class option { // General class for all the options
    private:
        double spot; //in $
        double strike; //in $
        double volatility; //in %
        double maturity; //in years
        double txinteret; //in %
    public:
        option(double S_t, double K, double vol, double maturity, double tx); //constructor
        double get_strike(); //get the strike price
        double get_spot(); //get the spot price
        double get_volatility(); //get the volatility
        double get_maturity(); //get the maturity
        double get_tx(); //get the interest rate
};        

// European option Pricing with analytic formula and Monte-Carlo---------------------------------------------------------------

class european_option : public option { // Heritage from the class option
    public:
        european_option(double S_t, double K, double vol, double maturity, double tx);
        
        double pricing_european_call_BS(); //pricing european call with BS
        double pricing_european_put_BS(); //pricing european put with BS
    
        double delta_european_call(); //Delta for a european call
        double delta_european_put(); //Delta for a european put
        
        double pricing_european_call_MC(int N); //pricing european call with MC, N the number of path
        double pricing_european_put_MC(int N); //pricing european put with MC, N the number of path
};

// For a lookback option with Monte-Carlo-----------------------------------------------------------------------------------------


class lookback_option : public option {
    public:
        lookback_option(double S_t, double K, double vol, double maturity, double tx);
        
        double lookback_pricingCallMC(int N); //pricing lookback call with MC, N the number of path
        double lookback_pricingPutMC(int N); //pricing lookback put with MC, N the number of path
        
};


// For an asian option with Monte-Carlo-------------------------------------------------------------------------------------

class asian_option : public option {
    public:
        asian_option(double S_t, double K, double vol, double maturity, double tx);
    
        double asian_pricingCallMC(int N); // Pricing asian call option with MC, N the number of path
        double asian_pricingPutMC(int N); // Pricing asian put option with MC, N the number of path
};


// For a barrier option------------------------------------------------------------------------------------------------

class barrier_option : public option { // The Barrier option class has two additional arguments : the Barrier B and the dividends q
    private:
        double q;
        double B;
    public:
        barrier_option(double S_t, double K, double vol, double maturity, double tx, double q, double B);
        double get_q();
        double get_B();
        double PriceCall();
        double PricePut();
        double Call_UpAndIn_UpAndOut();
        double Call_DownAndIn_DownAndOut();
        double Put_UpAndIn_UpAndOut();
        double Put_DownAndIn_DownAndOut();
};


#endif // option_pricing_h
