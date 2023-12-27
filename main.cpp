#include <iostream>
#include <math.h>
#include <random>
#include "mc.h"

using namespace std;

int main()
{   europeanOption option1(52.0, 50.0, 0.2, 2.0, 0.05);

    double K = option1.get_strike();
    double S0 = option1.get_spot();
    double sigma = option1.get_volatility();
    double T = option1.get_maturity();
    double r = option1.get_tx();

    
    cout << "strike :" << K << endl;
    cout << "spot :" << S0 << endl;
    cout << "volatility :" << sigma << endl;
    cout << "maturity :" << T << endl;
    cout << "taux d'intéret :" << r << endl;

    
    cout << option1.pricing_european_call_BS() << endl;
    cout << option1.pricing_european_put_BS() << endl;
    
    
    cout << option1.delta_european_call() << endl;
    cout << option1.delta_european_put() << endl;
    
    
    cout << option1.pricing_european_call_MC(10000) << endl; //convergence vers la valeur donnée par BS
    cout << option1.pricing_european_put_MC(10000) << endl; //convergence vers la valeur donnée par BS
    
    
    //----------------------------------------------------------------------------------------------------
    
    lookbackOption option2(50.0, 0.2, 0.25, 0.1);

    
    cout << option2.pricing_lookback_call_MC(1000) << endl;
    cout << option2.pricing_lookback_put_MC(1000) << endl;
    
    
    return 0;
}
