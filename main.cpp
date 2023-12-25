#include <iostream>
#include <math.h>
#include <random>
#include "mc.h"

using namespace std;

int main()
{   europeanOption option(52.0, 50.0, 0.2, 2.0, 0.05);

    double K = option.get_strike();
    double S0 = option.get_spot();
    double sigma = option.get_volatility();
    double T = option.get_maturity();
    double r = option.get_tx();

    
    cout << "strike :" << K << endl;
    cout << "spot :" << S0 << endl;
    cout << "volatility :" << sigma << endl;
    cout << "maturity :" << T << endl;
    cout << "taux d'intéret :" << r << endl;

    
    cout << option.pricing_european_call_BS() << endl;
    cout << option.pricing_european_put_BS() << endl;
    
    
    cout << option.delta_european_call() << endl;
    cout << option.delta_european_put() << endl;
    
    
    cout << option.pricing_european_call_MC(10000) << endl; //convergence vers la valeur donnée par BS
    cout << option.pricing_european_put_MC(10000) << endl; //convergence vers la valeur donnée par BS
    
    
    cout << option.pricing_american_call_MC(1000) << endl; //un american call est équivalent à un european call
    cout << option.pricing_american_put_MC(1000) << endl; //un american put est un peu plus cher que un european put 
    
    
    return 0;
}
