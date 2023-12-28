#include <iostream>
#include <math.h>
#include <random>
#include "mc.h"

using namespace std;

int main()
{   europeanOption option1(52.0, 50.0, 0.2, 2.0, 0.05);

    double K = option1.getstrike();
    double S0 = option1.getspot();
    double sigma = option1.getvolatility();
    double T = option1.getmaturity();
    double r = option1.gettx();

    
    cout << "strike :" << K << endl;
    cout << "spot :" << S0 << endl;
    cout << "volatility :" << sigma << endl;
    cout << "maturity :" << T << endl;
    cout << "taux d'intéret :" << r << endl;

    
    cout << option1.pricingCallBS() << endl;
    cout << option1.pricingPutBS() << endl;
    
    
    cout << option1.deltaCall() << endl;
    cout << option1.deltaPut() << endl;
    
    
    cout << option1.pricingCallMC(10000) << endl; //convergence vers la valeur donnée par BS
    cout << option1.pricingPutMC(10000) << endl; //convergence vers la valeur donnée par BS
    
    
    //----------------------------------------------------------------------------------------------------
    
    lookbackOption option2(50.0, 0.2, 0.25, 0.1);

    
    cout << option2.pricingCallMC(1000) << endl;
    cout << option2.pricingPutMC(1000) << endl;
    
    
    return 0;
}
