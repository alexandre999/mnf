#include <iostream>
#include <math.h>
#include <random>
#include <cmath>

double phi(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}


using namespace std;

double d1(double sigma, double T, double S0, double K, double r){
    return (log(S0/K) + (r + pow(sigma,2)/2)*T)/(sigma*sqrt(T));
}


double d2(double d1, double sigma, double T){
    return d1 - sigma*sqrt(T);
}



double C(double S0, double d1, double d2, double K, double r, double T){
    return S0*phi(d1) - K*exp(-r*T)*phi(d2);
}

random_device rd;
mt19937 gen(rd());
normal_distribution<> d(0,1);


double sigma = 0.14;
double T = 10;
double S0 = 205;
double K = 210;
double r = 0.0021;


double tab_S0[1000];
double tab_C[1000];
double tab_delta[1000];
double tab_gamma[1000];

float delta_eval(float S0, float d1, float d2, float K, float r, float T, float eps){
    float delta = (C(S0 + eps, d1, d2, K, r, T) - C(S0 - eps, d1, d2, K, r, T))/(2*eps);
    return delta;
}

float gamma_eval(double S0, double d1, double d2, double K, double r, double T, double eps){
    double a1 = C(S0 + eps, d1, d2, K, r, T);
    double a2 = C(S0 - eps, d1, d2, K, r, T);
    double a3 = C(S0, d1, d2, K, r, T);
    double a4 = (a1 + a2 - 2*a3);
    double a5 = eps*eps;
    double gamma = a4/a5;
    return gamma;
}


// array containing the draw of Xi
double tab_xi[1000];

void x_eval(double r, double T, double S0, double sigma, double K){
    for (int i = 0 ; i < 1000 ; i++){
        double xi = exp(-r*T)*fmax(0,S0*exp((r-pow(sigma,2)/2)*T + (sigma*sqrt(T)*d(gen)))-K);
        tab_xi[i] = xi;
    }
}


double monte_carlo[1000];
void monte_carlo_eval(){
    monte_carlo[0] = tab_xi[0];
    for (int i = 1 ; i < 1000 ; i++) {
        monte_carlo[i] = (monte_carlo[i-1]*(i-1) + tab_xi[i])/i;
    }
}



void fill_tab(int size){

    for (int i = 0 ; i < size ; i++)
    {
        tab_S0[i] = i;
        double d_1 = d1(sigma, T, tab_S0[i], K, r);
        double d_2 = d2(d_1, sigma, T);
        double C_ = C(tab_S0[i], d_1, d_2, K, r, T);
        //float delta = delta_eval(tab_S0[i], d_1, d_2, K, r, T, 0.01);
        double gamma = gamma_eval(tab_S0[i], d_1, d_2, K, r, T, 0.01);

        tab_C[i] = C_;
        //tab_delta[i] = delta;
        tab_gamma[i] = gamma;
    }

}

int main() {





    cout << d(gen) << endl;

    fill_tab(1000);
    cout << tab_C[10] << endl;

    x_eval(r, T, S0, sigma, K);
    monte_carlo_eval();
    for (int i = 0; i < 1000; i++)
        cout << monte_carlo[i] << endl;

    double d_1 = d1(sigma, T, S0, K, r);
    double d_2 = d2(d_1, sigma, T);
    double C_ = C(S0, d_1, d_2, K, r, T);
    cout << C_ << endl;
    cout << C_ << endl;
    return 0;




}