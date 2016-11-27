#include <iostream>
#include <math.h>
#include <random>
#include <cmath>
#include <vector>

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

int const T_max = 5000;

double tab_S0[T_max];
double tab_C[T_max];
double tab_delta[T_max];
double tab_gamma[T_max];

float delta_eval_talor(float S0, float d1, float d2, float K, float r, float T, float eps){
    float delta = (C(S0 + eps, d1, d2, K, r, T) - C(S0 - eps, d1, d2, K, r, T))/(2*eps);
    return delta;
}

float gamma_eval_taylor(double S0, double d1, double d2, double K, double r, double T, double eps){
    double a1 = C(S0 + eps, d1, d2, K, r, T);
    double a2 = C(S0 - eps, d1, d2, K, r, T);
    double a3 = C(S0, d1, d2, K, r, T);
    double a4 = (a1 + a2 - 2*a3);
    double a5 = eps*eps;
    double gamma = a4/a5;
    return gamma;
}



// array containing the draw of Xi
double tab_xi[T_max];

void x_eval(double r, double T, double S0, double sigma, double K){
    for (int i = 0 ; i < T_max ; i++){
        double xi = exp(-r*T)*fmax(0,S0*exp((r-pow(sigma,2)/2)*T + (sigma*sqrt(T)*d(gen)))-K);
        tab_xi[i] = xi;
    }
}


double monte_carlo[T_max];
void monte_carlo_eval(){
    // Evaluate the x_i using monte carlo simulation
    // Store the simu in monte_carlo array
    monte_carlo[0] = tab_xi[0];
    for (int i = 1 ; i < T_max ; i++) {
        monte_carlo[i] = (monte_carlo[i-1]*(i-1) + tab_xi[i])/i;
    }
}

double d_1 = d1(sigma, T, S0, K, r);
double d_2 = d2(d_1, sigma, T);
double C_ = C(S0, d_1, d_2, K, r, T);


double calculateSD(vector<double> data)
{
    double sum = 0.0, mean, std = 0.0;

    for(int i = 0; i < sizeof(data); ++i)
    {
        sum += data[i];
    }

    mean = sum/sizeof(data);
    for(int i = 0; i < sizeof(data); ++i)
        std += pow(data[i] - mean, 2);

    return sqrt(std / sizeof(data));
}

double monte_carlo_normalite[T_max];
void fill_monte_carlo_normalite(){
    monte_carlo_normalite[0] = monte_carlo[0];
    vector<double> sub_tab_xi;
    for (int i = 1 ; i < T_max ; i++) {
        sub_tab_xi.push_back(tab_xi[i]);
        double std = calculateSD(sub_tab_xi);
        monte_carlo_normalite[i] = (monte_carlo[i] - C_)/(std/sqrt(i));
    }
}


double gamma_eval(double sigma, double T, double S0, double K, double r){
    double d1_ = d1(sigma, T, S0, K, r);
    return exp(-pow(d1_,2))/(sigma*S0*sqrt(2*M_PI*T));
}


void fill_tab(int size){

    for (int i = 0 ; i < size ; i++)
    {
        tab_S0[i] = i;
        double d_1 = d1(sigma, T, tab_S0[i], K, r);
        double d_2 = d2(d_1, sigma, T);
        double C_ = C(tab_S0[i], d_1, d_2, K, r, T);

        //float delta = delta_eval(tab_S0[i], d_1, d_2, K, r, T, 0.01);
        double gamma = gamma_eval(sigma, T, tab_S0[i], K, r);

        tab_C[i] = C_;
        //tab_delta[i] = delta;
        tab_gamma[i] = gamma;
    }

}

int main() {

    cout << d(gen) << endl;

    fill_tab(1000);
    cout << tab_C[10] << endl;

    // drawing the xi to fill tab_xi
    x_eval(r, T, S0, sigma, K);

    monte_carlo_eval();

    fill_monte_carlo_normalite();

/*    for (int i = 0; i < 1000; i++)
        cout << monte_carlo_normalite[i] << endl;*/


    // we draw 1000 times from monte_carlo_normalite[10]
    double hist_n_10[1000];
    for (int j = 0; j < 1000; j++){
        x_eval(r, T, S0, sigma, K);
        monte_carlo_eval();
        fill_monte_carlo_normalite();
        hist_n_10[j] = monte_carlo_normalite[500];
    }
    for (int i = 0; i < 1000; i++)
        cout << hist_n_10[i] << endl;

    return 0;

}