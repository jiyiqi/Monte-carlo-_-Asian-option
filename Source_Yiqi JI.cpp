// Using windows 
// 17078089 Yiqi Ji
//Applied computational finance assignment 1

#include<iostream>
#include<cmath>
#include<ctime>
#include<fstream>
using namespace std;

double normal(double, double); // function prototype
double variance(double a[] , double, double);

int main() {

for (int m = 0; m < 4; m++){
  for (int n = 0; n < 4; n++){
		srand((unsigned)time(NULL));
		double NN[4] = { 100, 500, 1000, 10000 }; //100 and 500 are discrete sampling, 1000 and 10000 are continuous sampling
		double N_sim[4] = {  1, 100, 1000,10000 };//number of simulations
		double N = NN[m];
		double n_sim = N_sim[n];	
		double var[10001] = { 0 };
		//double N = 100;
		//double n_sim = 100;
		double S0 = 100, E = 100, IR = 0.05, sigma = 0.2; // variables & parameters 
		double mc_S0 = 100;
		double error = 0;
		double T = 1;
		double dt = T / N; // step-size for time  path-dependence
		double S = 0;
		double mc_S = 0;
		double ave_arithmetic = 100 / (N + 1), ave_geometric = pow(100, 1 / (N + 1));
		double arithmetic_fix_call = 0, arithmetic_fix_put = 0;
		double geometric_fix_call = 0, geometric_fix_put = 0;
		double arithmetic_float_call = 0, arithmetic_float_put = 0;
		double geometric_float_call = 0, geometric_float_put = 0;
		double arithmetic_fix_call_mean = 0, arithmetic_fix_put_mean = 0;
		double geometric_fix_call_mean = 0, geometric_fix_put_mean = 0;
		double arithmetic_float_call_mean = 0, arithmetic_float_put_mean = 0;
		double geometric_float_call_mean = 0, geometric_float_put_mean = 0;

		double arithmetic_fix_call_var = 0, arithmetic_fix_put_var = 0; //variance of payoff
		double geometric_fix_call_var = 0, geometric_fix_put_var = 0;
		double arithmetic_float_call_var = 0, arithmetic_float_put_var = 0;
		double geometric_float_call_var = 0, geometric_float_put_var = 0;

		double arithmetic_fix_call_var1[10000] = { 0 }, arithmetic_fix_put_var1[10000] = { 0 };
		double geometric_fix_call_var1[10000] = { 0 }, geometric_fix_put_var1[10000] = { 0 };
		double arithmetic_float_call_var1[10000] = { 0 }, arithmetic_float_put_var1[10000] = { 0 };
		double geometric_float_call_var1[10000] = { 0 }, geometric_float_put_var1[10000] = { 0 };

			for (int j = 1; j <= n_sim; j++){

				ave_arithmetic = 100 / (N + 1);
				ave_geometric = pow(100, 1 / (N + 1));
				S0 = 100;
				S = 0;
				//mc_S0 = 100;
				//mc_S = 0;
				for (unsigned short int i = 1; i <= N; i++){
					double time = i*dt;
					double dX = normal(1.0, 0.0)*sqrt(dt);
					//double S = S0*exp((IR - 0.5*sigma*sigma)*(dt)+sigma*dX);
					S = S0*(1 + IR*dt + sigma*dX); //ds = s*(1+r*dt+sigma*sqrt dt *normal)
					//mc_S = mc_S0*exp((IR - 0.5*sigma*sigma)*(dt)+sigma*dX);
					//mc_S0 = mc_S;
					S0 = S;
					ave_arithmetic += S / (N + 1);
					ave_geometric *= pow(S, 1 / (N + 1));
				}
				//error = error+abs(mc_S - S)/n_sim;
				arithmetic_fix_call = ave_arithmetic > E ? ave_arithmetic - E : 0; // payoff call	
				arithmetic_fix_put = E > ave_arithmetic ? E - ave_arithmetic : 0; // payoff put
				geometric_fix_call = ave_geometric > E ? ave_geometric - E : 0; // payoff call
				geometric_fix_put = E > ave_geometric ? E - ave_geometric : 0; // payoff put
				arithmetic_float_call = 1 * S > ave_arithmetic ? 1 * S - ave_arithmetic : 0; // payoff call
				arithmetic_float_put = ave_arithmetic > 1 * S ? ave_arithmetic - 1 * S : 0; // payoff put	
				geometric_float_call = 1 * S > ave_geometric ? 1 * S - ave_geometric : 0; // payoff call
				geometric_float_put = ave_geometric > 1 * S ? ave_geometric - 1 * S : 0; // payoff put

				arithmetic_fix_call_var1[j] = arithmetic_fix_call;
				arithmetic_fix_put_var1[j] = arithmetic_fix_put;
				geometric_fix_call_var1[j] = geometric_fix_call;
				geometric_fix_put_var1[j] = geometric_fix_put;
				arithmetic_float_call_var1[j] = arithmetic_float_call;
				arithmetic_float_put_var1[j] = arithmetic_float_put;
				geometric_float_call_var1[j] = geometric_float_call;
				geometric_float_put_var1[j] = geometric_float_put;

				arithmetic_fix_call_mean += arithmetic_fix_call / n_sim;
				arithmetic_fix_put_mean += arithmetic_fix_put / n_sim;
				geometric_fix_call_mean += geometric_fix_call / n_sim;
				geometric_fix_put_mean += geometric_fix_put / n_sim;
				arithmetic_float_call_mean += arithmetic_float_call / n_sim;
				arithmetic_float_put_mean += arithmetic_float_put / n_sim;
				geometric_float_call_mean += geometric_float_call / n_sim;
				geometric_float_put_mean += geometric_float_put / n_sim;
			}
			//cout << error;

			arithmetic_fix_call_var = variance(arithmetic_fix_call_var1, arithmetic_fix_call_mean,n_sim);
			arithmetic_fix_put_var = variance(arithmetic_fix_put_var1, arithmetic_fix_put_mean, n_sim);
			geometric_fix_call_var = variance(geometric_fix_call_var1, geometric_fix_call_mean, n_sim);
			geometric_fix_put_var = variance(geometric_fix_put_var1, geometric_fix_put_mean, n_sim);
			arithmetic_float_call_var = variance(arithmetic_float_call_var1, arithmetic_float_call_mean, n_sim);
			arithmetic_float_put_var = variance(arithmetic_float_put_var1, arithmetic_float_put_mean, n_sim);
			geometric_float_call_var = variance(geometric_float_call_var1, geometric_float_call_mean, n_sim);
			geometric_float_put_var = variance(geometric_float_put_var1, geometric_float_put_mean, n_sim);

			arithmetic_fix_call_mean = arithmetic_fix_call_mean*exp(-IR*T);  //discount
			arithmetic_fix_put_mean = arithmetic_fix_put_mean*exp(-IR*T);
			geometric_fix_call_mean = geometric_fix_call_mean*exp(-IR*T);
			geometric_fix_put_mean = geometric_fix_put_mean*exp(-IR*T);
			arithmetic_float_call_mean = arithmetic_float_call_mean*exp(-IR*T);
			arithmetic_float_put_mean = arithmetic_float_put_mean*exp(-IR*T);
			geometric_float_call_mean = geometric_float_call_mean*exp(-IR*T);
			geometric_float_put_mean = geometric_float_put_mean*exp(-IR*T);

			if (N == 100 || N == 500){
				cout << "It is discrete sampling and dt =" << dt << " , the number of simulations is " << n_sim << endl;
				cout << "The value of the discrete arithmetic call option (Fixed strike price) is " << endl << arithmetic_fix_call_mean << endl;
				cout << "The value of the discrete arithmetic put option (Fixed strike price) is " << endl << arithmetic_fix_put_mean << endl;
				cout << "The value of the discrete geometric call option (Fixed strike price) is " << endl << geometric_fix_call_mean << endl;
				cout << "The value of the discrete geometric put option (Fixed strike price) is " << endl << geometric_fix_put_mean << endl;
				cout << "The value of the discrete arithmetic call option (Floating strike price) is " << endl << arithmetic_float_call_mean << endl;
				cout << "The value of the discrete arithmetic put option (Floating strike price) is " << endl << arithmetic_float_put_mean << endl;
				cout << "The value of the discrete geometric call option (Floating strike price) is " << endl << geometric_float_call_mean << endl;
				cout << "The value of the discrete geometric put option (Floating strike price) is " << endl << geometric_float_put_mean << endl;
				/* print stand deviation to make error graph of discrete sampling
				cout << "Standard Deviation is " << arithmetic_fix_call_var << endl;
				cout << "Standard Deviation is " << arithmetic_fix_put_var << endl;
				cout << "Standard Deviation is " << geometric_fix_call_var << endl;
				cout << "Standard Deviation is " << geometric_fix_put_var << endl;
				cout << "Standard Deviation is " << arithmetic_float_call_var << endl;
				cout << "Standard Deviation is " << arithmetic_float_put_var << endl;
				cout << "Standard Deviation is " << geometric_float_call_var << endl;
				cout << "Standard Deviation is " << geometric_float_put_var << endl;
				*/
				cout << "-------------------------------------------------------------------------------" << endl;
			}

			if (N == 1000 || N == 10000){
				cout << "It is continuous sampling and dt =" << dt << " , the number of simulations is " << n_sim << endl;
				cout << "The value of the continuous arithmetic call option (Fixed strike price) is " << endl << arithmetic_fix_call_mean << endl;
				cout << "The value of the continuous arithmetic put option (Fixed strike price) is " << endl << arithmetic_fix_put_mean << endl;
				cout << "The value of the continuous geometric call option (Fixed strike price) is " << endl << geometric_fix_call_mean << endl;
				cout << "The value of the continuous geometric put option (Fixed strike price) is " << endl << geometric_fix_put_mean << endl;
				cout << "The value of the continuous arithmetic call option (Floating strike price) is " << endl << arithmetic_float_call_mean << endl;
				cout << "The value of the continuous arithmetic put option (Floating strike price) is " << endl << arithmetic_float_put_mean << endl;
				cout << "The value of the continuous geometric call option (Floating strike price) is " << endl << geometric_float_call_mean << endl;
				cout << "The value of the continuous geometric put option (Floating strike price) is " << endl << geometric_float_put_mean << endl;
				/* print stand deviation to make error graph of continuous sampling
				cout << "Standard Deviation is " << arithmetic_fix_call_var << endl;
				cout << "Standard Deviation is " << arithmetic_fix_put_var << endl;
				cout << "Standard Deviation is " << geometric_fix_call_var << endl;
				cout << "Standard Deviation is " << geometric_fix_put_var << endl;
				cout << "Standard Deviation is " << arithmetic_float_call_var << endl;
				cout << "Standard Deviation is " << arithmetic_float_put_var << endl;
				cout << "Standard Deviation is " << geometric_float_call_var << endl;
				cout << "Standard Deviation is " << geometric_float_put_var << endl;
				*/
				cout <<  "-------------------------------------------------------------------------------" << endl;
			}
		}
	}

	system("pause");
	return 0;
}

double variance(double a[], double mean, double n)
{
	double var = 0;
	for (int h = 1; h<=n; h++)
		var = (a[h] - mean)*(a[h] -mean);
	var = sqrt(var / n);
	return var;
}

double normal(double std, double mean)
{
	static int iset = 0;
	static double gset;
	double fac, r, v1, v2;
	// create two
	if (iset == 0)
	{
		r = 0;
		do
		{
			v1 = 2.0 *  rand() / RAND_MAX - 1.0;
			v2 = 2.0 *  rand() / RAND_MAX - 1.0;
			r = v1*v1 + v2*v2; // they define radius
		} while (r >= 1.0 || r == 0.0);
		// in unit circle? if not try again
		fac = sqrt((-2 * log(r)) / r); // Box-Muller transform 
		gset = (v1 * fac);
		iset = 1; // save one
		v2 = v2*fac*std + mean; // scale and return one
		return v2;
	}
	else{
		iset = 0;
		return (gset*std) + mean;
	}
//scale and return the saved one
}
 


