#include <iostream>
#include <cmath>
#include <cstdlib> //rand
#include <vector>
#include <time.h>
#include <fstream>


using namespace std;

// ================ tools =================
double max(double x, double y)
{
    return x > y ? x : y;
}

double uniform_draw()
{
    return (double) rand() / (RAND_MAX);
}

double gaussian_draw(double mu = 0, double sigma = 1)
{
    double u = ((double) rand() / (RAND_MAX));
    double v = ((double) rand() / (RAND_MAX));
    return mu + sigma * ( sqrt(-2*log(u)) * sin(2*M_PI*v) );
}







// ============ gaussian distribution =============
class Gaussian_distribution {
public:
    // Constructor, destructor
    Gaussian_distribution() {}
    // virtual ~Gaussian_distribution();
    
    // Distribution functions
    virtual double pdf(const double& x) const;
    virtual double cdf(const double& x) const;
    
    // Inverse cumulative distribution function (aka the probit function)
    virtual double inv_cdf(const double& quantile) const;
    
    // Descriptive statistics
    virtual double mean() const;   // equal to 0
    virtual double var() const;    // equal to 1
    virtual double stdev() const;  // equal to 1
    
    // Obtain a sequence of random draws from the standard normal distribution
    virtual void random_draws(const vector<double>& uniform_draws, vector<double>& dist_draws);
};

// Probability density function
double Gaussian_distribution::pdf(const double& x) const {
    return (1.0/sqrt(2.0 * M_PI)) * exp(-0.5*x*x);
}

// Cumulative density function
double Gaussian_distribution::cdf(const double& x) const {
    double k = 1.0/(1.0 + 0.2316419*x);
    double k_sum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));
    
    if (x >= 0.0) {
        return (1.0 - (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x) * k_sum);
    } else {
        return 1.0 - cdf(-x);
    }
}

// Inverse cumulative distribution function (probit function)
double Gaussian_distribution::inv_cdf(const double& quantile) const {
    // Beasley-Springer-Moro algorithm - Glasserman [2004]
    static double a[4] = {   2.50662823884, -18.61500062529, 41.39119773534, -25.44106049637};
    
    static double b[4] = {  -8.47351093090, 23.08336743743, -21.06224101826, 3.13082909833};
    
    static double c[9] = {0.3374754822726147, 0.9761690190917186, 0.1607979714918209,
                          0.0276438810333863, 0.0038405729373609, 0.0003951896511919,
                          0.0000321767881768, 0.0000002888167364, 0.0000003960315187};
    
    if (quantile >= 0.5 && quantile <= 0.92) {
        double num = 0.0;
        double denom = 1.0;
        
        for (int i=0; i<4; i++) {
            num += a[i] * pow((quantile - 0.5), 2*i + 1);
            denom += b[i] * pow((quantile - 0.5), 2*i);
        }
        return num/denom;
        
    } else if (quantile > 0.92 && quantile < 1) {
        double num = 0.0;
        
        for (int i=0; i<9; i++) {
            num += c[i] * pow((log(-log(1-quantile))), i);
        }
        return num;
        
    } else {
        return -1.0*inv_cdf(1-quantile);
    }
}

// Expectation/mean
double Gaussian_distribution::mean() const { return 0.0; }

// Variance
double Gaussian_distribution::var() const { return 1.0; }

// Standard Deviation
double Gaussian_distribution::stdev() const { return 1.0; }

void Gaussian_distribution::random_draws(const vector<double>& uniform_draws, vector<double>& dist_draws) {
    // The simplest method to calculate this is with the Box-Muller method
    if (uniform_draws.size() != dist_draws.size())
    {
        cout << "Draws vectors are of unequal size in standard normal dist." << endl;
    }
    if (uniform_draws.size() % 2 != 0)
    {
        cout << "Uniform draw vector size not an even number." << endl;
    }
    for (int i=0; i<uniform_draws.size() / 2; i++) {
        dist_draws[2*i] = sqrt(-2.0*log(uniform_draws[2*i])) * sin(2*M_PI*uniform_draws[2*i+1]);
        dist_draws[2*i+1] = sqrt(-2.0*log(uniform_draws[2*i])) * cos(2*M_PI*uniform_draws[2*i+1]);
    }
}











// ================ variance swap ==================
class Variance_swap {
public:
    Variance_swap(double maturity, double n = 0 , double sigma = 0):
    T_(maturity), notional_(n), variance_strike_(sigma) {}
    
    virtual ~Variance_swap();
    
private:
    double T_;
    double notional_;
    double variance_strike_;
};


// ================= heston model ==================
class Heston_model {
public:
    Heston_model(double r, double k, double t, double s, double rho):
    rate_(r), kappa_(k), theta_(t), sigma_(s), rho_(rho) {}
    
    // virtual ~Heston_model();
    
    // getters
    double get_rate()  {return rate_;}
    double get_kappa() {return kappa_;}
    double get_theta() {return theta_;}
    double get_sigma() {return sigma_;}
    double get_rho()   {return rho_;}
    
    // methods
    void correlated_draws(vector<double>& uncorr_draws, vector<double>& dist_draws);
    void vol_path_method(double T, const vector<double>& vol_draws, vector<double>& vol_path);
    void spot_path_method(double T, const vector<double>& spot_draws, const vector<double>& vol_path, vector<double>& spot_path);

private:
    double rate_;
    double kappa_;
    double theta_;
    double sigma_;
    double rho_;
};

void Heston_model::correlated_draws(vector<double>& uncorr_draws, vector<double>& corr_draws)
{
    for (int i(1); i<corr_draws.size(); i++)
    {
        uncorr_draws[i] = uncorr_draws[i-1] + gaussian_draw();
        corr_draws[i] = get_rho() * uncorr_draws[i] + corr_draws[i] * sqrt(1-get_rho()*get_rho());
    }
}

void Heston_model::vol_path_method(double T, const vector<double>& vol_draws, vector<double>& vol_path)
{
    size_t vec_size = vol_draws.size();
    double dt = T/static_cast<double>(vec_size);

    for (int i(1); i<vec_size; i++)
    {
        double v_max = max(vol_path[i-1], 0.0);
        vol_path[i] = vol_path[i-1] + get_kappa() * dt * (get_theta() - v_max) + get_sigma() * sqrt(v_max * dt) * vol_draws[i-1];
    }
}

void Heston_model::spot_path_method(double T, const vector<double>& spot_draws, const vector<double>& vol_path, vector<double>& spot_path)
{
    size_t vec_size = spot_draws.size();
    double dt = T/static_cast<double>(vec_size);
    
    for (int i(1); i<vec_size; i++)
    {
        double v_max = max(vol_path[i-1], 0.0);
        spot_path[i] = spot_path[i-1] * exp( (get_rate() - 0.5*v_max)*dt + sqrt(v_max*dt)*spot_draws[i-1]);
    }
}










// ============= monte carlo pricing ================
class Monte_carlo_pricing : public Heston_model {
public:
    Monte_carlo_pricing(double r, double k, double t, double s, double rho):
    Heston_model(r, k, t, s, rho) {}
    
    // virtual ~Monte_carlo_pricing();
    
    double m_function(double Vt_hat, double delta);
    double s_square_function(double Vt_hat, double delta);
    double psi_function(double s_square, double m);
    
    virtual double v_hat(double Vt_hat, double tau) = 0;
    double log_X_hat(double Xt_hat, double vt_hat, double tau);
    double X_hat(double Xt_hat, double vt_hat, double tau);
};

double Monte_carlo_pricing::m_function(double Vt_hat, double delta)
{
    return get_theta() + (Vt_hat - get_theta()) * exp(-get_kappa()*delta);
}

double Monte_carlo_pricing::s_square_function(double Vt_hat, double delta)
{
    double C = exp(-get_kappa()*delta);
    return (Vt_hat * pow(get_sigma(),2) * C * (1-C) ) /get_kappa() + get_theta() * pow(get_sigma()*(1-C), 2)/(2*get_kappa());
}

double Monte_carlo_pricing::psi_function(double s_square, double m)
{
    return s_square/(m*m);
}

double Monte_carlo_pricing::log_X_hat(double Xt_hat, double vt_hat, double tau)
{
    // first step - compute v_hat
    double vt_delta_hat = v_hat(vt_hat, tau);
    // second step - draw a uniform random number
    double u = uniform_draw();
    // third step - generate Z_v
    Gaussian_distribution gauss_dist;
    double Z = gauss_dist.inv_cdf(u);
    // final step - X_hat(t+delta)
    double cste = get_rho() * get_kappa() * get_theta() / get_sigma();
    double K0 = - cste * tau;
    double gamma1 = 0.5;
    double gamma2 = 0.5;
    double K1 = gamma1 * tau * (cste - 0.5) - get_rho()/get_sigma();
    double K2 = gamma2 * tau * (cste - 0.5) + get_rho()/get_sigma();
    double K3 = gamma1 * tau * (1 - get_rho()*get_rho());
    double K4 = gamma2 * tau * (1 - get_rho()*get_rho());
    
    return log(Xt_hat) + K0 + K1*vt_hat + K2*vt_delta_hat + sqrt(K3*vt_hat + K4*vt_delta_hat)*Z;
}

double Monte_carlo_pricing::X_hat(double Xt_hat, double vt_hat, double tau)
{
    // first step - compute v_hat
    double vt_delta_hat = v_hat(vt_hat, tau);
    // second step - draw a uniform random number
    double u = uniform_draw();
    // third step - generate Z_v
    Gaussian_distribution gauss_dist;
    double Z = gauss_dist.inv_cdf(u);
    // final step - X_hat(t+delta)
    double cste = get_rho() * get_kappa() * get_theta() / get_sigma();
    double K0 = - cste * tau;
    double gamma1 = 0.5;
    double gamma2 = 0.5;
    double K1 = gamma1 * tau * (cste - 0.5) - get_rho()/get_sigma();
    double K2 = gamma2 * tau * (cste - 0.5) + get_rho()/get_sigma();
    double K3 = gamma1 * tau * (1 - get_rho()*get_rho());
    double K4 = gamma2 * tau * (1 - get_rho()*get_rho());
    
    return Xt_hat * exp( K0 + K1*vt_hat + K2*vt_delta_hat + sqrt(K3*vt_hat + K4*vt_delta_hat)*Z );
}








// ============= TG scheme ==============
class TG_scheme : public Monte_carlo_pricing {
public:
    TG_scheme(double r, double k, double t, double s, double rho):
    Monte_carlo_pricing(r, k, t, s, rho) {}
    
    // virtual ~TG_scheme();
    
    // scheme for v
    double f_mu(double Psi);
    double f_sigma(double Psi);
    double mu(double Psi, double m);
    double sigma(double Psi, double s);
    double v_hat(double Vt_hat, double tau);
    
    // expectation
    double expectation(double A, double mu, double sigma);
};

double TG_scheme::f_mu(double Psi)
{
    return 1;
}

double TG_scheme::f_sigma(double Psi)
{
    return 1;
}

double TG_scheme::mu(double Psi, double m)
{
    return f_mu(Psi) * m;
}

double TG_scheme::sigma(double Psi, double s)
{
    return f_sigma(Psi) * s;
}

double TG_scheme::v_hat(double Vt_hat, double tau)
{
    // first step - compute m and s
    double m = m_function(Vt_hat, tau);
    double s_square = s_square_function(Vt_hat, tau);
    
    // second step - Psi
    double Psi = psi_function(s_square, m);
    
    // third step - compute mu and sigma
    double mu_ = mu(Psi, m);
    double sig_ = sigma(Psi, sqrt(s_square));
    
    // fourth step - draw a uniform random number
    double u = uniform_draw();
    
    // fifth step - generate Z_v
    Gaussian_distribution gauss_dist;
    double Z_v = gauss_dist.inv_cdf(u);
    
    // final step - V_hat(t+delta)
    return max(0, mu_ + sig_ * Z_v);
}

double TG_scheme::expectation(double A, double mu, double sigma)
{
    double d_plus = mu/sigma + A * sigma;
    double d_minus = mu/sigma;
    Gaussian_distribution gauss_dist;
    return exp(A*mu + pow(A*sigma, 2) / 2) * gauss_dist.cdf(d_plus) + gauss_dist.cdf( - d_minus);
}








// ============= QE scheme ==============
class QE_scheme : public Monte_carlo_pricing {
public:
    QE_scheme(double r, double k, double t, double s, double rho):
    Monte_carlo_pricing(r, k, t, s, rho) {}
    
    // virtual ~QE_scheme();
    
    double a_function(double m, double b_square);
    double b_square(double Psi);
    double p_function(double Psi);
    double beta_p(double p, double m);
    double beta_psi(double Psi, double m);
    double Psi_inversed(double u, double p, double beta);
    double v_hat_inf(double a, double b, double Z_v);
    double v_hat_sup(double p, double beta);
    double v_hat(double Vt_hat, double delta);
    
    // expectation
    double expectation(double A, double Vt_hat, double delta);
};


double QE_scheme::a_function(double m, double b_square)
{
    return m / (1 + b_square);
}

double QE_scheme::b_square(double Psi)
{
    double Psi_1(1.0 / Psi);
    double b_sq = 2 * Psi_1 - 1 + sqrt( 2 * Psi_1 * (2 * Psi_1 - 1) );
    if (b_sq > 0.0) { return b_sq; }
    else return 0.0;
}

double QE_scheme::p_function(double Psi)
{
    return (Psi - 1) / (Psi + 1);
}

double QE_scheme::beta_p(double p, double m)
{
    return (1 - p) / m;
}

double QE_scheme::beta_psi(double Psi, double m)
{
    return 2 / (m * (Psi + 1));
}

double QE_scheme::Psi_inversed(double u, double p, double beta)
{
    if (p < u <= 1.0)
    {
        return log( (1 - p)/(1 - u) ) / beta;
    }
    else return 0.0;
}

double QE_scheme::v_hat_inf(double a, double b_square, double Z_v)
{
    return a * pow(sqrt(b_square) + Z_v, 2);
}

double QE_scheme::v_hat_sup(double p, double beta)
{
    double u = uniform_draw();
    return Psi_inversed(u, p, beta);
}

double QE_scheme::v_hat(double Vt_hat, double delta)
{
    double Psi_c = 1.5;
    double m = m_function(Vt_hat, delta);
    double s_sq = s_square_function(Vt_hat, delta);
    double Psi = psi_function(s_sq, m);
   
    if (Psi <= Psi_c)
    {
        double b_sq = b_square(Psi);
        double a = a_function(m, b_sq);
        double u = uniform_draw();
        Gaussian_distribution gauss_dist;
        double Z_v = gauss_dist.inv_cdf(u);
        
        return v_hat_inf(a, b_sq, Z_v);
    }
    else // Psi > Psi_c
    {
        double p = p_function(Psi);
        double beta = beta_p(p, m);
        
        return v_hat_sup(p, beta);
    }
}

double QE_scheme::expectation(double A, double Vt_hat, double delta)
{
    double Psi_c = 1.5;
    double m = Monte_carlo_pricing::m_function(Vt_hat, delta);
    double s_sq = s_square_function(Vt_hat, delta);
    double Psi = psi_function(s_sq, m);
    
    if (Psi <= Psi_c)
    {
        double b_sq = b_square(Psi);
        double a = a_function(m, b_sq);
        
        return exp( (A*b_sq*a)/(1 - 2*A*a) ) / sqrt(1 - 2*A*a);
    }
    else // Psi > Psi_c
    {
        double p = p_function(Psi);
        double beta = beta_p(p, m);
        
        return beta * (1 - p) / (beta - A);
    }
}








// ============= analytical pricing ================
class Analytical_pricing : public Heston_model {
public:
    Analytical_pricing(double r, double k, double t, double s, double rho):
    Heston_model(r, k, t, s, rho) {}
    
    //Analytical_pricing(Heston_model heston_model):
    //Heston_model(heston_model) {}
    
    // virtual ~Analytical_pricing();
    
    // methods
    double C_function(double tau, double w);
    double D_function(double tau, double w);
    double C_derivative_1(double tau, double h);
    double D_derivative_1(double tau, double h);
    double C_derivative_2(double tau, double h);
    double D_derivative_2(double tau, double h);
    double G_function(double v, double tau, double h = 0.0001);
    double G_function_from_0(double v, double ti_1, double h = 0.0001);
    
    double a_function(double w);
    double b_function(double w);
};

// We should consider complex number if we strickly follow the analytical model,
// however because w is very small (it is h<<10^-3) we consider complex as a double
// Indeed, w is useful only when we consider derivatives with a w=h negligeable

double Analytical_pricing::a_function(double w)
{
    return get_kappa() - get_rho() * get_sigma() * w;
}

double Analytical_pricing::b_function(double w)
{
    double a(a_function(w));
    return sqrt( a*a + get_sigma()*get_sigma() * (w + w*w) );
}

double Analytical_pricing::C_function(double tau, double w)
{
    double a(a_function(w));
    double b(b_function(w));
    double cste(get_theta()*get_kappa()/pow(get_sigma(), 2));
    double g ( (a - b)/(a + b) );
    return (w - 1)*get_rate()*tau + cste*( (a - b)*tau - 2*log( (1 - g * exp(-b*tau) )/ (1 - g)) );
}

double Analytical_pricing::D_function(double tau, double w)
{
    double a(a_function(w));
    double b(b_function(w));
    double g( (a - b)/(a + b) );
    return ( (a - b) * (1 - exp(-b*tau)) ) / ( get_sigma()*get_sigma() * (1 - g * exp(-b*tau)) );
}

// Newton scheme first order
double Analytical_pricing::C_derivative_1(double tau, double h)
{
    return (C_function(tau, h) - C_function(tau, -h))/2*h;
}

double Analytical_pricing::D_derivative_1(double tau, double h)
{
    return (D_function(tau, h) - D_function(tau, -h))/2*h;
}

// Newton scheme second order
double Analytical_pricing::C_derivative_2(double tau, double h)
{
    return (C_function(tau, h) - 2*C_function(tau, 0) + C_function(tau, -h))/ (h*h);
}

double Analytical_pricing::D_derivative_2(double tau, double h)
{
    return (D_function(tau, h) - 2*D_function(tau, 0) + D_function(tau, -h))/ (h*h);
}

double Analytical_pricing::G_function(double v, double tau, double h)
{
    double C_1 = C_derivative_1(tau, h);
    double C_2 = C_derivative_2(tau, h);
    double D_1 = D_derivative_1(tau, h);
    double D_2 = D_derivative_2(tau, h);
    
    return - exp(-get_rate()*tau) * (pow(D_1*v,2) + (2*C_1*D_1 + D_2)*v + pow(C_1,2) + C_2);
}

double Analytical_pricing::G_function_from_0(double v, double ti_1, double h)
{
    double C_1 = C_derivative_1(ti_1, h);
    double C_2 = C_derivative_2(ti_1, h);
    double D_1 = D_derivative_1(ti_1, h);
    double D_2 = D_derivative_2(ti_1, h);
    
    double ci = 2 * get_kappa() / ( get_sigma()*get_sigma() * (1 - exp( - get_kappa() * ti_1)) );
    double q = 2 * get_kappa() * get_theta() / (get_sigma()*get_sigma());
    double wi = ci * v * exp(-get_kappa()*ti_1);
    
    return -pow(D_1,2)*( q+2*wi + pow(q+wi,2) ) /pow(ci,2) - (2*C_1*D_1 + D_2)*(q+wi)/ci - pow(C_1,2) + C_2;
}








// ======================= main ==========================

void write(string name, vector<double> data)
{
    fstream f ("/Users/user/Documents/Code/C++/Heston/_HestonModel/data/"+name+".txt", ios::out/*|ios::app*/);
    for (int i(0); i<data.size(); i++) { f << i << " " << data[i] << "\n";}
    f.close();
}


int main(int argc, char **argv)
{
    // === Initialization ===
    srand((unsigned)time(0));
    
    double T(3.0);
    double N(100.0);
    double S_0(100.0);
    double v_0(0.04);
    
    double rate(0.05);
    double kappa(0.3);
    double theta(0.04);
    double sigma(0.9);
    double rho(-0.5);
    
    double tau(1.0/8.0); // ti - ti-1

    
    // === Simulation ===
    
    Heston_model Hm(rate, kappa, theta, sigma, rho);
    vector<double> Brownian(N);
    vector<double> Brownian_2(N);
    Hm.correlated_draws(Brownian, Brownian_2);
    
    vector<double> Vol_path(N, v_0);
    Hm.vol_path_method(T, Brownian, Vol_path);
    
    vector<double> Spot_path(N, S_0);
    Hm.spot_path_method(T, Brownian_2, Vol_path, Spot_path);
    
    
    // === Analytical pricing ===
    Analytical_pricing Ap(rate, kappa, theta, sigma, rho);
    vector<double> A_Exceptation(N);
    for (int i(0); i<A_Exceptation.size(); i++) { A_Exceptation[i] = Ap.G_function(Vol_path[i], tau); }
    
    write("data_brownian", Brownian);
    write("data_brownian_2", Brownian_2);
    write("data_anal_vol", Vol_path);
    write("data_anal_spot", Spot_path);
    write("data_anal_exceptation", A_Exceptation);
    
    
    // === Truncated Gaussian scheme ===
    TG_scheme tgs (rate, kappa, theta, sigma, rho);
    vector<double> TG_v_hat(N, v_0);
    for (int i(1); i<TG_v_hat.size(); i++) { TG_v_hat[i] = tgs.v_hat(TG_v_hat[i-1], tau); }
    
    vector<double> TG_X_hat(N, S_0);
    for (int i(1); i<TG_X_hat.size(); i++) { TG_X_hat[i] = tgs.X_hat(TG_X_hat[i-1], TG_v_hat[i-1], tau); }
    
    vector<double> TG_X_hat_bis(N, S_0);
    for (int i(1); i<TG_X_hat_bis.size(); i++) { TG_X_hat_bis[i] = exp( tgs.log_X_hat(TG_X_hat_bis[i-1], TG_v_hat[i-1], tau) ); }
    
    write("data_TG_vol", TG_v_hat);
    write("data_TG_spot", TG_X_hat);
    write("data_TG_spot_bis", TG_X_hat_bis);
    // write("data_TG_exceptation", TG_Exceptation);
    
    
    // === Quadratic Exponential scheme ===
    QE_scheme qes (rate, kappa, theta, sigma, rho);
    vector<double> QE_v_hat(N, v_0);
    for (int i(1); i<QE_v_hat.size(); i++) { QE_v_hat[i] = qes.v_hat(QE_v_hat[i-1], tau); }
    
    vector<double> QE_X_hat(N, S_0);
    for (int i(1); i<QE_X_hat.size(); i++) { QE_X_hat[i] = qes.X_hat(QE_X_hat[i-1], QE_v_hat[i-1], tau); }
    
    write("data_QE_vol", QE_v_hat);
    write("data_QE_spot", QE_X_hat);
    // write("data_QE_exceptation", QE_Exceptation);
    
    cout << " analytic " << A_Exceptation[N] << endl;
    cout << " TG " << tgs.expectation(1.0, Spot_path[N], tau) << endl;
    cout << " QE " << qes.expectation(1.0, Spot_path[N], tau) << endl;
    cout << " QE bis " << qes.expectation(10.0, Spot_path[N], tau) << endl;
    
    return 0;
}
