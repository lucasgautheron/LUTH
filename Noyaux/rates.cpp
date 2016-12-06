#include <sstream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <array>
#include <algorithm>
#include <math.h>

#include <gsl/gsl_sf_fermi_dirac.h>

#include "TF1.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

#define M_ELECTRON 0.511 // MeV
#define HBARC 197.3269788  // MeV.fm 
#define FERMI_COUPLING 1.1663787e-11 // MeV^-2
#define CELERITY 2.99792458e23 // fm/s
#define MASS_UNIT 931.4940954 // MeV


enum { P_TEMPERATURE, P_DENSITY, P_FRACTION };

inline int Np(int Z)
{
    return min(8, max(0, Z-20));
}


inline int Nh(int N)
{
    return min(6, max(0, 40-N));
}

double fermi_dirac(double E, double T, double mu)
{
    return 1./(exp((E-mu)/T)+1);
}

/*double electron_capture_rate(int A, int Z, double E)
{
    const double Q = 1.29332 + 3;
    const double Vud = 0.97427;
    const double gA = 1.24;
    double rate = 1./M_PI;
    rate *= Vud * Vud * gA * gA * 2./7.;
    rate *= Np(Z) * Nh(A-Z);
    rate *= 1;
    return rate;
}

double integrated_electron_capture_rate(int A, int Z, double T, int N = 10000)
{
    double Emin = 0, Emax = 100*T;

    double rate = 0;
    const double dE = (Emax-Emin)/double(N);
    for(int i = 0; i < N; ++i)
    {
        const double E = Emin + (Emax-Emin)*(double(i)+0.5)/double(N);
        rate += E * E * fermi_dirac(E, T, M_ELECTRON) * electron_capture_rate(A, Z, E) * dE;
    }
    return rate;
}
*/

struct nuclear_data
{
    int A, Z;
    double m; // MeV
    double beta_q; // beta decay Q

    nuclear_data() {}
    nuclear_data(int _A, int _Z, double _m, double _bq)
    {
        A = _A; Z = _Z; m = _m; beta_q = _bq;
    }
};
typedef std::map< std::array<int, 2>, nuclear_data *> nuclear_array;


nuclear_array nuclear_table;

inline double nucleus_mass(int A, int Z)
{
    std::array<int, 2> elem = { A, Z };
    return (nuclear_table.count(elem)) ? nuclear_table[elem]->m : double(A)*MASS_UNIT;
}

inline double beta_decay_Q(int A, int Z)
{
    std::array<int, 2> mother = { A, Z }, daughter = { A, Z-1 };
    if(!nuclear_table.count(mother)) return -1e10;
    if(abs(nuclear_table[mother]->beta_q) >= 0.0001) return nuclear_table[mother]->beta_q;
    if(!nuclear_table.count(daughter)) return -1e10;
    return nuclear_table[mother]->m - nuclear_table[daughter]->m;
}

double integrated_electron_capture(int A, int Z, double T, double Q = 1e10, double mu = 0)
{
    if(Q>1e6) Q = beta_decay_Q(A, Z) + 3;
    const double Vud = 0.97427;
    const double gA = 1.24;

    int shell_factor = Np(Z) * Nh(A-Z);
    if(!shell_factor) return 0;
    double rate = 8*CELERITY/7. * gA*gA * Vud*Vud * pow(2*M_PI, -3.) * FERMI_COUPLING*FERMI_COUPLING * shell_factor;
    
    double integral = 0;
    const int N = 1000;
    double Emin = max(0, M_ELECTRON-Q), Emax = 30*T;
    const double dE = (Emax-Emin)/double(N);
    for(int i = 0; i < N; ++i)
    {
        const double E = Emin + (Emax-Emin)*(double(i)+0.5)/double(N);
        integral += E * E * fermi_dirac(E, T, mu) * (E+Q)*(E+Q) * sqrt(1-(M_ELECTRON/(E+Q))*(M_ELECTRON/(E+Q))) * dE;
    }
    rate *= integral;

    return rate;
}

double fast_electron_capture(int A, int Z, double T, double Q = 1e10, double mu_e = M_ELECTRON)
{
    const double beta = 4.6;
    const double K = 6146;
    if(Q>1e6) Q = beta_decay_Q(A, Z);
    if(Q<-1e6) return 0;
    const double dE = 2.5;
    const double chi = (Q-dE)/T;
    const double eta = chi+mu_e/T;

    double rate = 0.693 * beta/K * pow(T/M_ELECTRON, 5.);
    
    // gsl fermi integrals include a normalization factor 1/Gamma(j) hence the multiplication by (j-1)!
    rate *= 24*gsl_sf_fermi_dirac_int (4,eta) - 2 * chi * 6*gsl_sf_fermi_dirac_int (3,eta) + chi*chi * 2*gsl_sf_fermi_dirac_int (2, eta);
    return rate;
}

struct element
{
    int nucleus;
    int A,Z;
    double abundance;

    element() : nucleus(0), A(0), Z(0) { }
};

struct abundance_data
{
    int idx[3];
    double param[3];
    double np, nn;
    
    std::vector<element> elements;
};


typedef std::map< std::array<int, 3>, abundance_data *> abundance_array;

struct abundance_table
{
    std::vector<double> parameters[3];
    abundance_array abundances;
};

struct index_entry
{
    int idx;
    double val;
};

inline int nucleus_to_AZ(int nucleus, int &A, int &Z)
{
    A = nucleus/1000;
    Z = nucleus%1000;
}

int read_masses(const char *path, nuclear_array &table)
{
    std::ifstream infile(path);
    double val;
    int count = 0;
    std::string line;
    std::istringstream iss;
    while (std::getline(infile, line)) {
        iss.clear();
        iss.str(line);
        const char *str = line.c_str();
        int Z = atoi(str+11);
        int A = atoi(str+15);
        int base = atoi(str+96);
        double m = double(base) + 1e-6 * strtod(str+100, NULL);
        m *= MASS_UNIT;
        double beta_q = strtod(str+79, NULL)/1e3;
        printf("%d %d %e %.3f\n", Z, A, m, beta_q);
        std::array<int, 2> elem = { A, Z };
        table[elem] = new nuclear_data(A, Z, m, beta_q);
        ++count;
    }
    return count;
}

struct thermo_state
{
    double T, nb, Y, mu_e, lambda;

};

std::vector<thermo_state *> thermo_states[2];

int read_thermo_states(const char *path, std::vector<thermo_state *> &thermo_states)
{
    std::ifstream infile(path);
    double val;
    int count = 0;
    std::string line;
    std::istringstream iss;
    while (std::getline(infile, line)) {
        iss.clear();
        iss.str(line);
        thermo_state *ts = new thermo_state();
        iss >> ts->nb;
        ts->nb /= (1.674e15);
        iss >> ts->Y;
        iss >> val;
        iss >> ts->T;
        iss >> ts->mu_e;
        for(int k = 0; k < 8; ++k) iss >> val; 
        iss >> ts->lambda;
        ++count;
        thermo_states.push_back(ts);
    }
    return count;
}

struct nucleus_capture_data
{
    int A, Z;
    double Q;
    double T;
    double nb_Y; // nb*Y
    double mu_e;
    double lambda;
};

std::vector<nucleus_capture_data *> nucleus_capture;

int read_nucleus_capture(const char *path, std::vector<nucleus_capture_data *> &nucleus_capture)
{
    std::ifstream infile(path);
    double val;
    int count = 0;
    std::string line;
    std::istringstream iss;
    while (std::getline(infile, line)) {
        iss.clear();
        iss.str(line);
        nucleus_capture_data *nc = new nucleus_capture_data();
        iss >> nc->A;
        iss >> nc->Z;
        iss >> nc->Q;
        iss >> nc->T;
        nc->T /= 11.6594202899;
        iss >> nc->nb_Y;
        nc->nb_Y = pow(10., nc->nb_Y)/1.674e15;
        iss >> nc->mu_e;
        nc->mu_e += M_ELECTRON;
        for(int k = 0; k < 1; ++k) iss >> val; 
        iss >> nc->lambda;
        nc->lambda = pow(10., nc->lambda);
        ++count;
        nucleus_capture.push_back(nc);
    }
    return count;
}

int read_index(const char *path, std::vector<double> &index)
{
    std::ifstream infile(path);
    double val;
    int count = 0;
    std::string line;
    std::istringstream iss;
    while (std::getline(infile, line)) {
        iss.clear();
        iss.str(line);
        iss >> val;
        ++count;
        if(count <= 2) continue;
        index.push_back(val);
    }
    return count-2;
}

int read_compo_data(const char *path, abundance_table &table)
{
    read_index("EOS.t", table.parameters[0]);
    read_index("EOS.nb", table.parameters[1]);
    read_index("EOS.yq", table.parameters[2]);

    std::ifstream infile(path);
    
    int int_ph;
    double double_ph;
    int count = 0;
    std::string line;
    std::istringstream iss;
    while (std::getline(infile, line)) {
        abundance_data *ab = new abundance_data();
        iss.clear();
        iss.str(line);
        for(int k = 0; k < 3; ++k) iss >> ab->idx[k];
        iss >> int_ph;
        iss >> int_ph;
        iss >> int_ph;
        iss >> ab->nn;
        iss >> int_ph;
        iss >> ab->np;

        for(int k = 0; k < 3; ++k) ab->param[k] = table.parameters[k][ab->idx[k]];

        while(true)
        {
            element e;
            iss >> e.nucleus;
            iss >> e.abundance;
            if(!e.nucleus) break;

            nucleus_to_AZ(e.nucleus, e.A, e.Z);
            ab->elements.push_back(e);
            //printf("%d %d %f %e %e %e\n", e.A, e.Z, ab->param[0], integrated_electron_capture_rate(e.A, e.Z, ab->param[0]), integrated_electron_capture_rate(e.A, e.Z, ab->param[0], 100), integrated_electron_capture_rate(e.A, e.Z, ab->param[0], 10000));
        }

        element neutron, proton;
        neutron.A = 1;
        neutron.Z = 0;
        neutron.abundance = ab->nn;
        proton.A = 1;
        proton.Z = 1;
        proton.abundance = ab->np;
        ab->elements.push_back(neutron);
        ab->elements.push_back(proton);

        

        std::array<int, 3> conditions = {ab->idx[0], ab->idx[1], ab->idx[2]};
        table.abundances[conditions] = ab;
        ++count;
    }
    return count;
}

int get_lower_key(std::vector<double> &arr, double value)
{
    auto it = std::lower_bound(arr.begin(), arr.end(), value);
    if(it == arr.end()) return -1;
    return it-arr.begin()-1;
}

double element_abundance_interp(abundance_table &table, int A, int Z, double params[3], double *vlow = NULL, double *vhigh = NULL)
{
    int keys[3];
    for(int k = 0; k < 3; ++k) keys[k] = get_lower_key(table.parameters[k], params[k]);
 
    std::array<int, 3> lower = { keys[0], keys[1], keys[2] }, upper = { keys[0]+1, keys[1]+1, keys[2]+1 };
    std::array<int, 3> next[3] = { { keys[0]+1, keys[1], keys[2] }, { keys[0], keys[1]+1, keys[2] }, { keys[0], keys[1], keys[2]+1 } };
    double dv[3], dx[3], x[3];

    abundance_data *ptr = table.abundances[lower];
    double v = 0;
    if(ptr)
    {
        abundance_data ab = *ptr;
        for(int i = 0; i < ab.elements.size(); ++i) if(ab.elements[i].A == A && ab.elements[i].Z == Z) { v = ab.elements[i].abundance; break; }
    }

    for(int k = 0; k < 3; ++k)
    {
        dx[k] = keys[k]+1 < table.parameters[k].size() ? table.parameters[k][keys[k]+1]-table.parameters[k][keys[k]] : 1;
        dv[k] = 0;
        x[k] = table.parameters[k][keys[k]];

        if(!table.abundances.count(next[k])) continue;

        abundance_data ab_next = *table.abundances[next[k]];
        for(int i = 0; i < ab_next.elements.size(); ++i) if(ab_next.elements[i].A == A && ab_next.elements[i].Z == Z) { dv[k] = ab_next.elements[i].abundance - v; break; }
    }

    double v_next = 0;
    {
        abundance_data *ptr = table.abundances[upper];
        if(!ptr) v_next = v;
        else
        {
            abundance_data ab = *ptr;
            for(int i = 0; i < ab.elements.size(); ++i) if(ab.elements[i].A == A && ab.elements[i].Z == Z) { v_next = ab.elements[i].abundance; break; }
        }
    }

    double vals[] = { v+dv[0], v+dv[1], v+dv[2], v_next };

    if(vlow)
    {
        *vlow = *std::min_element(std::begin(vals), std::end(vals));
    }
    if(vhigh)
    {
        *vhigh = *std::max_element(std::begin(vals), std::end(vals));
    }

    for(int k = 0; k < 3; ++k) v += (params[k]-x[k]) * dv[k] / dx[k];

    return v;
}

double electron_potential(double nb, double Ye)
{
    //return M_ELECTRON*sqrt(1+pow(nb*Ye*1.143e9, 2./3.));
    return M_ELECTRON*sqrt(1+pow(nb*Ye*1.705199692e9, 2./3.));
}

int main(int argc, char *argv[])
{
    // read abundance table
    abundance_table table;
    read_compo_data("EOS.compo", table);
    
    // read mass table (for Q values)
    read_masses("mass.mas12", nuclear_table);

    // read CCSN trajectories
    read_thermo_states("trajectory_15", thermo_states[0]);
    read_thermo_states("trajectory_25", thermo_states[1]);

    // read single rates
    read_nucleus_capture("single_rates", nucleus_capture);
    
    FILE *fp = fopen("potential.res", "w+");
    for(int i = 0; i < 1000; ++i)
    {
        double n = pow(10., -8+7.*double(i)/1000.);
        fprintf(fp, "%e %e\n", n*1.674e-24*1e39, electron_potential(n, 1.));
    }
    fclose(fp);

    fp = fopen("single_rates.res", "w+");
    for(int i = 0; i < nucleus_capture.size(); ++i)
    {
        nucleus_capture_data *nc = nucleus_capture[i];
        if(nc->A < 2) continue;
        if(nc->T < 0.1) continue;
        if(nc->nb_Y < 1e-7) continue;

        fprintf(fp, "%d %d %e %e %e %e %e %e\n", nc->A, nc->Z, nc->T, nc->nb_Y, nc->lambda,
		fast_electron_capture(nc->A, nc->Z, nc->T, 1e10, electron_potential(nc->nb_Y, nc->mu_e)), fast_electron_capture(nc->A, nc->Z, nc->T, 1e10, electron_potential(nc->nb_Y, 1.)),
		integrated_electron_capture(nc->A, nc->Z, nc->T, 1e10, 0));
    }
    fclose(fp);

    fp = fopen("rates.res", "w+");

    const double Q[2] = { -18, 16 };
    for(int i = 0; i < 1000; ++i)
    {
        double _Q = Q[0]+(Q[1]-Q[0])*double(i)/1000.;
        double mu_e[2] = { electron_potential(1.32e-6, 0.447), electron_potential(1.12e-4, 0.361) };
        fprintf(fp, "%f %e %e %e %e\n", _Q, integrated_electron_capture(39, 21, 0.68, _Q, 0), integrated_electron_capture(39, 21, 1.3, _Q, 0), fast_electron_capture(39, 21, 0.68, _Q, mu_e[0]), fast_electron_capture(39, 21, 1.3, _Q, mu_e[1]));
    }
    fclose(fp);

    for(int k = 0; k < 2; ++k) {
            fp = fopen(k == 0 ? "total_rates_15.res" : "total_rates_25.res", "w+");
	    for(int i = 0; i < thermo_states[k].size(); ++i)
	    {
		thermo_state *ts = thermo_states[k][i];
		double rate = 0, fast_rate = 0, fast_rate_corr = 0;
		double mu_e = electron_potential(ts->nb, ts->Y);
		double total_abundance = 0, abundance_error = 0;
		double conditions[3] = {ts->T, ts->nb, ts->Y};

		for(int A = 2; A < 250; ++A) for(int Z = 0; Z <= A; ++Z) if(A >= 20)
		{
		    //if(abs(A-Z) > 50) continue;
		    double vl = 0, vh = 0;
		    double abundance = element_abundance_interp(table, A, Z, conditions, &vl, &vh);
		    total_abundance += abundance;
		    abundance_error += (vh-vl)*(vh-vl);
		    rate += abundance * integrated_electron_capture(A, Z, ts->T, 1e10, 0);
		    fast_rate += abundance * fast_electron_capture(A, Z, ts->T, 1e10, mu_e);
		    fast_rate_corr += abundance * fast_electron_capture(A, Z, ts->T, 1e10, ts->mu_e);
		}
	printf("%e %e %e %e %e %e %e %e\n", ts->T, ts->nb, ts->Y, ts->lambda, rate, fast_rate, rate/total_abundance, fast_rate/total_abundance);
		fprintf(fp, "%e %e %e %e %e %e %e %e %e\n", ts->T, ts->nb, ts->Y, ts->lambda, rate, fast_rate, rate/total_abundance, fast_rate/total_abundance, fast_rate * sqrt(abundance_error)/total_abundance / total_abundance);
	    }
        fclose(fp);
    }
    /*for(auto const& ab : table.abundances)
    {
         std::vector<element> &elements = ab.second->elements;
         if(ab.second->param[0] < 1 || ab.second->param[0] > 3) continue;
         if(ab.second->param[1] < 1e-6) continue;
         if(ab.second->param[2] < 0.25 || ab.second->param[2] > 0.5) continue;
         double rate = 0, fast_rate = 0;
         double total_abundance = 0;
         for(int i = 0; i < elements.size(); ++i)
         {
             double mu_e = electron_potential(ab.second->param[1], ab.second->param[2]);
             rate += elements[i].abundance * integrated_electron_capture(elements[i].A, elements[i].Z, ab.second->param[0], 1e10, 0);
             fast_rate += elements[i].abundance * fast_electron_capture(elements[i].A, elements[i].Z, ab.second->param[0], 1e10, mu_e);
             total_abundance += elements[i].abundance;
         }
         fprintf(fp, "%d %d %d %e %e %e %e %e %e %e\n", ab.first[0], ab.first[1], ab.first[2], ab.second->param[0], ab.second->param[1], ab.second->param[2], rate, fast_rate, rate/total_abundance, fast_rate/total_abundance);
    }
    fclose(fp);*/
    return 0;
}
