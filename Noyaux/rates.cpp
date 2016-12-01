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

#define M_ELECTRON 0.511 // MeV
#define HBARC 197.3269788  // MeV.fm 
#define FERMI_COUPLING 1.1663787e-11 // MeV^-2
#define CELERITY 2.99792458e23 // fm/s

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

enum { P_TEMPERATURE, P_DENSITY, P_FRACTION };

inline double Np(int Z)
{
    return min(8, max(0, Z-20));
}


inline double Nh(int N)
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

// from table or analytical expr. ?
double nucleus_mass(int A, int Z)
{
}

double integrated_electron_capture(int A, int Z, double T, int N = 10000)
{
    const double Q = 1.29332 + 3;
    const double Vud = 0.97427;
    const double gA = 1.24;

    double rate = 8*CELERITY/7. * gA*gA * Vud*Vud * pow(2*M_PI, -3.) * FERMI_COUPLING*FERMI_COUPLING 
                  * Np(Z) * Nh(A-Z);
    
    double integral = 0;
    double Emin = M_ELECTRON, Emax = 100*T;
    const double dE = (Emax-Emin)/double(N);
    for(int i = 0; i < N; ++i)
    {
        const double E = Emin + (Emax-Emin)*(double(i)+0.5)/double(N);
        integral += E * E * fermi_dirac(E, T, M_ELECTRON) * (E+Q)*(E+Q) * sqrt(1-(M_ELECTRON/(E+Q))*(M_ELECTRON/(E+Q))) * dE;
    }
    rate *= integral;

    return rate;
}

double fast_electron_capture(int A, int Z, double T, double mu_e = M_ELECTRON)
{
    const double beta = 4;
    const double K = 6146;
    const double Q = 4.6;
    const double dE = 2.5;
    const double chi = (Q-dE)/T;
    const double eta = chi+M_ELECTRON/T; // m -> mu

    double rate = 0.693 * beta/K * pow(T/M_ELECTRON, 5.);
    
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

int nucleus_to_AZ(int nucleus, int &A, int &Z)
{
    A = nucleus/1000;
    Z = nucleus%1000;
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
    if(!ptr) return 0;
    abundance_data ab = *ptr;

    double v = 0;
    for(int i = 0; i < ab.elements.size(); ++i) if(ab.elements[i].A == A && ab.elements[i].Z == Z) { v = ab.elements[i].abundance; break; }

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

int main(int argc, char *argv[])
{
    abundance_table table;
    read_compo_data("EOS.compo", table);

    FILE *fp = fopen("rates.res", "w+");
    for(auto const& ab : table.abundances)
    {
         std::vector<element> &elements = ab.second->elements;
         if(ab.second->param[0] < 1 || ab.second->param[0] > 3) continue;
         double rate = 0, fast_rate = 0;
         double total_abundance = 0;
         for(int i = 0; i < elements.size(); ++i)
         {
             rate += elements[i].abundance * integrated_electron_capture(elements[i].A, elements[i].Z, ab.second->param[0]);
             fast_rate += elements[i].abundance * fast_electron_capture(elements[i].A, elements[i].Z, ab.second->param[0]);
             total_abundance += elements[i].abundance;
         }
         fprintf(fp, "%d %d %d %e %e %e %e %e %e %e\n", ab.first[0], ab.first[1], ab.first[2], ab.second->param[0], ab.second->param[1], ab.second->param[2], rate, fast_rate, rate/total_abundance, fast_rate/total_abundance);
    }
    fclose(fp);
    return 0;
}
