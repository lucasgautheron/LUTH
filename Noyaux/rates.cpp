#include <sstream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <array>
#include <algorithm>

#include "TF1.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"

enum { P_TEMPERATURE, P_DENSITY, P_FRACTION };

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

typedef  std::map< std::array<int, 3>, abundance_data *> abundance_array;

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

        while(true)
        {
            element e;
            iss >> e.nucleus;
            iss >> e.abundance;
            if(!e.nucleus) break;

            nucleus_to_AZ(e.nucleus, e.A, e.Z);
            ab->elements.push_back(e);
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

        for(int k = 0; k < 3; ++k) ab->param[k] = table.parameters[k][ab->idx[k]];

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

    int N = 100;
    double Ymin = 0, Ymax = 0.6;
    double *x = new double[N], *y = new double[N], *z = new double[N], *exl = new double[N], *eyl = new double[N], *exh = new double[N], *eyh = new double[N];
    double *cmp = new double[N];
    
    const int A = 4;
    const int Z = 2;
    const double T = 1.1, nb = 1.08343634e-06;
    for(int i = 0; i < N; ++i)
    {
        x[i] = Ymin+(Ymax-Ymin)*((double(i)+0.5)/double(N));
        double params[3] = { T, nb, x[i] };
        double vlow, vhigh;
        exl[i] = exh[i] = 0;
        y[i] = element_abundance_interp(table, A, Z, params, &vlow, &vhigh );
        eyl[i] = y[i]-vlow;
        eyh[i] = vhigh-y[i];
        printf("%e %e %e\n", vlow, y[i], vhigh);
        z[i] = element_abundance(table, A, Z, T, nb, x[i]);
        printf("%e %e\n", x[i], y[i]);
        cmp[i] = NSE(T*11.5942028986, nb*1.67e-24*1e39, x[i], Z, A-Z);
        printf("%e %e %e\n", T*11.5942028986, nb*1.67e-24*1e39, cmp[i]);
    }

    TCanvas *c = new TCanvas("c","data",200, 10, 700, 500);
    TGraph *gr_interp, *gr_raw, *gr_cmp;

    gr_interp = new TGraphAsymmErrors(N,x,y,exl,exh,eyl,eyh);;
    gr_interp->SetTitle("interpolation;Y;n(78,28)");
    gr_interp->SetMarkerStyle(3);
    gr_interp->SetMarkerSize(0.3);
    gr_interp->SetMarkerColor(kRed);
    gr_interp->Draw("APZL");

    gr_raw = new TGraph(N,x,z);
    gr_raw->Draw("PZsame");
    gr_raw->SetMarkerColor(kBlue);

    gr_cmp = new TGraph(N,x,cmp);
    gr_cmp->SetMarkerColor(kGreen);
    gr_cmp->SetMarkerSize(1);
    gr_cmp->SetMarkerStyle(3);
    gr_cmp->Draw("PZsame");

    c->Print("n_78_28.pdf", "pdf");
    return 0;
}
