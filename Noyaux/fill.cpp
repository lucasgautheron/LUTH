#include <sstream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

struct element
{
    int nucleus;
    int A,Z;
    double abundance;

    element() : nucleus(0), A(0), Z(0) { }
};

struct abundance_data
{
    int Y_idx, nb_idx, T_idx;
    double Y, nb, T;
    double np, nn;
    
    std::vector<element> elements;
};

struct index_entry
{
    int idx;
    double val;
};

int nucleus_to_AZ(int nucleus, int &A, int &Z)
{
    Z = nucleus%1000;
    A = nucleus/1000;
}

int read_index(const char *path, std::map<int,double> &index)
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
        index[count-2] = val;
    }
    return count-2;
}

int read_compo_data(const char *path, std::vector<abundance_data *> &abundances)
{
    std::map<int,double> temperatures, densities, fractions; 
    read_index("EOS.t", temperatures);
    read_index("EOS.nb", densities);
    read_index("EOS.yq", fractions);

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
        iss >> ab->T_idx;
        iss >> ab->nb_idx;
        iss >> ab->Y_idx;
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

        ab->T = temperatures[ab->T_idx];
        ab->Y = fractions[ab->Y_idx];
        ab->nb = densities[ab->nb_idx];
        abundances.push_back(ab);
        ++count;
    }
    return count;
}

int main(int argc, char *argv[])
{
    std::vector<abundance_data *> abundances;
    read_compo_data("EOS.compo", abundances);

    FILE *fp = fopen("out.res", "w+");
    int T_idx = 36;
    int nb_idx = 87;
    int Y_idx = /*18*/27;
    for(int i = 0; i < abundances.size(); ++i) if(abundances[i]->T_idx == T_idx && abundances[i]->nb_idx == nb_idx && abundances[i]->Y_idx == Y_idx)
    {
        abundance_data ab = *abundances[i];
        for(int j = 0; j < ab.elements.size(); ++j)
        fprintf(fp, "%d %d %e\n", ab.elements[j].A-ab.elements[j].Z, ab.elements[j].Z, ab.elements[j].abundance);
    }
    fclose(fp);

    int max_nb = 0;
    std::vector<double> densities, mpvs[2], avgs[2], stdevs[2];
    
    for(int i = 0; i < abundances.size(); ++i) if(abundances[i]->T_idx == T_idx && abundances[i]->Y_idx == Y_idx)
    {
        abundance_data ab = *abundances[i];
        
        int mpv[2] = { 0, 0 };
        double mpv_val[2] = { 0, 0 };

        double x[2] = { 0, 0 }, xx[2] = { 0, 0 };

        double p = 0;
        for(int j = 0; j < ab.elements.size(); ++j)
        {
            element elem = ab.elements[j];
            if(elem.A >= 20) {
                p += elem.abundance;
                for(int k = 0; k < 2; ++k) 
                {
                    double v = k == 0 ? elem.A - elem.Z : elem.Z;

                    if(mpv_val[k] < elem.abundance) {
                        mpv_val[k] = elem.abundance;
                        mpv[k] = v;
                    }
                    x[k] += v * elem.abundance;
                    printf("%e %e\n", v, elem.abundance);
                    xx[k] += v * v * elem.abundance;
                }
            }
        }
        printf("%d %d %e\n", i, (int)ab.elements.size(), p);

        if(p)
        {
            for(int k = 0; k < 2; ++k) { x[k] /= p; xx[k] /= p; mpvs[k].push_back(mpv[k]); avgs[k].push_back(x[k]); stdevs[k].push_back(sqrt(xx[k]-x[k]*x[k])); }
            densities.push_back(ab.nb);
        }
    }

    int N = densities.size();
    double *x[2], *y[2], *ex[2], *ey[2];

    for(int k = 0; k < 2; ++k)
    {
        x[k] = new double[N];
        y[k] = new double[N];
        ex[k] = new double[N];
        ey[k] = new double[N];

       for(int i = 0; i < N; ++i)
       {
           x[k][i] = densities[i];
           y[k][i] = avgs[k][i];
           ex[k][i] = 0;
           ey[k][i] = stdevs[k][i];
       }
    }

    TCanvas *c = new TCanvas("c","data",200, 10, 700, 500);
    c->SetLogx();
    TGraphErrors *gr[2];

    for(int k = 0; k < 2; ++k)
    {
        gr[k] = new TGraphErrors(N,x[k],y[k],ex[k],ey[k]);
        gr[k]->SetMinimum(10);
        gr[k]->SetMaximum(85);
        gr[k]->SetMarkerStyle(k == 0 ? 8 : 21);
        gr[k]->SetLineColor(k == 0 ? kBlue : kRed);
        gr[k]->SetMarkerColor(k == 0 ? kBlue : kRed);
        gr[k]->SetFillColor(k == 0 ? kBlue : kRed);
        gr[k]->SetTitle(";n_B (fm^{-3});N,Z");
    }
    
   
    gr[0]->Draw("AP");

    gr[1]->Draw("Psame");
    
    c->Print("out.pdf", "pdf");
    return 0;
}
