#ifndef __SCANNING_TCHAIN_CHAIN_H
#define __SCANNING_TCHAIN_CHAIN_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"

#include "PDGData.h"

char scan_file_name[] = "filelist";
char scan_tree_name[] = "jam";

#ifdef MAXROOTFILES
const size_t MaxRootNum = (int) MAXROOTFILES;
#else
const size_t MaxRootNum = (int) 10000;
#endif

struct Scanner {
public:
    size_t fnum;
    size_t tentry;
    size_t f_head, f_idx;
    size_t iet, currentEntry;
    std::vector<std::string> flist;
    char fname[512], tname[512];
    TTree *ftree, *itree;
    TFile *rfile;

#ifdef B_MAXMUL
    static const size_t MaxMul = B_MAXMUL;
#else
    static const size_t MaxMul = 2000;
#endif
    float b;
    int mul, Npart;
    int pid[MaxMul];
    float px[MaxMul], py[MaxMul], pz[MaxMul];
    float E[MaxMul], mass[MaxMul], charge[MaxMul];

    int           ks[MaxMul];
    float         x[MaxMul];
    float         y[MaxMul];
    float         z[MaxMul];
    float         t[MaxMul];
    float         ft[MaxMul];

    Scanner(char *_fname = scan_file_name, char *_tname = scan_tree_name) {
        strcpy(fname, _fname);
        strcpy(tname, _tname);
        tentry = 0;
        f_idx = 0;
        ftree = NULL;
        itree = NULL;
        fnum = 0;
        iet = 0;
        currentEntry = 0;
        LoadList();
    }

    void LoadList() {
        ifstream fin(fname);
        std::string fentry;
        for (fnum = 0; fin>>fentry; ++fnum) {
            flist.push_back(fentry);
        }
    }

    void LoadTree() {
        //puts("Load Tree...");
        if (ftree != NULL) {
            delete ftree;
            ftree = NULL;
            gROOT->cd();
            rfile->Close();
            delete rfile;
            rfile = NULL;
        }
        rfile = TFile::Open(flist[f_idx].c_str());
        rfile->GetObject(tname, ftree);
        tentry = ftree->GetEntries();
        f_head = f_idx;

        ftree->SetBranchAddress("b", &b);
        ftree->SetBranchAddress("mul", &mul);
        ftree->SetBranchAddress("Npart", &Npart);
        ftree->SetBranchAddress("pid", pid);
        ftree->SetBranchAddress("px", px);
        ftree->SetBranchAddress("py", py);
        ftree->SetBranchAddress("pz", pz);
        //ftree->SetBranchAddress("E", E);
        //ftree->SetBranchAddress("mass", mass);
        //ftree->SetBranchAddress("charge", charge);
        ftree->SetBranchAddress("ks", ks);
        ftree->SetBranchAddress("x", x);
        ftree->SetBranchAddress("y", y);
        ftree->SetBranchAddress("z", z);
        ftree->SetBranchAddress("t", t);
        //ftree->SetBranchAddress("ft", ft);
    }

    bool hasNextEntry() {
        if (iet >= tentry) {
            if (f_idx >= fnum) return false;
            LoadTree();
            ++f_idx;
            iet = 0;
        }
        if (ftree->GetEntry(iet) == 0) return false;
        ++iet;
        ++currentEntry;

        {
            for (int ip = 0; ip < this->mul; ++ip) {
                int apid = abs(this->pid[ip]);
                const ParticleDataEntry* pdt = pdg[apid];

                if (pdt->pid == 0) continue;

                double m0 = pdt->m0;
                double px = this->px[ip];
                double py = this->py[ip];
                double pz = this->pz[ip];
                double p2 = px*px + py*py + pz*pz;

                this->E[ip] = sqrt(p2 + m0*m0);
                this->mass[ip] = m0;
            }
        }

        return true;
    }

private:
    PDGData pdg;
};
#endif
