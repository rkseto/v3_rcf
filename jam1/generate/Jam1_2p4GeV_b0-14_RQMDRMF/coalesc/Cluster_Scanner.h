////////////////////////////////////////////////////////
// 2020.11.26 modified by Xionghong He 
// save the freeze-out time  information and constituent
// nucleons momenta
////////////////////////////////////////////////////////



#ifndef __CLUSTER_SCANNING_TCHAIN_CHAIN_H
#define __CLUSTER_SCANNING_TCHAIN_CHAIN_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"

#include "inc/Scanner.h"
#include "inc/PDGData.h"

// Cluster Module
// Thanksgive to authors of JAM
#include "EParticle.h"
#include "NuclearCluster.h"
#define sign(a) ( ( (a) < 0 )  ?  -1   : ( (a) > 0 ) )

struct Cluster_Scanner {
#ifdef B_MAXMUL
    static const size_t MaxMul = B_MAXMUL;
    static const int cmaxv = B_MAXMUL/2;
#else
    static const size_t MaxMul = 2000;
    static const int cmaxv = 400;
#endif
    float b;
    int mul, Npart;
    int pid[MaxMul];
    float px[MaxMul], py[MaxMul], pz[MaxMul];
    float E[MaxMul], mass[MaxMul], charge[MaxMul];

    // Maybe later
    //int           ks[MaxMul];
    float         x[MaxMul];
    float         y[MaxMul];
    float         z[MaxMul];
    float         t[MaxMul];
    float         ft[MaxMul];

    // for clusters (for user also)
	int cmul, dmul, fmul; // mul for cluster, dead & free particle
    // the i-th cluster has npid[i] nucleons,
    // cpid[i][j] is the pid of the j-th nucleon in the i-th cluster
    int _cpid[cmaxv], *cpid[cmaxv], npid[cmaxv];
    float cpx[cmaxv], cpy[cmaxv], cpz[cmaxv], cE[cmaxv];
		float _cpx_p[cmaxv], _cpy_p[cmaxv], _cpz_p[cmaxv], _cmass_p[cmaxv];
		float _cx_p[cmaxv], _cy_p[cmaxv], _cz_p[cmaxv], _ct_p[cmaxv];
		float *cpx_p[cmaxv]; // momentum of the constituent nucleons
		float *cpy_p[cmaxv];
		float *cpz_p[cmaxv];
		float *cmass_p[cmaxv];
		float *cx_p[cmaxv]; // coordinates of the constituent nucleons
		float *cy_p[cmaxv];
		float *cz_p[cmaxv];
		float *ct_p[cmaxv];
    float cmass[cmaxv];
    // whether decay particles in cluster
    bool DepartCluster;
    // for cluster parameters
    int mstc106;
    float R0_d, P0_d; // deuteron coalescence parameters
    float R0_t, P0_t;	// triton and 3He coalescence parameters
    float R0_he4, P0_he4;	// 4He and more heavier coalescence parameters
    float R0_hl3, P0_hl3;	// hyper-triton coalescence parameters
    float R0_hl4, P0_hl4;	// hyper-h4 and more hlavier coalescence parameters

    Scanner &ev;

public:

    Cluster_Scanner(Scanner &_ev) : ev(_ev) {
        mstc106 = 14;
        R0_d = 3.575;
        P0_d = 0.285;
        R0_t = 3.575;
        P0_t = 0.285;
        R0_he4 = 3.575;
        P0_he4 = 0.285;
        R0_hl3 = 3.575;
        P0_hl3 = 0.285;
        R0_hl4 = 3.575;
        P0_hl4 = 0.285;
        DepartCluster = false;
    }

    bool hasNextEntry();
    void MakeCluster(Scanner &ev, std::vector<EParticle*> &part, std::vector<ENucleus*> &clust);

public:
    enum {
        ksDeuteron     = 0,
        ksTriton      = 1,
        ksHelium3      = 2,
        ksAlpha        = 3,
        ksAntiDeuteron = 4,
        ksAntiTriton  = 5,
        ksAntiHelium3  = 6,
        ksAntiAlpha    = 7,
        ksHyperTriton  = 11,
        ksHyperH4      = 12,
        ksUnknown      = 1024
    };

    int ClusterType(int idx);
};

bool Cluster_Scanner::hasNextEntry() {
    using namespace std;
    if (!ev.hasNextEntry()) return false;

    int i, oid, pr, clr, ncl, psize, msize, pidpos;

    vector<EParticle*> parts;
    vector<ENucleus*> clust; // Define vector of particles in an event;

    // Find Clusters
    parts.clear();
    clust.clear();
    MakeCluster(ev, parts, clust);

    // Record new event data;
    psize = parts.size();
    msize = clust.size();
    mul = psize;
    cmul = msize;
    fmul = mul;
    b = ev.b;
    Npart = ev.Npart;

    // Record particles
    for (pr = 0; pr < psize; ++pr) {
        oid = parts[pr]->oid;
        pid[pr] = ev.pid[oid];
        px[pr] = ev.px[oid];
        py[pr] = ev.py[oid];
        pz[pr] = ev.pz[oid];
        E[pr] = ev.E[oid];
        x[pr] = ev.x[oid];
        y[pr] = ev.y[oid];
        z[pr] = ev.z[oid];
        t[pr] = ev.t[oid];
        //ft[pr] = ev.ft[oid];
        mass[pr] = ev.mass[oid];
        delete parts[pr];
    }

    // Record clusters
    pidpos = 0;
    for (pr = 0; pr < msize; ++pr) {
        ncl = clust[pr]->pc.size();
        // Event.cpids[pr][clr] is the clr th particles in the pr th cluster
        // in EventScan.h of this data:
        // Event.cpids[pr] = &Event._cpid[Event.cpid[pr]];
        // where 0 <= clr < npid[pr]
        cpid[pr]  = &_cpid[pidpos]; // Where the pid series start for pr th cluster
				cpx_p[pr] = &_cpx_p[pidpos];
				cpy_p[pr] = &_cpy_p[pidpos];
				cpz_p[pr] = &_cpz_p[pidpos];
				cmass_p[pr] = &_cmass_p[pidpos];
				cx_p[pr] = &_cx_p[pidpos];
				cy_p[pr] = &_cy_p[pidpos];
				cz_p[pr] = &_cz_p[pidpos];
				ct_p[pr] = &_ct_p[pidpos];
        npid[pr] = ncl; // How many particles in pr th cluster
        for (clr = 0; clr < ncl; ++clr) {
            _cpid[pidpos]    = clust[pr]->pc[clr]->kf;
            _cpx_p[pidpos]   = clust[pr]->pc[clr]->p[1];
            _cpy_p[pidpos]   = clust[pr]->pc[clr]->p[2];
            _cpz_p[pidpos]   = clust[pr]->pc[clr]->p[3];
            _cmass_p[pidpos] = clust[pr]->pc[clr]->m;
            _cx_p[pidpos]    = clust[pr]->pc[clr]->r[1];
            _cy_p[pidpos]    = clust[pr]->pc[clr]->r[2];
            _cz_p[pidpos]    = clust[pr]->pc[clr]->r[3];
            _ct_p[pidpos]    = ev.t[clust[pr]->pc[clr]->oid];
            delete clust[pr]->pc[clr];
            ++pidpos;
        }
        cpx[pr] = clust[pr]->p[1];
        cpy[pr] = clust[pr]->p[2];
        cpz[pr] = clust[pr]->p[3];
        cmass[pr] = clust[pr]->m;
        //cE[pr] = clust[pr]->ex + cmass[pr]; //
        cE[pr] = clust[pr]->p[0]; //
        delete clust[pr];
    }

    dmul = pidpos - 1;

    int apid;
    double m0, m_clust, mfr;
    if (DepartCluster) { // decay the cluster
        for (pr = 0; pr < msize; ++pr) {
            m_clust = cmass[pr];
            for (i = 0; i < npid[pr]; ++i) {
                pid[mul] = cpid[pr][i];
                apid = abs(pid[mul]);
                switch (apid) {
                    case 2112: m0 = 0.93957; break;
                    case 2212: m0 = 0.93827; break;
                    case 3122: m0 = 1.11568; break; // lambda
                    case 3112: m0 = 1.19744; break; // Sigma-
                    case 3212: m0 = 1.19255; break; // Sigma0
                    case 3222: m0 = 1.18937; break; // Sigma+
                    case 3312: m0 = 1.32130; break; // Xi-
                    case 3322: m0 = 1.31490; break; // Xi0
                    case 3334: m0 = 1.67245; break; // Omega-
                    default:
                        cout << "NuclearCluster::getMass kf=?" << pid[mul] <<endl;
                        exit(1);
                }
                mfr = (m0 / m_clust);
                px[mul] = mfr * cpx[pr];
                py[mul] = mfr * cpy[pr];
                pz[mul] = mfr * cpz[pr];
                ++mul;
            }
        }
    }

    return true;
}

void Cluster_Scanner::MakeCluster(Scanner &ev, std::vector<EParticle*> &part, std::vector<ENucleus*> &clust) {
    int i, mul, pr, kf, n;
    //int kfa;
    int mstc106=this->mstc106;
    mul = ev.mul;

    using namespace std;
    vector<EParticle*> pbary, pmeson, precyle;
    vector<EParticle*>::iterator iter;
    EParticle *p;
    NuclearCluster* cluster=NULL;

    // From main() in jamcluster.cxx
    // TODO parameter
    // set parameer for nuclear potentials.
    cluster = new NuclearCluster(mstc106);
    // set parameter for coalecsce.
    float R0_d = this->R0_d, P0_d = this->P0_d;
    float R0_t = this->R0_t, P0_t = this->P0_t;
    float R0_he4 = this->R0_he4, P0_he4 = this->P0_he4;
    float R0_hl3 = this->R0_hl3, P0_hl3 = this->P0_hl3;
    float R0_hl4 = this->R0_hl4, P0_hl4 = this->P0_hl4;
    cluster->setRP0_d(R0_d, P0_d);
    cluster->setRP0_t(R0_t, P0_t);
    cluster->setRP0_he4(R0_he4, P0_he4);
    cluster->setRP0_hl3(R0_hl3, P0_hl3);
    cluster->setRP0_hl4(R0_hl4, P0_hl4);

    // from getParticles() from jamcluster.cxx
    bool skipBaryon;
    for (pr = 0; pr < mul; ++pr) {
        kf = ev.pid[pr];
        int kfa=abs(kf);
        p = new EParticle();
        p->ibar=sign(kf);
        if ((kfa%10000)/1000 == 0) p->ibar=0; // * A method to detact if baryon
        skipBaryon = false;
        if (p->ibar != 0) {
            switch (kfa) {
                case 2112: break;
                case 2212: break;
                case 3122: break; // lambda
                case 3112: break; // Sigma-
                case 3212: break; // Sigma0
                case 3222: break; // Sigma+
                case 3312: break; // Xi-
                case 3322: break; // Xi0
                case 3334: break; // Omega-
                default:
                    cout << "NuclearCluster::getMass kf=?" << kf<<endl;
                    skipBaryon = true;
            }
        }

        p->ks = ev.ks[pr];  // for JAM
        //p->ks = 1;   // for UrQMD
        p->kf = kf;
        p->m = ev.mass[pr];
        p->p[1] = ev.px[pr];
        p->p[2] = ev.py[pr];
        p->p[3] = ev.pz[pr];
        p->p[0] = ev.E[pr];
        p->r[1] = ev.x[pr];
        p->r[2] = ev.y[pr];
        p->r[3] = ev.z[pr];
        p->r[0] = ev.t[pr];
        p->oid = pr;

        if (skipBaryon) {// This particle is not interesting
            precyle.push_back(p);
            continue;
        }

        part.push_back(p);
    }

    // from nuclCluster of jamcluster.cxx
    pbary.clear();
    pmeson.clear();
    n = part.size();
    for(i = 0; i < n; i++) {
        if(part[i]->ibar !=0) {
            pbary.push_back(part[i]);
        } else {
            pmeson.push_back(part[i]);
        }
    }

    // findCluster has been redefined
    cluster->findCluster(pbary, clust); // Clusters will be inserted at the end of pbary
    part.clear();
    part.insert(part.end(), precyle.begin(), precyle.end()); // resume abandoned baryons at beginning
    for (iter = pbary.begin(); iter != pbary.end(); ++iter) // insert baryons
        if ((*iter)->ks != 30) { // The particle is not dead
            part.push_back(*iter);
        }

    // parts is the external vector
    part.insert(part.end(),pmeson.begin(),pmeson.end());

    delete cluster;
}

int Cluster_Scanner::ClusterType(int idx) {
    int csize = this->npid[idx]; // the size of the cluster
    int nz, nn, naz, nan, nl;
    nz = nn = naz = nan = nl = 0;

    for (int k = 0; k < csize; ++k) { // loop for all the nucleons in the cluster
        int pid = this->cpid[idx][k];
        switch (pid) {
            case 2212: // proton
                ++nz;
                break;
            case 2112: // neutron
                ++nn;
                break;
            case -2212: // anti-proton
                ++naz;
                break;
            case -2112: // anti-neutron
                ++nan;
                break;
            case 3122: // lambda
                ++nl;
                break;
        }
    } // end of loop for nucleons

    if (csize == 2 && nz == 1 && nn == 1) { // is deuteron
        return ksDeuteron;
    } else if (csize == 3 && nz == 1 && nn == 2) { // is triton
        return ksTriton;
    } else if (csize == 3 && nz == 2 && nn == 1) { // is helium-3
        return ksHelium3;
    } else if (csize == 4 && nz == 2 && nn == 2) { // is alpha
        return ksAlpha;
    } else if (csize == 2 && naz == 1 && nan == 1) { // anti-deuteron
        return ksAntiDeuteron;
    } else if (csize == 3 && naz == 1 && nan == 2) { // anti-triton
        return ksAntiTriton;
    } else if (csize == 3 && naz == 2 && nan == 1) { // anti-helium-3
        return ksAntiHelium3;
    } else if (csize == 4 && naz == 2 && nan == 2) { // anti-alpha
        return ksAntiAlpha;
    } else if (csize == 3 && nz == 1 && nn == 1 && nl == 1) { // is hyper-triton
        return ksHyperTriton;
    } else if (csize == 4 && nz == 1 && nn == 2 && nl == 1) { // is hyper-hydregon
        return ksHyperH4;
    } else { // is other
        return ksUnknown;
    }
} // int ClusterType

#endif
