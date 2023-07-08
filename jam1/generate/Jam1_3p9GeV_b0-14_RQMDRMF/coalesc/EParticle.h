#ifndef EParticle_h
#define EParticle_h

#include <vector>

class EParticle
{
public:
    int kf;        // PDG particle code, kf=91 for cluster
    int ks;        // status code ks>10 for dead particle
    int ibar;      // baryon number
    double p[4];   // momemtum p[0]=energy
    double m;      // mass
    double r[4];   // coordinate r[0]=time
    EParticle() {
        kf=ks=ibar=0; m=0.0;
        for(int i=0;i<4;i++) p[i]=r[i]=0.0;
    }
    virtual ~EParticle() {}

    // Polluted by Shu HE
    int oid;
};

class ENucleus : public EParticle
{
public:
    int iz;      // number of proton
    int in;      // number of neutron
    int iy;      // number of Lambda
    int is[3];   // number of Sigma is[0]=Sigma-, is[1]=Sigma0, is[2]=Sigma+
    int ix[2];   // number of Xi    ix[0]=Xi-, ix[1]=Xi0
    int ig;      // number of Omega-
    double ex;   // total energy of the nuclear cluster - mass (GeV)
    int j;       // angular momentum
    std::vector<EParticle*> pc;   // hadrons in the cluster

    ENucleus(): EParticle() {
        pc.clear();
        iz=in=iy=is[0]=is[1]=is[2]=ix[0]=ix[1]=ig=0;
        ex=0.0; j=0;
    }
    ~ENucleus() {
        pc.clear();
    }
    void add(EParticle* p)  {pc.push_back(p);}
    std::vector<EParticle*> getPc() {return pc;}
    int  getSize()          {return pc.size();}
    EParticle* getPc(int i) {return pc[i];}
};


#endif
