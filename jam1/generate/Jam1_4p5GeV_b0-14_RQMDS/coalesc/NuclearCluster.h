#ifndef NuclearCluster_h
#define NuclearCluster_h

#include <vector>
#include "EParticle.h"

struct qmdParam
{
    double pmu1,pmu2,vex1,vex2,clam;
    double rho0,t1,t3,gamma,gam1,width,fac32;
    int  optMCut,optP0,optPot, optPotDelta;
    int  isMomDep;
};

class NuclearCluster
{
private:
    int mstc106;
    double R0_d,P0_d;
    double R0_t,P0_t;
    double R0_he4,P0_he4;
    double R0_hl3,P0_hl3;
    double R0_hl4,P0_hl4;
    int optPotDelta;
    qmdParam qmdp;
		EParticle* combineParticle(EParticle* p1, EParticle* p2, int n);
		
public:
    NuclearCluster(int m) {
        mstc106=m;
        R0_d=3.575;  // fm
        P0_d=0.285; // GeV/c
        R0_t=3.575;  // fm
        P0_t=0.285; // GeV/c
        R0_he4=3.575;  // fm
        P0_he4=0.285; // GeV/c
        R0_hl3=3.575;  // fm
        P0_hl3=0.285; // GeV/c
        R0_hl4=3.575;  // fm
        P0_hl4=0.285; // GeV/c
        optPotDelta=0;
        setParam(m);
    }
    void setR0_d(double r) {R0_d=r;}
    void setP0_d(double p) {P0_d=p;}
    void setRP0_d(double r, double p) {R0_d=r,P0_d=p;}
    void setR0_t(double r) {R0_t=r;}
    void setP0_t(double p) {P0_t=p;}
    void setRP0_t(double r, double p) {R0_t=r,P0_t=p;}
    void setR0_he4(double r) {R0_he4=r;}
    void setP0_he4(double p) {P0_he4=p;}
    void setRP0_he4(double r, double p) {R0_he4=r,P0_he4=p;}
    void setR0_hl3(double r) {R0_hl3=r;}
    void setP0_hl3(double p) {P0_hl3=p;}
    void setRP0_hl3(double r, double p) {R0_hl3=r,P0_hl3=p;}
    void setR0_hl4(double r) {R0_hl4=r;}
    void setP0_hl4(double p) {P0_hl4=p;}
    void setRP0_hl4(double r, double p) {R0_hl4=r,P0_hl4=p;}
		
    double getMass(int kf);
    bool clust(EParticle* p1, EParticle* p2, double r0, double p0);
    void findCluster(std::vector<EParticle*>& p, std::vector<ENucleus*>& enclp);
    void rpCMsq(EParticle &p1, EParticle &p2,double& rr,double &pp);
    void setParam(int mstc106);
    double getEnergy(std::vector<EParticle*>& p);
    void boost(double* pcm, double p[4]);
    void  setClusterProperties(ENucleus* p);

};
#endif
