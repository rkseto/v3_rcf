#include "NuclearCluster.h"
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

#define sign(a) ( ( (a) < 0 )  ?  -1   : ( (a) > 0 ) )

inline double pow2(const double& x) {return x*x;}

const double hc=0.19732705;

double NuclearCluster::getMass(int kf)
{
	switch (abs(kf)) {
		case 2112: return 0.93957;
		case 2212: return 0.93827;
		case 3122: return 1.11568; // lambda
		case 3112: return 1.19744; // Sigma-
		case 3212: return 1.19255; // Sigma0
		case 3222: return 1.18937; // Sigma+
		case 3312: return 1.32130; // Xi-
		case 3322: return 1.31490; // Xi0
		case 3334: return 1.67245; // Omega-
		default:
							 cout << "NuclearCluster::getMass kf=?" << kf<<endl;
							 exit(1);
	}
}

// compute relative distance and momentum.
void NuclearCluster::rpCMsq(EParticle& p1, EParticle& p2,double& rr,double &pp)
{
	double pcm[4],dp[4],dx[4];
	for(int j=0;j<4;j++) {
		dx[j] = p1.r[j] - p2.r[j];
		dp[j] = p1.p[j] - p2.p[j];
		pcm[j]= p1.p[j] + p2.p[j];
	}
	double s=pcm[0]*pcm[0] - pcm[1]*pcm[1] - pcm[2]*pcm[2] - pcm[3]*pcm[3];
	double rsq=dx[1]*dx[1] + dx[2]*dx[2] + dx[3]*dx[3];
	double psq=dp[0]*dp[0] - dp[1]*dp[1] - dp[2]*dp[2] - dp[3]*dp[3];
	double rp=dx[1]*pcm[1] + dx[2]*pcm[2] + dx[3]*pcm[3];
	//double dm = p1.m*p1.m - p2.m*p2.m;
	double p1sq = pow2(p1.p[0])-pow2(p1.p[1])-pow2(p1.p[2])-pow2(p1.p[3]);
	double p2sq = pow2(p2.p[0])-pow2(p2.p[1])-pow2(p2.p[2])-pow2(p2.p[3]);
	double dm = p1sq - p2sq;

	rr =  rsq + rp*rp/s;
	pp = -psq + dm*dm/s;

}

void NuclearCluster::setParam(int mstc106)
{
	int mstd101=1;
	double pard101,pard102,pard103,pard104,pard105,pard106,pard107,pard108;

	if(mstc106 == 4) {
		pard101=-124.e-3; // Aich Hard
		pard102=  70.5e-3;
		pard103= 2.0;
		pard104=1.08;
		pard105=1.0;  // pmu1 dummy
		pard106=1.0;  // pmu2 dummy
		pard107= 0.0;   // vex1
		pard108= 0.0;   // vex2
		mstd101=0;  // momentum dependent potential is not used.
	} else if (mstc106 == 5) {
		pard101=-356e-3; // Aich Soft
		pard102= 303e-3;
		pard103=   1.16667;
		pard104=   1.08;
		pard105=1.0;  // pmu1 dummy
		pard106=1.0;  // pmu2 dummy
		pard107= 0.0;   // vex1
		pard108= 0.0;   // vex2
		mstd101=0; // momentum dependent potential is not used.

		//.....2015 A.Ohnishi new fit1 mstc(106)=11-14
	} else if (mstc106 == 11) {   // MEH1 K=436.5
		pard101= 0.00960346793263952; // alpha
		pard102= 0.0655471283733455;  // beta
		pard103= 2.0;               // gamma
		pard104= 2.05;              // The width of Gaussian L (fm^2)
		pard105= 2.02*hc;             // mu1
		pard106= 1.0*hc;              // mu2
		pard107= -0.38314;           // C1
		pard108=  0.33741;            // C2
	} else if (mstc106 == 12)  {     // MH1 K=370.919367532666
		pard101= -0.0122455748584756; // alpha
		pard102=  0.0873961711644606; // beta
		pard103= 5.0/3.0;             // gamma
		pard104= 2.05;              // The width of Gaussian L (fm^2)
		pard105= 2.02*hc;            // mu1
		pard106= 1.0*hc;             // mu2
		pard107= -0.38314;           // C1
		pard108=  0.33741;           // C2
	} else if (mstc106 == 13) {      // MM1 K= 305.37223915932
		pard101=-0.0777927032318212; // alpha
		pard102= 0.152943299537806;  // beta
		pard103= 4.0/3.0;            // gamma
		pard104= 2.05;               // The width of Gaussian L (fm^2)
		pard105= 2.1*hc;             // mu1
		pard106= 1.0*hc;             // mu2
		pard107= -0.38314;           // C1
		pard108=  0.33741;           // C2
	} else if (mstc106 == 14) {        // MS1 K=272.598674972647
		pard101= -0.208886959978512; // alpha
		pard102=  0.284037556284497; // beta
		pard103= 7.0/6.0;            // gamma
		pard104= 2.1;                // The width of Gaussian L (fm^2)
		pard105= 2.02*hc;            // mu1
		pard106= 1.0*hc;             // mu2
		pard107= -0.38314;           // C1
		pard108=  0.33741;           // C2
	} else {
		cout << " not implemented mstc106= " << mstc106 << endl;
		exit(1);
	}

	double pmu1=pard105;
	double pmu2=pard106;
	double vex1=pard107;
	double vex2=pard108;
	double clam=2.0;

	double rho0= 0.168;
	double t1=pard101/2.0/rho0;
	double t3=pard102/(pard103+1.0)/(pow(rho0,pard103));
	double width=pard104;
	double fac32=pow(4.0*M_PI*width,1.5);   //    [(4*pi*L)^3/2]
	double gamma = pard103;

	int mstc107=0;
	int mstc104=12;
	int mstc109=2;

	qmdp.pmu1=pmu1;
	qmdp.pmu2=pmu2;
	qmdp.vex1=vex1;
	qmdp.vex2=vex2;
	qmdp.clam=clam;
	qmdp.rho0=rho0;
	qmdp.t1=t1;
	qmdp.t3=t3;
	qmdp.width=width;
	qmdp.fac32=fac32;
	qmdp.gamma=gamma;
	qmdp.gam1=gamma-1.0;
	qmdp.optMCut=mstc107;
	qmdp.optP0=mstc109;
	qmdp.optPot=mstc104;
	qmdp.optPotDelta=optPotDelta;
	qmdp.isMomDep=mstd101;

}

double NuclearCluster::getEnergy(std::vector<EParticle*>& p)
{
	//...Purpose: to calculate single particle energy in RQMD/S
	// H_i =  1/(2*E_i) *[E_i^2 -vec{p}_i^2 - m_i^2 - 2m_i*V_i ]
	// here V_i = t1 * <rho_i> + t3 *<rho_i>^gamma : Skyrme Potential

	int n=p.size();
	double *rho = new double [n];
	std::vector< std::vector<double> > rhom, p22;
	rhom.resize(n);
	p22.resize(n);
	for(int i=0;i<n;i++) {
		rhom[i].resize(n);
		p22[i].resize(n);
	}

	for(int i=0;i<n;i++) {
		rho[i]=0.0;
		if(p[i]->ibar==0) continue;
		for(int j=i+1; j<n;j++) {
			if(p[j]->ibar==0) continue;
			double r22, p2;
			rpCMsq(*p[i],*p[j],r22,p2);
			rhom[i][j] = exp(-r22/(4.0*qmdp.width))/qmdp.fac32;
			rhom[j][i] = rhom[i][j];
			rho[i] += rhom[i][j];
			rho[j] += rhom[i][j];
			p22[i][j] = p2;
			p22[j][i] = p2;
		}
	}

	double etot = 0.0, etot2=0.0;
	for(int i=0;i<(int)p.size();i++) {
		double afac=1.0, bfac=1.0;
		if(qmdp.optPotDelta==1 && (p[i]->kf != 2212 || p[i]->kf != 2112)) {
			//if(optPotDelta==1 && (p[j].id == id_delt || p[j].id == id_delts)) {
			//afac =0.88;
			bfac =1.09;
		}

		// Skyrme type density dependent force.
		double vpot = afac*qmdp.t1*rho[i] + bfac*qmdp.t3*pow(rho[i],qmdp.gamma);

		// momentum dependent force.
		if(qmdp.isMomDep == 1) {

			double pmu1=qmdp.pmu1;
			double pmu2=qmdp.pmu2;
			double vex1=qmdp.vex1;
			double vex2=qmdp.vex2;
			double rho0=qmdp.rho0;

			for(int j=0;j<n;j++) {
				if(i==j) continue;
				double fac1 = 1.0 + p22[i][j]/(pmu1*pmu1);
				double fac2 = 1.0 + p22[i][j]/(pmu2*pmu2);
				double pmom2 = vex1/fac1 + vex2/fac2;
				vpot +=  pmom2/2.0/rho0*rhom[i][j];
			}
		}

		double pp = pow2(p[i]->p[1])+pow2(p[i]->p[2])+pow2(p[i]->p[3]);
		double m=sqrt(vpot*vpot + p[i]->m*p[i]->m);
		//double meff2=max(0.0,m*m* 2*m*vpot); 
		//etot += sqrt(meff2+pp);

		double em=getMass(p[i]->kf);
		etot += sqrt(em*em+pp + 2*em*vpot)-em;
		etot2 += sqrt(p[i]->m*p[i]->m +pp)-em;
		}

		/*
			 cout << " n= " << n << " etot= " << etot*1000
			 << " etot2= " << etot2*1000
			 <<endl;
			 cin.get();
			 */

		delete [] rho;
		rhom.clear();
		p22.clear();

		return etot;
		//return etot2;
	}

	void NuclearCluster::boost(double* pcm,double p[4])
	{
		// pcm[4]=srt, pcm[0]=energy
		double pcs=pcm[1]*p[1]+pcm[2]*p[2]+pcm[3]*p[3];
		double transf=(pcs/(pcm[4]+pcm[0])+p[0])/pcm[4];
		p[1] += pcm[1]*transf;
		p[2] += pcm[2]*transf;
		p[3] += pcm[3]*transf;
		p[0] += (pcm[4]*p[0] + pcs)/pcm[4];
	}

	EParticle* NuclearCluster::combineParticle(EParticle* p1, EParticle* p2, int n)
	{
		EParticle* mom = new EParticle();
		double pcm[5];
		pcm[0]=pcm[1]=pcm[2]=pcm[3]=0.0;
		double lcm[5];
		lcm[0]=lcm[1]=lcm[2]=lcm[3]=0.0;
		for(int i=0;i<4;i++) {
			pcm[i] = p1->p[i]+p2->p[i];
			if(i==0) lcm[i] = max(p1->r[0],p2->r[0]);
			else lcm[i] = p1->r[i]*n + p2->r[i];
		}

		mom->r[1] = lcm[1]/(n+1);
		mom->r[2] = lcm[2]/(n+1);
		mom->r[3] = lcm[3]/(n+1);
		mom->r[0] = lcm[0];
		mom->p[1] = pcm[1];
		mom->p[2] = pcm[2];
		mom->p[3] = pcm[3];
		mom->p[0] = pcm[0];
		mom->m = sqrt(pcm[0]*pcm[0] - pcm[1]*pcm[1] - pcm[2]*pcm[2] - pcm[3]*pcm[3]);

		mom->ibar = n+1;

		return mom;
	}

	void NuclearCluster::setClusterProperties(ENucleus* pc)
	{
		std::vector<EParticle*> p0 = pc->getPc();
		std::vector<EParticle*> p;
		p.clear();

		int n=p0.size();

		for(int i=0;i<n;i++) {
			EParticle* p1 = new EParticle();
			*p1=*p0[i];
			p.push_back(p1);
		}

		double pcm[5];
		pcm[0]=pcm[1]=pcm[2]=pcm[3]=0.0;
		double x=0.0, y=0.0, z=0.0, t=0.0;
		int iz=0, in=0,iy=0,ism=0, is0=0,isp=0,ixm=0,ix0=0,ig=0;
		int ib=0;
		for(int i=0;i<n;i++) {
			p[i]->ks=30;
			pcm[1] += p[i]->p[1];
			pcm[2] += p[i]->p[2];
			pcm[3] += p[i]->p[3];
			pcm[0] += p[i]->p[0];
			x += p[i]->r[1];
			y += p[i]->r[2];
			z += p[i]->r[3];
			t = max(t,p[i]->r[0]);
			ib += p[i]->ibar;
			int kfa=abs(p[i]->kf);
			if(kfa==2212) iz++;
			else if(kfa==2112) in++;
			else if(kfa==3122) iy++;
			else if(kfa==3112) ism++;
			else if(kfa==3212) is0++;
			else if(kfa==3222) isp++;
			else if(kfa==3312) ixm++;
			else if(kfa==3322) ix0++;
			else if(kfa==3334) ig++;
			else {
				cout << "setClusterProperties:: kfa? " << kfa <<endl;
			}
		}

		int kf0=p[0]->kf;
		pc->ibar=ib;
		pc->ks=2;
		pc->kf=91*sign(kf0);
		pc->iz=iz*sign(kf0);
		pc->in=in;
		pc->iy=iy*sign(kf0);
		pc->is[0]=ism*sign(kf0);
		pc->is[1]=is0*sign(kf0);
		pc->is[2]=isp*sign(kf0);
		pc->ix[0]=ixm*sign(kf0);
		pc->ix[1]=ix0*sign(kf0);
		pc->ig=ig*sign(kf0);
		pc->r[1] = x/n;
		pc->r[2] = y/n;
		pc->r[3] = z/n;
		pc->r[0] = t;
		pc->p[1] = pcm[1];
		pc->p[2] = pcm[2];
		pc->p[3] = pcm[3];
		pc->p[0] = pcm[0];
		pc->m = sqrt(pcm[0]*pcm[0] - pcm[1]*pcm[1] - pcm[2]*pcm[2] - pcm[3]*pcm[3]);

		pcm[1]=-pc->p[1];
		pcm[2]=-pc->p[2];
		pcm[3]=-pc->p[3];
		pcm[0]= pc->p[0];
		pcm[4]= pc->m;

		double pc0[4],cmr[4];
		pc0[1]=pc0[2]=pc0[3]=pc0[0]=0.0;
		cmr[1]=cmr[2]=cmr[3]=cmr[0]=0.0;

		for(int i=0;i< (int)p.size();i++) {

			//p[i]->r[0]=0.0;
			double pc[4],rc[4];
			for(int j=0; j<4;j++) {
				pc[j]=p[i]->p[j];
				rc[j]=p[i]->r[j];
			}
			double mfsq=pc[0]*pc[0]-pc[1]*pc[1]-pc[2]*pc[2]-pc[3]*pc[3];
			boost(pcm,pc);
			boost(pcm,rc);

			p[i]->p[1]=pc[1];
			p[i]->p[2]=pc[2];
			p[i]->p[3]=pc[3];
			p[i]->p[0]=sqrt(mfsq+pow2(pc[1])+pow2(pc[2])+pow2(pc[3]));

			p[i]->r[1]=rc[1];
			p[i]->r[2]=rc[2];
			p[i]->r[3]=rc[3];
			p[i]->r[0]=rc[0];

			pc0[1] += pc[1];
			pc0[2] += pc[2];
			pc0[3] += pc[3];
			pc0[0] += pc[0];

			for(int j=1;j<4;j++) cmr[j] += rc[j]*pc[0];
			cmr[0] += pc[0];

		}

		for(int j=1;j<4;j++) cmr[j] /= cmr[0];

		/*
			 cout << " pcmx= " << pc0[1]
			 << " pcmy= " << pc0[2]
			 << " pcmz= " << pc0[3]
			 << endl;
			 */


		for(int i=0;i< (int)p.size();i++) {
			for(int j=1;j<4;j++) p[i]->r[j] -= cmr[j];
		}


		// Compute angular momentum.
		double lx=0.0, ly=0.0, lz=0.0;
		for(int i=0;i< (int)p.size();i++) {
			lx += p[i]->r[2] * p[i]->p[3] - p[i]->r[3]*p[i]->p[2];
			ly += p[i]->r[3] * p[i]->p[1] - p[i]->r[1]*p[i]->p[3];
			lz += p[i]->r[1] * p[i]->p[2] - p[i]->r[2]*p[i]->p[1];

			//cout << " after getEnergy " << i << " px = " << p[i]->p[1] <<endl;
		}

		double l=static_cast<int>( sqrt(lx*lx + ly*ly + lz*lz )/hc + 0.5);
		pc->j=static_cast<int>(l);
		//pc.jz=int(lz/hc);

		// Compute the energy (minus mass) in GeV of the cluster. 
		// Note that this is not excitation energy
		// of the cluster, since binding energy is not subtracted. 
		pc->ex = getEnergy(p);

		for(int i=0;i<n;i++) {
			delete p[i];
		}
	}

	bool NuclearCluster::clust(EParticle* p1, EParticle* p2, double R0_nuclei, double P0_nuclei)
	{
		if(p1->ks >10) return false;
		if(p1->ibar==0) return false;

		if(p2->ks >10) return false;
		if(p2->ibar==0) return false;

		if(p1->kf*p2->kf<0) return false; // exclude B + antiB

		double rr, pp;
		rpCMsq(*p1,*p2,rr,pp);

		if(rr<= R0_nuclei*R0_nuclei && pp <= P0_nuclei*P0_nuclei) return true;

		return false;

	}

	// Shu HE dirted this work here, I redefined the parameters of this method
	void NuclearCluster::findCluster(std::vector<EParticle*>& p, std::vector<ENucleus*> &_clust)
	{
		int nv=p.size();

		double px0=0.0, py0=0.0, pz0=0.0, pe0=0.0;
		int ibar0=0;
		for(int i=0;i<nv;i++) {
			if(p[i]->ks>10) continue; // skip dead particles
			px0 += p[i]->p[1];
			py0 += p[i]->p[2];
			pz0 += p[i]->p[3];
			pe0 += p[i]->p[0];
			ibar0 += p[i]->ibar;
		} // calc the total momentum and baryon number conservation

		int *mscl =new int [nv] ; // count for particles in nucl[nv]
		int *num =new int [nv] ; // pointers to particle

		for(int i=0;i<nv;i++) {
			mscl[i]=1;
			num[i]=i; // set num[i] point the to particles[i] intially
		}

		// Algorithms (Depth-First Search) Find the connected components of a graph
		int  nclst=0; // count for nucleus formed
		int  icheck=0; // particles below it has been in a certain cluster
		int  nconstituent=1; // consitueont particles in the clust
		int  nused=0; // total particles used in clusters

		for(int i=0;i<nv;) { // for all particles in this searching
			icheck = i;

			int nProton(0), nNeutron(0);
			
			int i1=num[i];
			if(p[i1]->kf ==2212) nProton =1; 
			else if(p[i1]->kf == 2112)  nNeutron =1;
			else {i=i+1; continue;}

			EParticle *inP = p[i1]; // the intermediate particle in the coalescence

			for(int j=i+1;j<nv;j++) {
				int i2=num[j];
				
				if(p[i2]->kf != 2212 && p[i2]->kf != 2112) continue;
					
				if(clust(inP, p[i2], R0_he4, P0_he4)) { // if the particle i1, i2 are within the coalescence distance of 4He
					if(p[i2]->kf == 2212) nProton += 1; 
					else if(p[i2]->kf == 2112)  nNeutron += 1;
					
					int lp=num[icheck+1]; // swap particle num[icheck+1] with particle num[j]
					num[icheck+1]=i2; // thus the particle in a cluster adjoint
					num[j]=lp;
					
					inP = combineParticle(inP, p[i2], nconstituent);
					nconstituent++;
					icheck++; // push the particle i2 down below the limit
					if(nconstituent >= 4) break;
				}
			}
			
			// if all particles in a cluster been found
			if(nconstituent == 4 && nProton == 2 && nNeutron == 2) {
				for(int k=0; k<nconstituent; k++) {
						
					int lp1 = num[nused];
					num[nused] = num[icheck+1-nconstituent+k];
					num[icheck+1-nconstituent+k] = lp1;
					nused++;
				}
				mscl[nclst] = nconstituent; // particles in this cluster
				nclst++;
			} 
			
			icheck = icheck + 1;
		  i = icheck;
			nconstituent = 1;
		}

		if(nused<=icheck) {
			icheck = nused;
			int nused_he4 = nused;
			
			for(int i=nused_he4; i<nv;) { // for all particles in this searching

				int nProton(0), nNeutron(0);
				
				int i1=num[i];
				if(p[i1]->kf == 2212) nProton =1; 
				else if(p[i1]->kf == 2112)  nNeutron =1;
			  else {i=i+1; continue;}

				EParticle *inP = p[i1]; // the intermediate particle in the coalescence

				for(int j=i+1;j<nv;j++) {
					int i2=num[j];
				  if(p[i2]->kf != 2212 && p[i2]->kf != 2112) continue;
					
					if(clust(inP, p[i2], R0_t, P0_t)) { // if the particle i1, i2 are within the coalescence distance of 4He
						if(p[i2]->kf == 2212) nProton += 1; 
						else if(p[i2]->kf == 2112)  nNeutron += 1;
						
						int lp=num[icheck+1]; // swap particle num[icheck+1] with particle num[j]
						num[icheck+1]=i2; // thus the particle in a cluster adjoint
						num[j]=lp;
						
						inP = combineParticle(inP, p[i2], nconstituent);
						nconstituent++;
						icheck++; // push the particle i2 down below the limit
					  if(nconstituent >= 3) break;
					}
				}
				
				// if all particles in a cluster been found
				if(nconstituent == 3 && ((nProton==1 && nNeutron ==2) || (nProton==2 && nNeutron ==1)) ) {
					
				  for(int h=icheck+1; h<nv; h++) {  // coalescence of triton and a lambda
						
						int i3=num[h];
						if(p[i3]->kf != 3122) continue;
						
						if(clust(inP, p[i3], R0_hl4, P0_hl4)) { // if the particle i1, i2 are within the coalescence distance of H4L
							
							int lp=num[icheck+1]; 
							num[icheck+1]=i3; 
							num[h]=lp;
							
							nconstituent++;
							icheck++; // push the particle i2 down below the limit
							break;
						}
					}
				
					for(int k=0; k<nconstituent; k++) {
							
						int lp1 = num[nused];
						num[nused] = num[icheck-nconstituent+k+1];
						num[icheck-nconstituent+k+1] = lp1;
						nused++;
					}

					mscl[nclst] = nconstituent; // particles in this cluster
					nclst++;
				} 
			
				icheck = icheck + 1;
				i = icheck;
				nconstituent = 1;
			}
		}

		if(nused<=icheck) {
			icheck = nused;
			int nused_he3 = nused;
			
			for(int i=nused_he3;i<nv;) { // for all particles in this searching
				
				int nProton(0), nNeutron(0);
				
				int i1=num[i];
				if(p[i1]->kf ==2212) nProton = 1; 
				else if(p[i1]->kf == 2112)  nNeutron = 1;
				else {i = i+1; continue;}
				
				EParticle *inP = p[i1]; // the intermediate particle in the coalescence

				for(int j=i+1;j<nv;j++) {
					int i2=num[j];
				  if(p[i2]->kf != 2212 && p[i2]->kf != 2112) continue;
					
					if(clust(inP, p[i2], R0_d, P0_d)) { // if the particle i1, i2 are within the coalescence distance of 4He
						
						if(p[i2]->kf ==2212) nProton += 1; 
						else if(p[i2]->kf == 2112)  nNeutron += 1;
					
						int lp=num[icheck+1]; // swap particle num[icheck+1] with particle num[j]
						num[icheck+1]=i2; // thus the particle in a cluster adjoint
						num[j]=lp;
						
						inP = combineParticle(inP, p[i2], nconstituent);
						nconstituent++;
						icheck++; // push the particle i2 down below the limit
						if(nconstituent == 2) break;
					}
				}
				
				// if all particles in a deuteron been found
				if(nconstituent == 2 && nProton==1 && nNeutron ==1) {
					
					for(int h=icheck+1; h<nv; h++) {  // coalescence of deuteron and a lambda
						int i3=num[h];
						if(p[i3]->kf != 3122) continue;
						
						if(clust(inP, p[i3], R0_hl3, P0_hl3)) { // if the particle i1, i2 are within the coalescence distance of 4He
							
							int lp=num[icheck+1]; 
							num[icheck+1]=i3; 
							num[h]=lp;
							
							nconstituent++;
							icheck++; // push the particle i2 down below the limit
							break;
						}
					}
					
					for(int k=0; k<nconstituent; k++) {
						
						int lp1 = num[nused];
						num[nused] = num[icheck+1-nconstituent+k];
						num[icheck+1-nconstituent+k] = lp1;
						nused++;
					}

					mscl[nclst] = nconstituent; // particles in this cluster
					nclst++;
				} 
			
				icheck = icheck + 1;
				i = icheck;
				nconstituent = 1;
			}
		}


    std::vector<ENucleus*> pclust;
    pclust.clear();

    int ibar=0;
    int miss=0;
    int nn=0;
    int ii=0;
    for(int ic=0;ic<nclst;ic++) { // for all formed nucleus
        int nt=mscl[ic];
        nn += nt;
        int isave=ii;
        int iz=0, in=0,iy=0,ism=0, is0=0,isp=0,ixm=0,ix0=0,ig=0;
        for(int i=0;i<nt;i++) { // for all particles in the nucleus
            int j=num[ii];
            int kfa=abs(p[j]->kf);
            if(kfa==2212) iz++;
            else if(kfa==2112) in++;
            else if(kfa==3122) iy++;
            else if(kfa==3112) ism++;
            else if(kfa==3212) is0++;
            else if(kfa==3222) isp++;
            else if(kfa==3312) ixm++;
            else if(kfa==3322) ix0++;
            else if(kfa==3334) ig++;
            if(p[j]->ks<10) ibar += p[j]->ibar;
            ii++;
        }

        int inucl=1;
        if(nt==1) inucl=0; // see above, program has formed single particle as nucleus, so filter them now
        if(iz==0) inucl=0; // an nucleus must has at least one nuclear charge
        if(in==0) inucl=0; // and there must be neutron
        if(inucl==0) continue; // no nucleus formed, skip this ic

        ENucleus* pc = new ENucleus(); // Record a new nucleus
        int l=isave; // point to the start of this nucleus in particles vector
        int ib=0;
        for(int i=0;i<nt;i++) { // for all particles in the nucleus
            int j=num[l];
            p[j]->ks=30; // The particle formed the nucleus die
            ib += p[j]->ibar;
            pc->add(p[j]); // Record the particle into nucleus
            l++;
        }
        setClusterProperties(pc);
        pclust.push_back(pc);

        //int i0=num[isave];
        //p[i0]=pc;
        //ENucleus cp = static_cast<ENucleus>(p[i0]);
        //cout << " cluster! ibar= "<<  pc->ibar << " iz= " << pc->iz
        //     << " in= " << pc->in << " iy= " << pc->iy
        //     << " px= " << pc->p[1]
        //     << " py= " << pc->p[2]
        //     << " pz= " << pc->p[3]
        //     << endl;
        //cout << " ex= " << cp.ex << endl;


    } // end for ic in nclst

    // Stained by Shu He. I disable the following line and add the second line
    //p.insert(p.end(),pclust.begin(),pclust.end());
    _clust.insert(_clust.end(),pclust.begin(),pclust.end());

    delete [] mscl;
    delete [] num;

    if(ibar0 < ibar) {
	cout << " baryon number does not conserve? ibar0= "<< ibar0
	    << " ibar= " << ibar
	    << " miss= " << miss
	    << " nv " << nv
	    << " nn " << nn
	    << endl;
    }

    nv=p.size();
    double px=0.0, py=0.0,pz=0.0, pe=0.0;
    for(int i=0;i<nv;i++) {

        if(p[i]->ks>10) continue;
        px += p[i]->p[1];
        py += p[i]->p[2];
        pz += p[i]->p[3];
        pe += p[i]->p[0];

        /*
        cout << "ks= " << p[i].ks << " kf= " << p[i].kf
            << " p= " << p[i].p[1]
            << endl;
            */

    }
    // Stained by Shu HE, sorry, I need more learning
    nv = _clust.size();
    for(int i=0;i<nv;i++) {

        if(_clust[i]->ks>10) continue;
        px += _clust[i]->p[1];
        py += _clust[i]->p[2];
        pz += _clust[i]->p[3];
        pe += _clust[i]->p[0];

        /*
        cout << "ks= " << p[i].ks << " kf= " << p[i].kf
            << " p= " << p[i].p[1]
            << endl;
            */

    }

    if(abs(px-px0)>1e-2 || abs(py-py0)>1e-2 || abs(pz-pz0)>1e-2 || abs(pe-pe0)>1e-2) {
        cout << " total momentum does not conserve" << endl;
        cout << " nclust= " << nclst << " miss= " << miss << endl;
        cout << " nv= " << nv
             << " size= " << p.size()
             << " nn= " << nn
             << endl;
        cout << " px= " << px-px0 << " py= " << py-py0 << " pz= " << pz-pz0 
             << " pe= " << pe-pe0 <<endl;
        exit(1);
    }


}
