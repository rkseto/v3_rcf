#ifndef __PARTICLES_H_HEAD_PROTECT_
#define __PARTICLES_H_HEAD_PROTECT_

const int Proton = 2212; //p
const int Neutron = 2112; //n
const int Deltapp = 2224; //\Delta^{++}
const int Deltap = 2214; //\Delta^{+}
const int Delta0 = 2114; //\Delta^{0}
const int Deltam = 1114; //\Delta^(-)
const int Lambda = 3122; //\Lambda
const int Sigmap = 3222; //\Sigma^{+}
const int Sigma0 = 3212; //\Sigma^{0}
const int Sigmam = 3112; //\Sigma^{-}
const int Xi0 = 3322; //\Xi^{0}
const int Xim = 3312; //\Xi^{-}
const int Omegam = 3334; //\Omega^{-}
const int Pion = 211; //\pi
const int Kaon = 321; //K

bool isProton(int& pid){
	if(abs(pid)==Proton) return true;
	return false;
}
bool isProNeu(int& pid) {
	int apid = abs(pid);
	if (apid == Proton || apid == Neutron) return true;
	return false;
}
bool isProLam(int& pid) {
	int apid = abs(pid);
	if (apid == Proton || apid == Lambda) return true;
	return false;
}
bool isProDel(int& pid) {
	int apid = abs(pid);
	if (apid == Proton || apid == Deltapp || apid == Deltap || apid == Delta0 || apid == Deltam) return true;
	return false;
}
/*
bool isBaryon(int& pid){
	int apid = abs(pid);
	if(apid==2212||apid==2112||apid==2224||apid==2214||apid==2114||apid==1114||apid==3122||apid==3222||apid==3212||apid==3112||apid==3322||apid==3312||apid==3334)
		return true;
	return false;
    int apid = abs(pid);
    if ((apid % 10000)/1000 != 0) return true;
    return false;
}
*/
bool isBaryon(int& pid) {
    // From Pythia 8: ParticleData.cc
    int idSave = abs(pid);
    if (idSave <= 1000 || (idSave >= 1000000 && idSave <= 9000000)
        || idSave >= 9900000 ) return false;
    if (idSave%10 == 0 || (idSave/10)%10 == 0 || (idSave/100)%10 == 0
        || (idSave/1000)%10 == 0) return false;
    return true;
}
bool isMeson(int& pid) {
    return (pid / 1000 == 0 && pid % 1000 != 0);
}
bool isNeutron(int& pid) {
	int apid = abs(pid);
	if(apid==Neutron) return true;
	return false;
}
bool isDeltapp(int& pid) {
	int apid = abs(pid);
	if(apid==Deltapp) return true;
	return false;
}
bool isDeltap(int& pid) {
	int apid = abs(pid);
	if(apid==Deltap) return true;
	return false;
}
bool isDelta0(int& pid) {
	int apid = abs(pid);
	if(apid==Delta0) return true;
	return false;
}
bool isDeltam(int& pid) {
	int apid = abs(pid);
	if(apid==Deltam) return true;
	return false;
}
bool isLambda(int& pid) {
	int apid = abs(pid);
	if(apid==Lambda) return true;
	return false;
}
bool isSigmap(int& pid) {
	int apid = abs(pid);
	if(apid==Sigmap) return true;
	return false;
}
bool isSigma0(int& pid) {
	int apid = abs(pid);
	if(apid==Sigma0) return true;
	return false;
}
bool isSigmam(int& pid) {
	int apid = abs(pid);
	if(apid==Sigmam) return true;
	return false;
}
bool isXi0(int& pid) {
	int apid = abs(pid);
	if(apid==Xi0) return true;
	return false;
}
bool isXim(int& pid) {
	int apid = abs(pid);
	if(apid==Xim) return true;
	return false;
}
bool isOmegam(int& pid) {
	int apid = abs(pid);
	if(apid==Omegam) return true;
	return false;
}

#endif
