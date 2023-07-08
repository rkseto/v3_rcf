#ifndef TJAM_HD
#define TJAM_HD
#include <string>
#include <cstring>

namespace TJAM{

extern "C" {
	void c_jammain_();
	void c_jaminit_();
	void c_jamfin_();
	void c_jamevt_(int*);
	void c_jamcomp_(int*, int*);
	void c_p_(int*, int*, double*);
	void c_k_(int*, int*, double*);
	void c_v_(int*, int*, double*);
	void c_r_(int*, int*, double*);
	void c_wp_(int*, double*);
	void c_wk_(int*, double*);
	void c_wv_(int*, double*);
	void c_wr_(int*, double*);
	void c_cdata_(int*, int*, float*, float*, float*);
	void c_wdata_(int*, int*, float*, float*, float*, float*, float*, float*);
	void c_wedata_(int *cnpart, int *pid_ary, float *px_ary, float *py_ary, float *pz_ary, float *e_ary, float *m_ary, float *chg_ary, float *pt_ary, float *y_ary, float *eta_ary);
	void c_wxdata_(int *cnpart, int *ks_ary, int *pid_ary,
			float *x_ary, float *y_ary, float *z_ary, float *t_ary, float *ft_ary,
			float *px_ary, float *py_ary, float *pz_ary,
			float *e_ary, float *m_ary, float *chg_ary,
			float *pt_ary, float *ycm_ary, float *eta_ary);
	void c_cldata_(int *cnpart, int *ks_ary, int *pid_ary,
			float *x_ary, float *y_ary, float *z_ary, float *t_ary,
			float *fx_ary, float *fy_ary, float *fz_ary, float *ft_ary,
			float *px_ary, float *py_ary, float *pz_ary);
	void c_kchg_(int*, int*, double*);
	void c_pmas_(int*, int*, double*);
	void c_parf_(int*, double*);
	void c_vckm_(int*, int*, double*);
	void c_mdcy_(int*, int*, double*);
	void c_mdme_(int*, int*, double*);
	void c_brat_(int*, double*);
	void c_kfdp_(int*, int*, double*);
	void c_vq_(int*, int*, double*);
	void c_kq_(int*, int*, double*);
	void c_nv_(double*);
	void c_nbary_(double*);
	void c_nmeson_(double*);
	void c_mxcoll_(int*);
	void c_mentry_(double*);
	void c_icoll_(int*, int*, double*);
	void c_coll_(int*, int*, double*);
	void c_para_(int* mevent, double* bmin, double* bmax, char* cwin, char* frame, char* proj, char* targ, double* dt, int* nstep, double* nv, double* nbary, double* nmeson, int* mxcoll, double *mentry);
	void c_s_para_(int* mevent, double* bmin, double* bmax, char* cwin, char* frame, char* proj, char* targ, double* dt, int* nstep);
	void c_mstc_(int*, double*);
	void c_parc_(int*, double*);
	void c_mstd_(int*, double*);
	void c_pard_(int*, double*);
	void c_pare_(int*, double*);
	void c_mste_(int*, double*);
	void c_s_mstc_(int*, double*);
	void c_s_parc_(int*, double*);
	void c_s_mstd_(int*, double*);
	void c_s_pard_(int*, double*);
	void c_s_pare_(int*, double*);
	void c_s_mste_(int*, double*);
	void c_s_mdcy_(int*, int*, double*);
	void c_s_mdme_(int*, int*, double*);

	void c_on_timeana_();
	void c_off_timeana_();
	void c_tphase_(float* catime, int* Npart, int* cmul, int* ks_ary, int* pid_ary, float* rx_ary, float* ry_ary, float* rz_ary, float* rt_ary, float* rft_ary,
			float* px_ary, float* py_ary, float* pz_ary, float* e_ary, float* m_ary, float* chge_ary,
			float* pt_ary, float* ycm_ary, float* eta_ary);
	//void c_s_time_function_();
    void c_quiet_file_out_(int* id);
//////////////////////////////////////
// functions call by cjam.f
	void c_jam_time_record_(double*);
}

typedef void (*ValDeg1)(int*, double*);
typedef void (*ValDeg2)(int*, int*, double*);

class ValProxy{
	public:
	typedef void (ValProxy::*DegFunc)(double&) const;
	ValProxy(ValDeg1 sHandle, ValDeg1 gHandle) {
		sdeg1 = sHandle;
		gdeg1 = gHandle;
		sdeg = &ValProxy::_sdeg_1d;
		gdeg = &ValProxy::_gdeg_1d;
	}
	ValProxy(ValDeg2 sHandle, ValDeg2 gHandle) {
		sdeg2 = sHandle;
		gdeg2 = gHandle;
		sdeg = &ValProxy::_sdeg_2d;
		gdeg = &ValProxy::_gdeg_2d;
	}
	const ValProxy & operator = (double val) const {
		(this->*sdeg)(val);
		return *this;
	}
	ValProxy& at(int idx) {
		index = idx;
		return *this;
	}
	ValProxy& at(int idx, int idy) {
		index = idx;
		index2 = idy;
		return *this;
	}
	operator double() const {
		double ret;
		(this->*gdeg)(ret);
		return ret;
	}
	private:
	DegFunc sdeg, gdeg;
	ValDeg1 sdeg1, gdeg1;
	ValDeg2 sdeg2, gdeg2;
	int index, index2;
	void _sdeg_1d(double &val) const {
		int idx = index;
		double v = val;
		sdeg1(&idx, &v);
	}
	void _sdeg_2d(double &val) const {
		int idx1 = index;
		int idx2 = index2;
		double v = val;
		sdeg2(&idx1, &idx2, &v);
	}
	void _gdeg_1d(double &val) const {
		int idx = index;
		double ret;
		gdeg1(&idx, &ret);
		val = ret;
	}
	void _gdeg_2d(double &val) const {
		int idx1 = index;
		int idx2 = index2;
		double ret;
		gdeg2(&idx1, &idx2, &ret);
		val = ret;
	}
};

ValProxy vmstc(c_s_mstc_, c_mstc_), vparc(c_s_parc_, c_parc_),
	 vmstd(c_s_mstd_, c_mstd_), vpard(c_s_pard_, c_pard_),
	 vpare(c_s_pare_, c_pare_), vmste(c_s_mste_, c_mste_),
	 vmdcy(c_s_mdcy_, c_mdcy_), vmdme(c_s_mdme_, c_mdme_);

class TJam;

typedef void (*ctimeFuncT)(TJam&, float);
TJam* _TJam_tcallback_instance;
typedef void (*ValDeg)(int*, double*);

class TJam  {
	public:
	TJam() : record_time_evolution(false), record_time_function(NULL) {
		c_jammain_();
		get_para();
		c_off_timeana_();
		_TJam_tcallback_instance = NULL;
	}
	void jaminit() {
		set_para();
		if (record_time_evolution) {
			c_on_timeana_();
			_TJam_tcallback_instance=this;
                }
		c_jaminit_();
		get_para();
		_iev = 0;
	}
	void jamfin() {
		set_para();
		c_jamfin_();
		get_para();
	}
	void jamevt(const int &i) {
		int iv = i;
		set_para();
		c_jamevt_(&iv);
		get_para();
		_iev = i;
	}
    void quiet_file_out(int id) {
        c_quiet_file_out_(&id);
    }
	double p(const int& v, const int& i) const {
		double val;
		int iv = v, ii = i;
		c_p_(&iv,&ii,&val);
		return val;
	}
	double k(const int& v, const int& i) const {
		double val;
		int iv = v, ii = i;
		c_k_(&iv,&ii,&val);
		return val;
	}
	double v(const int& v, const int& i) const {
		double val;
		int iv = v, ii = i;
		c_v_(&iv,&ii,&val);
		return val;
	}
	double r(const int& v, const int& i) const {
		double val;
		int iv = v, ii = i;
		c_r_(&iv,&ii,&val);
		return val;
	}
	int mul() const {
		double val;
		c_nv_(&val);
		return int(val);
	}
	void write_p(const int& v, double *array) const {
		int iv = v;
		c_wp_(&iv, array);
	}
	void write_k(const int& v, double *array) const {
		int iv = v;
		c_wk_(&iv, array);
	}
	void write_v(const int& v, double *array) const {
		int iv = v;
		c_wv_(&iv, array);
	}
	void write_r(const int& v, double *array) const {
		int iv = v;
		c_wr_(&iv, array);
	}
	void write_cdata(int* Npart, int *pid, float *px, float *py, float *pz) {
		c_cdata_(Npart, pid, px, py, pz);
	}
	void write_pdata(int* Npart, int *pid, float *px, float *py, float *pz, float *e, float *mass, float *charge) {
		c_wdata_(Npart, pid, px, py, pz, e, mass, charge);
	}
    void write_pedata(int *Npart, int *pid, float *px, float *py, float *pz, float *e, float *mass, float *charge, float *pt, float *ycm, float *eta) {
		c_wedata_(Npart, pid, px, py, pz, e, mass, charge, pt, ycm, eta);
    }
    void write_pxdata(int *Npart, int *ks, int *pid, float *rx, float *ry, float *rz, float *rt, float *rft,
            float *px, float *py, float *pz, float *e, float *mass, float *charge,
            float *pt, float *ycm, float *eta) {
        c_wxdata_(Npart, ks, pid, rx, ry, rz, rt, rft, px, py, pz, e, mass, charge, pt, ycm, eta);
    }
    void write_cldata(int *Npart, int *ks, int *pid, float *rx, float *ry, float *rz, float *rt, 
		        float *frx, float *fry, float *frz, float *frt,
            float *px, float *py, float *pz) {
        c_cldata_(Npart, ks, pid, rx, ry, rz, rt, frx, fry, frz, frt, px, py, pz);
    }
	void write_time(float* catime, int* Npart, int* cmul, int *ks_ary, int* pid_ary,
			float* rx_ary, float* ry_ary, float* rz_ary, float *rt_ary, float *rft_ary,
			float* px_ary, float* py_ary, float* pz_ary,
			float* e_ary, float* m_ary, float* chge_ary,
			float* pt_ary, float* ycm_ary, float* eta_ary) {
		c_tphase_(catime, Npart, cmul, ks_ary, pid_ary,
			rx_ary, ry_ary, rz_ary, rt_ary, rft_ary,
			px_ary, py_ary, pz_ary,
			e_ary, m_ary, chge_ary,
			pt_ary, ycm_ary, eta_ary);
	}
	double kchg(const int &v, const int &i) const {
		double val;
		int iv = v, ii = i;
		c_kchg_(&iv,&ii,&val);
		return val;
	}
	double pmas(const int &v, const int &i) const {
		double val;
		int iv = v, ii = i;
		c_pmas_(&iv,&ii,&val);
		return val;
	}
	double parf(const int &v) const {
		double val;
		int iv = v;
		c_parf_(&iv, &val);
		return val;
	}
	double vckm(const int &v, const int &i) const {
		double val;
		int iv = v, ii = i;
		c_vckm_(&iv,&ii,&val);
		return val;
	}
	double brat(const int &v) const {
		double val;
		int iv = v;
		c_brat_(&iv, &val);
		return val;
	}
	double kfdp(const int &v, const int &i) const {
		double val;
		int iv = v, ii = i;
		c_kfdp_(&iv,&ii,&val);
		return val;
	}
	double vq(const int &v, const int &i) const {
		double val;
		int iv = v, ii = i;
		c_vq_(&iv,&ii,&val);
		return val;
	}
	double kq(const int &v, const int &i) const {
		double val;
		int iv = v, ii = i;
		c_kq_(&iv,&ii,&val);
		return val;
	}
	int jamcomp(const int& kf) const {
		int k = kf, ret;
		c_jamcomp_(&k, &ret);
		return ret;
	}
	double icoll(const int &v, const int &i) const {
		int iv =v, ii = i;
		double ret;
		c_icoll_(&iv, &ii, &ret);
		return ret;
	}
	double coll(const int &v, const int &i) const {
		int iv =v, ii = i;
		double ret;
		c_coll_(&iv, &ii, &ret);
		return ret;
	}
	ValProxy & mstc(const int& i) const {
		return vmstc.at(i);
	}
	ValProxy & parc(const int& i) const {
		return vparc.at(i);
	}
	ValProxy & mstd(const int& i) const {
		return vmstd.at(i);
	}
	ValProxy & pard(const int& i) const {
		return vpard.at(i);
	}
	ValProxy & mste(const int& i) const {
		return vmste.at(i);
	}
	ValProxy & pare(const int& i) const {
		return vpare.at(i);
	}
	ValProxy & mdcy(const int &v, const int &i) const {
		return vmdcy.at(v, i);
	}
	ValProxy & mdme(const int &v, const int &i) const {
		return vmdme.at(v, i);
	}
	int iev() const { return _iev; }
	int maxv() const { return 30000; }
	int mevent, nstep, mxcoll;
	double bmin, bmax, dt, nv, nbary, nmeson, mentry;
	std::string cwin, frame, proj, targ;
	bool record_time_evolution;
	ctimeFuncT record_time_function;
	private:
	char _cwin[20], _frame[20], _proj[20], _targ[20];
	int _iev;
	void cstr_ck() {
		sprintf(_cwin, "%-15s", cwin.c_str());
		sprintf(_frame, "%-8s", frame.c_str());
		sprintf(_proj, "%-8s", proj.c_str());
		sprintf(_targ, "%-8s", targ.c_str());
		_cwin[19]='\0';
		_frame[19]='\0';
		_proj[19]='\0';
		_targ[19]='\0';
	}
	void cstr_rt() {
		_cwin[19]='\0';
		_frame[19]='\0';
		_proj[19]='\0';
		_targ[19]='\0';
		cwin=_cwin;
		frame=_frame;
		proj=_proj;
		targ=_targ;
	}
	void set_para() {
		cstr_ck();
		c_s_para_(&mevent, &bmin, &bmax, _cwin, _frame, _proj, _targ, &dt, &nstep);
	}
	void get_para() {
		c_para_(&mevent, &bmin, &bmax, _cwin, _frame, _proj, _targ, &dt, &nstep, &nv, &nbary, &nmeson, &mxcoll, &mentry);
		cstr_rt();
	}
};

extern "C" {
void c_jam_time_record_(double *atime) {
	_TJam_tcallback_instance->record_time_function(*_TJam_tcallback_instance, float(*atime));
}
}

} // end of namespace TJAM

#endif
