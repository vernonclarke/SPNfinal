/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__glutsynapse
#define _nrn_initial _nrn_initial__glutsynapse
#define nrn_cur _nrn_cur__glutsynapse
#define _nrn_current _nrn_current__glutsynapse
#define nrn_jacob _nrn_jacob__glutsynapse
#define nrn_state _nrn_state__glutsynapse
#define _net_receive _net_receive__glutsynapse 
#define state state__glutsynapse 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define erev _p[0]
#define erev_columnindex 0
#define tau1_ampa _p[1]
#define tau1_ampa_columnindex 1
#define tau2_ampa _p[2]
#define tau2_ampa_columnindex 2
#define tau1_nmda _p[3]
#define tau1_nmda_columnindex 3
#define tau2_nmda _p[4]
#define tau2_nmda_columnindex 4
#define ratio _p[5]
#define ratio_columnindex 5
#define mg _p[6]
#define mg_columnindex 6
#define alpha _p[7]
#define alpha_columnindex 7
#define beta _p[8]
#define beta_columnindex 8
#define q _p[9]
#define q_columnindex 9
#define nmda_scale_factor _p[10]
#define nmda_scale_factor_columnindex 10
#define ampa_scale_factor _p[11]
#define ampa_scale_factor_columnindex 11
#define gmax_AMPA _p[12]
#define gmax_AMPA_columnindex 12
#define gmax_NMDA _p[13]
#define gmax_NMDA_columnindex 13
#define maxModNMDA _p[14]
#define maxModNMDA_columnindex 14
#define max2NMDA _p[15]
#define max2NMDA_columnindex 15
#define maxModAMPA _p[16]
#define maxModAMPA_columnindex 16
#define max2AMPA _p[17]
#define max2AMPA_columnindex 17
#define damod _p[18]
#define damod_columnindex 18
#define l1NMDA _p[19]
#define l1NMDA_columnindex 19
#define l2NMDA _p[20]
#define l2NMDA_columnindex 20
#define l1AMPA _p[21]
#define l1AMPA_columnindex 21
#define l2AMPA _p[22]
#define l2AMPA_columnindex 22
#define i _p[23]
#define i_columnindex 23
#define g _p[24]
#define g_columnindex 24
#define i_ampa _p[25]
#define i_ampa_columnindex 25
#define i_nmda _p[26]
#define i_nmda_columnindex 26
#define g_ampa _p[27]
#define g_ampa_columnindex 27
#define g_nmda _p[28]
#define g_nmda_columnindex 28
#define block _p[29]
#define block_columnindex 29
#define I _p[30]
#define I_columnindex 30
#define G _p[31]
#define G_columnindex 31
#define A _p[32]
#define A_columnindex 32
#define B _p[33]
#define B_columnindex 33
#define C _p[34]
#define C_columnindex 34
#define D _p[35]
#define D_columnindex 35
#define factor_nmda _p[36]
#define factor_nmda_columnindex 36
#define factor_ampa _p[37]
#define factor_ampa_columnindex 37
#define ical _p[38]
#define ical_columnindex 38
#define DA _p[39]
#define DA_columnindex 39
#define DB _p[40]
#define DB_columnindex 40
#define DC _p[41]
#define DC_columnindex 41
#define DD _p[42]
#define DD_columnindex 42
#define v _p[43]
#define v_columnindex 43
#define _g _p[44]
#define _g_columnindex 44
#define _tsav _p[45]
#define _tsav_columnindex 45
#define _nd_area  *_ppvar[0]._pval
#define _ion_ical	*_ppvar[2]._pval
#define _ion_dicaldv	*_ppvar[3]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static double _hoc_MgBlock(void*);
 static double _hoc_modulation(void*);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(Object* _ho) { void* create_point_process(int, Object*);
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt(void*);
 static double _hoc_loc_pnt(void* _vptr) {double loc_point_process(int, void*);
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(void* _vptr) {double has_loc_point(void*);
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(void* _vptr) {
 double get_loc_point_process(void*); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 "MgBlock", _hoc_MgBlock,
 "modulation", _hoc_modulation,
 0, 0
};
#define MgBlock MgBlock_glutsynapse
#define modulation modulation_glutsynapse
 extern double MgBlock( _threadargsproto_ );
 extern double modulation( _threadargsprotocomma_ double , double , double , double );
 /* declare global and static user variables */
#define ca_ratio_nmda ca_ratio_nmda_glutsynapse
 double ca_ratio_nmda = 0.1;
#define ca_ratio_ampa ca_ratio_ampa_glutsynapse
 double ca_ratio_ampa = 0.005;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "erev", "mV",
 "tau1_ampa", "ms",
 "tau2_ampa", "ms",
 "tau1_nmda", "ms",
 "tau2_nmda", "ms",
 "ratio", "1",
 "mg", "mM",
 "gmax_AMPA", "uS",
 "gmax_NMDA", "uS",
 "A", "uS",
 "B", "uS",
 "C", "uS",
 "D", "uS",
 "i", "nA",
 "g", "uS",
 0,0
};
 static double A0 = 0;
 static double B0 = 0;
 static double C0 = 0;
 static double D0 = 0;
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "ca_ratio_ampa_glutsynapse", &ca_ratio_ampa_glutsynapse,
 "ca_ratio_nmda_glutsynapse", &ca_ratio_nmda_glutsynapse,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(void* _vptr) {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[4]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"glutsynapse",
 "erev",
 "tau1_ampa",
 "tau2_ampa",
 "tau1_nmda",
 "tau2_nmda",
 "ratio",
 "mg",
 "alpha",
 "beta",
 "q",
 "nmda_scale_factor",
 "ampa_scale_factor",
 "gmax_AMPA",
 "gmax_NMDA",
 "maxModNMDA",
 "max2NMDA",
 "maxModAMPA",
 "max2AMPA",
 "damod",
 "l1NMDA",
 "l2NMDA",
 "l1AMPA",
 "l2AMPA",
 0,
 "i",
 "g",
 "i_ampa",
 "i_nmda",
 "g_ampa",
 "g_nmda",
 "block",
 "I",
 "G",
 0,
 "A",
 "B",
 "C",
 "D",
 0,
 0};
 static Symbol* _cal_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 46, _prop);
 	/*initialize range parameters*/
 	erev = 0;
 	tau1_ampa = 1.9;
 	tau2_ampa = 4.8;
 	tau1_nmda = 5.52;
 	tau2_nmda = 231;
 	ratio = 1;
 	mg = 1;
 	alpha = 0.096;
 	beta = 3.57;
 	q = 2;
 	nmda_scale_factor = 1;
 	ampa_scale_factor = 1;
 	gmax_AMPA = 0.001;
 	gmax_NMDA = 0.001;
 	maxModNMDA = 1;
 	max2NMDA = 1;
 	maxModAMPA = 1;
 	max2AMPA = 1;
 	damod = 0;
 	l1NMDA = 0;
 	l2NMDA = 0;
 	l1AMPA = 0;
 	l2AMPA = 0;
  }
 	_prop->param = _p;
 	_prop->param_size = 46;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_cal_sym);
 	_ppvar[2]._pval = &prop_ion->param[3]; /* ical */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_dicaldv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _net_receive(Point_process*, double*, double);
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _glutsynapse_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("cal", 2.0);
 	_cal_sym = hoc_lookup("cal_ion");
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 1,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 46, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "cal_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cal_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 glutsynapse /Users/euo9382/Documents/SPNmaster/mechanisms/single/glutsynapse.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[4], _dlist1[4];
 static int state(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset = 0; {
   DA = - A / tau1_nmda * q ;
   DB = - B / tau2_nmda * q ;
   DC = - C / tau1_ampa * q ;
   DD = - D / tau2_ampa * q ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
 DA = DA  / (1. - dt*( ( ( - 1.0 ) / tau1_nmda )*( q ) )) ;
 DB = DB  / (1. - dt*( ( ( - 1.0 ) / tau2_nmda )*( q ) )) ;
 DC = DC  / (1. - dt*( ( ( - 1.0 ) / tau1_ampa )*( q ) )) ;
 DD = DD  / (1. - dt*( ( ( - 1.0 ) / tau2_ampa )*( q ) )) ;
  return 0;
}
 /*END CVODE*/
 static int state (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) { {
    A = A + (1. - exp(dt*(( ( - 1.0 ) / tau1_nmda )*( q ))))*(- ( 0.0 ) / ( ( ( - 1.0 ) / tau1_nmda )*( q ) ) - A) ;
    B = B + (1. - exp(dt*(( ( - 1.0 ) / tau2_nmda )*( q ))))*(- ( 0.0 ) / ( ( ( - 1.0 ) / tau2_nmda )*( q ) ) - B) ;
    C = C + (1. - exp(dt*(( ( - 1.0 ) / tau1_ampa )*( q ))))*(- ( 0.0 ) / ( ( ( - 1.0 ) / tau1_ampa )*( q ) ) - C) ;
    D = D + (1. - exp(dt*(( ( - 1.0 ) / tau2_ampa )*( q ))))*(- ( 0.0 ) / ( ( ( - 1.0 ) / tau2_ampa )*( q ) ) - D) ;
   }
  return 0;
}
 
static void _net_receive (Point_process* _pnt, double* _args, double _lflag) 
{  double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _thread = (Datum*)0; _nt = (NrnThread*)_pnt->_vnt;   _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t; {
     if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = A;
    double __primary = (A + _args[0] * factor_nmda) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( ( - 1.0 ) / tau1_nmda )*( q ) ) ) )*( - ( 0.0 ) / ( ( ( - 1.0 ) / tau1_nmda )*( q ) ) - __primary );
    A += __primary;
  } else {
 A = A + _args[0] * factor_nmda ;
     }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = B;
    double __primary = (B + _args[0] * factor_nmda) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( ( - 1.0 ) / tau2_nmda )*( q ) ) ) )*( - ( 0.0 ) / ( ( ( - 1.0 ) / tau2_nmda )*( q ) ) - __primary );
    B += __primary;
  } else {
 B = B + _args[0] * factor_nmda ;
     }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = C;
    double __primary = (C + _args[0] * factor_ampa * ratio) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( ( - 1.0 ) / tau1_ampa )*( q ) ) ) )*( - ( 0.0 ) / ( ( ( - 1.0 ) / tau1_ampa )*( q ) ) - __primary );
    C += __primary;
  } else {
 C = C + _args[0] * factor_ampa * ratio ;
     }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = D;
    double __primary = (D + _args[0] * factor_ampa * ratio) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( ( - 1.0 ) / tau2_ampa )*( q ) ) ) )*( - ( 0.0 ) / ( ( ( - 1.0 ) / tau2_ampa )*( q ) ) - __primary );
    D += __primary;
  } else {
 D = D + _args[0] * factor_ampa * ratio ;
     }
 } }
 
double MgBlock ( _threadargsproto_ ) {
   double _lMgBlock;
 _lMgBlock = 1.0 / ( 1.0 + mg * exp ( - alpha * v ) / beta ) ;
   
return _lMgBlock;
 }
 
static double _hoc_MgBlock(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  MgBlock ( _p, _ppvar, _thread, _nt );
 return(_r);
}
 
double modulation ( _threadargsprotocomma_ double _lm1 , double _lm2 , double _ll1 , double _ll2 ) {
   double _lmodulation;
 _lmodulation = 1.0 + damod * ( ( _lm1 - 1.0 ) * _ll1 + ( _lm2 - 1.0 ) * _ll2 ) ;
   if ( _lmodulation < 0.0 ) {
     _lmodulation = 0.0 ;
     }
   
return _lmodulation;
 }
 
static double _hoc_modulation(void* _vptr) {
 double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   _p = ((Point_process*)_vptr)->_prop->param;
  _ppvar = ((Point_process*)_vptr)->_prop->dparam;
  _thread = _extcall_thread;
  _nt = (NrnThread*)((Point_process*)_vptr)->_vnt;
 _r =  modulation ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) , *getarg(3) , *getarg(4) );
 return(_r);
}
 
static int _ode_count(int _type){ return 4;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 4; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_cal_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_cal_sym, _ppvar, 3, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  A = A0;
  B = B0;
  C = C0;
  D = D0;
 {
   double _ltp ;
 if ( tau1_nmda / tau2_nmda > .9999 ) {
     tau1_nmda = .9999 * tau2_nmda ;
     }
   if ( tau1_ampa / tau2_ampa > .9999 ) {
     tau1_ampa = .9999 * tau2_ampa ;
     }
   A = 0.0 ;
   B = 0.0 ;
   _ltp = ( tau1_nmda * tau2_nmda ) / ( tau2_nmda - tau1_nmda ) * log ( tau2_nmda / tau1_nmda ) ;
   factor_nmda = - exp ( - _ltp / tau1_nmda ) + exp ( - _ltp / tau2_nmda ) ;
   factor_nmda = 1.0 / factor_nmda ;
   C = 0.0 ;
   D = 0.0 ;
   _ltp = ( tau1_ampa * tau2_ampa ) / ( tau2_ampa - tau1_ampa ) * log ( tau2_ampa / tau1_ampa ) ;
   factor_ampa = - exp ( - _ltp / tau1_ampa ) + exp ( - _ltp / tau2_ampa ) ;
   factor_ampa = 1.0 / factor_ampa ;
   }
 
}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _tsav = -1e20;
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   g_nmda = ( B - A ) * gmax_NMDA ;
   block = MgBlock ( _threadargs_ ) ;
   i_nmda = g_nmda * ( v - erev ) * block * nmda_scale_factor ;
   g_ampa = ( D - C ) * gmax_AMPA ;
   i_ampa = g_ampa * ( v - erev ) * ampa_scale_factor ;
   G = g_ampa + g_nmda ;
   I = i_ampa + i_nmda ;
   ical = i_ampa * ca_ratio_ampa + i_nmda * ca_ratio_nmda ;
   i = i_ampa * ( 1.0 - ca_ratio_ampa ) + i_nmda * ( 1.0 - ca_ratio_nmda ) ;
   }
 _current += i;
 _current += ical;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dical;
  _dical = ical;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dicaldv += (_dical - ical)/.001 * 1.e2/ (_nd_area);
 	}
 _g = (_g - _rhs)/.001;
  _ion_ical += ical * 1.e2/ (_nd_area);
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 {   state(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = A_columnindex;  _dlist1[0] = DA_columnindex;
 _slist1[1] = B_columnindex;  _dlist1[1] = DB_columnindex;
 _slist1[2] = C_columnindex;  _dlist1[2] = DC_columnindex;
 _slist1[3] = D_columnindex;  _dlist1[3] = DD_columnindex;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/euo9382/Documents/SPNmaster/mechanisms/single/glutsynapse.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "Updated Exp2Syn synapse with Mg-blocked nmda channel.\n"
  "\n"
  "Default values of parameters (time constants etc) set to match synaptic channels in \n"
  "striatal medium spiny neurons (Du et al., 2017; Chapman et al., 2003; Ding et al., 2008).\n"
  "\n"
  "Robert . Lindroos @ ki . se\n"
  "\n"
  "original comment:\n"
  "________________\n"
  "Two state kinetic scheme synapse described by rise time tau1,\n"
  "and decay time constant tau2. The normalized peak conductance is 1.\n"
  "Decay time MUST be greater than rise time.\n"
  "\n"
  "The solution of A->G->bath with rate constants 1/tau1 and 1/tau2 is\n"
  " A = a*exp(-t/tau1) and\n"
  " G = a*tau2/(tau2-tau1)*(-exp(-t/tau1) + exp(-t/tau2))\n"
  "	where tau1 < tau2\n"
  "\n"
  "If tau2-tau1 is very small compared to tau1, this is an alphasynapse with time constant tau2.\n"
  "If tau1/tau2 is very small, this is single exponential decay with time constant tau2.\n"
  "\n"
  "The factor is evaluated in the initial block \n"
  "such that an event of weight 1 generates a\n"
  "peak conductance of 1.\n"
  "\n"
  "Because the solution is a sum of exponentials, the\n"
  "coupled equations can be solved as a pair of independent equations\n"
  "by the more efficient cnexp method.\n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "NEURON {\n"
  "	POINT_PROCESS glutsynapse\n"
  "	RANGE tau1_ampa, tau2_ampa, tau1_nmda, tau2_nmda\n"
  "	RANGE erev, g, i\n"
  "	RANGE i_ampa, i_nmda, g_ampa, g_nmda, ratio, I, G, mg, q, block, alpha, beta\n"
  "	RANGE ampa_scale_factor, nmda_scale_factor\n"
  "    RANGE damod, maxModNMDA,max2NMDA,maxModAMPA,max2AMPA,l1NMDA,l2NMDA,l1AMPA,l2AMPA\n"
  "	RANGE gmax_AMPA, gmax_NMDA\n"
  "	\n"
  "	NONSPECIFIC_CURRENT i\n"
  "	USEION cal WRITE ical VALENCE 2\n"
  "}\n"
  "\n"
  "\n"
  "UNITS {\n"
  "	(nA) = (nanoamp)\n"
  "	(mV) = (millivolt)\n"
  "	(uS) = (microsiemens)\n"
  "}\n"
  "\n"
  "\n"
  "PARAMETER {\n"
  "	erev        = 0.0       (mV)\n"
  "	\n"
  "	tau1_ampa   = 1.9       (ms)\n"
  "    tau2_ampa   = 4.8       (ms)  : tau2 > tau1\n"
  "    tau1_nmda   = 5.52      (ms)  : Chapman et al 2003; table 1, adult rat (rise time, rt = 12.13. rt ~= 2.197*tau (wiki;rise time) -> tau = 12.13 / 2.197 ~= 5.52\n"
  "    tau2_nmda   = 231       (ms)  : Chapman et al 2003 (table 1; adult)\n"
  "    \n"
  "    ratio       = 1         (1)   : both components give same maximal amplitude of current\n"
  "    mg          = 1         (mM)\n"
  "    alpha       = 0.096			  : was 0.062\n"
  "    beta        = 3.57\n"
  "    q           = 2               : approx room temp -> \n"
  "    \n"
  "    nmda_scale_factor = 1\n"
  "    ampa_scale_factor = 1\n"
  "    \n"
  "    ca_ratio_ampa = 0.005\n"
  "    ca_ratio_nmda = 0.1\n"
  "\n"
  "	gmax_AMPA = 0.001		(uS)\n"
  "	gmax_NMDA = 0.001		(uS)\n"
  "    \n"
  "    maxModNMDA  = 1\n"
  "    max2NMDA    = 1\n"
  "    maxModAMPA  = 1\n"
  "    max2AMPA    = 1\n"
  "    damod       = 0\n"
  "    l1NMDA      = 0\n"
  "    l2NMDA      = 0\n"
  "    l1AMPA      = 0\n"
  "    l2AMPA      = 0\n"
  "}\n"
  "\n"
  "\n"
  "ASSIGNED {\n"
  "	v (mV)\n"
  "	i (nA)\n"
  "	g (uS)\n"
  "	factor_nmda\n"
  "	factor_ampa\n"
  "	i_ampa\n"
  "	i_nmda\n"
  "	g_ampa\n"
  "	g_nmda\n"
  "	block\n"
  "	I\n"
  "	G\n"
  "	ical (nA)\n"
  "}\n"
  "\n"
  "\n"
  "STATE {\n"
  "	A (uS)\n"
  "	B (uS)\n"
  "	C (uS)\n"
  "	D (uS)\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "INITIAL {\n"
  "	LOCAL tp\n"
  "	if (tau1_nmda/tau2_nmda > .9999) {\n"
  "		tau1_nmda = .9999*tau2_nmda\n"
  "	}\n"
  "	if (tau1_ampa/tau2_ampa > .9999) {\n"
  "		tau1_ampa = .9999*tau2_ampa\n"
  "	}\n"
  "	\n"
  "	: NMDA\n"
  "	A           = 0\n"
  "	B           = 0\n"
  "	tp          = (tau1_nmda*tau2_nmda)/(tau2_nmda - tau1_nmda) * log(tau2_nmda/tau1_nmda)\n"
  "	factor_nmda = -exp(-tp/tau1_nmda) + exp(-tp/tau2_nmda)\n"
  "	factor_nmda = 1/factor_nmda\n"
  "	\n"
  "	: AMPA\n"
  "	C           = 0\n"
  "	D           = 0\n"
  "	tp          = (tau1_ampa*tau2_ampa)/(tau2_ampa - tau1_ampa) * log(tau2_ampa/tau1_ampa)\n"
  "	factor_ampa = -exp(-tp/tau1_ampa) + exp(-tp/tau2_ampa)\n"
  "	factor_ampa = 1/factor_ampa\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE state METHOD cnexp\n"
  "	\n"
  "	: NMDA\n"
  "	: commenting out to input gmax_nmda\n"
  "	: g_nmda = (B - A) * modulation(maxModNMDA,max2NMDA,l1NMDA,l2NMDA)\n"
  "	g_nmda = (B - A) * gmax_NMDA\n"
  "	block  = MgBlock()\n"
  "	i_nmda = g_nmda * (v - erev) * block * nmda_scale_factor\n"
  "	\n"
  "	: AMPA\n"
  "	: commenting out to input gmax_ampa\n"
  "	: g_ampa = (D - C) * modulation(maxModAMPA,max2AMPA,l1AMPA,l2AMPA)\n"
  "	g_ampa = (D - C) * gmax_AMPA\n"
  "	i_ampa = g_ampa * (v - erev) * ampa_scale_factor\n"
  "	\n"
  "	: total current\n"
  "	G = g_ampa + g_nmda\n"
  "	I = i_ampa + i_nmda\n"
  "	\n"
  "	: splitting in ca and non ca currents\n"
  "	ical = i_ampa*ca_ratio_ampa  + i_nmda*ca_ratio_nmda\n"
  "    i = i_ampa*(1-ca_ratio_ampa) + i_nmda*(1-ca_ratio_nmda)\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "DERIVATIVE state {\n"
  "	A' = -A/tau1_nmda*q\n"
  "	B' = -B/tau2_nmda*q\n"
  "	C' = -C/tau1_ampa*q\n"
  "	D' = -D/tau2_ampa*q\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "NET_RECEIVE(weight (uS)) {\n"
  "	A = A + weight*factor_nmda\n"
  "	B = B + weight*factor_nmda\n"
  "	C = C + weight*factor_ampa*ratio\n"
  "	D = D + weight*factor_ampa*ratio\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "FUNCTION MgBlock() {\n"
  "    \n"
  "    MgBlock = 1 / (1 + mg * exp(-alpha * v) / beta )\n"
  "    \n"
  "}\n"
  "\n"
  "FUNCTION modulation(m1,m2,l1,l2) {\n"
  "    : returns modulation factor\n"
  "    \n"
  "    modulation = 1 + damod * ( (m1-1)*l1 + (m2-1)*l2 )\n"
  "    if (modulation < 0) {\n"
  "        modulation = 0\n"
  "    } \n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "\n"
  ;
#endif
