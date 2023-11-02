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
 
#define nrn_init _nrn_init__gaba2
#define _nrn_initial _nrn_initial__gaba2
#define nrn_cur _nrn_cur__gaba2
#define _nrn_current _nrn_current__gaba2
#define nrn_jacob _nrn_jacob__gaba2
#define nrn_state _nrn_state__gaba2
#define _net_receive _net_receive__gaba2 
 
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
#define gbar _p[0]
#define gbar_columnindex 0
#define i _p[1]
#define i_columnindex 1
#define minf _p[2]
#define minf_columnindex 2
#define e _p[3]
#define e_columnindex 3
#define v _p[4]
#define v_columnindex 4
#define _g _p[5]
#define _g_columnindex 5
 
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
 static void _hoc_a(void);
 static void _hoc_b(void);
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

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_gaba2", _hoc_setdata,
 "a_gaba2", _hoc_a,
 "b_gaba2", _hoc_b,
 0, 0
};
#define a a_gaba2
#define b b_gaba2
 extern double a( _threadargsprotocomma_ double );
 extern double b( _threadargsprotocomma_ double );
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gbar_gaba2", "siemens/cm2",
 "i_gaba2", "milliamp/cm2",
 "e_gaba2", "mV",
 0,0
};
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
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
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"gaba2",
 "gbar_gaba2",
 0,
 "i_gaba2",
 "minf_gaba2",
 "e_gaba2",
 0,
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 6, _prop);
 	/*initialize range parameters*/
 	gbar = 0;
 	_prop->param = _p;
 	_prop->param_size = 6;
 
}
 static void _initlists();
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _tonicgaba2_reg() {
	int _vectorized = 1;
  _initlists();
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 6, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 gaba2 /Users/euo9382/Documents/SPNmaster/mechanisms/single/tonicgaba2.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "gaba leak ";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
double a ( _threadargsprotocomma_ double _lv ) {
   double _la;
 double _lx ;
 if ( fabs ( _lv + 20.0 ) > 1e-5 ) {
     _lx = 0.1 * ( _lv + 20.0 ) ;
     }
   else {
     _lx = 0.1 ;
     }
   _la = ( 50.0 * _lx / ( 1.0 - exp ( - _lx ) ) ) ;
   
return _la;
 }
 
static void _hoc_a(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  a ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double b ( _threadargsprotocomma_ double _lv ) {
   double _lb;
 double _lx ;
 if ( fabs ( _lv - 10.0 ) > 1e-5 ) {
     _lx = - 0.08 * ( _lv - 10.0 ) ;
     }
   else {
     _lx = - 0.08 ;
     }
   _lb = ( 20.0 * _lx / ( 1.0 - exp ( - _lx ) ) ) ;
   
return _lb;
 }
 
static void _hoc_b(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  b ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{

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
   minf = a ( _threadargscomma_ v ) / ( a ( _threadargscomma_ v ) + b ( _threadargscomma_ v ) ) ;
   i = gbar * minf * ( v - e ) ;
   }
 _current += i;

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
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
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

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/euo9382/Documents/SPNmaster/mechanisms/single/tonicgaba2.mod";
static const char* nmodl_file_text = 
  "TITLE gaba leak \n"
  "\n"
  "COMMENT\n"
  "If want ohmic gaba use tonicgaba1\n"
  "From Pavlov et al., 2009\n"
  "Outwardly Rectifying Tonically Active GABAA Receptors in Pyramidal Cells Modulate Neuronal Offset, Not Gain \n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON{\n"
  "    SUFFIX gaba2\n"
  "    NONSPECIFIC_CURRENT i\n"
  "    RANGE  i, e, gbar, minf\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mA)  = (milliamp)\n"
  "	(mV)  =  (millivolt)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "    gbar = 0 (siemens/cm2)  \n"
  "}\n"
  "\n"
  "ASSIGNED{\n"
  "    v (mV)\n"
  "    i (milliamp/cm2)\n"
  "    minf     \n"
  "    e (mV)\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "    minf=a(v)/(a(v)+b(v))\n"
  "    i=gbar*minf*(v-e)\n"
  "}\n"
  "\n"
  "FUNCTION a(v(mV)) {\n"
  "     LOCAL x\n"
  "          if (fabs(v+20) > 1e-5) {\n"
  "               x = 0.1 * (v + 20)\n"
  "         }else{\n"
  "              x = 0.1\n"
  "         }\n"
  "    a = (50 * x / (1 - exp(-x)))\n"
  "}\n"
  "\n"
  "FUNCTION b(v(mV)) {\n"
  "    LOCAL x\n"
  "         if (fabs(v-10) > 1e-5) {\n"
  "              x = -0.08 * (v - 10)\n"
  "         } else {\n"
  "             x = -0.08\n"
  "         }\n"
  "    b = (20 * x / (1 - exp(-x)))\n"
  "}\n"
  ;
#endif
