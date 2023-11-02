#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _Im_reg(void);
extern void _bk_reg(void);
extern void _cadyn_reg(void);
extern void _cal12_reg(void);
extern void _cal13_reg(void);
extern void _caldyn_reg(void);
extern void _can_reg(void);
extern void _car_reg(void);
extern void _cav32_reg(void);
extern void _cav33_reg(void);
extern void _gabasynapse_reg(void);
extern void _glutsynapse_reg(void);
extern void _kaf_reg(void);
extern void _kas_reg(void);
extern void _kcnq_reg(void);
extern void _kdr_reg(void);
extern void _kir_reg(void);
extern void _naf_reg(void);
extern void _nap_reg(void);
extern void _sk_reg(void);
extern void _tonicgaba1_reg(void);
extern void _tonicgaba2_reg(void);
extern void _vecevent_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"Im.mod\"");
    fprintf(stderr, " \"bk.mod\"");
    fprintf(stderr, " \"cadyn.mod\"");
    fprintf(stderr, " \"cal12.mod\"");
    fprintf(stderr, " \"cal13.mod\"");
    fprintf(stderr, " \"caldyn.mod\"");
    fprintf(stderr, " \"can.mod\"");
    fprintf(stderr, " \"car.mod\"");
    fprintf(stderr, " \"cav32.mod\"");
    fprintf(stderr, " \"cav33.mod\"");
    fprintf(stderr, " \"gabasynapse.mod\"");
    fprintf(stderr, " \"glutsynapse.mod\"");
    fprintf(stderr, " \"kaf.mod\"");
    fprintf(stderr, " \"kas.mod\"");
    fprintf(stderr, " \"kcnq.mod\"");
    fprintf(stderr, " \"kdr.mod\"");
    fprintf(stderr, " \"kir.mod\"");
    fprintf(stderr, " \"naf.mod\"");
    fprintf(stderr, " \"nap.mod\"");
    fprintf(stderr, " \"sk.mod\"");
    fprintf(stderr, " \"tonicgaba1.mod\"");
    fprintf(stderr, " \"tonicgaba2.mod\"");
    fprintf(stderr, " \"vecevent.mod\"");
    fprintf(stderr, "\n");
  }
  _Im_reg();
  _bk_reg();
  _cadyn_reg();
  _cal12_reg();
  _cal13_reg();
  _caldyn_reg();
  _can_reg();
  _car_reg();
  _cav32_reg();
  _cav33_reg();
  _gabasynapse_reg();
  _glutsynapse_reg();
  _kaf_reg();
  _kas_reg();
  _kcnq_reg();
  _kdr_reg();
  _kir_reg();
  _naf_reg();
  _nap_reg();
  _sk_reg();
  _tonicgaba1_reg();
  _tonicgaba2_reg();
  _vecevent_reg();
}

#if defined(__cplusplus)
}
#endif
