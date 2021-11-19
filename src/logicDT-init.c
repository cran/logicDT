#include "sa-greedy.h"
#include <R_ext/Rdynload.h>

void R_init_logicDT(DllInfo *info) {
  R_RegisterCCallable("logicDT", "fitPETs_",  (DL_FUNC) &fitPETs_);
  R_RegisterCCallable("logicDT", "fitPET_",  (DL_FUNC) &fitPET_);
  R_RegisterCCallable("logicDT", "predict_",  (DL_FUNC) &predict_);
  R_RegisterCCallable("logicDT", "predictEnsemble_",  (DL_FUNC) &predictEnsemble_);
  R_RegisterCCallable("logicDT", "getDesignMatrix_",  (DL_FUNC) &getDesignMatrix_);
  R_RegisterCCallable("logicDT", "vim_permutation_",  (DL_FUNC) &vim_permutation_);
  R_RegisterCCallable("logicDT", "geneticProgramming_",  (DL_FUNC) &geneticProgramming_);
  R_RegisterCCallable("logicDT", "predictGP_",  (DL_FUNC) &predictGP_);
  R_RegisterCCallable("logicDT", "simulatedAnnealing_",  (DL_FUNC) &simulatedAnnealing_);
  R_RegisterCCallable("logicDT", "greedySearch_",  (DL_FUNC) &greedySearch_);
  R_RegisterCCallable("logicDT", "calcAUC_",  (DL_FUNC) &calcAUC_);
  R_RegisterCCallable("logicDT", "fit4plModel_",  (DL_FUNC) &fit4plModel_);

  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
