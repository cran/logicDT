#include "logicDT.h"

typedef struct _functional_grouped
{
  functional* fpl;
  double* beta;
  int n_groups;
} functional_grouped;
typedef struct _dataset_grouped
{
  dataset* data;
  int* G;
  int n_groups;
} dataset_grouped;

double binLogLikelihood2(int n, double* par, void* ex2);
void binLogLikelihoodGrad2(int n, double* par, double* gr, void* ex2);
void numericalGrad2(int n, double* par, double* gr, void* ex2);
double squaredError2(int n, double* par, void* ex2);
void squaredErrorGrad2(int n, double* par, double* gr, void* ex2);
SEXP fit4plModelWithGroups_(SEXP y, SEXP Z, SEXP G);
functional_grouped* fit4plModelWithGroups(int* bin_y, double* quant_y, int y_bin, double y_mean, double* Z, int N, int* obs_ind, int* G, int n_groups);

