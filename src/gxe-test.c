#include "gxe-test.h"

double binLogLikelihood2(int n, double* par, void* ex2) {
  double sum = 0;
  dataset_grouped* data2 = ex2;
  dataset* data = data2->data;
  int* y = data->bin_y;
  double* Z = data->Z;
  int* G = data2->G;
  int n_groups = data2->n_groups;
  int* obs_ind = data->obs_ind;
  int N = data->N;
  // Parameter scaling as in drc::drm
  double* par_scale = data->par_scale;
  double b = par[0] * par_scale[0];
  double c = par[1] * par_scale[1];
  double d = par[2] * par_scale[2];
  double e = par[3] * par_scale[3];
  double prob;
  int current_group;
  for(int i = 0; i < N; i++) {
    prob = c + (d-c)/(1+exp(b*(Z[obs_ind[i]] - e)));
    current_group = G[obs_ind[i]];
    if(current_group < n_groups - 1) prob += par[4 + current_group] * par_scale[4 + current_group];
    if(prob > UPPER_TOL_4PL)
      prob = UPPER_TOL_4PL;
    else if(prob < LOWER_TOL_4PL)
      prob = LOWER_TOL_4PL;
    if(y[obs_ind[i]]) {
      sum += log(prob);
    } else {
      sum += log(1 - prob);
    }
  }
  return -sum;
}

void binLogLikelihoodGrad2(int n, double* par, double* gr, void* ex2) {
  dataset_grouped* data2 = ex2;
  dataset* data = data2->data;
  int* y = data->bin_y;
  double* Z = data->Z;
  int* G = data2->G;
  int n_groups = data2->n_groups;
  int* obs_ind = data->obs_ind;
  int N = data->N;
  double* par_scale = data->par_scale;
  double b = par[0] * par_scale[0];
  double c = par[1] * par_scale[1];
  double d = par[2] * par_scale[2];
  double e = par[3] * par_scale[3];
  double buf, p;
  int current_group;
  double current_beta;

  int n_params = 4 + n_groups - 1;
  memset(gr, 0, n_params * sizeof(double));
  for(int i = 0; i < N; i++) {
    buf = exp(b*(Z[obs_ind[i]] - e));

    current_group = G[obs_ind[i]];
    if(current_group < n_groups - 1)
      current_beta = par[4 + current_group] * par_scale[4 + current_group];
    else
      current_beta = 0;

    p = c + (d-c)/(1+buf) + current_beta;
    if(p < LOWER_TOL_4PL || p > UPPER_TOL_4PL) continue;

    if(y[obs_ind[i]]) {
      gr[0] -= (d + current_beta)*(e - Z[obs_ind[i]])/((c+current_beta)*buf + d + current_beta) + (Z[obs_ind[i]]-e)/(buf + 1);
      // Utilize asymptotic behavior, if exponential explodes
      if(!R_FINITE(buf))
        gr[1] -= 1/(c + current_beta);
      else
        gr[1] -= buf/(buf*(c+current_beta) + d + current_beta);
      gr[2] -= 1/(d + current_beta + (c+current_beta)*buf);
      gr[3] -= b*(d+current_beta)/((c+current_beta)*buf + d + current_beta) - b/(buf + 1);
      if(current_group < n_groups - 1) {
        if(!R_FINITE(buf))
          gr[4 + current_group] -= 1/(c + current_beta);
        else
          gr[4 + current_group] -= 1/(c + (d-c)/(1+buf) + current_beta);
      }
    } else {
      gr[0] -= (d+current_beta-1)*(e - Z[obs_ind[i]])/((c+current_beta-1)*buf + d + current_beta - 1) + (Z[obs_ind[i]]-e)/(buf + 1);
      // Utilize asymptotic behavior, if exponential explodes
      if(!R_FINITE(buf))
        gr[1] -= 1/(c+current_beta-1);
      else
        gr[1] -= buf/(buf*(c+current_beta-1) + d + current_beta - 1);
      gr[2] -= 1/(d + current_beta - 1 + (c+current_beta-1)*buf);
      gr[3] -= b*(d+current_beta-1)/((c+current_beta-1)*buf + d + current_beta - 1) - b/(buf + 1);
      if(current_group < n_groups - 1) {
        if(!R_FINITE(buf))
          gr[4 + current_group] -= 1/(c + current_beta - 1);
        else
          gr[4 + current_group] -= 1/(c + (d-c)/(1+buf) + current_beta - 1);
      }
    }
  }
  for(int i = 0; i < n_params; i++) gr[i] *= par_scale[i];
}

void numericalGrad2(int n, double* par, double* gr, void* ex2) {
  dataset_grouped* data2 = ex2;
  dataset* data = data2->data;

  memset(gr, 0, n * sizeof(double));
  optimfn* fn = data->fn;
  double* par_copy = (double*) Calloc(n, double);
  memcpy(par_copy, par, n * sizeof(double));
  double h = 6.055454e-06; // Cube-root of machine epsilon 2.220446e-16
  double buf;
  double* par_scale = data->par_scale;
  for(int i = 0; i < n; i++) {
    par_copy[i] = par[i] + h/par_scale[i];
    buf = fn(n, par_copy, ex2);
    par_copy[i] = par[i] - h/par_scale[i];
    gr[i] = (buf - fn(n, par_copy, ex2))/(2.0 * h);
    par_copy[i] = par[i];
  }
  Free(par_copy);
  for(int i = 0; i < n; i++) gr[i] *= par_scale[i];
}

double squaredError2(int n, double* par, void* ex2) {
  double sum = 0;
  dataset_grouped* data2 = ex2;
  dataset* data = data2->data;
  double* y = data->quant_y;
  double* Z = data->Z;
  int* G = data2->G;
  int n_groups = data2->n_groups;

  int* obs_ind = data->obs_ind;
  int N = data->N;
  double* par_scale = data->par_scale;
  double b = par[0] * par_scale[0];
  double c = par[1] * par_scale[1];
  double d = par[2] * par_scale[2];
  double e = par[3] * par_scale[3];
  double buf;
  int current_group;
  for(int i = 0; i < N; i++) {
    buf = c + (d-c)/(1+exp(b*(Z[obs_ind[i]] - e)));
    current_group = G[obs_ind[i]];
    if(current_group < n_groups - 1) buf += par[4 + current_group] * par_scale[4 + current_group];
    sum += (y[obs_ind[i]] - buf) * (y[obs_ind[i]] - buf);
  }
  return sum/N;
}

void squaredErrorGrad2(int n, double* par, double* gr, void* ex2) {
  dataset_grouped* data2 = ex2;
  dataset* data = data2->data;
  double* y = data->quant_y;
  double* Z = data->Z;
  int* G = data2->G;
  int n_groups = data2->n_groups;
  int* obs_ind = data->obs_ind;
  int N = data->N;
  double* par_scale = data->par_scale;
  double b = par[0] * par_scale[0];
  double c = par[1] * par_scale[1];
  double d = par[2] * par_scale[2];
  double e = par[3] * par_scale[3];
  double buf;
  int current_group;
  double current_beta;

  int n_params = 4 + n_groups - 1;
  memset(gr, 0, n_params * sizeof(double));

  for(int i = 0; i < N; i++) {
    buf = exp(b*(Z[obs_ind[i]] - e));
    // Utilize asymptotic behavior, if exponential explodes
    current_group = G[obs_ind[i]];
    if(current_group < n_groups - 1)
      current_beta = par[4 + current_group] * par_scale[4 + current_group];
    else
      current_beta = 0;

    if(!R_FINITE(buf)) {
      gr[1] += 2 * (c + current_beta - y[obs_ind[i]]);
      if(current_group < n_groups - 1) gr[4 + current_group] += -2 * (y[obs_ind[i]] - (c + current_beta));
    } else {
      gr[0] += -2 * (c-d) * (e-Z[obs_ind[i]]) * buf * ((c-y[obs_ind[i]]+current_beta) * buf + d - y[obs_ind[i]] + current_beta) / ((buf+1)*(buf+1)*(buf+1));
      gr[1] += 2 * buf * ((c-y[obs_ind[i]]+current_beta) * buf + d - y[obs_ind[i]] + current_beta)/((buf+1)*(buf+1));
      gr[2] += 2*((c-y[obs_ind[i]]+current_beta) * buf + d - y[obs_ind[i]] + current_beta)/((buf+1)*(buf+1));
      gr[3] += -2 * b * (c-d) * buf * ((c-y[obs_ind[i]]+current_beta) * buf + d - y[obs_ind[i]] + current_beta) /((buf+1)*(buf+1)*(buf+1));
      if(current_group < n_groups - 1) gr[4 + current_group] += -2 * (y[obs_ind[i]] - (c + (d-c)/(1+buf) + current_beta));
    }
  }
  for(int i = 0; i < n_params; i++) gr[i] *= par_scale[i]/N;
}

SEXP fit4plModelWithGroups_(SEXP y, SEXP Z, SEXP G) {
  int* bin_y = NULL;
  double* quant_y = NULL;
  int y_bin = 0;
  int N = length(y);
  double* Z2 = REAL(Z);
  int* G2 = INTEGER(G);
  int n_groups = 0;
  for(int i = 0; i < N; i++) {
    if(G2[i] > n_groups) n_groups = G2[i];
  }
  n_groups++;
  int* obs_ind = (int*) Calloc(N, int);
  if(isInteger(y)) {
    bin_y = INTEGER(y);
    y_bin = 1;
  }
  else
    quant_y = REAL(y);
  double y_mean = 0;
  for(int i = 0; i < N; i++) obs_ind[i] = i;
  if(y_bin) {
    for(int i = 0; i < N; i++) y_mean += bin_y[i];
  } else {
    for(int i = 0; i < N; i++) y_mean += quant_y[i];
  }
  y_mean /= N;
  functional_grouped* model = fit4plModelWithGroups(bin_y, quant_y, y_bin, y_mean, Z2, N, obs_ind, G2, n_groups);
  // Correctly embed 4pL model
  SEXP model_R = PROTECT(allocVector(VECSXP, 4));
  SEXP bcde_R = allocVector(REALSXP, 4);
  SET_VECTOR_ELT(model_R, 0, bcde_R);
  SET_VECTOR_ELT(model_R, 1, ScalarLogical(y_bin));
  SEXP beta_R = allocVector(REALSXP, n_groups - 1);
  SET_VECTOR_ELT(model_R, 2, beta_R);
  SET_VECTOR_ELT(model_R, 3, ScalarInteger(n_groups));
  double* model_pointer = REAL(bcde_R);
  functional* fpl = model->fpl;
  model_pointer[0] = fpl->b;
  model_pointer[1] = fpl->c;
  model_pointer[2] = fpl->d;
  model_pointer[3] = fpl->e;
  double* beta = REAL(beta_R);
  for(int i = 0; i < n_groups - 1; i++) beta[i] = (model->beta)[i];
  Free(obs_ind);
  Free(model->fpl);
  Free(model->beta);
  Free(model);
  // Assign class "4pl.grouped"
  classgets(model_R, mkString("4pl.grouped"));
  UNPROTECT(1);
  return model_R;
}

functional_grouped* fit4plModelWithGroups(int* bin_y, double* quant_y, int y_bin, double y_mean, double* Z, int N, int* obs_ind, int* G, int n_groups) {
  optimfn* fn;
  optimgr* gr;

  functional_grouped* ret = (functional_grouped*) Calloc(1, functional_grouped);
  ret->n_groups = n_groups;
  double* beta = (double*) Calloc(n_groups, double);
  ret->beta = beta;
  memset(beta, 0, n_groups * sizeof(double));

  functional* initModel = fit4plModel(bin_y, quant_y, y_bin, y_mean, Z, N, obs_ind);
  if(doubleEquals(initModel->c, initModel->d)) {
    ret->fpl = initModel;
    return ret;
  }

  functional* fpl = (functional*) Calloc(1, functional);
  fpl->y_bin = y_bin;
  ret->fpl = fpl;

  // Get group means
  // Assume, G is coded as 0, ..., n_groups - 1
  int* group_counts = (int*) Calloc(n_groups, int);
  memset(group_counts, 0, n_groups * sizeof(int));
  int current_group;
  for(int i = 0; i < N; i++) {
    current_group = G[obs_ind[i]];
    group_counts[current_group]++;
    if(y_bin)
      beta[current_group] += bin_y[obs_ind[i]];
    else
      beta[current_group] += quant_y[obs_ind[i]];
  }
  for(int j = 0; j < n_groups; j++) beta[j] /= group_counts[j];

  // Calibrate the location paramters using the overall mean
  // and the last group
  for(int j = 0; j < n_groups; j++) beta[j] -= y_mean;
  for(int j = 0; j < n_groups - 1; j++) beta[j] -= beta[n_groups - 1];

  if(y_bin) {
    fn = &binLogLikelihood2;
    gr = &binLogLikelihoodGrad2;
    // gr = &numericalGrad2;
  } else {
    fn = &squaredError2;
    gr = &squaredErrorGrad2;
  }

  int n_params = 4 + n_groups - 1;
  // Scale parameters
  double* Pars = (double*) Calloc(n_params, double);
  Pars[0] = initModel->b;
  Pars[1] = initModel->c;
  Pars[2] = initModel->d;
  Pars[3] = initModel->e;
  for(int j = 0; j < n_groups - 1; j++) Pars[j + 4] = beta[j];
  double* par_scale = (double*) Calloc(n_params, double);
  memcpy(par_scale, Pars, n_params * sizeof(double));
  for(int i = 0; i < n_params; i++) {
    par_scale[i] = fabs(par_scale[i]);
    if(par_scale[i] < 1e-4) par_scale[i] = 1;
    Pars[i] /= par_scale[i];
  }
  double min_val = 0;
  int maxit = 500; // Standard value 100
  int trace = 0;
  int* mask = (int*) Calloc(n_params, int);
	for (int i = 0; i < n_params; i++) mask[i] = 1;
  double abstol = R_NegInf;
  double reltol = 1e-7; // Standard value 1.490116e-08
  int nREPORT = 10;
  dataset* ex = (dataset*) Calloc(1, dataset);
  ex->bin_y = bin_y;
  ex->quant_y = quant_y;
  ex->Z = Z;
  ex->obs_ind = obs_ind;
  ex->N = N;
  ex->par_scale = par_scale;
  ex->fn = fn; // For alternative numerical derivatives
  dataset_grouped* ex2 = (dataset_grouped*) Calloc(1, dataset_grouped);
  ex2->data = ex;
  ex2->G = G;
  ex2->n_groups = n_groups;
  int fncount = 0;
  int grcount = 0;
  int fail = 0;

  // Catch not finite error:
  double f = fn(n_params, Pars, ex2);
  if (!R_FINITE(f)) {
    fpl->c = y_mean;
    fpl->d = y_mean;
    fpl->b = 0;
    fpl->e = 0;
  } else {
    vmmin(n_params, Pars, &min_val, fn, gr, maxit, trace, mask, abstol, reltol, nREPORT, ex2, &fncount, &grcount, &fail);

    // Revert parameter scaling
    fpl->b = Pars[0] * par_scale[0];
    fpl->c = Pars[1] * par_scale[1];
    fpl->d = Pars[2] * par_scale[2];
    fpl->e = Pars[3] * par_scale[3];

    for(int i = 0; i < n_groups - 1; i++) beta[i] = Pars[4 + i] * par_scale[4 + i];
  }
  Free(Pars);
  Free(par_scale);
  Free(mask);
  Free(ex);
  Free(ex2);
  Free(initModel);
  Free(group_counts);
  return ret;
}

