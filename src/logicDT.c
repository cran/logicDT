#include "logicDT.h"

// extern float zeit = 0.0;

double calcDev(double* predictions, int* y, int N) {
  int i;
  double sum = 0.0;
  double buffer;

  for (i = 0; i < N; i++) {
    buffer = log(y[i] * predictions[i] + (1-y[i]) * (1-predictions[i]));
    if (isfinite(buffer)) {
      sum += buffer;
    } else {
      sum += log(0.001);
    }
  }

  return -2.0 * sum;
}

// Normalized cross entropy (deviance)
double calcNCE(double* predictions, int* y, int N) {
  int i;
  double sum = 0.0;
  double buffer;
  double p = 0.0;
  double lower_tol = 0.001;
  double upper_tol = 0.999;

  for (i = 0; i < N; i++) {
    buffer = log(y[i] * predictions[i] + (1-y[i]) * (1-predictions[i]));
    if (isfinite(buffer)) {
      sum += buffer;
    } else {
      sum += log(lower_tol);
    }
    p += y[i];
  }

  p /= N;
  if(p > upper_tol) p = upper_tol;
  else if(p < lower_tol) p = lower_tol;
  double denominator = p * log(p) + (1-p) * log(1-p);

  return (sum/N)/denominator;
}

double calcBrier(double* predictions, int* y, int N) {
  int i;
  double sum = 0.0;

  for (i = 0; i < N; i++) {
    sum += (predictions[i] - y[i]) * (predictions[i] - y[i]);
  }

  return sum/N;
}

double calcMis(int* predictions, int* y, int N) {
  int i;
  int sum = 0;

  for (i = 0; i < N; i++) {
    sum += (predictions[i] != y[i]);
  }

  return (double) sum/N;
}

double calcLikelihood(double* predictions, int* y, int N) {
  int i;
  double prod = 1.0;

  for (i = 0; i < N; i++) {
    prod *= pow(y[i] * predictions[i] + (1-y[i]) * (1-predictions[i]), 1.0/N);
  }

  return prod;
}

SEXP calcAUC_(SEXP probs, SEXP y, SEXP y_sorted_raw) {
  int y_sorted = asLogical(y_sorted_raw);
  double auc;
  if(y_sorted)
    auc = calcAUCSorted(REAL(probs), INTEGER(y), length(y));
  else
    auc = calcAUCUnsorted(REAL(probs), INTEGER(y), length(y));
  return ScalarReal(auc);
}

double calcAUCSorted(double* predictions, int* y, int N) {
  // Requires sorted y (and predictions)

  int N_0;
  for(N_0 = 0; N_0 < N; N_0++) {
    if(y[N_0])
      break;
  }
  int N_1 = N - N_0;

  double* preds_0 = (double*) Calloc(N_0, double);
  double* preds_1 = (double*) Calloc(N_1, double);

  memcpy(preds_0, predictions, N_0 * sizeof(double));
  memcpy(preds_1, predictions + N_0, N_1 * sizeof(double));

  qsort(preds_0, N_0, sizeof(double), cmp_double);
  qsort(preds_1, N_1, sizeof(double), cmp_double);

  double ranksum = 0.0;
  int i, j;
  int ind_buffer_old = 0;
  int ind_buffer_new = 0;
  double rank_buffer = 0.0;
  for(i = 0; i < N_1; i++) {
    if(i > 0) {
      if(doubleEquals(preds_1[i], preds_1[i-1])) {
        ranksum += ind_buffer_old + rank_buffer;
        continue;
      } else {
        ranksum += ind_buffer_new;
      }
    }
    rank_buffer = 0.0;
    for(j = ind_buffer_new; j < N_0; j++) {
      if(doubleEquals(preds_1[i], preds_0[j])) {
        ranksum += 0.5;
        rank_buffer += 0.5;
      }
      else if(preds_1[i] > preds_0[j]) {
        ranksum++;
        rank_buffer++;
      }
      else
        break;
    }
    ind_buffer_old = ind_buffer_new;
    ind_buffer_new = j;
  }

  Free(preds_0);
  Free(preds_1);

  return ranksum/(N_0 * N_1);
}

double calcAUCUnsorted(double* predictions, int* y, int N) {
  y_prob_pair_t* y_probs = (y_prob_pair_t*) Calloc(N, y_prob_pair_t);

  int i, j;
  int N_0 = 0;
  for(i = 0; i < N; i++) {
    if(!y[i])
      N_0++;
    y_probs[i].y = y[i];
    y_probs[i].prob = predictions[i];
  }

  qsort(y_probs, N, sizeof(y_prob_pair_t), cmp_y_probs_int);

  int N_1 = N - N_0;

  y_prob_pair_t* preds_0 = y_probs;
  y_prob_pair_t* preds_1 = y_probs + N_0;

  qsort(preds_0, N_0, sizeof(y_prob_pair_t), cmp_y_probs_double);
  qsort(preds_1, N_1, sizeof(y_prob_pair_t), cmp_y_probs_double);

  double ranksum = 0.0;
  int ind_buffer_old = 0;
  int ind_buffer_new = 0;
  double rank_buffer = 0.0;
  for(i = 0; i < N_1; i++) {
    if(i > 0) {
      if(doubleEquals(preds_1[i].prob, preds_1[i-1].prob)) {
        ranksum += ind_buffer_old + rank_buffer;
        continue;
      } else {
        ranksum += ind_buffer_new;
      }
    }
    rank_buffer = 0.0;
    for(j = ind_buffer_new; j < N_0; j++) {
      if(doubleEquals(preds_1[i].prob, preds_0[j].prob)) {
        ranksum += 0.5;
        rank_buffer += 0.5;
      }
      else if(preds_1[i].prob > preds_0[j].prob) {
        ranksum++;
        rank_buffer++;
      }
      else
        break;
    }
    ind_buffer_old = ind_buffer_new;
    ind_buffer_new = j;
  }

  Free(y_probs);

  return ranksum/(N_0 * N_1);
}

/*
double calcAUC2(double* predictions, int* y, int N) {

  double sum = 0.0;
  double o;
  int counter = 0;

  for (int i = 0; i < N; i++) {
    if(!y[i])
      continue;
    for(int j = 0; j < N; j++) {
      if(y[j])
        continue;
      o = predictions[i] - predictions[j];
      if(doubleEquals(o, 0.0))
        sum += 0.5;
      else
        sum += (o > 0);
      counter++;
    }
  }

  if(!counter)
    return 0.0;
  else
    return sum/counter;
}*/

double calcMSE(double* predictions, double* y, int N) {
  int i;
  double sum = 0.0;

  for (i = 0; i < N; i++) {
    sum += (predictions[i] - y[i]) * (predictions[i] - y[i]);
  }

  return sum/N;
}

SEXP vim_permutation_(SEXP ensemble, SEXP X_val_raw, SEXP y_val_raw, SEXP Z_val_raw, SEXP permutation_raw, SEXP disj_raw, SEXP real_n_conj_raw, SEXP scoring_rule_raw, SEXP y_bin_raw, SEXP leaves_raw) {
  int* permutation;

  int i, j, var_index;
  int N;
  int n_folds = length(ensemble);
  int leaves = asInteger(leaves_raw);

  double* Z_temp = NULL;

  int n_conj = nrows(disj_raw);
  int n_vars = ncols(disj_raw);

  int real_n_conj = asInteger(real_n_conj_raw);

  int scoring_rule = asInteger(scoring_rule_raw);
  int y_bin = asLogical(y_bin_raw);

  int pred_type;
  if(scoring_rule != 2 || !y_bin) {
    pred_type = 0;
  } else {
    pred_type = 1;
  }

  int* dm_val;
  int* dm_val_perm;
  node* tree;

  SEXP current_pet;
  double orig_score;
  double perm_score;

  SEXP vim_raw = PROTECT(allocMatrix(REALSXP, real_n_conj, n_folds));
  double* vim = REAL(vim_raw);

  SEXP orig_scores_raw = PROTECT(allocMatrix(REALSXP, real_n_conj, n_folds));
  SEXP perm_scores_raw = PROTECT(allocMatrix(REALSXP, real_n_conj, n_folds));
  double* orig_scores = REAL(orig_scores_raw);
  double* perm_scores = REAL(perm_scores_raw);

  pet_preds_t* orig_preds_buffer;
  pet_preds_t* perm_preds_buffer;

  for (i = 0; i < n_folds; i++) {
    permutation = INTEGER(VECTOR_ELT(permutation_raw, i));
    N = length(VECTOR_ELT(permutation_raw, i));

    dm_val = getDesignMatrixIntern(INTEGER(VECTOR_ELT(X_val_raw, i)), N, INTEGER(disj_raw), n_conj, n_vars, real_n_conj);
    dm_val_perm = (int*) Calloc(N * real_n_conj, int);

    current_pet = VECTOR_ELT(ensemble, i);
    rebuild_tree(current_pet);
    tree = (node*) R_ExternalPtrAddr(VECTOR_ELT(current_pet, 5));
    if(!isNull(Z_val_raw))
      Z_temp = REAL(VECTOR_ELT(Z_val_raw, i));
    orig_preds_buffer = predictIntern(tree, dm_val, Z_temp, N, pred_type, leaves);

    if(!y_bin) {
      orig_score = calcMSE(orig_preds_buffer->prob_preds, REAL(VECTOR_ELT(y_val_raw, i)), N);
    } else {
      if(scoring_rule == 0) {
        orig_score = calcDev(orig_preds_buffer->prob_preds, INTEGER(VECTOR_ELT(y_val_raw, i)), N);
      } else if (scoring_rule == 1) {
        orig_score = calcBrier(orig_preds_buffer->prob_preds, INTEGER(VECTOR_ELT(y_val_raw, i)), N);
      } else if (scoring_rule == 2) {
        orig_score = calcMis(orig_preds_buffer->class_preds, INTEGER(VECTOR_ELT(y_val_raw, i)), N);
      } else if (scoring_rule == 3) {
        orig_score = -calcLikelihood(orig_preds_buffer->prob_preds, INTEGER(VECTOR_ELT(y_val_raw, i)), N);
      } else if (scoring_rule == 5) {
        orig_score = calcNCE(orig_preds_buffer->prob_preds, INTEGER(VECTOR_ELT(y_val_raw, i)), N);
      } else {
        orig_score = 1 - calcAUCUnsorted(orig_preds_buffer->prob_preds, INTEGER(VECTOR_ELT(y_val_raw, i)), N);
      }
    }

    for (var_index = 0; var_index < real_n_conj; var_index++) {
      memcpy(dm_val_perm, dm_val, N * real_n_conj * sizeof(int));

      for (j = 0; j < N; j++) {
        dm_val_perm[var_index*N + j] = dm_val[var_index*N + permutation[j] - 1];
      }

      perm_preds_buffer = predictIntern(tree, dm_val_perm, Z_temp, N, pred_type, leaves);

      if(!y_bin) {
        perm_score = calcMSE(perm_preds_buffer->prob_preds, REAL(VECTOR_ELT(y_val_raw, i)), N);
      } else {
        if(scoring_rule == 0) {
          perm_score = calcDev(perm_preds_buffer->prob_preds, INTEGER(VECTOR_ELT(y_val_raw, i)), N);
        } else if (scoring_rule == 1) {
          perm_score = calcBrier(perm_preds_buffer->prob_preds, INTEGER(VECTOR_ELT(y_val_raw, i)), N);
        } else if (scoring_rule == 2) {
          perm_score = calcMis(perm_preds_buffer->class_preds, INTEGER(VECTOR_ELT(y_val_raw, i)), N);
        } else if (scoring_rule == 3) {
          perm_score = -calcLikelihood(perm_preds_buffer->prob_preds, INTEGER(VECTOR_ELT(y_val_raw, i)), N);
        } else if (scoring_rule == 5) {
          perm_score = calcNCE(perm_preds_buffer->prob_preds, INTEGER(VECTOR_ELT(y_val_raw, i)), N);
        } else {
          perm_score = 1 - calcAUCUnsorted(perm_preds_buffer->prob_preds, INTEGER(VECTOR_ELT(y_val_raw, i)), N);
        }
      }
      vim[i*real_n_conj + var_index] = perm_score - orig_score;
      orig_scores[i*real_n_conj + var_index] = orig_score;
      perm_scores[i*real_n_conj + var_index] = perm_score;

      Free(perm_preds_buffer->prob_preds);
      if(pred_type)
        Free(perm_preds_buffer->class_preds);
      Free(perm_preds_buffer);
    }

    Free(dm_val);
    Free(dm_val_perm);
    Free(orig_preds_buffer->prob_preds);
    if(pred_type)
      Free(orig_preds_buffer->class_preds);
    Free(orig_preds_buffer);
  }

  SEXP ret = PROTECT(allocVector(VECSXP, 3));
  SET_VECTOR_ELT(ret, 0, vim_raw);
  SET_VECTOR_ELT(ret, 1, orig_scores_raw);
  SET_VECTOR_ELT(ret, 2, perm_scores_raw);

  UNPROTECT(4);
  return ret;
}

SEXP predictEnsemble_(SEXP ensemble, SEXP X_raw, SEXP Z_raw, SEXP type_raw, SEXP leaves_raw) {
  int n_folds = length(ensemble);
  int i, j;

  int N = nrows(X_raw);
  int type = asLogical(type_raw);

  SEXP prob_preds_raw = PROTECT(allocVector(REALSXP, N));
  double* prob_preds = REAL(prob_preds_raw);
  memset(prob_preds, 0, N * sizeof(double));

  SEXP prob_type = PROTECT(ScalarLogical(0));

  double* current_pred;

  for(i = 0; i < n_folds; i++) {
    current_pred = REAL(PROTECT(predict_(VECTOR_ELT(ensemble, i), X_raw, Z_raw, prob_type, leaves_raw)));
    for(j = 0; j < N; j++) {
      prob_preds[j] += current_pred[j];
    }
    UNPROTECT(1);
  }
  for(j = 0; j < N; j++) {
    prob_preds[j] /= n_folds;
  }

  if(!type) {
    UNPROTECT(2);
    return prob_preds_raw;
  } else {
    SEXP class_preds_raw = PROTECT(allocVector(INTSXP, N));
    int* class_preds = INTEGER(class_preds_raw);
    for(j = 0; j < N; j++) {
      if (prob_preds[j] > 0.5)
        class_preds[j] = 1;
      else
        class_preds[j] = 0;
    }
    UNPROTECT(3);
    return class_preds_raw;
  }
}

static void _finalizer(SEXP tree) {
  if (R_ExternalPtrAddr(tree) == NULL)
    return;
  // Rprintf("finalizing\n");
  node* ptr = (node*) R_ExternalPtrAddr(tree);
  tree_destroy(ptr);
  R_ClearExternalPtr(tree);
}

SEXP C_PET_TO_R_PET(pet_t* pet, int N) {
  SEXP train_preds_R = PROTECT(allocVector(REALSXP, N));
  memcpy(REAL(train_preds_R), pet->train_preds, N * sizeof(double));
  SEXP splits_R = PROTECT(allocVector(INTSXP, pet->number_of_nodes));
  memcpy(INTEGER(splits_R), pet->splits, pet->number_of_nodes * sizeof(int));
  SEXP splits_bin_or_cont_R = PROTECT(allocVector(INTSXP, pet->number_of_nodes));
  memcpy(INTEGER(splits_bin_or_cont_R), pet->splits_bin_or_cont, pet->number_of_nodes * sizeof(int));
  SEXP split_points_R = PROTECT(allocVector(REALSXP, pet->number_of_nodes));
  memcpy(REAL(split_points_R), pet->split_points, pet->number_of_nodes * sizeof(double));
  SEXP preds_R = PROTECT(allocVector(REALSXP, pet->number_of_nodes));
  memcpy(REAL(preds_R), pet->preds, pet->number_of_nodes * sizeof(double));

  SEXP bin_tree = PROTECT(R_MakeExternalPtr(pet->tree, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(bin_tree, _finalizer, TRUE);

  SEXP pet_R = PROTECT(allocVector(VECSXP, 8));
  SET_VECTOR_ELT(pet_R, 0, splits_R);
  SET_VECTOR_ELT(pet_R, 1, splits_bin_or_cont_R);
  SET_VECTOR_ELT(pet_R, 2, split_points_R);
  SET_VECTOR_ELT(pet_R, 3, preds_R);
  SET_VECTOR_ELT(pet_R, 4, train_preds_R);
  SET_VECTOR_ELT(pet_R, 5, bin_tree);
  if(pet->model_list != NULL) {
    SEXP model_list_R = PROTECT(allocVector(VECSXP, pet->number_of_nodes));
    double* bcde;
    SEXP func_pred_R, bcde_R;
    functional* current_func;
    for(int i = 0; i < pet->number_of_nodes; i++) {
      if((pet->splits)[i] == 0) {
        func_pred_R = allocVector(VECSXP, 2);
        SET_VECTOR_ELT(model_list_R, i, func_pred_R);
        bcde_R = allocVector(REALSXP, 4);
        SET_VECTOR_ELT(func_pred_R, 0, bcde_R);
        SET_VECTOR_ELT(func_pred_R, 1, ScalarLogical(pet->y_bin));
        // Assign class "4pl"
        classgets(func_pred_R, mkString("4pl"));
        bcde = REAL(bcde_R);
        current_func = (pet->model_list)[i];
        bcde[0] = current_func->b;
        bcde[1] = current_func->c;
        bcde[2] = current_func->d;
        bcde[3] = current_func->e;
      } else {
        SET_VECTOR_ELT(model_list_R, i, R_NilValue);
      }
    }
    SET_VECTOR_ELT(pet_R, 6, model_list_R);
    UNPROTECT(1);
  } else {
    SET_VECTOR_ELT(pet_R, 6, R_NilValue);
  }
  SET_VECTOR_ELT(pet_R, 7, ScalarLogical(pet->y_bin));

  pet_destroy(pet, 0);

  UNPROTECT(7);
  return pet_R;
}

SEXP fitPETs_(SEXP X_train_raw, SEXP y_train_raw, SEXP X_val_raw, SEXP y_val_raw, SEXP Z_train_raw, SEXP Z_val_raw, SEXP use_validation_raw, SEXP y_bin_raw, SEXP nodesize_raw, SEXP cp_raw, SEXP smoothing_raw, SEXP mtry_raw, SEXP covariable_mode_raw, SEXP disj_raw, SEXP real_n_conj_raw, SEXP scoring_rule_raw, SEXP return_full_model_raw) {
  int use_validation = asLogical(use_validation_raw);
  int y_bin = asLogical(y_bin_raw);
  int scoring_rule = INTEGER(scoring_rule_raw)[0];
  int return_full_model = asLogical(return_full_model_raw);
  int nodesize = INTEGER(nodesize_raw)[0];
  double cp = REAL(cp_raw)[0];
  int smoothing = INTEGER(smoothing_raw)[0];
  int mtry = asInteger(mtry_raw);
  int covariable_mode = asInteger(covariable_mode_raw);
  int* disj = INTEGER(disj_raw);
  int n_conj = nrows(disj_raw);
  int n_vars = ncols(disj_raw);
  int real_n_conj = asInteger(real_n_conj_raw);
  int n_folds = length(X_train_raw);

  SEXP return_obj;
  if(return_full_model) {
    return_obj = PROTECT(allocVector(VECSXP, n_folds));
  }

  pet_ensemble_t* pets_intern = fitPETsIntern(X_train_raw, y_train_raw, X_val_raw, y_val_raw, Z_train_raw, Z_val_raw, use_validation, y_bin, nodesize, cp, smoothing, mtry, covariable_mode, disj, n_conj, n_vars, real_n_conj, scoring_rule, return_full_model);
  pet_t** petsss = pets_intern->pets;

  if(return_full_model) {
    for(int i = 0; i < n_folds; i++) {
      SET_VECTOR_ELT(return_obj, i, C_PET_TO_R_PET(petsss[i], length(VECTOR_ELT(y_train_raw, i))));
    }
    Free(petsss);
  } else {
    return_obj = PROTECT(ScalarReal(pets_intern->score));
  }

  Free(pets_intern);
  UNPROTECT(1);
  return return_obj;
}

pet_ensemble_t* fitPETsIntern(SEXP X_train_raw, SEXP y_train_raw, SEXP X_val_raw, SEXP y_val_raw, SEXP Z_train_raw, SEXP Z_val_raw, int use_validation, int y_bin, int nodesize, double cp, int smoothing, int mtry, int covariable_mode, int* disj, int n_conj, int n_vars, int real_n_conj, int scoring_rule, int return_full_model) {
  int n_folds = length(X_train_raw);
  int i;
  int N;
  int N_val;
  SEXP current_X_train, current_X_val;
  int* X_train;
  int* bin_y_train = NULL;
  double* quant_y_train = NULL;
  int* X_val;
  int* bin_y_val = NULL;
  double* quant_y_val = NULL;

  double* Z_train = NULL;
  double* Z_val = NULL;
  int pZ = 0;
  int use_Z = !isNull(Z_train_raw);
  if(use_Z)
    pZ = ncols(VECTOR_ELT(Z_train_raw, 0));

  pet_t* current_pet;
  pet_preds_t* predictions_raw;
  double* predictions;
  //SEXP X_train, y_train, X_val, y_val;
  //SEXP current_pet, predictions;

  int pred_type = 0;

  double scores = 0.0;

  pet_ensemble_t* return_obj = (pet_ensemble_t*) Calloc(1, pet_ensemble_t);
  return_obj->n_pets = n_folds;

  if(return_full_model)
    return_obj->pets = (pet_t**) Calloc(n_folds, pet_t*);

  for (i = 0; i < n_folds; i++) {
    /*X_train = PROTECT(getDesignMatrix_(VECTOR_ELT(X_train_raw, i), disj_raw, real_n_conj_raw));*/
    current_X_train = VECTOR_ELT(X_train_raw, i);
    N = nrows(current_X_train);
    X_train = getDesignMatrixIntern(INTEGER(current_X_train), N, disj, n_conj, n_vars, real_n_conj);
    if(y_bin)
      bin_y_train = INTEGER(VECTOR_ELT(y_train_raw, i));
    else
      quant_y_train = REAL(VECTOR_ELT(y_train_raw, i));

    if(use_Z)
      Z_train = REAL(VECTOR_ELT(Z_train_raw, i));

    current_pet = fitPETIntern(X_train, bin_y_train, quant_y_train, y_bin, Z_train, N, real_n_conj, pZ, nodesize, cp, smoothing, mtry, covariable_mode);
    Free(X_train);

    if(use_validation) {
      current_X_val = VECTOR_ELT(X_val_raw, i);
      N_val = nrows(current_X_val);
      X_val = getDesignMatrixIntern(INTEGER(current_X_val), N_val, disj, n_conj, n_vars, real_n_conj);
      if(use_Z)
        Z_val = REAL(VECTOR_ELT(Z_val_raw, i));
      predictions_raw = predictIntern(current_pet->tree, X_val, Z_val, N_val, pred_type, 1);
      predictions = predictions_raw->prob_preds;
      Free(predictions_raw);
      Free(X_val);
    } else {
      N_val = N;
      predictions = current_pet->train_preds;
    }
    if(y_bin) {
      bin_y_val = INTEGER(VECTOR_ELT(y_val_raw, i));
      if (!scoring_rule)
        scores += calcDev(predictions, bin_y_val, N_val);
      else if (scoring_rule == 1)
        scores += calcBrier(predictions, bin_y_val, N_val);
      else if (scoring_rule == 5)
        scores += calcNCE(predictions, bin_y_val, N_val);
      else
        scores += 1 - calcAUCSorted(predictions, bin_y_val, N_val);
    }
    else {
      quant_y_val = REAL(VECTOR_ELT(y_val_raw, i));
      scores += calcMSE(predictions, quant_y_val, N_val);
    }

    if(use_validation)
      Free(predictions);

    if(return_full_model) {
      (return_obj->pets)[i] = current_pet;
    } else {
      pet_destroy(current_pet, 1);
    }
  }

  return_obj->score = scores/n_folds;
  return return_obj;
}

SEXP fitPET_(SEXP X_raw, SEXP y_raw, SEXP Z_raw, SEXP nodesize_raw, SEXP cp_raw, SEXP smoothing_raw, SEXP mtry_raw, SEXP covariable_mode_raw) {
  int* X = INTEGER(X_raw);

  int* bin_y = NULL;
  double* quant_y = NULL;
  if(isInteger(y_raw))
    bin_y = INTEGER(y_raw);
  else
    quant_y = REAL(y_raw);

  double* Z = NULL;
  int pZ = 0;
  if(!isNull(Z_raw)) {
    Z = REAL(Z_raw);
    pZ = ncols(Z_raw);
  }
  int nodesize = INTEGER(nodesize_raw)[0];
  double cp = REAL(cp_raw)[0];
  int smoothing = INTEGER(smoothing_raw)[0];
  int mtry = asInteger(mtry_raw);
  int covariable_mode = asInteger(covariable_mode_raw);

  int N = length(y_raw);
  int p = ncols(X_raw);
  pet_t* pet = fitPETIntern(X, bin_y, quant_y, isInteger(y_raw), Z, N, p, pZ, nodesize, cp, smoothing, mtry, covariable_mode);

  SEXP pet_R = C_PET_TO_R_PET(pet, N);
  return pet_R;
}

pet_t* fitPETIntern(int* X, int* bin_y, double* quant_y, int y_bin, double* Z, int N, int p, int pZ, int nodesize, double cp, int smoothing, int mtry, int covariable_mode) {
  int i, j;

//  linked_list *splits_preds = malloc(sizeof(linked_list));
  linked_list *splits_preds = (linked_list*) Calloc(1, linked_list);
  splits_preds->next = NULL;
  linked_list *current_split_pred = splits_preds;

  /*SEXP train_preds_R = PROTECT(allocVector(REALSXP, N));
  double* train_preds = REAL(train_preds_R);*/
  double* train_preds = (double*) Calloc(N, double);

  logic_stack_t *stack = stack_new();

  node* tree;
  tree = (node*) Calloc(1, node);
  tree->left = NULL;
  tree->right = NULL;
  tree->obs_ind = (int*) Calloc(N, int);
  for(i = 0; i < N; i++) {
    (tree->obs_ind)[i] = i;
  }
  tree->N_k = N;
  tree->func_pred = NULL;

  stack_push(stack, tree);

  int best_index = -1;
  double max_imp_decrease;
  double imp_decrease;
  int bin_or_cont = 0;
  double best_split_point = 0;

  int N_k, N_L, N_R, N_k_1;
  double p_k;

  int N_L_1, N_R_1;

  double N_k_sum, N_L_sum, N_R_sum;
  double N_k_sum_2, N_L_sum_2, N_R_sum_2;

  int best_N_L;

  int buffer_L, buffer_R;

  int number_of_nodes = 0;

  int* obs_ind;

  node* left_child;
  node* right_child;

  double current_leaf_pred;

  // mtry coding: -1: none, 0: sqrt(p), 1-p: 1-p
  if(mtry == 0)
    mtry = (int) sqrt(p);
  else if(mtry >= p)
    mtry = -1;
  int n_split_cands = p;
  int m, which_var;
  int* available_vars_raw = NULL;
  int* available_vars = NULL;
  if(mtry > 0) {
    available_vars_raw = (int*) Calloc(p, int);
    for(i = 0; i < p; i++)
      available_vars_raw[i] = i;
    available_vars = (int*) Calloc(p, int);
    GetRNGstate();
    n_split_cands = mtry;
  }

  while(stack->top != NULL) {
    node* knot = stack_pop(stack);
    //node* knot = queue_pop(stack);

    number_of_nodes++;

    best_index = -1;
    max_imp_decrease = cp;
    bin_or_cont = 0;
    best_split_point = 0;

    N_k = knot->N_k;
    N_k_1 = 0;
    obs_ind = knot->obs_ind;

    N_k_sum = 0;
    N_k_sum_2 = 0;

    p_k = (double) N_k/N;

    if(y_bin) {
      for (i = 0; i < N_k; i++) {
        if (bin_y[obs_ind[i]] == 1)
          N_k_1++;
      }
    } else {
      for (i = 0; i < N_k; i++) {
        N_k_sum += quant_y[obs_ind[i]];
        N_k_sum_2 += quant_y[obs_ind[i]] * quant_y[obs_ind[i]];
      }
    }

    // Stopping criterion #1
    if (y_bin && (N_k_1 > N_k - nodesize || N_k_1 < nodesize)) {
    // if (N_k <= nodesize || (y_bin && (N_k_1 == N_k || N_k_1 == 0))) { // nodesize as in ranger or randomForest
    /* ranger and randomForest aim at controlling the minimum node size before splitting
       We aim at controllin the minimum node size after splitting such that
       each terminal node holds sufficient examples for an appropriate risk prediction */
      if(y_bin)
        current_leaf_pred = calcLeafProb(N_k_1, N_k, smoothing);
      else
        current_leaf_pred = N_k_sum/N_k;
      make_leaf(knot, current_leaf_pred, train_preds);
      current_split_pred = set_values_and_next(current_split_pred, 0, 0, 0, current_leaf_pred);
      // preds[number_of_nodes - 1] = (double) N_k_1/N_k;
      continue;
    }

    // Evaluate splits
    if(mtry > 0)
      memcpy(available_vars, available_vars_raw, p * sizeof(int));

    for (m = 0; m < n_split_cands; m++) {
      if(mtry > 0) {
        which_var = unif_rand() * (p - m);
        j = available_vars[which_var];
        for(i = which_var; i < p - m - 1; i++) {
          available_vars[i] = available_vars[i+1];
        }
      } else j = m;

      N_L = 0;
      N_L_1 = 0;

      N_L_sum = 0;
      N_L_sum_2 = 0;

      for (i = 0; i < N_k; i++) {
        if (X[j*N + obs_ind[i]] == 0) {
          N_L++;
          if (y_bin && bin_y[obs_ind[i]] == 1)
            N_L_1++;
          if (!y_bin) {
            N_L_sum += quant_y[obs_ind[i]];
            N_L_sum_2 += quant_y[obs_ind[i]] * quant_y[obs_ind[i]];
          }
        }
      }

      N_R = N_k - N_L;
      N_R_1 = N_k_1 - N_L_1;

      N_R_sum = N_k_sum - N_L_sum;
      N_R_sum_2 = N_k_sum_2 - N_L_sum_2;

      // Stopping criterion #2
      if (N_L < nodesize || N_R < nodesize) { // NOT as in ranger or randomForest
        continue;
      }

      if(y_bin)
        imp_decrease = gini_decrease((double) N_k_1/N_k, (double) N_L/N_k, (double) N_L_1/N_L, (double) N_R_1/N_R);
      else
        imp_decrease = mse_decrease(N_k, N_L, N_R, N_k_sum, N_L_sum, N_R_sum, N_k_sum_2, N_L_sum_2, N_R_sum_2);

      if(p_k * imp_decrease > max_imp_decrease) {
      // if(imp_decrease > max_imp_decrease) {
        max_imp_decrease = p_k * imp_decrease;
        best_index = j;
        best_N_L = N_L;
      }
    }

    y_Z_pair_t* Z_sorted;
    if(pZ > 0 && covariable_mode == 1)
      Z_sorted = (y_Z_pair_t*) Calloc(N_k, y_Z_pair_t);

    for (j = 0; j < pZ && covariable_mode == 1; j++) {
      N_L = 0;
      N_L_1 = 0;

      N_L_sum = 0;
      N_L_sum_2 = 0;

      for(i = 0; i < N_k; i++) {
        (Z_sorted[i]).Z = Z[j*N + obs_ind[i]];
        if(y_bin)
          (Z_sorted[i]).bin_y = bin_y[obs_ind[i]];
        else
          (Z_sorted[i]).quant_y = quant_y[obs_ind[i]];
      }
      qsort(Z_sorted, N_k, sizeof(y_Z_pair_t), cmp_y_Z_pair);

      for(i = 0; i < N_k - 1; i++) {
        N_L++;
        if (y_bin && Z_sorted[i].bin_y)
          N_L_1++;
        if (!y_bin) {
          N_L_sum += Z_sorted[i].quant_y;
          N_L_sum_2 += Z_sorted[i].quant_y * Z_sorted[i].quant_y;
        }

        // if (doubleEquals(Z_sorted[i+1].Z, Z_sorted[i].Z))
        if (Z_sorted[i+1].Z == Z_sorted[i].Z)
          continue;

        N_R = N_k - N_L;
        N_R_1 = N_k_1 - N_L_1;

        N_R_sum = N_k_sum - N_L_sum;
        N_R_sum_2 = N_k_sum_2 - N_L_sum_2;

        // Stopping criterion #2
        if (N_L < nodesize || N_R < nodesize) { // NOT as in ranger or randomForest
          continue;
        }

        if(y_bin)
          imp_decrease = gini_decrease((double) N_k_1/N_k, (double) N_L/N_k, (double) N_L_1/N_L, (double) N_R_1/N_R);
        else
          imp_decrease = mse_decrease(N_k, N_L, N_R, N_k_sum, N_L_sum, N_R_sum, N_k_sum_2, N_L_sum_2, N_R_sum_2);

        if(p_k * imp_decrease > max_imp_decrease) {
        // if(imp_decrease > max_imp_decrease) {
          max_imp_decrease = p_k * imp_decrease;
          best_index = j;
          bin_or_cont = 1;
          best_split_point = (Z_sorted[i].Z + Z_sorted[i+1].Z)/2;
          best_N_L = N_L;
        }
      }
    }

    if(pZ > 0 && covariable_mode == 1)
      Free(Z_sorted);

    // Stopping criterion #3
    if (best_index == -1) {
      if(y_bin)
        current_leaf_pred = calcLeafProb(N_k_1, N_k, smoothing);
      else
        current_leaf_pred = N_k_sum/N_k;
      make_leaf(knot, current_leaf_pred, train_preds);
      current_split_pred = set_values_and_next(current_split_pred, 0, 0, 0, current_leaf_pred);
      continue;
    }

    left_child = (node*) Calloc(1, node);
    right_child = (node*) Calloc(1, node);
    left_child->left = NULL;
    right_child->left = NULL;
    left_child->right = NULL;
    right_child->right = NULL;
    left_child->N_k = best_N_L;
    right_child->N_k = N_k - best_N_L;
    left_child->obs_ind = (int*) Calloc(best_N_L, int);
    right_child->obs_ind = (int*) Calloc((N_k - best_N_L), int);
    left_child->func_pred = NULL;
    right_child->func_pred = NULL;

    buffer_L = 0;
    buffer_R = 0;

    for (i = 0; i < N_k; i++) {
      if(!bin_or_cont) {
        if (X[best_index*N + obs_ind[i]] == 0)
          (left_child->obs_ind)[buffer_L++] = obs_ind[i];
        else
          (right_child->obs_ind)[buffer_R++] = obs_ind[i];
      } else {
        if (Z[best_index*N + obs_ind[i]] <= best_split_point)
          (left_child->obs_ind)[buffer_L++] = obs_ind[i];
        else
          (right_child->obs_ind)[buffer_R++] = obs_ind[i];
      }
    }

    current_split_pred = set_values_and_next(current_split_pred, best_index + 1, bin_or_cont, best_split_point, 0.0);
    // splits[number_of_nodes - 1] = best_index + 1;

    knot->leaf = 0;
    knot->split = best_index;
    knot->split_bin_or_cont = bin_or_cont;
    knot->split_point = best_split_point;
    knot->left = left_child;
    knot->right = right_child;

    stack_push(stack, right_child);
    stack_push(stack, left_child);
  }

  if(mtry > 0) {
    PutRNGstate();
    Free(available_vars_raw);
    Free(available_vars);
  }

  stack_destroy(stack);

  functional** model_list = NULL;
  // If there is at least one quantitative covariable and 4pl fitting is desired,
  // use the first available covariable. (Additional ones are ignored.)
  if(pZ > 0 && covariable_mode == 2) {
    model_list = functionalLeaves(tree, number_of_nodes, bin_y, quant_y, y_bin, Z);
    // Update train_preds:
    pet_preds_t* func_preds = predictIntern(tree, X, Z, N, 0, 1);
    Free(train_preds);
    train_preds = func_preds->prob_preds;
    Free(func_preds);
  }

  int* splits_pointer = (int*) Calloc(number_of_nodes, int);
  int* splits_bin_or_cont_pointer = (int*) Calloc(number_of_nodes, int);
  double* split_points_pointer = (double*) Calloc(number_of_nodes, double);
  double* preds_pointer = (double*) Calloc(number_of_nodes, double);

  linked_list* old_element;
  for(i = 0; i < number_of_nodes; i++) {
    splits_pointer[i] = splits_preds->split;
    splits_bin_or_cont_pointer[i] = splits_preds->split_bin_or_cont;
    split_points_pointer[i] = splits_preds->split_point;
    preds_pointer[i] = splits_preds->pred;

    old_element = splits_preds;
    splits_preds = splits_preds->next;
    Free(old_element);
  }
  Free(splits_preds);

  pet_t* pet = (pet_t*) Calloc(1, pet_t);
  pet->splits = splits_pointer;
  pet->splits_bin_or_cont = splits_bin_or_cont_pointer;
  pet->split_points = split_points_pointer;
  pet->preds = preds_pointer;
  pet->train_preds = train_preds;
  pet->tree = tree;
  pet->number_of_nodes = number_of_nodes;
  pet->model_list = model_list;
  pet->y_bin = y_bin;

  return pet;
}

SEXP predict_(SEXP pet, SEXP X_raw, SEXP Z_raw, SEXP type_raw, SEXP leaves_raw) {
  rebuild_tree(pet);
  node* tree = (node*) R_ExternalPtrAddr(VECTOR_ELT(pet, 5));

  int* X = INTEGER(X_raw);
  double* Z = NULL;
  if(!isNull(Z_raw))
    Z = REAL(Z_raw);
  int type = LOGICAL(type_raw)[0];
  int leaves = asInteger(leaves_raw);
  int N = nrows(X_raw);

  pet_preds_t* pet_preds = predictIntern(tree, X, Z, N, type, leaves);
  SEXP ret;
  if(type == 1) {
    ret = PROTECT(allocVector(INTSXP, N));
    memcpy(INTEGER(ret), pet_preds->class_preds, N * sizeof(int));
    Free(pet_preds->class_preds);
  } else {
    ret = PROTECT(allocVector(REALSXP, N));
    memcpy(REAL(ret), pet_preds->prob_preds, N * sizeof(double));
  }

  Free(pet_preds->prob_preds);
  Free(pet_preds);

  UNPROTECT(1);
  return ret;
}

pet_preds_t* predictIntern(node* tree, int* X, double* Z, int N, int type, int leaves) {
  node* current_node;

  double* prob_preds_pointer = (double*) Calloc(N, double);
  int* class_preds_pointer = NULL;

  int i;
  for (i = 0; i < N; i++) {
    current_node = tree;

    while(!(current_node->leaf)) {
      if(!(current_node->split_bin_or_cont)) {
        if (X[(current_node->split)*N + i] == 0)
          current_node = current_node->left;
        else
          current_node = current_node->right;
      } else {
        if (Z[(current_node->split)*N + i] <= current_node->split_point)
          current_node = current_node->left;
        else
          current_node = current_node->right;
      }
    }

    if(current_node->func_pred == NULL || !leaves)
      prob_preds_pointer[i] = current_node->pred;
    else
      prob_preds_pointer[i] = eval4plModel(current_node->func_pred, Z[i]);
  }

  if (type == 1) {
    class_preds_pointer = (int*) Calloc(N, int);
    for (i = 0; i < N; i++) {
      if (prob_preds_pointer[i] > 0.5)
        class_preds_pointer[i] = 1;
      else
        class_preds_pointer[i] = 0;
    }
  }

  pet_preds_t* ret_pointer = (pet_preds_t*) Calloc(1, pet_preds_t);
  ret_pointer->class_preds = class_preds_pointer;
  ret_pointer->prob_preds = prob_preds_pointer;
  return ret_pointer;
}

SEXP getDesignMatrix_(SEXP X_raw, SEXP disj_raw, SEXP real_n_conj_raw) {
  int* X = INTEGER(X_raw);
  int p = ncols(X_raw);
  int N = nrows(X_raw);
  int* disj = INTEGER(disj_raw);
  int n_conj = nrows(disj_raw);
  int n_vars = ncols(disj_raw);
  int real_n_conj = INTEGER(real_n_conj_raw)[0];

  int* dm = getDesignMatrixIntern(X, N, disj, n_conj, n_vars, real_n_conj);

  SEXP design_matrix_raw = PROTECT(allocMatrix(INTSXP, N, real_n_conj));
  memcpy(INTEGER(design_matrix_raw), dm, N * real_n_conj * sizeof(int));
  Free(dm);

  int max_var_length;
  if (p < 10)
    max_var_length = 2;
  else if (p < 100)
    max_var_length = 3;
  else
    max_var_length = 4;
  int conj_length = n_vars * (max_var_length + 1);
  /*char* str_buffer = malloc(conj_length * real_n_conj * sizeof(char));*/
  char* str_buffer = (char*) Calloc(conj_length * real_n_conj, char);

  int written_chars;

  SEXP dimnames = PROTECT(allocVector(VECSXP, 2));
  SEXP colnames = PROTECT(allocVector(STRSXP, real_n_conj));

  for (int i = 0; i < real_n_conj; i++) {
    written_chars = sprintf(str_buffer + i * conj_length, "%d", disj[i]);
    for (int j = 1; j < n_vars && disj[j*n_conj + i] != NA_INTEGER; j++) {
      // strcat(str_buffer + i * conj_length, )
      written_chars += sprintf(str_buffer + i * conj_length + written_chars, "^%d", disj[j*n_conj + i]);
    }
    SET_STRING_ELT(colnames, i, mkChar(str_buffer + i * conj_length));
  }

  SET_VECTOR_ELT(dimnames, 0, getAttrib(design_matrix_raw, R_RowNamesSymbol));
  // SET_VECTOR_ELT(dimnames, 0, getAttrib(X_raw, R_RowNamesSymbol));
  SET_VECTOR_ELT(dimnames, 1, colnames);
  setAttrib(design_matrix_raw, R_DimNamesSymbol, dimnames);

  /*free(str_buffer);*/
  Free(str_buffer);
  UNPROTECT(3);
  return design_matrix_raw;
}

int* getDesignMatrixIntern(int* X, int N, int* disj, int n_conj, int n_vars, int real_n_conj) {
  // int real_n_conj = arrangeNAs(disj, n_conj, n_vars);

  int* design_matrix = (int*) Calloc(N * real_n_conj, int);

  int i, j, k;
  for (k = 0; k < N; k++) {
    for (i = 0; i < real_n_conj; i++) {
      design_matrix[i*N + k] = 1;
      for (j = 0; j < n_vars && disj[j*n_conj + i] != NA_INTEGER; j++) {
        if (disj[j*n_conj + i] < 0) {
          if (X[(-disj[j*n_conj + i]-1)*N + k] == 1) {
            design_matrix[i*N + k] = 0;
            break;
          }
        } else {
          if (X[(disj[j*n_conj + i]-1)*N + k] == 0) {
            design_matrix[i*N + k] = 0;
            break;
          }
        }
      }
    }
  }
  return design_matrix;
}

/*
int arrangeNAs(int* disj, int n_conj, int n_vars) {
  // Move NA rows to the bottom
  int i, j, k;
  int real_n_conj = 0;
  int n_NA = 0;
  bool swap;
  for (i = 0; i < n_conj;) {
    swap = true;
    for (j = 0; j < n_vars; j++) {
      if (disj[j*n_conj + i] != NA_INTEGER) {
        swap = false;
      }
    }
    if (swap) {
      // Only NA rows below
      if (n_conj - (i+1) <= n_NA) {
        break;
      }
      for (k = i; k < n_conj - 1 - n_NA; k++) {
        for (j = 0; j < n_vars; j++) {
          disj[j*n_conj + k] = disj[j*n_conj + k + 1];
          disj[j*n_conj + k + 1] = NA_INTEGER;
        }
      }
      n_NA++;
    } else {
      real_n_conj++;
      i++;
    }
  }

  // Move NA entries to the right
  for (i = 0; i < real_n_conj; i++) {
    n_NA = 0;
    // Break if only NAs on the right
    for (j = 0; j < n_vars - 1 - n_NA;) {
      if (disj[j*n_conj + i] == NA_INTEGER) {
        for (k = j; k < n_vars - 1 - n_NA; k++) {
          disj[k*n_conj + i] = disj[(k+1)*n_conj + i];
          disj[(k+1)*n_conj + i] = NA_INTEGER;
        }
        n_NA++;
      } else {
        j++;
      }
    }
  }

  return real_n_conj;
}
*/

double gini_decrease(double p_k_1, double p_L, double p_L_1, double p_R_1) {
  return 2*(p_k_1 * (1-p_k_1) - (p_L * (1-p_L_1) * p_L_1 + (1-p_L) * (1-p_R_1) * p_R_1));
}

double mse_impurity(int N_k, double y_sum, double y_sum_2) {
  double y_pred = y_sum/N_k;
  return y_sum_2 - 2 * y_pred * y_sum + y_pred * y_pred * N_k;
}

double mse_decrease(int N_k, int N_k_L, int N_k_R, double N_k_sum, double N_L_sum, double N_R_sum, double N_k_sum_2, double N_L_sum_2, double N_R_sum_2) {
  double decr = mse_impurity(N_k, N_k_sum, N_k_sum_2) - mse_impurity(N_k_L, N_L_sum, N_L_sum_2) - mse_impurity(N_k_R, N_R_sum, N_R_sum_2);
  return decr/N_k;
}

int stack_destroy(logic_stack_t *stack) {
  if (stack == NULL) {
    return ERR_INVAL;
  }
  while (stack->top != NULL) {
    struct stack_frame_s *frame = stack->top;
    stack->top = frame->next;
    /*free(frame);*/
    Free(frame);
  }
  /*free(stack);*/
  Free(stack);
  return SUCCESS;
}

int stack_empty(logic_stack_t *stack) {
  if (stack == NULL || stack->top == NULL) {
    return TRUE;
  } else {
    return FALSE;
  }
}

logic_stack_t *stack_new(void) {
  /*logic_stack_t *stack = malloc(sizeof(*stack));*/
  logic_stack_t *stack = (logic_stack_t*) Calloc(1, logic_stack_t);
  if (stack == NULL) {
    return NULL;
  }
  stack->top = NULL;
  return stack;
}

void *stack_pop(logic_stack_t *stack) {
  if (stack == NULL || stack->top == NULL) {
    return NULL;
  }
  struct stack_frame_s *frame = stack->top;
  void *data = frame->data;
  stack->top = frame->next;
  /*free(frame);*/
  Free(frame);
  return data;
}

void *queue_pop(logic_stack_t *stack) {
  if (stack == NULL || stack->top == NULL) {
    return NULL;
  }

  struct stack_frame_s *frame = stack->top;
  struct stack_frame_s *frame2 = frame->next;

  if (frame2 == NULL) {
    void *data = frame->data;
    stack->top = NULL;
    /*free(frame);*/
    Free(frame);
    return data;
  }

  while (frame2->next != NULL) {
    frame = frame->next;
    frame2 = frame->next;
  }

  void *data = frame2->data;
  frame->next = NULL;
  /*free(frame2);*/
  Free(frame2);
  return data;
}

int stack_push(logic_stack_t *stack, void *data) {
  if (stack == NULL) {
    return ERR_INVAL;
  }
  /*struct stack_frame_s *frame = malloc(sizeof(*frame));*/
  struct stack_frame_s *frame = (stack_frame_t*) Calloc(1, stack_frame_t);
  if (frame == NULL) {
    return ERR_NOMEM;
  }
  frame->data = data;
  frame->next = stack->top;
  stack->top = frame;
  return SUCCESS;
}

/*
void split_node(node* knot, int split) {
  knot->left = (node*) malloc(sizeof(node));
  knot->right = (node*) malloc(sizeof(node));
  knot->leaf = false;
  knot->split = split;
}
*/

double calcLeafProb(int N_k_1, int N_k, int smoothing) {
  if(!smoothing)
    return (double) N_k_1/N_k;
  else
    return (double) (N_k_1 + 1)/(N_k + 2);
}

void make_leaf(node* knot, double p_k_1, double* train_preds) {
  knot->leaf = 1;
  knot->pred = p_k_1;
  knot->split = -1;

  int* obs_ind = knot->obs_ind;
  for (int i = 0; i < knot->N_k; i++)
    train_preds[obs_ind[i]] = p_k_1;
}

void tree_destroy(node* tree) {
  if (tree == NULL)
    return;
  if (tree->left != NULL)
    tree_destroy(tree->left);
  if (tree->right != NULL)
    tree_destroy(tree->right);
  if (tree->obs_ind != NULL)
    Free(tree->obs_ind);
  if (tree->func_pred != NULL)
    Free(tree->func_pred);
  Free(tree);
}

void pet_destroy(pet_t* pet, int destroy_tree) {
  if(destroy_tree)
    tree_destroy(pet->tree);
  Free(pet->splits);
  Free(pet->splits_bin_or_cont);
  Free(pet->split_points);
  Free(pet->preds);
  Free(pet->train_preds);
  if(pet->model_list != NULL)
    Free(pet->model_list);
  Free(pet);
}

void rebuild_tree(SEXP pet) {
  node* c_pet = (node*) R_ExternalPtrAddr(VECTOR_ELT(pet, 5));
  if(c_pet != NULL) return;

  int number_of_nodes = length(VECTOR_ELT(pet, 0));
  int* splits = INTEGER(VECTOR_ELT(pet, 0));
  int* splits_bin_or_cont = INTEGER(VECTOR_ELT(pet, 1));
  double* split_points = REAL(VECTOR_ELT(pet, 2));
  double* preds = REAL(VECTOR_ELT(pet, 3));
  SEXP model_list_R = VECTOR_ELT(pet, 6);
  double* func_buffer = NULL;
  int covariable_mode = !isNull(model_list_R);
  if(covariable_mode && splits[0] == 0)
    func_buffer = REAL(VECTOR_ELT(VECTOR_ELT(model_list_R, 0), 0));
  int y_bin = asLogical(VECTOR_ELT(pet, 7));


  node* tree;
  tree = (node*) Calloc(1, node);
  node* current_node = tree;
  tree->left = NULL;
  tree->right = NULL;
  tree->split = splits[0] - 1;
  tree->split_bin_or_cont = splits_bin_or_cont[0];
  tree->split_point = split_points[0];
  tree->pred = preds[0];
  if(splits[0] != 0) tree->leaf = 0; else tree->leaf = 1;
  if(func_buffer != NULL) {
    tree->func_pred = (functional*) Calloc(1, functional);
    tree->func_pred->y_bin = y_bin;
    tree->func_pred->b = func_buffer[0];
    tree->func_pred->c = func_buffer[1];
    tree->func_pred->d = func_buffer[2];
    tree->func_pred->e = func_buffer[3];
  } else
    tree->func_pred = NULL;

  logic_stack_t *stack = stack_new();

  for(int i = 0; i < number_of_nodes; i++) {
    if(splits[i] != 0) {
      current_node->left = (node*) Calloc(1, node);
      stack_push(stack, current_node);
      current_node = current_node->left;
    } else if(i+1 < number_of_nodes) {
      current_node = stack_pop(stack);
      current_node->right = (node*) Calloc(1, node);
      current_node = current_node->right;
    }

    if(splits[i] != 0 || i+1 < number_of_nodes) {
      current_node->split = splits[i+1] - 1;
      current_node->split_bin_or_cont = splits_bin_or_cont[i+1];
      current_node->split_point = split_points[i+1];
      current_node->pred = preds[i+1];
      if(splits[i+1] != 0) current_node->leaf = 0; else current_node->leaf = 1;
      if(covariable_mode && splits[i+1] == 0) {
        func_buffer = REAL(VECTOR_ELT(VECTOR_ELT(model_list_R, i+1), 0));
        current_node->func_pred = (functional*) Calloc(1, functional);
        current_node->func_pred->y_bin = y_bin;
        current_node->func_pred->b = func_buffer[0];
        current_node->func_pred->c = func_buffer[1];
        current_node->func_pred->d = func_buffer[2];
        current_node->func_pred->e = func_buffer[3];
      } else
        current_node->func_pred = NULL;
    }
  }

  stack_destroy(stack);
  SEXP bin_tree = PROTECT(R_MakeExternalPtr(tree, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(bin_tree, _finalizer, TRUE);
  SET_VECTOR_ELT(pet, 5, bin_tree);
}

linked_list* set_values_and_next(linked_list* l, int split, int split_bin_or_cont, double split_point, double pred) {
  l->split = split;
  l->split_bin_or_cont = split_bin_or_cont;
  l->split_point = split_point;
  l->pred = pred;
  /*l->next = (linked_list*) malloc(sizeof(linked_list));*/
  l->next = (linked_list*) Calloc(1, linked_list);
  return l->next;
}

// 4pl model inside terminal nodes
functional** functionalLeaves(node* tree, int number_of_nodes, int* bin_y, double* quant_y, int y_bin, double* Z) {
  node* current_node;
  functional** model_list = (functional**) Calloc(number_of_nodes, functional*);
  logic_stack_t *stack = stack_new();
  stack_push(stack, tree);
  int i = 0;

  while(stack->top != NULL) {
    current_node = stack_pop(stack);
    if(current_node->leaf) {
      current_node->func_pred = fit4plModel(bin_y, quant_y, y_bin, current_node->pred, Z, current_node->N_k, current_node->obs_ind);
      model_list[i] = current_node->func_pred;
    } else {
      model_list[i] = NULL;
      stack_push(stack, current_node->right);
      stack_push(stack, current_node->left);
    }
    i++;
  }
  stack_destroy(stack);
  return model_list;
}

double binLogLikelihood(int n, double* par, void* ex) {
  double sum = 0;
  dataset* data = ex;
  int* y = data->bin_y;
  double* Z = data->Z;
  int* obs_ind = data->obs_ind;
  int N = data->N;
  // Parameter scaling as in drc::drm
  double* par_scale = data->par_scale;
  double b = par[0] * par_scale[0];
  double c = par[1] * par_scale[1];
  double d = par[2] * par_scale[2];
  double e = par[3] * par_scale[3];
  double prob;
  double lower_tol = 1e-12;
  double upper_tol = 1-1e-12;
  for(int i = 0; i < N; i++) {
    prob = c + (d-c)/(1+exp(b*(Z[obs_ind[i]] - e)));
    if(prob > upper_tol)
      prob = upper_tol;
    else if(prob < lower_tol)
      prob = lower_tol;
    if(y[obs_ind[i]]) {
      sum += log(prob);
    } else {
      sum += log(1 - prob);
    }
  }
  return -sum;
}

void binLogLikelihoodGrad(int n, double* par, double* gr, void* ex) {
  memset(gr, 0, 4 * sizeof(double));
  dataset* data = ex;
  int* y = data->bin_y;
  double* Z = data->Z;
  int* obs_ind = data->obs_ind;
  int N = data->N;
  double* par_scale = data->par_scale;
  double b = par[0] * par_scale[0];
  double c = par[1] * par_scale[1];
  double d = par[2] * par_scale[2];
  double e = par[3] * par_scale[3];
  double buf;
  for(int i = 0; i < N; i++) {
    buf = exp(b*(Z[obs_ind[i]] - e));
    if(y[obs_ind[i]]) {
      gr[0] -= d*(e - Z[obs_ind[i]])/(c*buf + d) + (Z[obs_ind[i]]-e)/(buf + 1);
      // Utilize asymptotic behavior, if exponential explodes
      if(!R_FINITE(buf))
        gr[1] -= 1/c;
      else
        gr[1] -= buf/(buf*c + d);
      gr[2] -= 1/(d + c*buf);
      gr[3] -= b*d/(c*buf + d) - b/(buf + 1);
    } else {
      gr[0] -= (d-1)*(e - Z[obs_ind[i]])/((c-1)*buf + d - 1) + (Z[obs_ind[i]]-e)/(buf + 1);
      // Utilize asymptotic behavior, if exponential explodes
      if(!R_FINITE(buf))
        gr[1] -= 1/(c-1);
      else
        gr[1] -= buf/(buf*(c-1) + d - 1);
      gr[2] -= 1/(d - 1 + (c-1)*buf);
      gr[3] -= b*(d-1)/((c-1)*buf + d - 1) - b/(buf + 1);
    }
  }
  gr[0] *= par_scale[0];
  gr[1] *= par_scale[1];
  gr[2] *= par_scale[2];
  gr[3] *= par_scale[3];
}

double squaredError(int n, double* par, void* ex) {
  double sum = 0;
  dataset* data = ex;
  double* y = data->quant_y;
  double* Z = data->Z;
  int* obs_ind = data->obs_ind;
  int N = data->N;
  double* par_scale = data->par_scale;
  double b = par[0] * par_scale[0];
  double c = par[1] * par_scale[1];
  double d = par[2] * par_scale[2];
  double e = par[3] * par_scale[3];
  double buf;
  for(int i = 0; i < N; i++) {
    buf = c + (d-c)/(1+exp(b*(Z[obs_ind[i]] - e)));
    sum += (y[obs_ind[i]] - buf) * (y[obs_ind[i]] - buf);
  }
  return sum/N;
}

void squaredErrorGrad(int n, double* par, double* gr, void* ex) {
  memset(gr, 0, 4 * sizeof(double));
  dataset* data = ex;
  double* y = data->quant_y;
  double* Z = data->Z;
  int* obs_ind = data->obs_ind;
  int N = data->N;
  double* par_scale = data->par_scale;
  double b = par[0] * par_scale[0];
  double c = par[1] * par_scale[1];
  double d = par[2] * par_scale[2];
  double e = par[3] * par_scale[3];
  double buf;
  for(int i = 0; i < N; i++) {
    buf = exp(b*(Z[obs_ind[i]] - e));
    // Utilize asymptotic behavior, if exponential explodes
    if(!R_FINITE(buf)) {
      gr[1] += 2 * (c-y[obs_ind[i]]);
    } else {
      gr[0] += -2 * (c-d) * (e-Z[obs_ind[i]]) * buf * ((c-y[obs_ind[i]]) * buf + d - y[obs_ind[i]]) / ((buf+1)*(buf+1)*(buf+1));
      gr[1] += 2 * buf * ((c-y[obs_ind[i]]) * buf + d - y[obs_ind[i]])/((buf+1)*(buf+1));
      gr[2] += 2*((c-y[obs_ind[i]]) * buf + d - y[obs_ind[i]])/((buf+1)*(buf+1));
      gr[3] += -2 * b * (c-d) * buf * ((c-y[obs_ind[i]]) * buf + d - y[obs_ind[i]]) /((buf+1)*(buf+1)*(buf+1));
    }
  }
  gr[0] *= par_scale[0]/N;
  gr[1] *= par_scale[1]/N;
  gr[2] *= par_scale[2]/N;
  gr[3] *= par_scale[3]/N;
}

SEXP fit4plModel_(SEXP y, SEXP Z) {
  int* bin_y = NULL;
  double* quant_y = NULL;
  int y_bin = 0;
  int N = length(y);
  double* Z2 = REAL(Z);
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
  functional* model = fit4plModel(bin_y, quant_y, y_bin, y_mean, Z2, N, obs_ind);
  // Correctly embed 4pL model
  SEXP model_R = PROTECT(allocVector(VECSXP, 2));
  SEXP bcde_R = allocVector(REALSXP, 4);
  SET_VECTOR_ELT(model_R, 0, bcde_R);
  SET_VECTOR_ELT(model_R, 1, ScalarLogical(y_bin));
  double* model_pointer = REAL(bcde_R);
  model_pointer[0] = model->b;
  model_pointer[1] = model->c;
  model_pointer[2] = model->d;
  model_pointer[3] = model->e;
  Free(obs_ind);
  Free(model);
  // Assign class "4pl"
  classgets(model_R, mkString("4pl"));
  UNPROTECT(1);
  return model_R;
}

functional* fit4plModel(int* bin_y, double* quant_y, int y_bin, double y_mean, double* Z, int N, int* obs_ind) {
  optimfn* fn;
  optimgr* gr;

  double b, c, d, e;
  double scaler = 0.001;

  functional* ret = (functional*) Calloc(1, functional);
  ret->y_bin = y_bin;

  if(y_bin) {
    if(doubleEquals(y_mean, 0) || doubleEquals(y_mean, 1)) {
      ret->c = y_mean;
      ret->d = y_mean;
      ret->b = 0;
      ret->e = 0;
      return ret;
    }
    c = -scaler;
    d = 1 + scaler;
    fn = &binLogLikelihood;
    gr = &binLogLikelihoodGrad;
  } else {
    double min_y = R_PosInf;
    double max_y = R_NegInf;
    for(int i = 0; i < N; i++) {
      if(quant_y[obs_ind[i]] < min_y) min_y = quant_y[obs_ind[i]];
      if(quant_y[obs_ind[i]] > max_y) max_y = quant_y[obs_ind[i]];
    }
    if(doubleEquals(min_y, max_y)) {
      ret->c = y_mean;
      ret->d = y_mean;
      ret->b = 0;
      ret->e = 0;
      return ret;
    }
    double diff = max_y - min_y;
    c = min_y - scaler * diff;
    d = max_y + scaler * diff;
    fn = &squaredError;
    gr = &squaredErrorGrad;
  }
  // Response for linear model finding initial values for b and e
  double* initResp = (double*) Calloc(N, double);
  // Z2 to ensure that each entry in Z2 corresponds to each entry in initResp
  double* Z2 = (double*) Calloc(N, double);
  if(y_bin) {
    for(int i = 0; i < N; i++)
      initResp[i] = log((d-bin_y[obs_ind[i]])/(bin_y[obs_ind[i]]-c));
  } else {
    for(int i = 0; i < N; i++)
      initResp[i] = log((d-quant_y[obs_ind[i]])/(quant_y[obs_ind[i]]-c));
  }
  for(int i = 0; i < N; i++)
    Z2[i] = Z[obs_ind[i]];
  double* initialPars = fitLinearModel(Z2, initResp, N);
  b = initialPars[1];
  e = -initialPars[0]/b;
  Free(initResp);
  Free(Z2);
  Free(initialPars);
  // Scale parameters
  double* Pars = (double*) Calloc(4, double);
  Pars[0] = b;
  Pars[1] = c;
  Pars[2] = d;
  Pars[3] = e;
  double* par_scale = (double*) Calloc(4, double);
  memcpy(par_scale, Pars, 4 * sizeof(double));
  for(int i = 0; i < 4; i++) {
    par_scale[i] = fabs(par_scale[i]);
    if(par_scale[i] < 1e-4) par_scale[i] = 1;
    Pars[i] /= par_scale[i];
  }
  double min_val = 0;
  int maxit = 500; // Standard value 100
  int trace = 0;
  int* mask = (int*) Calloc(4, int);
	for (int i = 0; i < 4; i++) mask[i] = 1;
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
  int fncount = 0;
  int grcount = 0;
  int fail = 0;

  // Catch not finite error:
  double f = fn(4, Pars, ex);
  if (!R_FINITE(f)) {
    ret->c = y_mean;
    ret->d = y_mean;
    ret->b = 0;
    ret->e = 0;
  } else {
    vmmin(4, Pars, &min_val, fn, gr, maxit, trace, mask, abstol, reltol, nREPORT, ex, &fncount, &grcount, &fail);
    /* if(fail) {
      ret->c = y_mean;
      ret->d = y_mean;
      ret->b = 0;
      ret->e = 0;
    } else {
      // Revert parameter scaling
      ret->b = Pars[0] * par_scale[0];
      ret->c = Pars[1] * par_scale[1];
      ret->d = Pars[2] * par_scale[2];
      ret->e = Pars[3] * par_scale[3];
    } */

    // Ignore fail, since it only captures whether the maximum number of iterations was reached
    // Revert parameter scaling
    ret->b = Pars[0] * par_scale[0];
    ret->c = Pars[1] * par_scale[1];
    ret->d = Pars[2] * par_scale[2];
    ret->e = Pars[3] * par_scale[3];
  }
  Free(Pars);
  Free(par_scale);
  Free(mask);
  Free(ex);
  return ret;
}

double eval4plModel(functional* func_pred, double Z) {
  double ret;
  if(doubleEquals(func_pred->c, func_pred->d))
    ret = func_pred->c;
  else
    ret = func_pred->c + (func_pred->d - func_pred->c)/(1+exp(func_pred->b*(Z - func_pred->e)));
  if(func_pred->y_bin) {
    if(ret > 1)
      ret = 1;
    else if(ret < 0)
      ret = 0;
  }
  return ret;
}

double* fitLinearModel(double* x, double* y, int N) {
  double x_mean = 0;
  double y_mean = 0;
  for(int i = 0; i < N; i++) {
    x_mean += x[i];
    y_mean += y[i];
  }
  x_mean /= N;
  y_mean /= N;
  double numerator = 0;
  double denominator = 0;
  for(int i = 0; i < N; i++) {
    numerator += x[i] * y[i];
    denominator += x[i] * x[i];
  }
  double* beta = (double*) Calloc(2, double);
  beta[1] = (numerator - N * x_mean * y_mean)/(denominator - N * x_mean * x_mean);
  beta[0] = y_mean - beta[1] * x_mean;
  return beta;
}

int doubleEquals(double a, double b) {
  return (fabs(a - b) <= 0.000001 * fabs(a));
}

int cmp_integer(const void* value1, const void* value2) {
  return (*(int*)value1 - *(int*)value2);
}

int cmp_double(const void* value1, const void* value2) {
  return (*(double*)value1 > *(double*)value2) ? 1 : (*(double*)value1 < *(double*)value2) ? -1:0 ;
}

int cmp_integer_direct(int value1, int value2) {
  return (value1 - value2);
}

int cmp_double_direct(double value1, double value2) {
  return (value1 > value2) ? 1 : (value1 < value2) ? -1:0 ;
}

int cmp_y_probs_int(const void* value1, const void* value2) {
  return cmp_integer_direct((*(y_prob_pair_t*)value1).y, (*(y_prob_pair_t*)value2).y);
}

int cmp_y_probs_double(const void* value1, const void* value2) {
  return cmp_double_direct((*(y_prob_pair_t*)value1).prob, (*(y_prob_pair_t*)value2).prob);
}

int cmp_y_Z_pair(const void* value1, const void* value2) {
  return cmp_double_direct((*(y_Z_pair_t*)value1).Z, (*(y_Z_pair_t*)value2).Z);
}

