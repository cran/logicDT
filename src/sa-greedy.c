#include "sa-greedy.h"

int prevented_evals;
int* mtry_vars_raw;
int* mtry_vars;
int total_length_global;

double buildModel(SEXP X_train, SEXP y_train, SEXP Z_train, SEXP Z_val, int* disj, int n_conj, int n_vars, int nodesize, double cp, int smoothing, int mtry, int covariable_mode, int scoring_rule, SEXP X_val, SEXP y_val, int use_validation, int y_bin, eval_models_list** models) {
  int i;
  for(i = 0; i < n_conj; i++) {
    if(disj[i] == NA_INTEGER)
      break;
  }
  int real_n_conj = i;

  if(models != NULL) {
    int total_length = n_conj * n_vars;
    int hash = calcDisjHash2(disj, total_length);
    eval_models_list* model_entry = models[hash];
    eval_models_list* old_model_entry = model_entry;
    int cmp_val;
    int change_start_point = 1;
    while(model_entry != NULL) {
      cmp_val = cmp_disj(disj, model_entry->disj, total_length);
      if(cmp_val == 0) {
        prevented_evals++;
        return model_entry->score;
      }
      else if(cmp_val < 0) {
        old_model_entry = model_entry;
        model_entry = model_entry->next;
        change_start_point = 0;
      }
      else break;
    }

    pet_ensemble_t* eval = fitPETsIntern(X_train, y_train, X_val, y_val, Z_train, Z_val, use_validation, y_bin, nodesize, cp, smoothing, mtry, covariable_mode, disj, n_conj, n_vars, real_n_conj, scoring_rule, 0);
    double score = eval->score;
    Free(eval);
    eval_models_list* new_model_entry = (eval_models_list*) Calloc(1, eval_models_list);
    new_model_entry->score = score;
    new_model_entry->disj = (int*) Calloc(total_length, int);
    memcpy(new_model_entry->disj, disj, total_length * sizeof(int));
    new_model_entry->next = model_entry;
    if(change_start_point) models[hash] = new_model_entry;
    else old_model_entry->next = new_model_entry;
    return score;
  }

  pet_ensemble_t* eval = fitPETsIntern(X_train, y_train, X_val, y_val, Z_train, Z_val, use_validation, y_bin, nodesize, cp, smoothing, mtry, covariable_mode, disj, n_conj, n_vars, real_n_conj, scoring_rule, 0);
  double score = eval->score;
  Free(eval);
  return score;
}

sa_eval_t* evaluateModel(SEXP X_train, SEXP y_train, SEXP Z_train, SEXP Z_val, double t, int acc_type, int* disj, int n_conj, int n_vars, double old_score, int nodesize, double cp, int smoothing, int mtry, int covariable_mode, int scoring_rule, SEXP X_val, SEXP y_val, int use_validation, int y_bin, eval_models_list** models) {
  double new_score = buildModel(X_train, y_train, Z_train, Z_val, disj, n_conj, n_vars, nodesize, cp, smoothing, mtry, covariable_mode, scoring_rule, X_val, y_val, use_validation, y_bin, models);

  double acc, rnd;
  if(!acc_type) {
    acc = exp((old_score - new_score)/t);
    GetRNGstate();
    rnd = unif_rand();
    PutRNGstate();
  } else {
    acc = old_score - new_score;
    rnd = -t;
  }

  sa_eval_t* ret_obj = (sa_eval_t*) Calloc(1, sa_eval_t);

  if(acc > rnd) {
    ret_obj->score = new_score;
    ret_obj->acc = 1;
  }
  else {
    ret_obj->score = old_score;
    ret_obj->acc = 0;
  }

  return ret_obj;
}

double getGPScore(double score, int y_bin, int scoring_rule) {
  // Handle deviance, brier, and MSE as in [0, infinity)
  // Use that AUC and NCE are bounded
  if(y_bin) {
    if(scoring_rule <= 1) {
      if(score == 0)
        score = 1e30;
      else
        score = 1/score;
    }
    else if(scoring_rule == 5)
      score = 1 - score;
    else
      score = 1 - score - 0.5;
  } else {
    if(score == 0)
      score = 1e30;
    else
      score = 1/score;
  }
  score = (score > 0) ? (score) : (0);
  return score;
}

void calcSharedFitness(gen_t* generation, int n_ind, int max_conj, int max_vars, double sigma) {
  if(sigma == 0)
    return;
  double d;
  for(int i = 0; i < n_ind; i++)
    generation[i].fitness = 0;
  for(int i = 0; i < n_ind; i++) {
    for(int j = i; j < n_ind; j++) {
      d = 0;
      for(int k = 0; k < max_conj; k++) {
        // d += generation[i].disj[k] != generation[j].disj[k];
        for(int l = 0; l < max_vars; l++) {
          d += generation[i].disj[l * max_conj + k] != generation[j].disj[l * max_conj + k];
          if(generation[i].disj[l * max_conj + k] == NA_INTEGER || generation[i].disj[l * max_conj + k] == NA_INTEGER)
            break;
        }
      }
      d /= max_vars;
      // d += 5 * (generation[i].score - generation[j].score) * (generation[i].score - generation[j].score);
      if(d < sigma) {
        generation[i].fitness += 1 - d/sigma;
        if(i != j)
          generation[j].fitness += 1 - d/sigma;
      }
    }
    // generation[i].fitness /= n_ind - 1;
    generation[i].fitness = generation[i].score/generation[i].fitness;
  }
}

int* tournamentSelection(gen_t* generation, int n_ind, int k) {
  int* ret = (int*) Calloc(k, int);
  for(int i = 0; i < k; i++) {
    ret[i] = unif_rand() * n_ind;
  }
  qsort(ret, k, sizeof(int), cmp_integer);
  return ret;
}

int* mutateGeneration(gen_t* generation, int n_ind, int max_vars, int max_conj, int p, int type, int which_ind, int allow_conj_removal) {
  // Type coding
  // 0: Replace variable, 1: Add variable, 2: Remove variable, 3: Add conjunction, 4: Remove conjunction

  int i, j, n_conj;
  double rnd;
  int* disj_temp = generation[which_ind].disj;
  for(n_conj = 0; n_conj < max_conj; n_conj++) {
    if(disj_temp[n_conj] == NA_INTEGER)
      break;
  }

  int n_vars_total = 0;
  for(i = 0; i < n_conj; i++) {
    for(j = 0; j < max_vars; j++) {
      if(disj_temp[j * max_conj + i] == NA_INTEGER)
        break;
      n_vars_total++;
    }
  }

  if((type == 3 && (n_conj == max_conj || max_vars <= n_vars_total)) || (type == 4 && (n_conj == 1 || !allow_conj_removal)))
    return NULL;

  int* disj2 = (int*) Calloc(max_vars * max_conj, int);
  memcpy(disj2, disj_temp, max_vars * max_conj * sizeof(int));

  if(type < 3) {
    // Conjunction modifications
    rnd = unif_rand();
    int which_conj = rnd * n_conj;
    int n_vars_here;
    for(n_vars_here = 0; n_vars_here < max_vars; n_vars_here++) {
      if(disj_temp[n_vars_here * max_conj + which_conj] == NA_INTEGER)
        break;
    }

    if((type == 1 && max_vars <= n_vars_total) || (type == 2 && n_vars_here == 1)) {
      Free(disj2);
      return NULL;
    }

    rnd = unif_rand();
    int which_var = rnd * n_vars_here;

    if(type < 2) {
      int* unused_vars = (int*) Calloc(p, int);
      memset(unused_vars, 0, p * sizeof(int));
      int* available_vars = (int*) Calloc(2 * p, int);
      for(i = 0; i < n_vars_here; i++) {
        unused_vars[abs(disj_temp[i * max_conj + which_conj]) - 1] = 1;
      }
      j = 0;
      for(i = 0; i < p; i++) {
        if(!unused_vars[i]) {
          available_vars[j] = i+1;
          available_vars[j+1] = -(i+1);
          j += 2;
        }
      }

      if(type == 0) {
        rnd = unif_rand();
        int swap_var = rnd * (j+1);

        if(swap_var == j) {
          disj2[which_var * max_conj + which_conj] = -disj2[which_var * max_conj + which_conj];
        } else {
          disj2[which_var * max_conj + which_conj] = available_vars[swap_var];
        }
      } else {
        rnd = unif_rand();
        int add_var = rnd * j;
        disj2[n_vars_here * max_conj + which_conj] = available_vars[add_var];
      }
      Free(unused_vars);
      Free(available_vars);
    } else {
      disj2[which_var * max_conj + which_conj] = disj2[(n_vars_here-1) * max_conj + which_conj];
      disj2[(n_vars_here-1) * max_conj + which_conj] = NA_INTEGER;
    }
  } else if (type == 3) {
    rnd = unif_rand();
    int add_var = rnd * 2 * p;
    disj2[n_conj] = 1 + (add_var/2);
    if(add_var % 2 == 1)
      disj2[n_conj] = -disj2[n_conj];
  } else {
    rnd = unif_rand();
    int which_rem = rnd * n_conj;

    for(i = 0; i < max_vars; i++) {
      disj2[i * max_conj + which_rem] = disj2[i * max_conj + (n_conj-1)];
      disj2[i * max_conj + (n_conj-1)] = NA_INTEGER;
    }
  }
  return disj2;
}

gp_eval_t* geneticProgrammingStep(SEXP X_train, SEXP y_train, int max_vars, int max_conj, SEXP Z_train, SEXP Z_val, gen_t* generation, int n_ind, double best_score, int nodesize, double cp, int smoothing, int mtry, int covariable_mode, int scoring_rule, SEXP X_val, SEXP y_val, int use_validation, int y_bin, int allow_conj_removal, int conjsize, SEXP X, eval_models_list** models) {

  int p = ncols(VECTOR_ELT(X_train, 0));
  double rnd;

  // Crossover
  int cross_ind1, cross_ind2;
  double prob;
  if(best_score == 0) {
    rnd = unif_rand();
    cross_ind1 = rnd * n_ind;
    cross_ind2 = cross_ind1;
    while(cross_ind2 == cross_ind1) {
      rnd = unif_rand();
      cross_ind2 = rnd * n_ind;
    }
  } else {
    cross_ind1 = unif_rand() * n_ind;
    prob = generation[cross_ind1].score/best_score;
    rnd = unif_rand();
    while(prob < rnd) {
      cross_ind1 = unif_rand() * n_ind;
      prob = generation[cross_ind1].score/best_score;
      rnd = unif_rand();
    }
    cross_ind2 = unif_rand() * n_ind;
    prob = generation[cross_ind2].score/best_score;
    rnd = unif_rand();
    while(cross_ind2 == cross_ind1 || prob < rnd) {
      cross_ind2 = unif_rand() * n_ind;
      prob = generation[cross_ind2].score/best_score;
      rnd = unif_rand();
    }
  }

  // Random conjunction in the first individual is going to be replaced
  // by a random conjunction of the second individual
  int n_conj1;
  for(n_conj1 = 0; n_conj1 < max_conj; n_conj1++) {
    if(generation[cross_ind1].disj[n_conj1] == NA_INTEGER)
      break;
  }
  int n_conj2;
  for(n_conj2 = 0; n_conj2 < max_conj; n_conj2++) {
    if(generation[cross_ind2].disj[n_conj2] == NA_INTEGER)
      break;
  }

  /*rnd = unif_rand();
  int conj_ind1 = rnd * n_conj1;*/
  int conj_ind1, conj_ind11;
  rnd = unif_rand();
  int conj_ind2 = rnd * n_conj2;

  if(conj_ind2 >= n_conj1) {
    conj_ind1 = unif_rand() * (n_conj1 + 1); // Replace or add
    conj_ind11 = unif_rand() * n_conj1; // Source conjunction
  }
  else {
    conj_ind1 = conj_ind2;
    conj_ind11 = conj_ind2;
  }

  generation[n_ind].disj = (int*) Calloc(max_vars * max_conj, int);
  generation[n_ind + 1].disj = (int*) Calloc(max_vars * max_conj, int);
  memcpy(generation[n_ind].disj, generation[cross_ind1].disj, max_vars * max_conj * sizeof(int));
  memcpy(generation[n_ind + 1].disj, generation[cross_ind2].disj, max_vars * max_conj * sizeof(int));
  int i;
  for(i = 0; i < max_vars; i++) {
    generation[n_ind].disj[i * max_conj + conj_ind1] = generation[cross_ind2].disj[i * max_conj + conj_ind2];
    generation[n_ind + 1].disj[i * max_conj + conj_ind2] = generation[cross_ind1].disj[i * max_conj + conj_ind11];
  }
  generation[n_ind].score = getGPScore(buildModel(X_train, y_train, Z_train, Z_val, generation[n_ind].disj, max_conj, max_vars, nodesize, cp, smoothing, mtry, covariable_mode, scoring_rule, X_val, y_val, use_validation, y_bin, models), y_bin, scoring_rule);
  generation[n_ind + 1].score = getGPScore(buildModel(X_train, y_train, Z_train, Z_val, generation[n_ind + 1].disj, max_conj, max_vars, nodesize, cp, smoothing, mtry, covariable_mode, scoring_rule, X_val, y_val, use_validation, y_bin, models), y_bin, scoring_rule);
  generation[n_ind].fitness = generation[n_ind].score;
  generation[n_ind + 1].fitness = generation[n_ind + 1].score;
  n_ind += 2;

  // Mutations (5)
  int* disj2;
  int iter = 2;
  int old_n_ind = n_ind - 2;
  for(i = 0; i < 5; i++) {
    if(best_score == 0) {
      cross_ind1 = unif_rand() * old_n_ind;
    } else {
      cross_ind1 = unif_rand() * old_n_ind;
      prob = generation[cross_ind1].score/best_score;
      rnd = unif_rand();
      while(prob < rnd) {
        cross_ind1 = unif_rand() * old_n_ind;
        prob = generation[cross_ind1].score/best_score;
        rnd = unif_rand();
      }
    }
    disj2 = mutateGeneration(generation, n_ind, max_vars, max_conj, p, i, cross_ind1, allow_conj_removal);
    if(disj2 != NULL) {
      generation[n_ind].disj = disj2;
      generation[n_ind].score = getGPScore(buildModel(X_train, y_train, Z_train, Z_val, disj2, max_conj, max_vars, nodesize, cp, smoothing, mtry, covariable_mode, scoring_rule, X_val, y_val, use_validation, y_bin, models), y_bin, scoring_rule);
      generation[n_ind].fitness = generation[n_ind].score;
      n_ind++;
      iter++;
    }
  }

  gp_eval_t* ret = (gp_eval_t*) Calloc(1, gp_eval_t);
  ret->generation = generation;
  ret->n_ind = n_ind;
  ret->iter = iter;
  return ret;
}

SEXP geneticProgramming_(SEXP X_train, SEXP y_train, SEXP max_vars_raw, SEXP max_conj_raw, SEXP max_gen_raw, SEXP gp_sigma_raw, SEXP gp_fs_interval_raw, SEXP Z_train, SEXP Z_val, SEXP nodesize_raw, SEXP cp_raw, SEXP smoothing_raw, SEXP mtry_raw, SEXP covariable_mode_raw, SEXP scoring_rule_raw, SEXP X_val, SEXP y_val, SEXP use_validation_raw, SEXP y_bin_raw, SEXP allow_conj_removal_raw, SEXP conjsize_raw, SEXP X) {
  int n_vars = asInteger(max_vars_raw);
  int n_conj = asInteger(max_conj_raw);
  int max_gen = asInteger(max_gen_raw);
  int p = ncols(VECTOR_ELT(X_train, 0));
  double gp_sigma = asReal(gp_sigma_raw);
  int gp_fs_interval = asInteger(gp_fs_interval_raw);

  int nodesize = asInteger(nodesize_raw);
  double cp = asReal(cp_raw);
  int smoothing = asInteger(smoothing_raw);
  int mtry = asInteger(mtry_raw);
  int covariable_mode = asInteger(covariable_mode_raw);
  int scoring_rule = asInteger(scoring_rule_raw);
  int use_validation = asLogical(use_validation_raw);
  int y_bin = asLogical(y_bin_raw);
  int allow_conj_removal = asLogical(allow_conj_removal_raw);
  int conjsize = asInteger(conjsize_raw);

  int pop_size = p;
  int n_ind = pop_size;
  int total_iter = pop_size;
  // Ensure that 6 slots are always available
  int reserved_slots = 7;
  int allocd = pop_size + reserved_slots;

  double best_score = 0;
  double best_fitness = 0;

  eval_models_list** models = (eval_models_list**) Calloc(PRIME, eval_models_list*);
  prevented_evals = 0;

  GetRNGstate();

  gen_t* generation = (gen_t*) Calloc(allocd, gen_t);
  for(int i = 0; i < n_ind; i++) {
    generation[i].disj = (int*) Calloc(n_vars * n_conj, int);
    for(int j = 0; j < n_vars * n_conj; j++) {
      generation[i].disj[j] = NA_INTEGER;
    }
    // generation[i].disj[0] = (unif_rand() * p) + 1;
    generation[i].disj[0] = i + 1;
    generation[i].score = getGPScore(buildModel(X_train, y_train, Z_train, Z_val, generation[i].disj, n_conj, n_vars, nodesize, cp, smoothing, mtry, covariable_mode, scoring_rule, X_val, y_val, use_validation, y_bin, models), y_bin, scoring_rule);
    generation[i].fitness = generation[i].score;
  }
  qsort(generation, n_ind, sizeof(gen_t), cmp_gen_score);
  best_score = generation[0].score;
  best_fitness = generation[0].fitness;

  gp_eval_t* current_eval;
  double rnd, prob;
  int l, m, old_n_ind;
  for(int i = 0; i < max_gen; i++) {
    if(allocd < n_ind + reserved_slots) {
      allocd += 1000;
      generation = (gen_t*) Realloc(generation, allocd, gen_t);
    }

    current_eval = geneticProgrammingStep(X_train, y_train, n_vars, n_conj, Z_train, Z_val, generation, n_ind, best_score, nodesize, cp, smoothing, mtry, covariable_mode, scoring_rule, X_val, y_val, use_validation, y_bin, allow_conj_removal, conjsize, X, models);
    total_iter += current_eval->iter;

    if(i % gp_fs_interval == 0)
      calcSharedFitness(generation, current_eval->n_ind, n_conj, n_vars, gp_sigma);
    qsort(generation, current_eval->n_ind, sizeof(gen_t), cmp_gen_fitness);
    best_fitness = generation[0].fitness;
    best_score = 0;
    for(int j = 0; j < current_eval->n_ind; j++) {
      if(best_score < generation[j].score)
        best_score = generation[j].score;
    }

    if(best_fitness != 0) {
      // Stochastic Acceptance I
      old_n_ind = n_ind;
      n_ind = 2;
      l = 2;
      m = 0;
      for(int j = 2; j < current_eval->n_ind; j++) {
        rnd = unif_rand();
        prob = generation[l].fitness/best_fitness;
        // Always keep newly bred children and mutated individuals
        // But only if it's not the last generation
        if(prob > rnd || (j >= old_n_ind && i != max_gen - 1)) {
          n_ind++;
          l++;
        } else {
          Free(generation[l].disj);
          for(int k = l; k < current_eval->n_ind - 1 - m; k++) {
            generation[k].disj = generation[k+1].disj;
            generation[k].score = generation[k+1].score;
            generation[k].fitness = generation[k+1].fitness;
          }
          m++;
        }
      }
    } else {
      n_ind = current_eval->n_ind;
    }
    Free(current_eval);

    Rprintf("\r Generation %d/%d (%.0f%%) | Number of Individuals: %d", i+1, max_gen, (double) (i+1)/max_gen * 100, n_ind);
  }
  Rprintf("\n");

  PutRNGstate();

  destroy_eval_models(models);

  // Remove duplicate disj's
  total_length_global = n_conj * n_vars;
  qsort(generation, n_ind, sizeof(gen_t), cmp_gen_conj);
  for(int i = 0; i < n_ind - 1; i++) {
    while((generation[i+1].disj != NULL) && (cmp_disj_fixed(generation[i].disj, generation[i+1].disj) == 0)) {
      Free(generation[i+1].disj);
      for(int k = i+1; k < n_ind - 1; k++) {
        generation[k].disj = generation[k+1].disj;
        generation[k].score = generation[k+1].score;
      }
      generation[n_ind - 1].disj = NULL;
      n_ind--;
    }
  }
  qsort(generation, n_ind, sizeof(gen_t), cmp_gen_score);

  SEXP ret_obj = PROTECT(allocVector(VECSXP, 5));
  SEXP disjs = PROTECT(allocVector(VECSXP, n_ind));
  SEXP scores = PROTECT(allocVector(REALSXP, n_ind));

  SEXP disj_R;
  double* scores_pointer = REAL(scores);

  for(int i = 0; i < n_ind; i++) {
    disj_R = PROTECT(allocMatrix(INTSXP, n_conj, n_vars));
    memcpy(INTEGER(disj_R), generation[i].disj, n_conj * n_vars * sizeof(int));
    SET_VECTOR_ELT(disjs, i, disj_R);
    scores_pointer[i] = generation[i].score;
    Free(generation[i].disj);
  }

  SET_VECTOR_ELT(ret_obj, 0, disjs);
  SET_VECTOR_ELT(ret_obj, 1, scores);
  SET_VECTOR_ELT(ret_obj, 2, ScalarInteger(total_iter));
  SET_VECTOR_ELT(ret_obj, 3, ScalarInteger(prevented_evals));
  SET_VECTOR_ELT(ret_obj, 4, ScalarReal(best_score));

  Free(generation);
  UNPROTECT(3 + n_ind);
  return ret_obj;
}

SEXP predictGP_(SEXP model, SEXP X_raw, SEXP Z_raw, SEXP type_raw, SEXP n_models_raw, SEXP leaves_raw) {
  // Type coding
  // 0: Best, 1: All, 2: n_models best

  SEXP disjs = getListElement(model, "disj");
  SEXP ensemble = getListElement(model, "ensemble");
  int type = asInteger(type_raw);
  int n_models = asInteger(n_models_raw);
  int leaves = asInteger(leaves_raw);
  int n_ind = length(disjs);
  n_models = (n_models > n_ind) ? (n_ind) : (n_models);

  SEXP pet;
  SEXP inner_ensemble;
  SEXP disj_raw;
  node* tree;
  int N = nrows(X_raw);
  int* X = INTEGER(X_raw);
  int real_n_conj;
  int* disj;
  int* dm;
  pet_preds_t* pet_preds;
  double* Z = NULL;
  if(!isNull(Z_raw))
    Z = REAL(Z_raw);

  int max_ind;
  if(type == 0)
    max_ind = 1;
  else if(type == 1)
    max_ind = n_ind;
  else
    max_ind = n_models;

  SEXP ret = PROTECT(allocVector(REALSXP, N));
  double* prob_preds = REAL(ret);
  memset(prob_preds, 0, N * sizeof(double));

  int n_val = length(VECTOR_ELT(ensemble, 0));

  for(int i = 0; i < max_ind; i++) {
    inner_ensemble = VECTOR_ELT(ensemble, i);
    disj_raw = VECTOR_ELT(disjs, i);
    disj = INTEGER(disj_raw);
    for(real_n_conj = 0; real_n_conj < nrows(disj_raw); real_n_conj++) {
      if(disj[real_n_conj] == NA_INTEGER)
        break;
    }

    dm = getDesignMatrixIntern(X, N, disj, nrows(disj_raw), ncols(disj_raw), real_n_conj);

    for(int j = 0; j < n_val; j++) {
      pet = VECTOR_ELT(inner_ensemble, j);
      rebuild_tree(pet);
      tree = (node*) R_ExternalPtrAddr(VECTOR_ELT(pet, 5));
      pet_preds = predictIntern(tree, dm, Z, N, 0, leaves);
      for(int k = 0; k < N; k++) {
        prob_preds[k] += (pet_preds->prob_preds)[k];
      }
      Free(pet_preds->prob_preds);
      Free(pet_preds);
    }
    Free(dm);
  }

  for(int k = 0; k < N; k++) {
    prob_preds[k] /= max_ind * n_val;
  }

  UNPROTECT(1);
  return ret;
}

sa_eval_t* simulatedAnnealingStep(SEXP X_train, SEXP y_train, int max_vars, int max_conj, SEXP Z_train, SEXP Z_val, int* disj, int n_conj_raw, int n_vars_raw, double t, int acc_type, double score, int nodesize, double cp, int smoothing, int mtry, int covariable_mode, int scoring_rule, SEXP X_val, SEXP y_val, int use_validation, int y_bin, int allow_conj_removal, int conjsize, SEXP X, eval_models_list** models) {
  int p = ncols(VECTOR_ELT(X_train, 0));

  int* disj2 = (int*) Calloc(n_conj_raw * n_vars_raw, int);
  memcpy(disj2, disj, n_conj_raw * n_vars_raw * sizeof(int));

  int n_conj;
  for(n_conj = 0; n_conj < n_conj_raw; n_conj++) {
    if(disj2[n_conj] == NA_INTEGER)
      break;
  }
  int* n_vars = (int*) Calloc(n_conj, int);
  int n_vars_total = 0;
  int i, j;
  for(i = 0; i < n_conj; i++) {
    n_vars[i] = 0;
    for(j = 0; j < n_vars_raw; j++) {
      if(disj2[j * n_conj_raw + i] == NA_INTEGER)
        break;
      n_vars[i]++;
    }
    n_vars_total += n_vars[i];
  }

  GetRNGstate();

  int poss_add_moves = 2 * p * (n_conj < max_conj) * (n_vars_total < max_vars);
  int poss_rem_moves = n_conj * (n_conj > 1) * allow_conj_removal;

  int* poss_mod_moves = (int*) Calloc(n_conj, int);
  int poss_mod_moves_count = 0;

  for(i = 0; i < n_conj; i++) {
    poss_mod_moves[i] = 2 * (p - n_vars[i]) * (n_vars_total < max_vars) + n_vars[i] * (n_vars[i] > 1) + n_vars[i] * (2 * (p - n_vars[i]) + 1);
    poss_mod_moves_count += poss_mod_moves[i];
  }

  int total_moves = poss_add_moves + poss_rem_moves + poss_mod_moves_count;

  double add_prob = (double) poss_add_moves/total_moves;
  double rem_prob = (double) poss_rem_moves/total_moves;
  // double mod_prob = (double) poss_mod_moves_count/total_moves;

  int main_move;
  int mod_move;
  double rnd = unif_rand();
  if(rnd < add_prob)
    main_move = 0;
  else if(rnd < add_prob + rem_prob)
    main_move = 1;
  else
    main_move = 2;

  int which_conj = 0;

  if(main_move == 2) {
    double current_prob = 0.0;
    rnd = unif_rand();
    for(i = 0; i < n_conj; i++) {
      current_prob += (double) poss_mod_moves[i]/poss_mod_moves_count;
      if(rnd < current_prob) {
        which_conj = i;
        break;
      }
    }

    int n_vars_here = n_vars[which_conj];
    int poss_mod_add_moves = 2 * (p - n_vars_here) * (n_vars_total < max_vars);
    int poss_mod_rem_moves = n_vars_here * (n_vars_here > 1);
    int poss_mod_mod_moves = n_vars_here * (2 * (p - n_vars_here) + 1);
    int total_mod_moves = poss_mod_add_moves + poss_mod_rem_moves + poss_mod_mod_moves;

    double mod_add_prob = (double) poss_mod_add_moves/total_mod_moves;
    double mod_rem_prob = (double) poss_mod_rem_moves/total_mod_moves;
    // double mod_mod_prob = (double) poss_mod_mod_moves/total_mod_moves;

    rnd = unif_rand();
    if(rnd < mod_add_prob)
      mod_move = 0;
    else if(rnd < mod_add_prob + mod_rem_prob)
      mod_move = 1;
    else
      mod_move = 2;

    int* unused_vars = (int*) Calloc(p, int);
    memset(unused_vars, 0, p * sizeof(int));
    int* available_vars = (int*) Calloc(2 * p, int);
    for(i = 0; i < n_vars_here; i++) {
      unused_vars[abs(disj2[i * n_conj_raw + which_conj]) - 1] = 1;
    }
    j = 0;
    for(i = 0; i < p; i++) {
      if(!unused_vars[i]) {
        available_vars[j] = i+1;
        available_vars[j+1] = -(i+1);
        j += 2;
      }
    }

    if(mod_move == 2) {
      rnd = unif_rand();
      int which_var = rnd * n_vars_here;
      if(which_var == n_vars_here)
        which_var--;

      rnd = unif_rand();
      int swap_var = rnd * (j+1);
      if(swap_var == j+1)
        swap_var--;

      if(swap_var == j) {
        disj2[which_var * n_conj_raw + which_conj] = -disj2[which_var * n_conj_raw + which_conj];
      } else {
        disj2[which_var * n_conj_raw + which_conj] = available_vars[swap_var];
      }
    } else if(mod_move == 0) {
      rnd = unif_rand();
      int swap_var = rnd * j;
      if(swap_var == j)
        swap_var--;

      disj2[n_vars_here * n_conj_raw + which_conj] = available_vars[swap_var];
    } else {
      rnd = unif_rand();
      int which_var = rnd * n_vars_here;
      if(which_var == n_vars_here)
        which_var--;

      disj2[which_var * n_conj_raw + which_conj] = disj2[(n_vars_here-1) * n_conj_raw + which_conj];
      disj2[(n_vars_here-1) * n_conj_raw + which_conj] = NA_INTEGER;
    }

    Free(unused_vars);
    Free(available_vars);

  } else if (main_move == 0) {
    rnd = unif_rand();
    int swap_var = rnd * (2*p);
    if(swap_var == 2*p)
      swap_var--;

    disj2[n_conj] = 1 + (swap_var/2);
    if(swap_var % 2 == 1)
      disj2[n_conj] = -disj2[n_conj];
  } else {
    rnd = unif_rand();
    int which_rem = rnd * n_conj;
    if(which_rem == n_conj)
      which_rem--;

    for(i = 0; i < n_vars_raw; i++) {
      disj2[i * n_conj_raw + which_rem] = disj2[i * n_conj_raw + (n_conj-1)];
      disj2[i * n_conj_raw + (n_conj-1)] = NA_INTEGER;
    }
  }

  PutRNGstate();

  Free(n_vars);
  Free(poss_mod_moves);

  if(main_move == 2 && (mod_move == 2 || mod_move == 0)) {
    int* sub_disj = (int*) Calloc(n_vars_raw, int);
    for(i = 0; i < n_vars_raw; i++) {
      sub_disj[i] = disj2[i * n_conj_raw + which_conj];
    }
    int* dm = getDesignMatrixIntern(INTEGER(X), nrows(X), sub_disj, 1, n_vars_raw, 1);
    int conjsum = 0;
    for(i = 0; i < nrows(X); i++) {
      conjsum += dm[i];
    }
    Free(dm);
    Free(sub_disj);

    if (conjsum < conjsize || conjsum > nrows(X) - conjsize) {
      Free(disj2);
      return simulatedAnnealingStep(X_train, y_train, max_vars, max_conj, Z_train, Z_val, disj, n_conj_raw, n_vars_raw, t, acc_type, score, nodesize, cp, smoothing, mtry, covariable_mode, scoring_rule, X_val, y_val, use_validation, y_bin, allow_conj_removal, conjsize, X, models);
    }
  }

  sa_eval_t* eval = evaluateModel(X_train, y_train, Z_train, Z_val, t, acc_type, disj2, n_conj_raw, n_vars_raw, score, nodesize, cp, smoothing, mtry, covariable_mode, scoring_rule,
                                   X_val, y_val, use_validation, y_bin, models);
  if(!(eval->acc))
    memcpy(disj2, disj, n_conj_raw * n_vars_raw * sizeof(int));
  eval->disj = disj2;
  return eval;
}

SEXP simulatedAnnealing_(SEXP X_train, SEXP y_train, SEXP max_vars_raw, SEXP max_conj_raw, SEXP Z_train, SEXP Z_val, SEXP disj, SEXP t_raw, SEXP score, SEXP nodesize_raw, SEXP cp_raw, SEXP smoothing_raw, SEXP mtry_raw, SEXP covariable_mode_raw, SEXP scoring_rule_raw, SEXP X_val, SEXP y_val, SEXP use_validation_raw, SEXP y_bin_raw, SEXP allow_conj_removal_raw, SEXP conjsize_raw, SEXP X, SEXP cooling_schedule) {
  int lower_temp;
  double t = asReal(t_raw);
  double end_temp = asReal(getListElement(cooling_schedule, "real_end_temp"));
  int markov_iter = asInteger(getListElement(cooling_schedule, "markov_iter"));
  double markov_leave_frac = asReal(getListElement(cooling_schedule, "markov_leave_frac"));
  int acc_type = asInteger(getListElement(cooling_schedule, "acc_type2"));
  int frozen_def = asInteger(getListElement(cooling_schedule, "frozen_def2"));
  double frozen_acc_frac = asReal(getListElement(cooling_schedule, "frozen_acc_frac"));
  int frozen_markov_count = asInteger(getListElement(cooling_schedule, "frozen_markov_count"));
  int frozen_markov_mode = asInteger(getListElement(cooling_schedule, "frozen_markov_mode2"));
  int cooling_type = asInteger(getListElement(cooling_schedule, "type2"));
  double lambda = asReal(getListElement(cooling_schedule, "lambda"));
  double step_temp = asReal(getListElement(cooling_schedule, "step_temp"));
  int remember_models = asLogical(getListElement(cooling_schedule, "remember_models"));
  int print_iter = asInteger(getListElement(cooling_schedule, "print_iter"));
  //SEXP eval, disj2, min_conj, t_package;
  double min_score = asReal(score);
  int* current_acc = (int*) Calloc(markov_iter, int);
  double* current_scores = (double*) Calloc(markov_iter, double);
  int i, j, acc_sum, real_acc_sum;
  int total_iter = 0;
  int frozen = 0;
  //SEXP current_score_package;
  double current_score, score_sum, score_mean;
  current_score = min_score;
  double score_sd = 0.0;
  double real_acc_ratio = 0.0;
  int n_conj = nrows(disj);
  int n_vars = ncols(disj);
  int max_vars = asInteger(max_vars_raw);
  int max_conj = asInteger(max_conj_raw);
  int nodesize = asInteger(nodesize_raw);
  double cp = asReal(cp_raw);
  int smoothing = asInteger(smoothing_raw);
  int mtry = asInteger(mtry_raw);
  int covariable_mode = asInteger(covariable_mode_raw);
  int scoring_rule = asInteger(scoring_rule_raw);
  int use_validation = asLogical(use_validation_raw);
  int y_bin = asLogical(y_bin_raw);
  int allow_conj_removal = asLogical(allow_conj_removal_raw);
  int conjsize = asInteger(conjsize_raw);

  int* disj2 = (int*) Calloc(n_conj * n_vars, int);
  memcpy(disj2, INTEGER(disj), n_conj * n_vars * sizeof(int));
  int* min_conj = disj2;

  // int* disj_buffer = disj;
  int protect_min_conj = 1;

  sa_eval_t* eval;

  eval_models_list** models = NULL;
  if(remember_models) models = (eval_models_list**) Calloc(PRIME, eval_models_list*);
  prevented_evals = 0;

  while(t >= end_temp) {
    lower_temp = 0;
    acc_sum = 0;
    real_acc_sum = 0;
    score_sum = 0.0;
    for(i = 0; i < markov_iter; i++) {
      eval = simulatedAnnealingStep(X_train, y_train, max_vars, max_conj, Z_train, Z_val, disj2, n_conj, n_vars, t, acc_type, current_score, nodesize, cp, smoothing, mtry, covariable_mode, scoring_rule, X_val, y_val, use_validation, y_bin, allow_conj_removal, conjsize, X, models);
      if(!protect_min_conj)
        Free(disj2);
      else
        protect_min_conj = 0;
      disj2 = eval->disj;
      current_score = eval->score;
      current_acc[i] = eval->acc;
      Free(eval);

      current_scores[i] = current_score;

      if (current_score <= min_score) {
        Free(min_conj);
        min_score = current_score;
        min_conj = disj2;
        protect_min_conj = 1;
      }

      acc_sum += current_acc[i];
      if(i > 0 && current_acc[i])
        real_acc_sum += !doubleEquals(current_scores[i] - current_scores[i-1], 0);
      score_sum += current_score;

      if(acc_sum > markov_leave_frac * markov_iter) {
        lower_temp = 1;
      }

      if (i == markov_iter - 1 || lower_temp) {
        score_mean = score_sum/(i+1);
        score_sd = 0.0;
        for(j = 0; j < i+1; j++) {
          score_sd += (current_scores[j] - score_mean) * (current_scores[j] - score_mean);
        }
        score_sd = sqrt(score_sd/i);
        real_acc_ratio = (double) real_acc_sum/(i+1);
        if(print_iter > 0) {
          Rprintf("log10(t)=%5.2f, i=%4d, acc/ratio=%.2f, real acc/ratio=%.2f, score/sd=%10.6f\n", log10(t), i + 1, (double) acc_sum/(i+1),
                  real_acc_ratio, score_sd);
        }
        break;
      }
    }

    total_iter += i + 1;

    if ((frozen_def && score_sd < 10e-10) || (!frozen_def && real_acc_ratio < frozen_acc_frac)) {
      frozen++;
      if (frozen >= frozen_markov_count)
        break;
    }
    else if (frozen_markov_mode == 1)
      frozen = 0;

    if(cooling_type == 0) {
      // Adaptive cooling
      if(score_sd >= 10e-10) {
        t *= exp(-lambda * t/score_sd);
      }
    } else {
      // Geometric cooling
      t /= pow(10, step_temp);
    }
  }

  destroy_eval_models(models);

  if(!protect_min_conj)
    Free(disj2);

  Free(current_acc);
  Free(current_scores);

  SEXP ret_obj = PROTECT(allocVector(VECSXP, 4));
  SEXP min_conj_R = PROTECT(allocMatrix(INTSXP, n_conj, n_vars));
  memcpy(INTEGER(min_conj_R), min_conj, n_conj * n_vars * sizeof(int));
  SET_VECTOR_ELT(ret_obj, 0, min_conj_R);
  SET_VECTOR_ELT(ret_obj, 1, ScalarReal(min_score));
  SET_VECTOR_ELT(ret_obj, 2, ScalarInteger(total_iter));
  SET_VECTOR_ELT(ret_obj, 3, ScalarInteger(prevented_evals));

  Free(min_conj);
  UNPROTECT(2);
  return ret_obj;
}

gs_eval_t* greedyStep(SEXP X_train, SEXP y_train, int max_vars, int max_conj, SEXP Z_train, SEXP Z_val, int* disj, int n_conj_raw, int n_vars_raw, double score, int mtry_greedy, int greedy_mod, int greedy_rem, int nodesize, double cp, int smoothing, int mtry, int covariable_mode, int scoring_rule, SEXP X_val, SEXP y_val, int use_validation, int y_bin, int allow_conj_removal, int conjsize, SEXP X) {
  int p = ncols(VECTOR_ELT(X_train, 0));

  int* disj2 = (int*) Calloc(n_conj_raw * n_vars_raw, int);
  memcpy(disj2, disj, n_conj_raw * n_vars_raw * sizeof(int));
  int* min_conj = (int*) Calloc(n_conj_raw * n_vars_raw, int);
  memcpy(min_conj, disj, n_conj_raw * n_vars_raw * sizeof(int));

  int n_conj;
  for(n_conj = 0; n_conj < n_conj_raw; n_conj++) {
    if(disj2[n_conj] == NA_INTEGER)
      break;
  }
  int* n_vars = (int*) Calloc(n_conj, int);
  int n_vars_total = 0;
  int i, j, k, l;
  for(i = 0; i < n_conj; i++) {
    n_vars[i] = 0;
    for(j = 0; j < n_vars_raw; j++) {
      if(disj2[j * n_conj_raw + i] == NA_INTEGER)
        break;
      n_vars[i]++;
    }
    n_vars_total += n_vars[i];
  }

  double current_score;
  double min_score = score;
  int iter = 0;
  int move_type = -1;

  double mtry_factor = 1.0;
  if(mtry_greedy) {
    int poss_moves = p * (n_conj < max_conj) * (n_vars_total < max_vars);
    poss_moves += n_conj * (n_conj > 1) * allow_conj_removal * greedy_rem;
    for(i = 0; i < n_conj; i++)
      poss_moves += 2 * (p - n_vars[i]) * (n_vars_total < max_vars) + n_vars[i] * (n_vars[i] > 1) * greedy_rem + n_vars[i] * (2 * (p - n_vars[i]) + 1) * greedy_mod;
    mtry_factor = sqrt(poss_moves)/poss_moves;
    GetRNGstate();
  }
  int mtry_result;

  // Add
  if (n_conj < max_conj && n_vars_total < max_vars) {
    for(j = 0; j < p * mtry_factor; j++) {
      mtry_result = drawNumberWithReplacement(p, j, mtry_greedy);

      disj2[n_conj] = mtry_result + 1;
      iter++;
      current_score = buildModel(X_train, y_train, Z_train, Z_val, disj2, n_conj_raw, n_vars_raw, nodesize, cp, smoothing, mtry, covariable_mode, scoring_rule, X_val, y_val, use_validation, y_bin, NULL);
      if(current_score < min_score) {
        min_score = current_score;
        memcpy(min_conj, disj2, n_conj_raw * n_vars_raw * sizeof(int));
        move_type = 0;
      }
    }
    disj2[n_conj] = NA_INTEGER;
  }

  // Remove
  if (n_conj > 1 && allow_conj_removal && greedy_rem) {
    for(j = 0; j < n_conj * mtry_factor; j++) {
      mtry_result = drawNumberWithReplacement(n_conj, j, mtry_greedy);

      memcpy(disj2, disj, n_conj_raw * n_vars_raw * sizeof(int));
      for(i = 0; i < n_vars_raw; i++) {
        disj2[i * n_conj_raw + mtry_result] = disj2[i * n_conj_raw + (n_conj-1)];
        disj2[i * n_conj_raw + (n_conj-1)] = NA_INTEGER;
      }
      iter++;
      current_score = buildModel(X_train, y_train, Z_train, Z_val, disj2, n_conj_raw, n_vars_raw, nodesize, cp, smoothing, mtry, covariable_mode, scoring_rule, X_val, y_val, use_validation, y_bin, NULL);
      if(current_score < min_score) {
        min_score = current_score;
        memcpy(min_conj, disj2, n_conj_raw * n_vars_raw * sizeof(int));
        move_type = 1;
      }
    }
  }

  int* unused_vars = (int*) Calloc(p, int);
  int* available_vars = (int*) Calloc(2 * p, int);
  int n_free_vars;

  int* sub_disj = (int*) Calloc(n_vars_raw, int);
  int* dm;
  int conjsum;

  memcpy(disj2, disj, n_conj_raw * n_vars_raw * sizeof(int));
  int conj_buffer, n_vars_here;

  // Modify
  for(i = 0; i < n_conj; i++) {
    n_vars_here = n_vars[i];

    memset(unused_vars, 0, p * sizeof(int));
    for(j = 0; j < n_vars_here; j++) {
      unused_vars[abs(disj[j * n_conj_raw + i]) - 1] = 1;
    }
    n_free_vars = 0;
    for(j = 0; j < p; j++) {
      if(!unused_vars[j]) {
        available_vars[n_free_vars] = j+1;
        available_vars[n_free_vars+1] = -(j+1);
        n_free_vars += 2;
      }
    }

    // Add
    if (n_vars_total < max_vars) {
      for(j = 0; j < n_free_vars * mtry_factor; j++) {
        mtry_result = drawNumberWithReplacement(n_free_vars, j, mtry_greedy);

        disj2[n_vars_here * n_conj_raw + i] = available_vars[mtry_result];

        for(k = 0; k < n_vars_raw; k++) {
          sub_disj[k] = disj2[k * n_conj_raw + i];
        }
        dm = getDesignMatrixIntern(INTEGER(X), nrows(X), sub_disj, 1, n_vars_raw, 1);
        conjsum = 0;
        for(k = 0; k < nrows(X); k++) {
          conjsum += dm[k];
        }
        Free(dm);
        if (conjsum < conjsize || conjsum > nrows(X) - conjsize) {
          continue;
        }

        iter++;
        current_score = buildModel(X_train, y_train, Z_train, Z_val, disj2, n_conj_raw, n_vars_raw, nodesize, cp, smoothing, mtry, covariable_mode, scoring_rule, X_val, y_val, use_validation, y_bin, NULL);
        if(current_score < min_score) {
          min_score = current_score;
          memcpy(min_conj, disj2, n_conj_raw * n_vars_raw * sizeof(int));
          move_type = 2;
        }
      }
      disj2[n_vars_here * n_conj_raw + i] = NA_INTEGER;
    }

    // Remove
    if (n_vars_here > 1 && greedy_rem) {
      for(j = 0; j < n_vars_here * mtry_factor; j++) {
        mtry_result = drawNumberWithReplacement(n_vars_here, j, mtry_greedy);

        conj_buffer = disj2[mtry_result * n_conj_raw + i];
        disj2[mtry_result * n_conj_raw + i] = disj2[(n_vars_here-1) * n_conj_raw + i];
        disj2[(n_vars_here-1) * n_conj_raw + i] = NA_INTEGER;

        iter++;
        current_score = buildModel(X_train, y_train, Z_train, Z_val, disj2, n_conj_raw, n_vars_raw, nodesize, cp, smoothing, mtry, covariable_mode, scoring_rule, X_val, y_val, use_validation, y_bin, NULL);
        if(current_score < min_score) {
          min_score = current_score;
          memcpy(min_conj, disj2, n_conj_raw * n_vars_raw * sizeof(int));
          move_type = 3;
        }

        disj2[(n_vars_here-1) * n_conj_raw + i] = disj2[mtry_result * n_conj_raw + i];
        disj2[mtry_result * n_conj_raw + i] = conj_buffer;
      }
    }

    // Modify
    if(greedy_mod) {
      for(j = 0; j < n_vars_here; j++) {
        conj_buffer = disj2[j * n_conj_raw + i];

        for(l = 0; l < (n_free_vars+1) * mtry_factor; l++) {
          mtry_result = drawNumberWithReplacement(n_free_vars+1, l, mtry_greedy);

          if(mtry_result == n_free_vars) {
            disj2[j * n_conj_raw + i] = -disj2[j * n_conj_raw + i];
          } else {
            disj2[j * n_conj_raw + i] = available_vars[mtry_result];
          }

          for(k = 0; k < n_vars_raw; k++) {
            sub_disj[k] = disj2[k * n_conj_raw + i];
          }
          dm = getDesignMatrixIntern(INTEGER(X), nrows(X), sub_disj, 1, n_vars_raw, 1);
          conjsum = 0;
          for(k = 0; k < nrows(X); k++) {
            conjsum += dm[k];
          }
          Free(dm);
          if (conjsum < conjsize || conjsum > nrows(X) - conjsize) {
            continue;
          }

          iter++;
          current_score = buildModel(X_train, y_train, Z_train, Z_val, disj2, n_conj_raw, n_vars_raw, nodesize, cp, smoothing, mtry, covariable_mode, scoring_rule, X_val, y_val, use_validation, y_bin, NULL);
          if(current_score < min_score) {
            min_score = current_score;
            memcpy(min_conj, disj2, n_conj_raw * n_vars_raw * sizeof(int));
            move_type = 4;
          }
        }
        disj2[j * n_conj_raw + i] = conj_buffer;
      }
    }
  }

  if(mtry_greedy)
    PutRNGstate();

  Free(disj2);
  Free(n_vars);
  Free(unused_vars);
  Free(available_vars);
  Free(sub_disj);

  gs_eval_t* eval = (gs_eval_t*) Calloc(1, gs_eval_t);
  eval->disj = min_conj;
  eval->score = min_score;
  eval->iter = iter;
  eval->move_type = move_type;
  return eval;
}

SEXP greedySearch_(SEXP X_train, SEXP y_train, SEXP max_vars_raw, SEXP max_conj_raw, SEXP Z_train, SEXP Z_val, SEXP disj_raw, SEXP score, SEXP mtry_greedy_raw, SEXP greedy_mod_raw, SEXP greedy_rem_raw, SEXP nodesize_raw, SEXP cp_raw, SEXP smoothing_raw, SEXP mtry_raw, SEXP covariable_mode_raw, SEXP scoring_rule_raw, SEXP X_val, SEXP y_val, SEXP use_validation_raw, SEXP y_bin_raw, SEXP allow_conj_removal_raw, SEXP conjsize_raw, SEXP X) {
  double min_score = asReal(score);
  int total_iter = 0;
  int n_conj = nrows(disj_raw);
  int n_vars = ncols(disj_raw);
  int max_vars = asInteger(max_vars_raw);
  int max_conj = asInteger(max_conj_raw);
  int nodesize = asInteger(nodesize_raw);
  double cp = asReal(cp_raw);
  int smoothing = asInteger(smoothing_raw);
  int mtry = asInteger(mtry_raw);
  int covariable_mode = asInteger(covariable_mode_raw);
  int scoring_rule = asInteger(scoring_rule_raw);
  int use_validation = asLogical(use_validation_raw);
  int y_bin = asLogical(y_bin_raw);
  int allow_conj_removal = asLogical(allow_conj_removal_raw);
  int conjsize = asInteger(conjsize_raw);
  int* disj = INTEGER(disj_raw);
  int* min_conj = (int*) Calloc(n_conj * n_vars, int);
  memcpy(min_conj, disj, n_conj * n_vars * sizeof(int));

  int mtry_greedy = asLogical(mtry_greedy_raw);
  if(mtry_greedy) {
    int p = ncols(VECTOR_ELT(X_train, 0));
    int max_mtry_vars = max_conj * 2 * max_vars * p;
    mtry_vars_raw = (int*) Calloc(max_mtry_vars, int);
    mtry_vars = (int*) Calloc(max_mtry_vars, int);
    for(int i = 0; i < max_mtry_vars; i++)
      mtry_vars_raw[i] = i;
  }

  int move_counts[5] = { 0 };
  int greedy_mod = asLogical(greedy_mod_raw);
  int greedy_rem = asLogical(greedy_rem_raw);

  while(1) {
    gs_eval_t* eval = greedyStep(X_train, y_train, max_vars, max_conj, Z_train, Z_val, min_conj, n_conj, n_vars, min_score, mtry_greedy, greedy_mod, greedy_rem, nodesize, cp, smoothing, mtry, covariable_mode, scoring_rule, X_val, y_val, use_validation, y_bin, allow_conj_removal, conjsize, X);
    total_iter += eval->iter;
    // if (eval->score >= min_score) {
    // if (doubleEquals(min_score, eval->score) || eval->score > min_score) {
    if (eval->move_type == -1) {
      Free(eval->disj);
      Free(eval);
      break;
    }
    Free(min_conj);
    min_conj = eval->disj;
    min_score = eval->score;
    move_counts[eval->move_type]++;
    Free(eval);
  }

  if(mtry_greedy) {
    Free(mtry_vars_raw);
    Free(mtry_vars);
  }

  SEXP ret_obj = PROTECT(allocVector(VECSXP, 4));
  SEXP min_conj_R = PROTECT(allocMatrix(INTSXP, n_conj, n_vars));
  memcpy(INTEGER(min_conj_R), min_conj, n_conj * n_vars * sizeof(int));
  SET_VECTOR_ELT(ret_obj, 0, min_conj_R);
  SET_VECTOR_ELT(ret_obj, 1, ScalarReal(min_score));
  SET_VECTOR_ELT(ret_obj, 2, ScalarInteger(total_iter));
  SEXP move_counts_R = PROTECT(allocVector(INTSXP, 5));
  memcpy(INTEGER(move_counts_R), move_counts, 5 * sizeof(int));
  SET_VECTOR_ELT(ret_obj, 3, move_counts_R);

  Free(min_conj);
  UNPROTECT(3);
  return ret_obj;
}

SEXP getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);

  for (int i = 0; i < length(list); i++) {
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  }
  return elmt;
}

int calcDisjHash(int* disj, int n_conj, int n_vars, int real_n_conj) {
  int hash_sum = 0;
  for(int i = 0; i < real_n_conj; i++) {
    for(int j = 0; j < n_vars; j++) {
      if(disj[j * n_conj + i] == NA_INTEGER)
        break;
      hash_sum += disj[j * n_conj + i];
    }
  }
  return abs(hash_sum) % PRIME;
}

int calcDisjHash2(int* disj, int total_length) {
  char str[STRLENGTH];
  int index = 0;
  int chars_left;
  for(int i = 0; i < total_length; i++) {
    if(disj[i] == NA_INTEGER) continue;
    chars_left = STRLENGTH - index;
    if(chars_left < 2) break;
    index += snprintf(str + index, chars_left, "%d", abs(disj[i]));
  }
  str[DIGITS] = '\0';
  return atoi(str);
}

int cmp_disj(int* disj1, int* disj2, int total_length) {
  for(int i = 0; i < total_length; i++) {
    if(disj1[i] < disj2[i]) return -1;
    else if(disj1[i] > disj2[i]) return 1;
  }
  return 0;
}

int cmp_disj_fixed(int* disj1, int* disj2) {
  return cmp_disj(disj1, disj2, total_length_global);
}

void destroy_eval_models(eval_models_list** models) {
  if(models == NULL) return;
  eval_models_list* current;
  eval_models_list* next;
  for(int i = 0; i < PRIME; i++) {
    next = models[i];
    while(next != NULL) {
      current = next;
      next = next->next;
      Free(current->disj);
      Free(current);
    }
  }
  Free(models);
}

int drawNumberWithReplacement(int total, int iter, int random) {
  if(!random)
    return iter;
  if(!iter)
    memcpy(mtry_vars, mtry_vars_raw, total * sizeof(int));
  double rnd = unif_rand();
  int mtry_which = rnd * (total - iter);
  int mtry_result = mtry_vars[mtry_which];
  for(int mtry_i = mtry_which; mtry_i < total - iter - 1; mtry_i++)
    mtry_vars[mtry_i] = mtry_vars[mtry_i+1];
  return mtry_result;
}

int cmp_gen_score(const void* value1, const void* value2) {
  // Decreasing
  return -cmp_double_direct((*(gen_t*)value1).score, (*(gen_t*)value2).score);
}

int cmp_gen_fitness(const void* value1, const void* value2) {
  // Decreasing
  return -cmp_double_direct((*(gen_t*)value1).fitness, (*(gen_t*)value2).fitness);
}

int cmp_gen_conj(const void* value1, const void* value2) {
  return cmp_disj_fixed((*(gen_t*)value1).disj, (*(gen_t*)value2).disj);
}



