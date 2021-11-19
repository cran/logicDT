#include "logicDT.h"
#include <Rmath.h>

#define PRIME 10000000
#define STRLENGTH 8
#define DIGITS 7

extern int prevented_evals;
extern int* mtry_vars_raw;
extern int* mtry_vars;
extern int total_length_global;

typedef struct _sa_eval
{
  int* disj;
  double score;
  int acc;
} sa_eval_t;

typedef struct _gs_eval
{
  int* disj;
  double score;
  int iter;
  int move_type;
} gs_eval_t;

typedef struct _eval_models_list
{
  int* disj;
  double score;
  struct _eval_models_list *next;
} eval_models_list;

typedef struct _gen
{
  int* disj;
  double score;
  double fitness;
} gen_t;

typedef struct _gp_eval
{
  gen_t* generation;
  int n_ind;
  int iter;
} gp_eval_t;

double buildModel(SEXP X_train, SEXP y_train, SEXP Z_train, SEXP Z_val, int* disj, int n_conj, int n_vars, int nodesize, double cp, int smoothing, int mtry, int covariable_mode, int scoring_rule, SEXP X_val, SEXP y_val, int use_validation, int y_bin, eval_models_list** models);
sa_eval_t* evaluateModel(SEXP X_train, SEXP y_train, SEXP Z_train, SEXP Z_val, double t, int acc_type, int* disj, int n_conj, int n_vars, double old_score, int nodesize, double cp, int smoothing, int mtry, int covariable_mode, int scoring_rule, SEXP X_val, SEXP y_val, int use_validation, int y_bin, eval_models_list** models);
double getGPScore(double score, int y_bin, int scoring_rule);
void calcSharedFitness(gen_t* generation, int n_ind, int max_conj, int max_vars, double sigma);
int* tournamentSelection(gen_t* generation, int n_ind, int k);
int* mutateGeneration(gen_t* generation, int n_ind, int max_vars, int max_conj, int p, int type, int which_ind, int allow_conj_removal);
gp_eval_t* geneticProgrammingStep(SEXP X_train, SEXP y_train, int max_vars, int max_conj, SEXP Z_train, SEXP Z_val, gen_t* generation, int n_ind, double best_score, int nodesize, double cp, int smoothing, int mtry, int covariable_mode, int scoring_rule, SEXP X_val, SEXP y_val, int use_validation, int y_bin, int allow_conj_removal, int conjsize, SEXP X, eval_models_list** models);
SEXP geneticProgramming_(SEXP X_train, SEXP y_train, SEXP max_vars_raw, SEXP max_conj_raw, SEXP max_gen_raw, SEXP gp_sigma_raw, SEXP gp_fs_interval_raw, SEXP Z_train, SEXP Z_val, SEXP nodesize_raw, SEXP cp_raw, SEXP smoothing_raw, SEXP mtry_raw, SEXP covariable_mode_raw, SEXP scoring_rule_raw, SEXP X_val, SEXP y_val, SEXP use_validation_raw, SEXP y_bin_raw, SEXP allow_conj_removal_raw, SEXP conjsize_raw, SEXP X);
SEXP predictGP_(SEXP model, SEXP X_raw, SEXP Z_raw, SEXP type_raw, SEXP n_models_raw, SEXP leaves_raw);
sa_eval_t* simulatedAnnealingStep(SEXP X_train, SEXP y_train, int max_vars, int max_conj, SEXP Z_train, SEXP Z_val, int* disj, int n_conj_raw, int n_vars_raw, double t, int acc_type, double score, int nodesize, double cp, int smoothing, int mtry, int covariable_mode, int scoring_rule, SEXP X_val, SEXP y_val, int use_validation, int y_bin, int allow_conj_removal, int conjsize, SEXP X, eval_models_list** models);
SEXP simulatedAnnealing_(SEXP X_train, SEXP y_train, SEXP max_vars_raw, SEXP max_conj_raw, SEXP Z_train, SEXP Z_val, SEXP disj, SEXP t_raw, SEXP score, SEXP nodesize_raw, SEXP cp_raw, SEXP smoothing_raw, SEXP mtry_raw, SEXP covariable_mode_raw, SEXP scoring_rule_raw, SEXP X_val, SEXP y_val, SEXP use_validation_raw, SEXP y_bin_raw, SEXP allow_conj_removal_raw, SEXP conjsize_raw, SEXP X, SEXP cooling_schedule);
gs_eval_t* greedyStep(SEXP X_train, SEXP y_train, int max_vars, int max_conj, SEXP Z_train, SEXP Z_val, int* disj, int n_conj_raw, int n_vars_raw, double score, int mtry_greedy, int greedy_mod, int greedy_rem, int nodesize, double cp, int smoothing, int mtry, int covariable_mode, int scoring_rule, SEXP X_val, SEXP y_val, int use_validation, int y_bin, int allow_conj_removal, int conjsize, SEXP X);
SEXP greedySearch_(SEXP X_train, SEXP y_train, SEXP max_vars_raw, SEXP max_conj_raw, SEXP Z_train, SEXP Z_val, SEXP disj_raw, SEXP score, SEXP mtry_greedy_raw, SEXP greedy_mod_raw, SEXP greedy_rem_raw, SEXP nodesize_raw, SEXP cp_raw, SEXP smoothing_raw, SEXP mtry_raw, SEXP covariable_mode_raw, SEXP scoring_rule_raw, SEXP X_val, SEXP y_val, SEXP use_validation_raw, SEXP y_bin_raw, SEXP allow_conj_removal_raw, SEXP conjsize_raw, SEXP X);
SEXP getListElement(SEXP list, const char *str);
int calcDisjHash(int* disj, int n_conj, int n_vars, int real_n_conj);
int calcDisjHash2(int* disj, int total_length);
int cmp_disj(int* disj1, int* disj2, int total_length);
int cmp_disj_fixed(int* disj1, int* disj2);
void destroy_eval_models(eval_models_list** models);
int drawNumberWithReplacement(int total, int iter, int random);
int cmp_gen_score(const void* value1, const void* value2);
int cmp_gen_fitness(const void* value1, const void* value2);
int cmp_gen_conj(const void* value1, const void* value2);

