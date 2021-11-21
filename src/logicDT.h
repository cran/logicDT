#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <stdlib.h>

#include <time.h>

#include <R_ext/Applic.h>

#define SUCCESS 0
#define ERR_INVAL 1
#define ERR_NOMEM 2
#define FALSE 0
#define TRUE 1

typedef struct logic_stack_s logic_stack_t;
typedef struct stack_frame_s stack_frame_t;
typedef struct _functional
{
  double b,c,d,e;
  int y_bin;
} functional;
typedef struct _node
{
  struct _node *left, *right;
  int leaf;
  int split;
  int split_bin_or_cont;
  double split_point;
  int* obs_ind;
  int N_k;
  double pred;
  functional* func_pred;
} node;
typedef struct _dataset
{
  int* bin_y;
  double* quant_y;
  double* Z;
  int* obs_ind;
  int N;
  double* par_scale;
} dataset;
typedef struct _linked_list
{
  int split;
  int split_bin_or_cont; // 0: Binary split, 1: Continuous split
  double split_point;
  double pred;
  struct _linked_list *next;
} linked_list;
struct stack_frame_s {
  struct stack_frame_s *next;
  void *data;
};
struct logic_stack_s {
  struct stack_frame_s *top;
};
typedef struct _pet
{
  int* splits;
  int* splits_bin_or_cont;
  double* split_points;
  double* preds;
  double* train_preds;
  node* tree;
  int number_of_nodes;
  functional** model_list;
  int y_bin;
} pet_t;
typedef struct _pet_preds
{
  int* class_preds;
  double* prob_preds;
} pet_preds_t;
typedef struct _pet_ensemble
{
  pet_t** pets;
  int n_pets;
  double score;
} pet_ensemble_t;
typedef struct _y_prob_pair
{
  int y;
  double prob;
} y_prob_pair_t;
typedef struct _y_Z_pair
{
  int bin_y;
  double quant_y;
  double Z;
} y_Z_pair_t;

int stack_destroy(logic_stack_t *stack);
int stack_empty(logic_stack_t *stack);
logic_stack_t *stack_new(void);
void *stack_pop(logic_stack_t *stack);
void *queue_pop(logic_stack_t *stack);
int stack_push(logic_stack_t *stack, void *data);

/*void split_node(node* knot, int split);*/
double calcLeafProb(int N_k_1, int N_k, int smoothing);
void make_leaf(node* knot, double p_k_1, double* train_preds);
void tree_destroy(node* tree);
// static void _finalizer(SEXP tree);
void pet_destroy(pet_t* pet, int destroy_tree);
void rebuild_tree(SEXP pet);

linked_list* set_values_and_next(linked_list* l, int split, int split_bin_or_cont, double split_point, double pred);

double calcDev(double* predictions, int* y, int N);
double calcNCE(double* predictions, int* y, int N);
double calcBrier(double* predictions, int* y, int N);
double calcMis(int* predictions, int* y, int N);
double calcLikelihood(double* predictions, int* y, int N);
SEXP calcAUC_(SEXP probs, SEXP y, SEXP y_sorted_raw);
double calcAUCSorted(double* predictions, int* y, int N);
double calcAUCUnsorted(double* predictions, int* y, int N);
double calcMSE(double* predictions, double* y, int N);
SEXP C_PET_TO_R_PET(pet_t* pet, int N);
SEXP fitPETs_(SEXP X_train_raw, SEXP y_train_raw, SEXP X_val_raw, SEXP y_val_raw, SEXP Z_train_raw, SEXP Z_val_raw, SEXP use_validation_raw, SEXP y_bin_raw, SEXP nodesize_raw, SEXP cp_raw, SEXP smoothing_raw, SEXP mtry_raw, SEXP covariable_mode_raw, SEXP disj_raw, SEXP real_n_conj_raw, SEXP scoring_rule_raw, SEXP return_full_model_raw);
pet_ensemble_t* fitPETsIntern(SEXP X_train_raw, SEXP y_train_raw, SEXP X_val_raw, SEXP y_val_raw, SEXP Z_train_raw, SEXP Z_val_raw, int use_validation, int y_bin, int nodesize, double cp, int smoothing, int mtry, int covariable_mode, int* disj, int n_conj, int n_vars, int real_n_conj, int scoring_rule, int return_full_model);
SEXP fitPET_(SEXP X_raw, SEXP y_raw, SEXP Z_raw, SEXP nodesize_raw, SEXP cp_raw, SEXP smoothing_raw, SEXP mtry_raw, SEXP covariable_mode_raw);
pet_t* fitPETIntern(int* X, int* bin_y, double* quant_y, int y_bin, double* Z, int N, int p, int pZ, int nodesize, double cp, int smoothing, int mtry, int covariable_mode);
SEXP predict_(SEXP pet, SEXP X_raw, SEXP Z_raw, SEXP type_raw, SEXP leaves_raw);
pet_preds_t* predictIntern(node* tree, int* X, double* Z, int N, int type, int leaves);
SEXP predictEnsemble_(SEXP ensemble, SEXP X_raw, SEXP Z_raw, SEXP type_raw, SEXP leaves_raw);
SEXP getDesignMatrix_(SEXP X_raw, SEXP disj_raw, SEXP real_n_conj_raw);
int* getDesignMatrixIntern(int* X, int N, int* disj, int n_conj, int n_vars, int real_n_conj);
/*int arrangeNAs(int* disj, int n_conj, int n_vars);*/
SEXP vim_permutation_(SEXP ensemble, SEXP X_val_raw, SEXP y_val_raw, SEXP Z_val_raw, SEXP permutation_raw, SEXP disj_raw, SEXP real_n_conj_raw, SEXP scoring_rule_raw, SEXP y_bin_raw, SEXP leaves_raw);

// #define [i,j]  [i*p+j]
double gini_decrease(double p_k_1, double p_L, double p_L_1, double p_R_1);
double mse_impurity(int N_k, double y_sum, double y_sum_2);
double mse_decrease(int N_k, int N_k_L, int N_k_R, double N_k_sum, double N_L_sum, double N_R_sum, double N_k_sum_2, double N_L_sum_2, double N_R_sum_2);

int doubleEquals(double a, double b);
int cmp_integer(const void* value1, const void* value2);
int cmp_double(const void* value1, const void* value2);
int cmp_integer_direct(int value1, int value2);
int cmp_double_direct(double value1, double value2);
int cmp_y_probs_int(const void* value1, const void* value2);
int cmp_y_probs_double(const void* value1, const void* value2);
int cmp_y_Z_pair(const void* value1, const void* value2);

functional** functionalLeaves(node* tree, int number_of_nodes, int* bin_y, double* quant_y, int y_bin, double* Z);
double binLogLikelihood(int n, double* par, void* ex);
void binLogLikelihoodGrad(int n, double* par, double* gr, void* ex);
double squaredError(int n, double* par, void* ex);
void squaredErrorGrad(int n, double* par, double* gr, void* ex);
functional* fit4plModel(int* bin_y, double* quant_y, int y_bin, double y_mean, double* Z, int N, int* obs_ind);
SEXP fit4plModel_(SEXP y, SEXP Z);
double eval4plModel(functional* func_pred, double Z);
double* fitLinearModel(double* x, double* y, int N);


