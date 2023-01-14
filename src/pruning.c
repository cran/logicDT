#include "pruning.h"

static void _finalizer(SEXP tree) {
  if (R_ExternalPtrAddr(tree) == NULL)
    return;
  // Rprintf("finalizing\n");
  node* ptr = (node*) R_ExternalPtrAddr(tree);
  tree_destroy(ptr);
  R_ClearExternalPtr(tree);
}

SEXP prune_(SEXP pet, SEXP y_raw, SEXP Z_raw) {
  node* tree = (node*) R_ExternalPtrAddr(VECTOR_ELT(pet, 5));

  int* bin_y = NULL; double* quant_y = NULL;
  int y_bin;
  if(isInteger(y_raw)) {
    y_bin = 1;
    bin_y = INTEGER(y_raw);
  } else {
    y_bin = 0;
    quant_y = REAL(y_raw);
  }
  double* Z = NULL;
  if(!isNull(Z_raw)) {
    Z = REAL(Z_raw);
  }

  int number_of_nodes = length(VECTOR_ELT(pet, 0));
  int covariable_mode = asInteger(VECTOR_ELT(pet, 8));
  if(covariable_mode >= 2) {
    functional** model_list = functionalLeaves(tree, number_of_nodes, bin_y, quant_y, y_bin, Z, covariable_mode, 0, 1);
    Free(model_list);
  }

  linked_list2* prune_list = prune(tree);

  int n_prunes = 1;
  linked_list2* current_list_item = prune_list;
  while(current_list_item->next != NULL) {
    n_prunes++;
    current_list_item = current_list_item->next;
  }

  SEXP ret = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(ret, 0, allocVector(REALSXP, n_prunes));
  SET_VECTOR_ELT(ret, 1, allocVector(VECSXP, n_prunes));
  double* alphas = REAL(VECTOR_ELT(ret, 0));
  SEXP pets = VECTOR_ELT(ret, 1);

  SEXP current_pet;
  current_list_item = prune_list;
  linked_list2* buffer;

  node* current_node;
  logic_stack_t *stack = stack_new();
  SEXP model_list_R = NULL;
  int* splits_pointer; int* splits_bin_or_cont_pointer; double* split_points_pointer; double* preds_pointer;
  double* bcde;
  SEXP func_pred_R, bcde_R;
  functional* current_func;
  int j;

  for(int i = 0; i < n_prunes; i++) {
    alphas[i] = current_list_item->alpha;
    SET_VECTOR_ELT(pets, i, shallow_duplicate(pet));
    current_pet = VECTOR_ELT(pets, i);

    SET_VECTOR_ELT(current_pet, 5, R_MakeExternalPtr(current_list_item->tree, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(VECTOR_ELT(current_pet, 5), _finalizer, TRUE);

    stack_push(stack, current_list_item->tree);
    number_of_nodes = getNumberOfNodes(current_list_item->tree);

    SET_VECTOR_ELT(current_pet, 0, allocVector(INTSXP, number_of_nodes));
    SET_VECTOR_ELT(current_pet, 1, allocVector(INTSXP, number_of_nodes));
    SET_VECTOR_ELT(current_pet, 2, allocVector(REALSXP, number_of_nodes));
    SET_VECTOR_ELT(current_pet, 3, allocVector(REALSXP, number_of_nodes));

    splits_pointer = INTEGER(VECTOR_ELT(current_pet, 0));
    splits_bin_or_cont_pointer = INTEGER(VECTOR_ELT(current_pet, 1));
    split_points_pointer = REAL(VECTOR_ELT(current_pet, 2));
    preds_pointer = REAL(VECTOR_ELT(current_pet, 3));

    if(covariable_mode >= 2) {
      SET_VECTOR_ELT(current_pet, 6, allocVector(VECSXP, number_of_nodes));
      model_list_R = VECTOR_ELT(current_pet, 6);
    } else {
      SET_VECTOR_ELT(current_pet, 6, R_NilValue);
    }

    j = 0;

    while(stack->top != NULL) {
      current_node = stack_pop(stack);

      splits_pointer[j] = current_node->split + 1;
      splits_bin_or_cont_pointer[j] = current_node->split_bin_or_cont;
      split_points_pointer[j] = current_node->split_point;
      preds_pointer[j] = current_node->pred;

      if(covariable_mode >= 2) {
        current_func = current_node->func_pred;
        if(current_func != NULL) {
          func_pred_R = allocVector(VECSXP, 3);
          SET_VECTOR_ELT(model_list_R, j, func_pred_R);
          bcde_R = allocVector(REALSXP, 4);
          SET_VECTOR_ELT(func_pred_R, 0, bcde_R);
          SET_VECTOR_ELT(func_pred_R, 1, ScalarLogical(y_bin));
          SET_VECTOR_ELT(func_pred_R, 2, ScalarInteger(current_func->func_type));
          // Assign class
          if(current_func->func_type == 0)
            classgets(func_pred_R, mkString("4pl"));
          else
            classgets(func_pred_R, mkString("linear"));
          bcde = REAL(bcde_R);
          bcde[0] = current_func->b;
          bcde[1] = current_func->c;
          bcde[2] = current_func->d;
          bcde[3] = current_func->e;
        } else {
          SET_VECTOR_ELT(model_list_R, j, R_NilValue);
        }
      }

      if(!(current_node->leaf)) {
        stack_push(stack, current_node->right);
        stack_push(stack, current_node->left);
      }
      j++;
    }
    buffer = current_list_item;
    current_list_item = current_list_item->next;
    Free(buffer);
  }
  stack_destroy(stack);
  UNPROTECT(1);
  return ret;
}

linked_list2* prune(node* tree) {
  node* current_node;
  logic_stack_t *stack = stack_new();
  stack_push(stack, tree);

  // Pre-prune tree by cutting off leaves with an impurity decrease of 0
  linked_list2* prune_list = pre_prune(tree);
  linked_list2* current_list_item = prune_list;

  node* removed_inner_node = NULL;
  double min_g;
  double g, leaf_imp_sum, leaf_imp;
  double p_k;
  int n_leaves;
  int N = tree->N_k;
  node* new_tree = tree;

  while(!(new_tree->leaf)) {
    new_tree = copyTree(new_tree);
    stack_push(stack, new_tree);
    min_g = R_PosInf;
    while(stack->top != NULL) {
      current_node = stack_pop(stack);
      if(!(current_node->leaf)) {
        leaf_imp_sum = calcWeightedLeafImpurities(current_node, N);
        p_k = (double) (current_node->N_k)/N;
        leaf_imp = p_k * (current_node->impurity);
        n_leaves = getNumberOfLeaves(current_node);
        g = (leaf_imp - leaf_imp_sum)/(n_leaves - 1);
        if(g < min_g) {
          min_g = g;
          removed_inner_node = current_node;
        }

        stack_push(stack, current_node->right);
        stack_push(stack, current_node->left);
      }
    }
    makeInnerNode(removed_inner_node);
    current_list_item->next = (linked_list2*) Calloc(1, linked_list2);
    current_list_item = current_list_item->next;
    current_list_item->alpha = min_g;
    current_list_item->tree = new_tree;
  }

  stack_destroy(stack);
  return prune_list;
}

int getNumberOfNodes(node* tree) {
  node* current_node;
  logic_stack_t *stack = stack_new();
  stack_push(stack, tree);
  int number_of_nodes = 0;

  while(stack->top != NULL) {
    current_node = stack_pop(stack);
    number_of_nodes++;
    if(!(current_node->leaf)) {
      stack_push(stack, current_node->right);
      stack_push(stack, current_node->left);
    }
  }
  stack_destroy(stack);
  return number_of_nodes;
}

int getNumberOfLeaves(node* tree) {
  node* current_node;
  logic_stack_t *stack = stack_new();
  stack_push(stack, tree);
  int n_leaves = 0;

  while(stack->top != NULL) {
    current_node = stack_pop(stack);
    if(current_node->leaf) {
      n_leaves++;
    } else {
      stack_push(stack, current_node->right);
      stack_push(stack, current_node->left);
    }
  }
  stack_destroy(stack);
  return n_leaves;
}

double calcWeightedLeafImpurities(node* tree, int N) {
  node* current_node;
  logic_stack_t *stack = stack_new();
  stack_push(stack, tree);
  double w_imp = 0;
  double p_k;

  while(stack->top != NULL) {
    current_node = stack_pop(stack);
    if(current_node->leaf) {
      p_k = (double) (current_node->N_k)/N;
      w_imp += p_k * (current_node->impurity);
    } else {
      stack_push(stack, current_node->right);
      stack_push(stack, current_node->left);
    }
  }
  stack_destroy(stack);
  return w_imp;
}

linked_list2* pre_prune(node* tree) {
  node* current_node;
  logic_stack_t *stack = stack_new();
  double imp_decrease;
  int pruned_nodes = 1;

  linked_list2* prune_list = (linked_list2*) Calloc(1, linked_list2);
  linked_list2* current_list_item = prune_list;
  current_list_item->alpha = 0;
  current_list_item->tree = copyTree(tree);
  node* new_tree = tree;

  while(pruned_nodes != 0) {
    new_tree = copyTree(new_tree);
    stack_push(stack, new_tree);
    pruned_nodes = 0;
    while(stack->top != NULL) {
      current_node = stack_pop(stack);
      if(current_node->leaf) continue;
      if(current_node->left->leaf && current_node->right->leaf) {
        imp_decrease = impurity_decrease(current_node->impurity, current_node->left->impurity, current_node->right->impurity, (double) (current_node->left->N_k)/(current_node->N_k));
        if(fabs(imp_decrease) <= 1e-7) {
          current_list_item->next = (linked_list2*) Calloc(1, linked_list2);
          current_list_item = current_list_item->next;
          current_list_item->alpha = 0;
          current_list_item->tree = new_tree;
          makeInnerNode(current_node);
          pruned_nodes++;
          break;
        }
      } else {
        stack_push(stack, current_node->right);
        stack_push(stack, current_node->left);
      }
    }
  }
  tree_destroy(new_tree);
  stack_destroy(stack);
  return prune_list;
}

void makeInnerNode(node* current_node) {
  current_node->leaf = 1;
  current_node->split = -1;
  tree_destroy(current_node->left); tree_destroy(current_node->right);
  current_node->left = NULL; current_node->right = NULL;
}

node* copyTree(node* tree) {
  node* new_tree = Calloc(1, node);
  node* new_node = new_tree;
  node* current_node;

  logic_stack_t *stack_old = stack_new();
  logic_stack_t *stack_new2 = stack_new();
  stack_push(stack_old, tree);
  stack_push(stack_new2, new_node);
  while(stack_old->top != NULL) {
    current_node = stack_pop(stack_old);
    new_node = stack_pop(stack_new2);

    new_node->leaf = current_node->leaf;
    new_node->split = current_node->split;
    new_node->split_bin_or_cont = current_node->split_bin_or_cont;
    new_node->split_point = current_node->split_point;
    new_node->N_k = current_node->N_k;
    new_node->pred = current_node->pred;
    new_node->ll = current_node->ll;
    new_node->impurity = current_node->impurity;
    new_node->obs_ind = Calloc(new_node->N_k, int);
    memcpy(new_node->obs_ind, current_node->obs_ind, current_node->N_k * sizeof(int));
    if(current_node->func_pred != NULL) {
      new_node->func_pred = Calloc(1, functional);
      new_node->func_pred->b = current_node->func_pred->b;
      new_node->func_pred->c = current_node->func_pred->c;
      new_node->func_pred->d = current_node->func_pred->d;
      new_node->func_pred->e = current_node->func_pred->e;
      new_node->func_pred->y_bin = current_node->func_pred->y_bin;
      new_node->func_pred->func_type = current_node->func_pred->func_type;
    }

    if(!(current_node->leaf)) {
      new_node->left = Calloc(1, node); new_node->right = Calloc(1, node);
      stack_push(stack_old, current_node->right);
      stack_push(stack_old, current_node->left);
      stack_push(stack_new2, new_node->right);
      stack_push(stack_new2, new_node->left);
    } else {
      new_node->left = NULL; new_node->right = NULL;
    }
  }
  stack_destroy(stack_old);
  stack_destroy(stack_new2);
  return new_tree;
}

SEXP simplifyTree_(SEXP pet, SEXP transform_R) {
  node* tree = (node*) R_ExternalPtrAddr(VECTOR_ELT(pet, 5));
  node* current_node;
  int* splits_pointer = INTEGER(VECTOR_ELT(pet, 0));
  int* transform = INTEGER(transform_R);

  logic_stack_t *stack = stack_new();
  stack_push(stack, tree);

  int j = 0;
  int old_split, new_split;

  while(stack->top != NULL) {
    current_node = stack_pop(stack);

    if(!(current_node->leaf)) {
      old_split = current_node->split;
      new_split = transform[old_split];
      current_node->split = new_split;
      splits_pointer[j] = new_split + 1;

      stack_push(stack, current_node->right);
      stack_push(stack, current_node->left);
    }
    j++;
  }

  stack_destroy(stack);
  return pet;
}


