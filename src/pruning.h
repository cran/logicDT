#include "logicDT.h"

typedef struct _linked_list2
{
  double alpha;
  node* tree;
  struct _linked_list2 *next;
} linked_list2;



SEXP prune_(SEXP pet, SEXP y_raw, SEXP Z_raw);
linked_list2* prune(node* tree);
int getNumberOfNodes(node* tree);
int getNumberOfLeaves(node* tree);
double calcWeightedLeafImpurities(node* tree, int N);
linked_list2* pre_prune(node* tree);
void makeInnerNode(node* current_node);
node* copyTree(node* tree);
SEXP simplifyTree_(SEXP pet, SEXP transform_R);




