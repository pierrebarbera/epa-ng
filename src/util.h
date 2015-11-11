#define DEFAULT_BRANCH_LENGTH 0.000001

// attribution!

static void set_missing_branch_length_recursive(pll_utree_t * tree, 
                                                double length)
{
  if (tree)
  {
    /* set branch length to length if not set */
    if (!tree->length)
      tree->length = length;

    if (tree->next)
    {
      if (!tree->next->length)
        tree->next->length = length;

      if (!tree->next->next->length)
        tree->next->next->length = length;

      set_missing_branch_length_recursive(tree->next->back, length);
      set_missing_branch_length_recursive(tree->next->next->back, length);
    }
  }
};

/* set missing branch lengths to length */
static void set_missing_branch_length(pll_utree_t * tree, double length)
{
  set_missing_branch_length_recursive(tree, length);
  set_missing_branch_length_recursive(tree->back, length);
};

static void set_unique_clv_indices_recursive(pll_utree_t * tree, const int num_tip_nodes)
{
  if (tree && tree->next)
  {
    int idx = tree->clv_index;
    /* new index is in principle old index * 3 + 0 for the first traversed, + 1 for the
      second etc., however we need to account for the first num_tip_nodes entries, as 
      the tip nodes only have one clv  */
    idx = (idx - num_tip_nodes) * 3 + num_tip_nodes;
    tree->clv_index = idx;
    //tree->scaler_index = idx;
    tree->next->clv_index = ++idx;
    //tree->next->scaler_index = idx;
    tree->next->next->clv_index = ++idx;
    //tree->next->next->scaler_index = idx;

    // recurse
    set_unique_clv_indices_recursive(tree->next->back, num_tip_nodes);
    set_unique_clv_indices_recursive(tree->next->next->back, num_tip_nodes);
  }
};

static void set_unique_clv_indices(pll_utree_t * tree, const int num_tip_nodes)
{
  set_unique_clv_indices_recursive(tree, num_tip_nodes);
  set_unique_clv_indices_recursive(tree->back, num_tip_nodes);
};

/* a callback function for performing a full traversal */
static int cb_full_traversal(pll_utree_t * node)
{
  return 1;
}