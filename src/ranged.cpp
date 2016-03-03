#include "ranged.hpp"

#include <algorithm>

#include "pll_util.hpp"

static void core_update_partial_tiptip(unsigned int states,
                                      unsigned int sites,
                                      unsigned int rate_cats,
                                      unsigned int log2_maxstates,
                                      unsigned int log2_states,
                                      unsigned int log2_rates,
                                      double * lh_statepair,
                                      double * parent_clv,
                                      const char * left_tipchars,
                                      const char * right_tipchars,
                                      unsigned int attrib)
{
  unsigned int j,k,n;

  unsigned int span = states * rate_cats;

  double * local_lh_statepair;

  (void)attrib;

  for (n = 0; n < sites; ++n)
  {
    j = (unsigned int)(left_tipchars[n]);
    k = (unsigned int)(right_tipchars[n]);

    local_lh_statepair = lh_statepair;
    local_lh_statepair += ((j << log2_maxstates) + k) << (log2_states+log2_rates);

    memcpy(parent_clv, local_lh_statepair, span*sizeof(double));

    parent_clv += span;
  }
}

static void core_update_partial_tipinner(unsigned int states,
                                      unsigned int sites,
                                      unsigned int rate_cats,
                                      unsigned int * revmap,
                                      double * parent_clv,
                                      unsigned int * parent_scaler,
                                      const char * tip_tipchars,
                                      const double * inner_clv,
                                      const double * tip_matrix,
                                      const double * inner_matrix,
                                      unsigned int attrib)
{
  unsigned int i,j,k,n;
  unsigned int scaling;

  const double * lmat;
  const double * rmat;

  unsigned int span = states * rate_cats;

  (void)attrib;

  for (n = 0; n < sites; ++n)
  {
    lmat = tip_matrix;
    rmat = inner_matrix;

    scaling = (parent_scaler) ? 1 : 0;

    for (k = 0; k < rate_cats; ++k)
    {
      for (i = 0; i < states; ++i)
      {
        double terma = 0;
        double termb = 0;
        unsigned int lstate = revmap[(int)tip_tipchars[n]];
        for (j = 0; j < states; ++j)
        {
          if (lstate & 1)
            terma += lmat[j];

          termb += rmat[j] * inner_clv[j];

          lstate >>= 1;
        }
        parent_clv[i] = terma*termb;
        lmat += states;
        rmat += states;

        scaling = scaling && (parent_clv[i] < PLL_SCALE_THRESHOLD);
      }
      parent_clv += states;
      inner_clv  += states;
    }
    /* if *all* entries of the site CLV were below the threshold then scale
       (all) entries by PLL_SCALE_FACTOR */
    if (scaling)
    {
      parent_clv -= span;
      for (i = 0; i < span; ++i)
        parent_clv[i] *= PLL_SCALE_FACTOR;
      parent_clv += span;
      parent_scaler[n] += 1;
    }
  }
}

/*  Update a single partial likelihood vector (CLV) based on a range. Requires the
    CLV to be appropriately initialized (0.0 within the range, 1.0 without)*/
void update_partial_ranged(pll_partition_t * partition, pll_utree_t * tree, Range& range)
{
  auto begin = range.begin;
  auto span = range.span;

  auto parent = tree;
  auto child1 = parent->next->back;
  auto child2 = parent->next->next->back;

  // if this is a tip-inner case, ensure that child 1 is the tip
  if (is_char_tip(child1) != is_char_tip(child2))
    if (is_char_tip(child2))
    {
      auto tmp = child2;
      child2 = child1;
      child1 = tmp;
    }

  auto parent_clv_index = parent->clv_index;
  int parent_scaler_index = parent->clv_index;
  auto child1_clv_index = child1->clv_index;
  auto child1_scaler_index = child1->scaler_index;
  auto child1_matrix_index = child1->pmatrix_index;
  auto child2_clv_index = child2->clv_index;
  auto child2_scaler_index = child2->scaler_index;
  auto child2_matrix_index = child2->pmatrix_index;

#ifdef __AVX
  auto attrib = PLL_ATTRIB_ARCH_AVX;
#else
  auto attrib = PLL_ATTRIB_ARCH_SSE;
#endif

  // determine how many entries we need to skip for clvs and tipchars
  const auto states = partition->states;
  auto clv_size = states * partition->rate_cats;
  auto skip_clv = begin * clv_size;
  auto tipchars_size = states;
  auto skip_tipchars = begin * tipchars_size;

  // quick probe to check if properly initialized
  if (begin > 0)
    assert(partition->clv[parent_clv_index][0] != 0.0);
  else if (span < partition->sites)
    assert(partition->clv[parent_clv_index][span] != 0.0);

  // scalers are per site, so we just skip <begin>
  // check if scalers exist, and if so shift them accordingly
  unsigned int * parent_scaler = (parent_scaler_index == PLL_SCALE_BUFFER_NONE) ?
      nullptr : partition->scale_buffer[parent_scaler_index] + begin;
  unsigned int * child1_scaler = (child1_scaler_index == PLL_SCALE_BUFFER_NONE) ?
      nullptr : partition->scale_buffer[child1_scaler_index] + begin;
  unsigned int * child2_scaler = (child2_scaler_index == PLL_SCALE_BUFFER_NONE) ?
      nullptr : partition->scale_buffer[child2_scaler_index] + begin;

  if (is_char_tip(child1) && is_char_tip(child2))
  {
    core_update_partial_tiptip(
      states,
      span, // begin + span = first CLV entry not used in computation
      partition->rate_cats,
      partition->log2_maxstates,
      partition->log2_states,
      partition->log2_rates,
      partition->lh_statepair,
      partition->clv[parent_clv_index] + skip_clv,
      // parent_scaler,
      partition->tipchars[child1_clv_index] + skip_tipchars,
      partition->tipchars[child2_clv_index] + skip_tipchars,
      // partition->pmatrix[child1_matrix_index],
      // partition->pmatrix[child2_matrix_index],
      // child1_scaler,
      // child2_scaler,
      attrib);
  }
  else if (is_char_tip(child1) or is_char_tip(child2))
  {
    core_update_partial_tipinner(
      states,
      span, // begin + span = first CLV entry not used in computation
      partition->rate_cats,
      partition->revmap,
      partition->clv[parent_clv_index] + skip_clv,
      parent_scaler,
      partition->tipchars[child1_clv_index] + skip_tipchars,
      partition->clv[child2_clv_index] + skip_clv,
      partition->pmatrix[child1_matrix_index],
      partition->pmatrix[child2_matrix_index],
      // child1_scaler,
      // child2_scaler,
      attrib);
  }
  else
    pll_core_update_partial(partition->states,
      span, // begin + span = first CLV entry not used in computation
      partition->rate_cats,
      partition->clv[parent_clv_index] + skip_clv,
      parent_scaler,
      partition->clv[child1_clv_index] + skip_clv,
      partition->clv[child2_clv_index] + skip_clv,
      partition->pmatrix[child1_matrix_index],
      partition->pmatrix[child2_matrix_index],
      child1_scaler,
      child2_scaler,
      attrib);
}

static double core_compute_edge_logl(const unsigned int states,
                                    const unsigned int states_padded,
                                    const unsigned int sites,
                                    const unsigned int rate_cats,
                                    const double prop_invar,
                                    const unsigned int mixture,
                                    const double * parent_clv,
                                    const unsigned int * parent_scaler,
                                    const double * child_clv,
                                    const unsigned int * child_scaler,
                                    const double * pmatrix,
                                    const double * frequencies,
                                    const double * rate_weights,
                                    const unsigned int * pattern_weights,
                                    const int * invariant)
{
  unsigned int n,i,j,k;
  double logl = 0.0;
  double terma, terma_r, termb;
  double site_lk, inv_site_lk;

  const double * local_pmatrix;
  const double * local_frequencies;

  unsigned int scale_factors;

  for (n = 0; n < sites; ++n)
  {
    local_pmatrix = pmatrix;
    local_frequencies = frequencies;
    terma = 0;
    for (i = 0; i < rate_cats; ++i)
    {
      terma_r = 0;
      for (j = 0; j < states; ++j)
      {
        termb = 0;
        for (k = 0; k < states; ++k)
        {
          termb += pmatrix[k] * child_clv[k];
        }
        terma_r += parent_clv[j] * local_frequencies[j] * termb;
        pmatrix += states;
      }
      terma += terma_r * rate_weights[i];
      parent_clv += states;
      child_clv += states;
      if (mixture > 1)
        local_frequencies += states_padded;
    }

    site_lk = terma;

    /* account for invariant sites */
    if (prop_invar > 0)
    {
      inv_site_lk = (invariant[n] == -1) ?
                        0 : local_frequencies[invariant[n]];

      site_lk = site_lk * (1. - prop_invar) +
                inv_site_lk * prop_invar;
    }

    scale_factors = (parent_scaler) ? parent_scaler[n] : 0;
    scale_factors += (child_scaler) ? child_scaler[n] : 0;

    logl += log(site_lk) * pattern_weights[n];
    if (scale_factors)
      logl += scale_factors * log(PLL_SCALE_THRESHOLD);
  }

  return logl;
}

static double core_compute_edge_logl_tipinner(const unsigned int states,
                                    const unsigned int states_padded,
                                    const unsigned int sites,
                                    const unsigned int rate_cats,
                                    const double prop_invar,
                                    const unsigned int mixture,
                                    const double * parent_clv,
                                    const unsigned int * parent_scaler,
                                    const char * child_tipchars,
                                    const double * pmatrix,
                                    const double * frequencies,
                                    const double * rate_weights,
                                    const unsigned int * pattern_weights,
                                    const unsigned int * revmap,
                                    const int * invariant)
{
  unsigned int n,i,j,k;
  double logl = 0;
  double terma, terma_r, termb;
  double site_lk, inv_site_lk;

  unsigned int scale_factors;

  const double * local_pmatrix;
  const double * local_frequencies;

  /* child is the tip sequence, gets its tipchar */
  const char * tipchar = child_tipchars;
  unsigned int cstate;

  for (n = 0; n < sites; ++n)
  {
    local_pmatrix = pmatrix;
    local_frequencies = frequencies;
    terma = 0;
    for (i = 0; i < rate_cats; ++i)
    {
      terma_r = 0;
      for (j = 0; j < states; ++j)
      {
        termb = 0;
        cstate = revmap[(int)(*tipchar)];
        for (k = 0; k < states; ++k)
        {
          if (cstate & 1)
            termb += local_pmatrix[k];
          cstate >>= 1;
        }
        terma_r += parent_clv[j] * local_frequencies[j] * termb;
        local_pmatrix += states;
      }
      terma += terma_r * rate_weights[i];
      parent_clv += states;
      if (mixture > 1)
        local_frequencies += states_padded;
    }

    site_lk = terma;

    /* account for invariant sites */
    if (prop_invar > 0)
    {
      inv_site_lk = (invariant[n] == -1) ?
                        0 : local_frequencies[invariant[n]];

      site_lk = site_lk * (1. - prop_invar) +
                inv_site_lk * prop_invar;
    }

    scale_factors = (parent_scaler) ? parent_scaler[n] : 0;

    logl += log(site_lk) * pattern_weights[n];
    if (scale_factors)
      logl += scale_factors * log(PLL_SCALE_THRESHOLD);

    tipchar++;
  }

  return logl;
}

double compute_edge_logl_ranged(pll_partition_t * partition, pll_utree_t * tree, Range& range)
{
  const unsigned int freqs_index = 0; // TODO change if we ever have more than one
  auto begin = range.begin;
  auto span = range.span;

  auto parent = tree;
  auto child = parent->back;

  bool tip_inner = false;

  // if this is a tip-inner case, ensure that child 1 is the tip
  if (is_char_tip(child) != is_char_tip(parent))
  {
    if (is_char_tip(parent))
    {
      auto tmp = parent;
      parent = child;
      child = tmp;
    }
    tip_inner = true;
  }

  auto parent_clv_index = parent->clv_index;
  int parent_scaler_index = parent->clv_index;
  auto child_clv_index = child->clv_index;
  auto child_scaler_index = child->scaler_index;
  auto child_matrix_index = child->pmatrix_index;

  // determine how many entries we need to skip for clvs and tipchars
  const auto states = partition->states;
  auto clv_size = states * partition->rate_cats;
  auto skip_clv = begin * clv_size;
  auto tipchars_size = states;
  auto skip_tipchars = begin * tipchars_size;

  const double * pmatrix = partition->pmatrix[child_matrix_index];

  // scalers are per site, so we just skip <begin>
  // check if scalers exist, and if so shift them accordingly
  unsigned int * parent_scaler = (parent_scaler_index == PLL_SCALE_BUFFER_NONE) ?
      nullptr : partition->scale_buffer[parent_scaler_index] + begin;
  unsigned int * child_scaler = (child_scaler_index == PLL_SCALE_BUFFER_NONE) ?
      nullptr : partition->scale_buffer[child_scaler_index] + begin;

  if (tip_inner) // TODO always do inner inner if pattern tip isn't set
    return core_compute_edge_logl_tipinner(
      states,
      partition->states_padded,
      span,
      partition->rate_cats,
      partition->prop_invar[freqs_index],
      partition->mixture,
      partition->clv[parent_clv_index] + skip_clv,
      parent_scaler,
      partition->tipchars[child_clv_index] + skip_tipchars,
      pmatrix,
      partition->frequencies[freqs_index],
      partition->rate_weights,
      partition->pattern_weights,
      partition->revmap,
      partition->invariant

    );
  else
    return core_compute_edge_logl(
      states,
      partition->states_padded,
      span,
      partition->rate_cats,
      partition->prop_invar[freqs_index],
      partition->mixture,
      partition->clv[parent_clv_index] + skip_clv,
      parent_scaler,
      partition->clv[child_clv_index] + skip_clv,
      child_scaler,
      pmatrix,
      partition->frequencies[freqs_index],
      partition->rate_weights,
      partition->pattern_weights,
      partition->invariant);

}

void fill_without(pll_partition_t * partition, unsigned int clv_index, Range& range, double value)
{
  auto size = partition->states * partition->rate_cats;
  auto skip_clv = range.begin * size;
  auto clv = partition->clv[clv_index];
  std::fill_n(clv, skip_clv, value);
  auto first_after_skip = (range.begin + range.span) * size;
  std::fill_n(clv + first_after_skip,
    (partition->sites * size) - first_after_skip, value);
}
