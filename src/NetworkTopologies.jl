# These are taken from: Cotterell, J., & Sharpe, J. (2010). An atlas of gene regulatory networks reveals multiple three‐gene mechanisms for interpreting morphogen gradients. Molecular systems biology, 6(1), 425.

const w_feed_forward = [0 0 0 1 ; 1 0 0 0 ; 1 -1 1 0];
const w_mutual_inh = [0 0 0 1 ; 1 0 -1 0 ; 1 -1 0 0];
const w_frozen_osc = [1 -1 -1 1 ; 0 1 0 0 ; 1 -1 0 0];
const w_overlap_dom = [0 0 -1 1 ; 1 0 0 0 ; -1 1 0 0];
const w_bistable = [0 0 0 1 ; -1 0 1 0 ; 0 -1 1 0];
const w_classical = [0 0 0 1 ; -1 1 0 0 ; -1 -1 1 0];