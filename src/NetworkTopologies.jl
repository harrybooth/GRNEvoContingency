const w_feed_forward = [0 0 0 1 ; 1 0 0 0 ; 1 -1 1 0]
const w_mutual_inh = [0 0 0 1 ; 1 0 -1 0 ; 1 -1 0 0]
const w_frozen_osc = [1 -1 -1 1 ; 0 1 0 0 ; 1 -1 0 0]
const w_overlap_dom = [0 0 -1 1 ; 1 0 0 0 ; -1 1 0 0]
const w_bistable = [0 0 0 1 ; -1 0 1 0 ; 0 -1 1 0]
const w_classical = [0 0 0 1 ; -1 1 0 0 ; -1 -1 1 0];

const w_ex1 =  [0.0 0.0 -0.0588219 2.08643;
            0.388038 0.0 0.0 0.0;
            -1.47171 0.0672493 0.0 0.0];

;

function random_w()
    0.9995 .^ rand(0:10000,Ng,Ng+1) .* 10 .* rand(Ng,Ng+1)
end