#pragma once

#include "hyperparams.hpp"
#include "params.hpp"
#include "clusters.hpp"


double calc_elbo(const HyperParams& HP, const Params& params, const Clusters& clusters);

