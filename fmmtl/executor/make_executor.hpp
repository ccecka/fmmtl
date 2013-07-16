#pragma once

#include "fmmtl/executor/EvalSimple.hpp"

template <typename Context, typename Options>
EvaluatorBase<Context>* make_evaluator(Context& c, Options& opts) {
  // Determine thetype of Evaluator
  // Statically from the type of Options
  // Dynamically from the Options input

  // For now
  return make_eval_simple(c, opts);
}
