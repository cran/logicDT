#' Define the cooling schedule for simulated annealing
#'
#' This function should be used to configure a search
#' with simulated annealing.
#'
#' \code{type = "adapative"} (default)
#' automatically choses the temperature steps by using the
#' standard deviation of the scores in a Markov chain
#' together with the current temperature to
#' evaluate if equilibrium is achieved. If the standard
#' deviation is small or the temperature is high,
#' equilibrium can be assumed leading
#' to a strong temperature reduction. Otherwise, the
#' temperature is only merely lowered.
#' The parameter \code{lambda} is essential to control
#' how fast the schedule will be executed and, thus,
#' how many total iterations will be performed.
#'
#' \code{type = "geometric"} is the conventional
#' approach which requires more finetuning. Here,
#' temperatures are uniformly lowered on a log10 scale.
#' Thus, a start and an end temperature have to be
#' supplied.
#'
#' @param type Type of cooling schedule. "adaptive"
#'   (default) or "geometric"
#' @param start_temp Start temperature on a log10 scale.
#'   Only used if \code{auto_start_temp = FALSE}.
#' @param end_temp End temperature on a log10 scale.
#'   Only used if \code{type = "geometric"}.
#' @param lambda Cooling parameter for the adaptive
#'   schedule. Values between 0.01 and 0.1 are recommended
#'   such that in total, several hundred thousand
#'   iterations are performed. Lower values lead to a more
#'   fine search with more iterations while higher values
#'   lead to a more rigorous search with less total
#'   iterations.
#' @param total_iter Total number of iterations that
#'   should be performed. Only used for the geometric
#'   cooling schedule.
#' @param markov_iter Number of iterations for each
#'   Markov chain. The standard value does not need to be
#'   tuned, since the temperature steps and number of
#'   iterations per chain act complementary to each other,
#'   i.e., less iterations can be compensated by smaller
#'   temperature steps.
#' @param markov_leave_frac Fraction of accepted moves
#'   leading to an early temperature reduction. This is
#'   primarily used at (too) high temperatures lowering
#'   the temperature if essentially a random walk is
#'   performed. E.g., a value of 0.5 together with
#'   \code{markov_iter = 1000} means that the chain will
#'   be left if \eqn{0.5 \cdot 1000 = 500} states were
#'   accepted in a single chain.
#' @param acc_type Type of acceptance function. The
#'   standard "probabilistic" uses the conventional
#'   function
#'   \eqn{\exp((\mathrm{Score}_\mathrm{old} -
#'   \mathrm{Score}_\mathrm{new})/t)}
#'   for calculating the acceptance probability.
#'   "deterministc" accepts the new state, if and only
#'   if
#'   \eqn{\mathrm{Score}_\mathrm{new} -
#'   \mathrm{Score}_\mathrm{old} < t}.
#' @param frozen_def How to define a frozen chain.
#'   "acc" means that if less than
#'   \eqn{\texttt{frozen\_acc\_frac} \cdot
#'   \texttt{markov\_iter}} states with different scores
#'   were accepted in a single chain, this chain is
#'   marked as frozen. Several frozen chains indicate
#'   that the search is finished.
#' @param frozen_acc_frac If \code{frozen_def =
#'   frozen_acc_frac}, this parameter determines the
#'   fraction of iterations which defines a frozen
#'   chain.
#' @param frozen_markov_count How many frozen chains
#'   need to be observed for finishing the search?
#' @param frozen_markov_mode Do the frozen chains
#'   have to occur consecutively ("consecutive")
#'   or is the total number of frozen chains
#'   relevant ("total")?
#' @param start_temp_steps If \code{auto_start_temp =
#'   TRUE}, how many iterations should be used for
#'   estimating the ideal start temperature?
#' @param start_acc_ratio Which acceptance ratio
#'   should be achieved with the automatically
#'   configured start temperature?
#' @param auto_start_temp Should the start
#'   temperature be configured automatically?
#'   TRUE or FALSE
#' @param remember_models Should already evaluated
#'   models be saved in a 2-dimensional hash table
#'   to prevent fitting the same trees multiple
#'   times?
#' @param print_iter After how many iterations
#'   shall a progress report be printed?
#' @return An object of class \code{cooling.schedule}
#'   which is a list of all necessary cooling parameters.
#'
#' @export
cooling.schedule <- function(type = "adaptive", start_temp = 1, end_temp = -1, lambda = 0.01,
                             total_iter = 200000, markov_iter = 1000,
                             markov_leave_frac = 1.0,
                             acc_type = "probabilistic",
                             frozen_def = "acc", frozen_acc_frac = 0.01,
                             frozen_markov_count = 5, frozen_markov_mode = "total",
                             start_temp_steps = 10000, start_acc_ratio = 0.95,
                             auto_start_temp = TRUE, remember_models = TRUE,
                             print_iter = 1000) {
  cs <- list(type = type, start_temp = as.numeric(start_temp), end_temp = as.numeric(end_temp), lambda = as.numeric(lambda),
             total_iter = as.integer(total_iter), markov_iter = as.integer(markov_iter),
             markov_leave_frac = as.numeric(markov_leave_frac),
             acc_type = acc_type, acc_type2 = as.integer(ifelse(acc_type == "deterministic", 1, 0)),
             frozen_def = frozen_def, frozen_def2 = as.integer(ifelse(frozen_def == "sd", 1, 0)),
             frozen_acc_frac = as.numeric(frozen_acc_frac),
             frozen_markov_count = as.integer(frozen_markov_count),
             frozen_markov_mode = frozen_markov_mode, frozen_markov_mode2 = as.integer(ifelse(frozen_markov_mode == "consecutive", 1, 0)),
             start_temp_steps = as.integer(start_temp_steps), start_acc_ratio = as.numeric(start_acc_ratio),
             auto_start_temp = auto_start_temp, remember_models = remember_models,
             print_iter = as.integer(print_iter),
             type2 = as.integer(ifelse(type == "geometric", 1, 0)))
  class(cs) <- "cooling.schedule"

  cs$real_start_temp <- as.numeric(10^(start_temp))

  if(type == "adaptive") {
    cs$real_end_temp <- as.numeric(0)
    cs$temp_func <- function(t, sd) ifelse(sd > 0, t * exp(-cs$lambda * t/sd), t)
  }
  else if(type == "geometric") {
    cs$real_end_temp <- as.numeric(10^(end_temp))
    cs$n_steps <- as.integer(floor(total_iter/markov_iter))
    cs$step_temp <- as.numeric((start_temp - end_temp)/cs$n_steps)
    cs$temp_func <- function(t, sd) t/10^(cs$step_temp)
  }
  return(cs)
}

