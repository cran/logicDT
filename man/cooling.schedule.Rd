% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cooling.schedule.R
\name{cooling.schedule}
\alias{cooling.schedule}
\title{Define the cooling schedule for simulated annealing}
\usage{
cooling.schedule(
  type = "adaptive",
  start_temp = 1,
  end_temp = -1,
  lambda = 0.01,
  total_iter = 2e+05,
  markov_iter = 1000,
  markov_leave_frac = 1,
  acc_type = "probabilistic",
  frozen_def = "acc",
  frozen_acc_frac = 0.01,
  frozen_markov_count = 5,
  frozen_markov_mode = "total",
  start_temp_steps = 10000,
  start_acc_ratio = 0.95,
  auto_start_temp = TRUE,
  remember_models = TRUE,
  print_iter = 1000
)
}
\arguments{
\item{type}{Type of cooling schedule. \code{"adaptive"}
(default) or \code{"geometric"}}

\item{start_temp}{Start temperature on a log10 scale.
Only used if \code{auto_start_temp = FALSE}.}

\item{end_temp}{End temperature on a log10 scale.
Only used if \code{type = "geometric"}.}

\item{lambda}{Cooling parameter for the adaptive
schedule. Values between 0.01 and 0.1 are recommended
such that in total, several hundred thousand
iterations are performed. Lower values lead to a more
fine search with more iterations while higher values
lead to a more coarse search with less total
iterations.}

\item{total_iter}{Total number of iterations that
should be performed. Only used for the geometric
cooling schedule.}

\item{markov_iter}{Number of iterations for each
Markov chain. The standard value does not need to be
tuned, since the temperature steps and number of
iterations per chain act complementary to each other,
i.e., less iterations can be compensated by smaller
temperature steps.}

\item{markov_leave_frac}{Fraction of accepted moves
leading to an early temperature reduction. This is
primarily used at (too) high temperatures lowering
the temperature if essentially a random walk is
performed. E.g., a value of 0.5 together with
\code{markov_iter = 1000} means that the chain will
be left if \eqn{0.5 \cdot 1000 = 500} states were
accepted in a single chain.}

\item{acc_type}{Type of acceptance function. The
standard \code{"probabilistic"} uses the conventional
function
\eqn{\exp((\mathrm{Score}_\mathrm{old} -
\mathrm{Score}_\mathrm{new})/t)}
for calculating the acceptance probability.
\code{"deterministic"} accepts the new state, if and only
if
\eqn{\mathrm{Score}_\mathrm{new} -
\mathrm{Score}_\mathrm{old} < t}.}

\item{frozen_def}{How to define a frozen chain.
\code{"acc"} means that if less than
\eqn{\texttt{frozen\_acc\_frac} \cdot
\texttt{markov\_iter}} states with different scores
were accepted in a single chain, this chain is
marked as frozen. \code{"sd"} declares a chain as frozen if the
corresponding score standard deviation is zero.
Several frozen chains indicate that the search is finished.}

\item{frozen_acc_frac}{If \code{frozen_def = "acc"}, this parameter
determines the fraction of iterations that define a frozen chain.}

\item{frozen_markov_count}{Number of frozen chains that
need to be observed for finishing the search.}

\item{frozen_markov_mode}{Do the frozen chains
have to occur consecutively (\code{"consecutive"})
or is the total number of frozen chains
relevant (\code{"total"})?}

\item{start_temp_steps}{Number of iterations that should be used for
estimating the ideal start temperature if \code{auto_start_temp =
TRUE} is set.}

\item{start_acc_ratio}{Acceptance ratio that
should be achieved with the automatically
configured start temperature.}

\item{auto_start_temp}{Should the start
temperature be configured automatically?
\code{TRUE} or \code{FALSE}}

\item{remember_models}{Should already evaluated
models be saved in a 2-dimensional hash table
to prevent fitting the same trees multiple
times?}

\item{print_iter}{Number of iterations after which
a progress report shall be printed.}
}
\value{
An object of class \code{cooling.schedule}
  which is a list of all necessary cooling parameters.
}
\description{
This function should be used to configure a search
with simulated annealing.
}
\details{
\code{type = "adapative"} (default)
automatically choses the temperature steps by using the
standard deviation of the scores in a Markov chain
together with the current temperature to
evaluate if equilibrium is achieved. If the standard
deviation is small or the temperature is high,
equilibrium can be assumed leading
to a strong temperature reduction. Otherwise, the
temperature is only merely lowered.
The parameter \code{lambda} is essential to control
how fast the schedule will be executed and, thus,
how many total iterations will be performed.

\code{type = "geometric"} is the conventional
approach which requires more finetuning. Here,
temperatures are uniformly lowered on a log10 scale.
Thus, a start and an end temperature have to be
supplied.
}
