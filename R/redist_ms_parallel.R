#####################################################
# Author: Cory McCartan
# Institution: Harvard University
# Date Created: 2021/01/31
# Purpose: parallel merge-split
####################################################

#' Parallel Merge-Split/Recombination MCMC Redistricting Sampler
#'
#' `redist_mergesplit_parallel()` runs [redist_mergesplit()] on several
#' chains in parallel.
#'
#' @inherit redist_mergesplit details
#'
#' @inheritParams redist_mergesplit
#' @param chains the number of parallel chains to run. Each chain will have
#' `nsims` draws. If `init_plan` is sampled, each chain will be initialized
#' with its own sampled plan.
#' @param init_plan The initial state of the map, provided as a single vector
#' to be shared across all chains, or a matrix with `chains` columns.
#' If not provided, will default to the reference map of the map object, or if
#' none exists, will sample a random initial state using redist_smc. You can
#' also request a random initial state for each chain by setting
#' init_plan="sample".
#' @param ncores the number of parallel processes to run. Defaults to the
#' maximum available.
#' @param cl_type the cluster type (see [makeCluster()]). Safest is `"PSOCK"`,
#' but `"FORK"` may be appropriate in some settings.
#' @param return_all if `TRUE` return all sampled plans; otherwise, just return
#' the final plan from each chain.
#'
#' @returns A [`redist_plans`] object with all of the simulated plans, and an
#' additional `chain` column indicating the chain the plan was drawn from.
#'
#' @inherit redist_mergesplit references
#'
#' @examples \dontrun{
#' data(fl25)
#' fl_map <- redist_map(fl25, ndists = 3, pop_tol = 0.1)
#' sampled <- redist_mergesplit_parallel(fl_map, nsims = 100, chains = 100)
#' }
#'
#' @concept simulate
#' @md
#' @export
redist_mergesplit_parallel <- function(map, nsims, chains = 1,
                                       warmup = if (is.null(init_plan)) 10 else max(100, nsims %/% 5),
                                       thin = 1L, init_plan = NULL, counties = NULL, compactness = 1,
                                       constraints = list(), constraint_fn = function(m) rep(0, ncol(m)),
                                       adapt_k_thresh = 0.99, k = NULL, ncores = NULL,
                                       cl_type = "PSOCK", return_all = TRUE, init_name = NULL,
                                       silly_adj_fix = FALSE,
                                       verbose = FALSE, silent = FALSE) {
    if (!missing(constraint_fn)) cli_warn("{.arg constraint_fn} is deprecated.")

    map <- validate_redist_map(map)
    V <- nrow(map)
    adj <- get_adj(map)
    ndists <- attr(map, "ndists")
    thin <- as.integer(thin)

    chains <- as.integer(chains)
    stopifnot(chains > 1)

    if (compactness < 0)
        cli_abort("{.arg compactness} must be non-negative.")
    if (adapt_k_thresh < 0 | adapt_k_thresh > 1)
        cli_abort("{.arg adapt_k_thresh} must lie in [0, 1].")
    if (nsims <= warmup)
        cli_abort("{.arg nsims} must be greater than {.arg warmup}.")
    if (thin < 1 || thin > nsims - warmup)
        cli_abort("{.arg thin} must be a positive integer, and no larger than {.arg nsims - warmup}.")
    if (nsims < 1)
        cli_abort("{.arg nsims} must be positive.")

    exist_name <- attr(map, "existing_col")
    counties <- rlang::eval_tidy(rlang::enquo(counties), map)
    if (is.null(init_plan)) {
        if (!is.null(exist_name)) {
            init_plans <- matrix(rep(vctrs::vec_group_id(get_existing(map)), chains), ncol = chains)
            if (is.null(init_name)) {
                init_names <- rep(exist_name, chains)
            } else {
                init_names <- rep(init_name, chains)
            }
        } else {
            init_plan <- "sample"
        }
    } else if (!is.null(init_plan)) {
        if (is.matrix(init_plan)) {
            stopifnot(ncol(init_plan) == chains)
            init_plans <- init_plan
        } else {
            init_plans <- matrix(rep(as.integer(init_plan), chains), ncol = chains)
        }

        if (is.null(init_name))
            init_names <- paste0("<init> ", seq_len(chains))
        else
            init_names <- rep(init_name, chains)
    }
    if (isTRUE(init_plan == "sample")) {
        if (!silent) cat("Sampling initial plans with SMC")
        init_plans <- get_plans_matrix(
            redist_smc(map, chains, counties, compactness, constraints,
                       resample = TRUE, adapt_k_thresh = adapt_k_thresh,
                       ref_name = FALSE, verbose = verbose, silent = silent, ncores = 1))
        if (is.null(init_name))
            init_names <- paste0("<init> ", seq_len(chains))
        else
            init_names <- paste(init_name, seq_len(chains))
    }

    # check init
    if (nrow(init_plans) != V)
        cli_abort("{.arg init_plan} must be as long as the number of units as `map`.")
    if (max(init_plans) != ndists)
        cli_abort("{.arg init_plan} must have the same number of districts as `map`.")
    if (any(apply(init_plans, 2, function(x) contiguity(adj, x)) != 1))
        cli_warn("{.arg init_plan} should have contiguous districts.")

    if (is.null(counties)) {
        counties <- rep(1, V)
    } else {
        if (any(is.na(counties)))
            cli_abort("{.arg counties} must not contain missing values.")

        if (silly_adj_fix) {
            for (j in seq_len(ndists)) {
                idx_distr <- which(init_plans[, 1] == j)
                adj_distr <- redist.reduce.adjacency(adj, idx_distr)
                component <- contiguity(adj_distr, vctrs::vec_group_id(counties[idx_distr]))
                counties[idx_distr] <- paste0(
                    j, ":",
                    dplyr::if_else(component > 1,
                                   paste0(as.character(counties[idx_distr]), "-", component),
                                   as.character(counties[idx_distr]))
                )
            }
            counties <- vctrs::vec_group_id(counties)
        } else {
            counties = vctrs::vec_group_id(counties)
            # handle discontinuous counties
            component <- contiguity(adj, vctrs::vec_group_id(counties))
            counties <- dplyr::if_else(component > 1,
                                       paste0(as.character(counties), "-", component),
                                       as.character(counties)) %>%
                as.factor() %>%
                as.integer()
        }
    }

    # Other constraints
    if (!inherits(constraints, "redist_constr")) {
        constraints <- new_redist_constr(eval_tidy(enquo(constraints), map))
    }
    if (any(c("edges_removed", "log_st") %in% names(constraints))) {
        cli_warn(c("{.var edges_removed} or {.var log_st} constraint found in
           {.arg constraints} and will be ignored.",
            ">" = "Adjust using {.arg compactness} instead."))
    }
    if (any(c("poslby", "fry_hold") %in% names(constraints)) && compactness == 1) {
        cli_warn("{.var polsby} or {.var fry_hold} constraint found in {.arg constraints}
                 with {.arg compactness == 1). This may disrupt efficient sampling.")
    }
    constraints <- as.list(constraints) # drop data attribute

    verbosity <- 1
    if (verbose) verbosity <- 3
    if (silent) verbosity <- 0
    if (is.null(k)) k <- 0

    pop_bounds <- attr(map, "pop_bounds")
    pop <- map[[attr(map, "pop_col")]]
    init_pop <- pop_tally(init_plans, pop, ndists)
    if (any(init_pop < pop_bounds[1]) | any(init_pop > pop_bounds[3]))
        cli_abort("Provided initialization does not meet population bounds.")
    if (any(pop >= pop_bounds[3])) {
        too_big <- as.character(which(pop >= pop_bounds[3]))
        cli_abort(c("Unit{?s} {too_big} ha{?ve/s/ve}
                    population larger than the maximum district size.",
                    "x" = "Redistricting impossible."))
    }

    control = list(adapt_k_thresh=adapt_k_thresh, do_mh=TRUE)
    x <- ms_plans(1, adj, init_plans[, 1], counties, pop, ndists, pop_bounds[2],
                  pop_bounds[1], pop_bounds[3], compactness,
                  list(), control, 0L, 1L, verbosity = 0)
    k <- x$est_k
    rm(x)

    # set up parallel
    if (is.null(ncores)) ncores <- parallel::detectCores()
    ncores <- min(ncores, chains)
    of <- ifelse(Sys.info()[['sysname']] == 'Windows',
                 tempfile(pattern = paste0('ms_', substr(Sys.time(), 1, 10)), fileext = '.txt'),
                 '')
    if (!silent) {
        cl <- parallel::makeCluster(ncores, outfile = of, methods = FALSE,
                                    useXDR = .Platform$endian != "little")
    } else {
        cl <- parallel::makeCluster(ncores, methods = FALSE,
                                    useXDR = .Platform$endian != "little")
    }

    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))

    out_par <- foreach::foreach(chain = seq_len(chains), .inorder = FALSE, .packages="redist") %dorng% {
        if (!silent) cat("Starting chain ", chain, "\n", sep = "")
        run_verbosity <- if (chain == 1 || verbosity == 3) verbosity else 0
        t1_run <- Sys.time()
        algout <- ms_plans(nsims, adj, init_plans[, chain], counties, pop,
                           ndists, pop_bounds[2], pop_bounds[1], pop_bounds[3],
                           compactness, constraints, control, k, thin, run_verbosity)
        t2_run <- Sys.time()

        algout$l_diag <- list(
            runtime = as.numeric(t2_run - t1_run, units = "secs"),
            prethinned_sims = nsims,
            thin = thin,
            warmup = warmup
        )

        algout$mh <- mean(as.logical(algout$mhdecisions))

        warmup_idx <- c(seq_len(1 + warmup %/% thin), nsims %/% thin + 2L)
        if (return_all) {
            algout$plans <- algout$plans[, -warmup_idx, drop = FALSE]
        } else {
            algout$plans <- algout$plans[, nsims + 1L, drop = FALSE]
        }
        storage.mode(algout$plans) <- "integer"

        warmup_idx <- c(seq_len(warmup %/% thin), nsims %/% thin + 1L)
        if (!return_all) {
            algout$mhdecisions <- as.logical(algout$mhdecisions[nsims])
        } else {
            algout$mhdecisions <- as.logical(algout$mhdecisions[-warmup_idx])
        }

        algout
    }


    plans <- lapply(out_par, function(algout) {
        algout$plans
    })
    each_len <- ncol(plans[[1]])
    plans <- do.call(cbind, plans)
    storage.mode(plans) <- "integer"

    mh <- sapply(out_par, function(algout) {
        algout$mh
    })
    l_diag <- lapply(out_par, function(algout) algout$l_diag)

    acceptances <- sapply(out_par, function(algout) {
        algout$mhdecisions
    })


    out <- new_redist_plans(plans = plans, map = map, algorithm = "mergesplit",
                            wgt = NULL, resampled = FALSE,
                            compactness = compactness,
                            constraints = constraints,
                            ndists = ndists,
                            adapt_k_thresh = adapt_k_thresh,
                            mh_acceptance = mh,
                            version = packageVersion("redist"),
                            diagnostics = l_diag) %>%
        mutate(chain = rep(seq_len(chains), each = each_len*ndists),
               mcmc_accept = rep(acceptances, each = ndists))

    if (!is.null(init_names) && !isFALSE(init_name)) {
        if (all(init_names[1] == init_names)) {
            out <- add_reference(out, init_plans[, 1], init_names[1])
        } else {
            out <- Reduce(function(cur, idx) {
                add_reference(cur, init_plans[, idx], init_names[idx]) %>%
                    mutate(chain = dplyr::coalesce(chain, idx))
            }, rev(seq_len(chains)), init = out)
        }
    }

    dplyr::relocate(out, chain, .after = "draw")
}

utils::globalVariables("chain")
