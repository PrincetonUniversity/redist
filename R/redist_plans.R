##############################################
## Author: Cory McCartan
## Institution: Harvard University
## Date Created: 2021/01/28
## Purpose: redist functions for a tidy workflow
##############################################


# constructors and reconstructors -----------------------------------------


# plans has n_precinct columns and n_sims rows
# map is a redist_map
# algorithm is one of "smc" or "mcmc"
# wgt is the weights before any resampling or truncation
# ... will depend on the algorithm
new_redist_plans <- function(plans, map, algorithm, wgt, resampled = TRUE, ndists = attr(map, "ndists"), ...) {
    n_sims <- ncol(plans)
    if (n_sims < 1) cli_abort("Need at least one simulation draw.")

    n_prec <- nrow(plans)
    map_dists <- attr(map, "ndists")
    partial <- ndists < map_dists

    prec_pop <- map[[attr(map, "pop_col")]]
    if (!partial) {
        distr_range <- 1:ndists
        distr_pop <- pop_tally(plans, prec_pop, ndists)
    } else {
        distr_range <- 1:ndists - 1L
        distr_pop <- pop_tally(plans, prec_pop, ndists)
        pl_tmp <- plans + 1L
        distr_pop <- pop_tally(pl_tmp, prec_pop, ndists)
    }

    attr_names <- c("redist_attr", "plans", "ndists", "algorithm", "wgt",
        "resampled", "ndists", "merge_idx", "prec_pop",
        names(list(...)))

    if (is.null(colnames(plans))) {
        draw_fac = as.factor(1:n_sims)
    } else {
        draw_fac = character(n_sims)
        ref_idx = which(nchar(colnames(plans)) > 0)
        draw_fac[ref_idx] = colnames(plans)[ref_idx]
        draw_fac[-ref_idx] = as.character(seq_len(n_sims - length(ref_idx)))
        draw_fac = factor(draw_fac, levels=draw_fac)
    }

    structure(tibble(draw = rep(draw_fac, each = ndists),
        district = rep(distr_range, n_sims),
        total_pop = as.numeric(distr_pop)),
    plans = plans, ndists = ndists, algorithm = algorithm, wgt = wgt,
    resampled = resampled, merge_idx = attr(map, "merge_idx"),
    prec_pop = prec_pop, redist_attr = attr_names, ...,
    class = c("redist_plans", "tbl_df", "tbl", "data.frame"))
}

validate_redist_plans <- function(x) {
    if (names(x)[1] != "draw") cli_abort("First column must be named \"{.field draw}\"")
    if (!is.factor(x$draw)) cli_abort("{.field draw} column must be a factor")

    plan_m <- attr(x, "plans")
    if (is.null(plan_m)) cli_abort("Missing plans matrix")

    min_distr <- colmin(plan_m)
    max_distr <- colmax(plan_m)
    if (any(min_distr != 1) || any(diff(max_distr) != 0))
        cli_abort("District numbers must start at 1 and run sequentially to the number of districts.")

    x
}

reconstruct.redist_plans <- function(data, old) {
    if (!("draw" %in% names(data)))
        return(data)

    if (!missing(old)) {
        for (name in attr(old, "redist_attr")) {
            if (is.null(attr(data, name)))
                attr(data, name) <- attr(old, name)
        }
    }

    classes <- c("tbl_df", "tbl", "data.frame")
    if (inherits(data, "grouped_df"))
        classes <- c("grouped_df", classes)

    if (inherits(data, "rowwise_df"))
        classes <- c("rowwise_df", classes)

    class(data) <- c("redist_plans", classes)

    data
}

#' A set of redistricting plans
#'
#' A \code{redist_plans} object is essentially a data frame of summary
#' information on each district and each plan, along with the matrix of district
#' assignments and information about the simulation process used to generate the
#' plans.
#'
#' The first two columns of the data frame will be \code{draw}, a factor indexing
#' the simulation draw, and \code{district}, an integer indexing the districts
#' within a plan. The data frame will therefore have \code{n_sims*ndists} rows.
#' As a data frame, the usual \code{\link{dplyr}} methods will work.
#'
#' Other useful methods for \code{redist_plans} objects:
#' * \code{\link{summary.redist_plans}}
#' * \code{\link{add_reference}}
#' * \code{\link{subset_sampled}}
#' * \code{\link{subset_ref}}
#' * \code{\link{pullback}}
#' * \code{\link{number_by}}
#' * \code{\link{match_numbers}}
#' * \code{\link{is_county_split}}
#' * \code{\link{prec_assignment}}
#' * \code{\link{plan_distances}}
#' * \code{\link{get_plans_matrix}}
#' * \code{\link{get_plans_weights}}
#' * \code{\link{get_sampling_info}}
#' * \code{\link{as.matrix.redist_plans}}
#' * \code{\link{plot.redist_plans}}
#'
#' @param plans a matrix with \code{n_precinct} columns and \code{n_sims} rows,
#' or a single vector of precinct assignments.
#' @param map a \code{\link{redist_map}} object
#' @param algorithm the algorithm used to generate the plans (usually "smc" or "mcmc")
#' @param wgt the weights to use, if any.
#' @param ... Other named attributes to set
#'
#' @returns a new \code{redist_plans} object.
#'
#' @examples
#' data(iowa)
#'
#' iowa <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05, total_pop = pop)
#' rsg_plan <- redist.rsg(iowa$adj, iowa$pop, ndists = 4, pop_tol = 0.05)$plan
#' redist_plans(rsg_plan, iowa, "rsg")
#'
#' @md
#' @concept analyze
#' @export
redist_plans <- function(plans, map, algorithm, wgt = NULL, ...) {
    if (is.numeric(plans) && length(plans) == nrow(map)) {
        plans <- matrix(as.integer(plans), ncol = 1)
    }
    if (!is.matrix(plans)) cli_abort("{.arg plans} must be a matrix.")
    if (nrow(plans) != nrow(map)) cli_abort("{.arg plans} matrix must have as many rows as {.arg map} has precincts.")
    if (!inherits(map, "redist_map")) cli_abort("{.arg map} must be a {.cls redist_map}")

    if (min(plans) == 0L) plans <- plans + 1L
    storage.mode(plans) <- "integer"

    obj <- new_redist_plans(plans, map, algorithm, wgt = wgt,
        resampled = FALSE, ...)
    validate_redist_plans(obj)
}


# getters / setters ------------------------------------------------------------


#' Extract the matrix of district assignments from a redistricting simulation
#'
#' @param x the \code{redist_plans} object
#' @param ... ignored
#' @return matrix
#' @concept analyze
#' @export
get_plans_matrix <- function(x) {
    if (!inherits(x, "redist_plans")) cli_abort("Not a {.cls redist_plans}")
    attr(x, "plans")
}
#' @rdname get_plans_matrix
#' @method as.matrix redist_plans
#' @return matrix
#' @export
as.matrix.redist_plans <- function(x, ...) get_plans_matrix(x)

# internal -- no check performed!
set_plan_matrix <- function(x, mat) {
    attr(x, "plans") <- mat
    x
}

#' Extract the sampling weights from a redistricting simulation.
#'
#' May be \code{NULL} if no weights exist (MCMC or optimization methods).
#'
#' @param plans,object the \code{redist_plans} object
#'
#' @returns A numeric vector of weights, with an additional attribute
#' \code{resampled} indicating whether the plans have been resampled according
#' to these weights. If weights have been resampled, this returns the weights
#' before resampling (i.e., they do not correspond to the resampled plans).
#'
#' @concept analyze
#' @export
get_plans_weights <- function(plans) {
    if (!inherits(plans, "redist_plans")) cli_abort("Not a {.cls redist_plans}")
    wgt <- attr(plans, "wgt")
    if (!is.null(wgt))
        attr(wgt, "resampled") <- attr(plans, "resampled")
    wgt
}

#' @rdname get_plans_weights
#' @param ... Ignored.
#' @importFrom stats weights
#' @method weights redist_plans
#' @return numeric vector
#' @export
weights.redist_plans <- function(object, ...) {
    get_plans_weights(object)
}

get_n_ref <- function(x) {
    if (!inherits(x, "redist_plans")) cli_abort("Not a {.cls redist_plans}")
    plans_m <- get_plans_matrix(x)
    if (is.null(colnames(plans_m))) return(0)
    sum(nchar(colnames(plans_m)) > 0)
}

#' Extract the sampling information from a redistricting simulation
#'
#' @param plans the \code{redist_plans} object
#'
#' @returns a list of parameters and information about the sampling problem.
#'
#' @concept analysis
#' @export
get_sampling_info <- function(plans) {
    if (!inherits(plans, "redist_plans")) cli_abort("Not a {.cls redist_plans}")
    all_attr <- attributes(plans)

    all_attr$names <- NULL
    all_attr$row.names <- NULL
    all_attr$class <- NULL
    all_attr$plans <- NULL
    all_attr$redist_attr <- NULL

    all_attr
}


#' Add a reference plan to a set of plans
#'
#' This function facilitates comparing an existing (i.e., non-simulated)
#' redistricting plan to a set of simulated plans.
#'
#' @param plans a \code{redist_plans} object
#' @param ref_plan an integer vector containing the reference plan. It will be
#' renumbered to 1..\code{ndists}.
#' @param name a human-readable name for the reference plan. Defaults to the
#' name of \code{ref_plan}.
#'
#' @returns a modified \code{redist_plans} object containing the reference plan
#' @concept analyze
#' @export
add_reference <- function(plans, ref_plan, name = NULL) {
    if (!inherits(plans, "redist_plans")) cli_abort("{.arg plans} must be a {.cls redist_plans}")
    if (isTRUE(attr(plans, "partial")))
        cli_abort("Reference plans not supported for partial plans objects.")

    plan_m <- get_plans_matrix(plans)
    if (!is.numeric(ref_plan)) cli_abort("{.arg ref_plan} must be numeric")
    if (length(ref_plan) != nrow(plan_m))
        cli_abort("{.arg ref_plan} must have the same number of precincts as {.arg plans}")

    if (is.null(name)) {
        ref_str <- deparse(substitute(ref_plan))
        if (stringr::str_detect(ref_str, stringr::fixed("$")))
            name <- strsplit(ref_str, "$", fixed = TRUE)[[1]][2]
        else
            name <- ref_str
    } else {
        if (!is.character(name)) cli_abort("{.arg name} must be a {.cls chr}")
    }

    ref_plan <- vctrs::vec_group_id(ref_plan)
    ndists <- max(ref_plan)
    if (ndists != max(plan_m[, 1]))
        cli_abort("{.arg ref_plan} has a different number of districts than {.arg plans}")

    # first the matrix
    plan_m <- cbind(ref_plan, plan_m)
    colnames(plan_m)[1] <- name

    # then the dataframe
    prec_pop <- attr(plans, "prec_pop")
    if (!is.null(prec_pop))
        distr_pop <- pop_tally(matrix(ref_plan, ncol = 1), prec_pop, ndists)
    else
        distr_pop <- rep(NA_real_, ndists)

    if (is.ordered(plans$district)) {
        rg_labels = range(as.integer(as.character(levels(plans$district))))
        if (any(rg_labels != c(1L, attr(plans, "ndists")))) {
            cli_abort(c("Cannot add a reference plan to a set of plans which
                        have relabeled district numbers that don't start at 1.",
                        ">"="Match the district labels on the unmatched plans with
                            {.fn match_numbers}")
            )
        }

        # good to go
        plans$district = as.integer(plans$district)
        cli_inform(c("Coercing {.val district} column to integers.",
                     "i"="You may want to run {.fn match_numbers} again to fix district labels.\n"))
    }

    if (name %in% levels(plans$draw)) cli_abort("Reference plan name already exists")
    fct_levels <- c(name, levels(plans$draw))
    new_draw <- rep(factor(fct_levels, levels = fct_levels), each = ndists)
    x <- dplyr::bind_rows(
        tibble(district = 1:ndists,
               total_pop = as.numeric(distr_pop)),
        plans[, -match("draw", names(plans))]
    ) %>%
        dplyr::mutate(draw = new_draw, .before = "district")

    exist_wgts <- get_plans_weights(plans)
    if (!is.null(exist_wgts))
        attr(plans, "wgt") <- c(0, exist_wgts)

    reconstruct.redist_plans(x, set_plan_matrix(plans, plan_m))
}

#' Subset to sampled or reference draws
#'
#' @param plans the \code{redist_plans} object
#' @param matrix if \code{TRUE}, the default, also subset the plans matrix. If
#' the plans matrix is not needed, turning this off may save some time.
#'
#' @returns a \code{redist_plans} object, with only rows corresponding to
#' simulated (or reference) draws remaining.
#'
#' @concept analyze
#' @export
subset_sampled <- function(plans, matrix = TRUE) {
    plans_m <- get_plans_matrix(plans)
    if (is.null(colnames(plans_m))) return(plans)

    nm_lengths <- nchar(colnames(plans_m))
    draw_ints <- as.integer(plans$draw)
    idxs <- which(nm_lengths[draw_ints] == 0)

    out <- vctrs::vec_slice(plans, idxs)
    out$draw <- droplevels(out$draw)

    idxs <- which(nm_lengths[unique(draw_ints)] == 0)
    attr(out, "wgt") <- attr(out, "wgt")[idxs]
    if (isTRUE(matrix)) {
        out <- set_plan_matrix(out, plans_m[, idxs, drop = FALSE])
    }

    out
}

#' @rdname subset_sampled
#' @export
subset_ref <- function(plans, matrix = TRUE) {
    plans_m <- get_plans_matrix(plans)
    if (is.null(colnames(plans_m))) return(plans)

    nm_lengths <- nchar(colnames(plans_m))
    draw_ints <- as.integer(plans$draw)
    idxs <- which(nm_lengths[draw_ints] > 0)

    out <- vctrs::vec_slice(plans, idxs)
    out$draw <- droplevels(out$draw)

    idxs <- which(nm_lengths[unique(draw_ints)] > 0)
    attr(out, "wgt") <- attr(out, "wgt")[idxs]
    out <- set_plan_matrix(out, plans_m[, idxs, drop = FALSE])

    out
}

#' Extract the Metropolis Hastings Acceptance Rate
#'
#' @param plans the \code{redist_plans} object
#'
#' @returns a numeric acceptance rate
#'
#' @concept analysis
#' @export
get_mh_acceptance_rate <- function(plans) {
    if (!inherits(plans, "redist_plans")) cli_abort("Not a {.cls redist_plans}")
    alg <- attr(plans, "algorithm")

    if (alg %in% c("flip", "mergesplit")) {
        attr(plans, "mh_acceptance")
    } else {
        NA_real_
    }
}

# generics ----------------------------------------------------------------


#' @method dplyr_row_slice redist_plans
#' @export
dplyr_row_slice.redist_plans <- function(data, i, ...) {
    if (is.logical(i)) i <- which(i)

    draws <- rle(as.integer(data$draw))
    draws$values <- seq_along(draws$values)
    draws_left <- unique(inverse.rle(draws)[i])
    y <- vctrs::vec_slice(data, i)
    plans_m <- get_plans_matrix(data)

    # if we don't have every district present in every row
    # this check is necessary but not sufficient for what we want
    if (length(i) > 0 && "district" %in% colnames(y)) {
        distrs <- table(as.integer(y$district))
        ndists <- max(plans_m[, 1])
        if (any(distrs != distrs[1]) || length(distrs) != ndists)
            cli_warn(c("Some districts may have been dropped. This will prevent summary statistics from working correctly.",
                ">" = "To avoid this message, coerce using {.fun as_tibble}."))
    }

    if (length(draws_left) != ncol(plans_m)) {
        attr(y, "wgt") <- attr(y, "wgt")[draws_left]
        y <- set_plan_matrix(y, plans_m[, draws_left, drop = FALSE])
    }

    if (is.factor(y$draw)) {
        y$draw <- droplevels(y$draw)
    }

    y
}

# 'template' is the old df
#' @method dplyr_reconstruct redist_plans
#' @export
dplyr_reconstruct.redist_plans <- function(data, template) {
    reconstruct.redist_plans(data, template)
}

#' Combine multiple sets of redistricting plans
#'
#' Only works when all the sets are compatible---generated from the same map,
#' with the same number of districts.  Sets of plans will be indexed by the
#' `chain` column.
#'
#' @param ... The [`redist_plans`] objects to combine.  If named arguments are
#' provided, the names will be used in the `chain` column; otherwise, numbers
#' will be used for the `chain` column.
#' @param deparse.level Ignored.
#'
#' @return A new [`redist_plans`] object.
#'
#' @md
#' @concept analyze
#' @export
rbind.redist_plans <- function(..., deparse.level = 1) {
    objs <- rlang::list2(...)
    n_obj <- length(objs)
    if (n_obj == 1) return(objs[[1]])

    # check types
    n_prec <- nrow(get_plans_matrix(objs[[1]]))
    prec_pop <- attr(objs[[1]], "prec_pop")
    ndists <- attr(objs[[1]], "ndists")
    constr <- attr(objs[[1]], "constraints")
    resamp <- attr(objs[[1]], "resampled")
    comp <- attr(objs[[1]], "compactness")
    distr_ord <- is.ordered(objs[[1]]$district)
    for (i in 2:n_obj) {
        if (nrow(get_plans_matrix(objs[[i]])) != n_prec)
            cli_abort("Number of precincts must match for all sets of plans.")
        if (!identical(attr(objs[[i]], "prec_pop"), prec_pop))
            cli_abort("Precinct populations must match for all sets of plans.")
        if (!identical(attr(objs[[i]], "ndists"), ndists))
            cli_abort("Number of districts must match for all sets of plans.")
        if (attr(objs[[i]], "resampled") != resamp)
            cli_abort("Some sets of plans are resampled while others are not.")
        if (!is.null(comp)) {
            if (attr(objs[[i]], "compactness") != comp) {
                cli_warn("Compactness values differ across sets of plans.")
                comp <- NA
            }
        } else {
            if (!is.null(attr(objs[[i]], "compactness"))) {
                cli_warn("Some compactness values were non-NULL. Set to {.val NA}.")
                comp <- NA
            }
        }
        if (!identical(attr(objs[[i]], "constraints"), constr)) {
            cli_inform("Constraints may not match for all sets of plans.")
            constr <- NA
        }
        if (is.ordered(objs[[i]]$district) != distr_ord) {
            cli_abort(c("Some sets of plans have had district numbers matched to a reference plan,
                         while others have not. This may cause problems in analysis.",
                        "i"="Do one of the following:",
                        ">"="Match the district labels on the unmatched plans with
                            {.fn match_numbers} [recommended]",
                        ">"="Convert the matched plans district labels to integers with
                            {.code as.integer(district)}")
            )
        }
    }

    # split by chain and then recombine
    split_if <- function(pl) {
        if ("chain" %in% names(pl)) {
            pl$chain[is.na(pl$chain)] = -1L
            split(pl, pl$chain)
        } else {
            list(pl)
        }
    }
    ret <- objs |>
        lapply(split_if) |>
        unlist(recursive = FALSE) |>
        unname() |>
        lapply(dplyr::as_tibble) |>
        dplyr::bind_rows(.id = "chain")

    ret <- reconstruct.redist_plans(ret, objs[[1]])
    attr(ret, "compactness") <- comp
    attr(ret, "constraints") <- constr
    attr(ret, "ndists") <- ndists
    attr(ret, "prec_pop") <- prec_pop
    attr(ret, "diagnostics") <- do.call(c, lapply(objs, function(x) attr(x, "diagnostics")))
    attr(ret, "plans") <- do.call(cbind, lapply(objs, function(x) get_plans_matrix(x)))
    attr(ret, "wgt") <- do.call(c, lapply(objs, function(x) get_plans_weights(x)))
    attr(ret, "n_eff") <- sum(do.call(c, lapply(objs, function(x) attr(x, "n_eff"))))

    ret
}


#' Print method for \code{redist_plans}
#'
#' @param x a [redist_plans] object
#' @param \dots additional arguments (ignored)
#'
#' @method print redist_plans
#' @return The original object, invisibly.
#' @export
print.redist_plans <- function(x, ...) {
    plans_m <- get_plans_matrix(x)
    n_ref <- get_n_ref(x)
    n_samp <- ncol(plans_m) - n_ref
    nd <- attr(x, "ndists")
    if (is.null(nd)) nd <- max(plans_m[, 1])

    fmt_comma <- function(x) format(x, nsmall = 0, digits = 1, big.mark = ",")
    if (n_ref > 0)
        cli_text("A {.cls redist_plans} containing {fmt_comma(n_samp)}{cli::qty(n_samp)}
                 sampled plan{?s} and {n_ref} reference plan{?s}")
    else
        cli_text("A {.cls redist_plans} containing
                 {fmt_comma(n_samp)}{cli::qty(n_samp)} sampled plan{?s}")

    if (ncol(plans_m) == 0) return(invisible(x))

    alg_name <- c(mcmc = "Flip Markov chain Monte Carlo",
        smc = "Sequential Monte Carlo",
        mergesplit = "Merge-split Markov chain Monte Carlo",
        rsg = "random seed-and-grow",
        crsg = "compact random seed-and-grow",
        enumpart = "Enumpart",
        shortburst = "short bursts",
        none = "a custom collection")[attr(x, "algorithm")]
    if (is.na(alg_name)) alg_name <- "an unknown algorithm"

    cli_text("Plans have {nd} district{?s} from a
             {fmt_comma(nrow(plans_m))}{cli::qty(nrow(plans_m))}-unit map,
             and were drawn using {alg_name}.")

    merge_idx <- attr(x, "merge_idx")
    if (!is.null(merge_idx))
        cat("Merged from another map with reindexing:",
            utils::capture.output(str(merge_idx, vec.len = 2)), "\n", sep = "")

    if (!is.null(attr(x, "wgt"))) {
        if (attr(x, "resampled"))
            cat("With plans resampled from weights\n")
        else
            cat("With weighted plans not resampled\n")
    }

    cat("Plans matrix:", utils::capture.output(str(plans_m, give.attr = FALSE)),
        "\n", sep = "")

    utils::getS3method("print", "tbl")(x)

    invisible(x)
}
