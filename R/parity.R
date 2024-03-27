#' Calculates Maximum Deviation from Population Parity
#'
#' Computes the deviation from population parity from a plan.
#' Higher values indicate that (at least) a single district in the map deviates
#' from population parity. See Details.
#'
#' @details With a map with \code{pop} representing the populations of each district,
#' the deviation from population parity is given as \code{max(abs(pop - parity) / parity)}
#' where \code{parity = sum(pop)/length(pop)} is the population size for the
#' average district.
#' Therefore, the metric can be thought of as the maximum percent deviation from
#' equal population. For example, a value of 0.03 in this metric indicates that
#' all districts are within 3 percent of population parity.
#'
#' @param plans A matrix with one row for each precinct and one column for each
#' map. Required.
#' @param total_pop A numeric vector with the population for every precinct.
#'
#' @return numeric vector with the population parity for each column
#'
#' @concept analyze
#' @export
redist.parity <- function(plans, total_pop) {
    if (!is.numeric(total_pop)) {
        cli_abort("{.arg total_pop} must be a numeric vector")
    }
    if (!is.matrix(plans)) {
        plans <- matrix(plans, ncol = 1)
    }
    if (!is.matrix(plans)) {
        cli_abort("{.arg plans} must be a matrix")
    }

    if (length(total_pop) != nrow(plans)) {
        cli_abort(".arg plans} and {.arg total_pop} must have same number of precincts.")
    }

    rg <- range(plans[, 1])
    if (rg[1] == 0) {
        plans <- plans + 1
        n_distr <- rg[2] + 1
    } else {
        n_distr <- rg[2]
    }

    max_dev(plans, total_pop, n_distr)
}


#' Calculates Sparse Population Moves to Minimize Population Deviation
#'
#' This function computes a minimal set of population moves (e.g., 5 people from
#' district 1 to district 3) to maximally balance the population between
#' districts. The moves are only allowed between districts that share the
#' territory of a county, so that any boundary adjustments are guaranteed to
#' preserve all unbroken county boundaries.
#'
#' @param map a [redist_map]
#' @param plan an integer vector containing the plan to be balanced.
#' Tidy-evaluated.
#' @param counties an optional vector of counties, whose boundaries will be
#' preserved. Tidy-evaluated.
#' @param penalty the larger this value, the more to encourage sparsity.
#'
#' @returns a list with components:
#' \describe{
#' \item{`moves`}{A tibble describing the population moves}
#' \item{`pop_old`}{The current district populations}
#' \item{`pop_new`}{The district populations after the moves}
#' }
#'
#' @examples
#' data(iowa)
#' iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01)
#' min_move_parity(iowa_map, cd_2010)
#'
#' @concept analyze
#' @md
#' @export
min_move_parity <- function(map, plan, counties = NULL, penalty = 0.2) {
    adj <- get_adj(map)
    V <- length(adj)
    nd <- attr(map, "ndists")

    plan <- eval_tidy(enquo(plan), map)
    if (!is.numeric(plan) && all(plan > 0) && length(plan) == V)
        cli_abort("{.arg plan} must be a positive integer vector with one entry per precinct.")

    if (missing(counties)) {
        counties <- rep(1L, V)
    } else {
        counties <- vctrs::vec_group_id(eval_tidy(enquo(counties), map))
    }

    distr_adj <- get_plan_graph(adj, length(adj), plan, nd)
    edges <- do.call(rbind, lapply(seq_along(distr_adj), function(i) {
        tibble(from = i, to = distr_adj[[i]] + 1L)
    })) %>%
        rowwise() %>%
        filter(.data$from < .data$to,
            any(unique(counties[plan == .data$from]) %in% counties[plan == .data$to])) %>%
        ungroup()

    n_edge <- nrow(edges)
    e_idx <- as.matrix(mutate(edges, i = row_number()))
    diff_mat <- matrix(0, nrow = nd, ncol = n_edge)
    diff_mat[e_idx[, -2]] <- -1
    diff_mat[e_idx[, -1]] <- 1

    pops <- pop_tally(matrix(plan, ncol = 1), map[[attr(map, "pop_col")]], nd)
    discrep <- pops - mean(pops)

    fn_balance <- function(move, alpha = 0.1) {
        sum(abs(discrep + diff_mat %*% move)) + alpha*sum(abs(move))
    }
    gr_balance <- function(move, alpha = 0.1) {
        t(sign(discrep + diff_mat %*% move)) %*% diff_mat + alpha*sign(move)
    }

    res <- optim(rep(0, n_edge), fn_balance, gr_balance, alpha = penalty,
        method = "BFGS", control = list(maxit = 1e3, reltol = 1e-9, abstol = 1))

    move <- round(res$par)
    pops_new <- pops + diff_mat %*% move
    from_old <- edges$from
    edges <- mutate(edges,
        from = if_else(move < 0, .data$to, .data$from),
        to = if_else(move < 0, from_old, .data$to),
        move = abs(move)) %>%
        filter(.data$move > 0)

    list(moves = edges,
        pop_old = pops[, 1],
        pop_new = pops_new[, 1])
}

#' Calculates Sparse Population Moves to Reduce Population Deviation to 1
#'
#' This function computes a minimal set of population moves (e.g., 5 people from
#' district 1 to district 3) to balance the population between districts such
#' that the deviation is at most 1 between districts.
#'
#' @param map a [redist_map]
#' @param plan an integer vector containing the plan to be balanced.
#' @param county_splits a boolean value that indicates whether to recommend a
#' set of moves in which no additional counties are split (in testing)
#' @param no_transfer a dataframe consisting of a from, to, and pop column that
#' lists the districts between which no transfer should occur
#'
#' @returns a dataframe consisting of a from, to, and pop column that describes
#' how much population should be moved from/to each district
#'
#' @examples
#' data(iowa)
#' optimal_transfer(redist_map(iowa, existing_plan = cd_2010), iowa$cd_2010)
#'
#' @concept analyze
#' @md
#' @export
optimal_transfer <- function(map, plan, county_splits = FALSE, no_transfer = data.frame()){

    map <- validate_redist_map(map)

    if (!is.numeric(plan) && all(plan > 0) && length(plan) == nrow(map))
        cli_abort("{.arg plan} must be a positive integer vector with one entry per precinct.")

    n_dists = length(unique(plan))
    move_list = gen_move_list(map, plan, n_dists)
    adj_matrix = gen_adj_matrix(map, plan, county_splits, no_transfer)

    # Check that some population can be moved in and out of each district
    adj_mat_check <- rowSums(adj_matrix) + colSums(adj_matrix)
    no_pop_move <- c()
    for (i in 1:n_dists){
        if (adj_mat_check[i] == 0){
            no_pop_move <- c(no_pop_move, i)
        }
    }

    # Remove districts for which no population can be moved
    if (length(no_pop_move) != 0){
        for (i in no_pop_move){
            # Decrease districts by 1
            n_dists <- n_dists - 1
            # Remove the district from the adjacency matrix
            adj_matrix <- adj_matrix[-c(i), -c(i)]
            # Generate a new move list
            move_list = move_list[-c(i)]
        }
    }

    mat = solve_transfer(move_list, adj_matrix)

    # Replace districts for which no population can be moved
    for (i in no_pop_move){
        # Increase districts by 1
        n_dists <- n_dists + 1
        # Add district back into the matrix with all 0s
        if (n_dists != i){
            mat <- rbind(mat[1:i, ], rep(0, n_dists-1), mat[-(1:i), ])
            mat <- cbind(mat[,1:i], rep(0, n_dists), mat[,-(1:i)])
        }
        else {
            mat <- rbind(mat[1:(i-1), ], rep(0, n_dists-1))
            mat <- cbind(mat[,1:(i-1)], rep(0, n_dists))
        }

    }

    # Turn final matrix into a dataframe
    out = data.frame()
    for (i in 1:n_dists){
        for (j in 1:n_dists){
            if (mat[i,j] != 0){
                out <- rbind(out, c(i, j, mat[i,j]))
            }
        }
    }
    if(nrow(out) > 0){
        colnames(out) <- c("from", "to", "pop")
    }
    return(out)
}

# Generate a list of all of the district population deviations
gen_move_list <- function(redist_map, plan, n_dists){
    map <- tibble::as_tibble(redist_map)
    map$plan <- plan
    pop_ideal <- floor(sum(redist_map$pop)/n_dists)
    move_list <- map %>%
        group_by(plan) %>%
        summarize(pop_dev = pop_ideal - sum(pop)) %>%
        arrange(plan) %>%
        pull(pop_dev)

    # Adjust the population of the districts to equal to the total population
    # of the state
    counter = 1;
    while(sum(move_list) < 0){
        move_list[counter] = move_list[counter] + 1
        counter = counter + 1
    }

    return(move_list)
}

# Generate an adjacency matrix for the districts
gen_adj_matrix <- function(redist_map, plan, county_splits = FALSE, no_transfer = data.frame()){
    # Number of districts
    n_dists = length(unique(plan))

    # Join plan to redist_map
    redist_map$plan <- plan

    # Create adjacency list representation for districts
    adj_list <- get_plan_graph(redist_map$adj, nrow(redist_map), plan, n_dists)

    if (county_splits){
        # Create list of counties in each district
        counties <- tibble::as_tibble(redist_map) %>%
            group_by(plan) %>%
            summarize(county = list(unique(county)))

        # Convert adjacency list into adjacency matrix
        adj_matrix = matrix(rep(0, n_dists^2), nrow = n_dists)
        for (i in 1:n_dists){
            for (num in adj_list[[i]]){
                # Check for an overlapping county
                if (any(counties$county[[i]] %in% counties$county[[num + 1]])){
                    adj_matrix[i, num + 1] = 1
                }
            }
        }
    }
    else {
        # Convert adjacency list into adjacency matrix
        adj_matrix = matrix(rep(0, n_dists^2), nrow = n_dists)
        for (i in 1:n_dists){
            for (num in adj_list[[i]]){
                adj_matrix[i, num + 1] = 1
            }
        }
    }

    if (nrow(no_transfer) != 0){
        for (i in 1:nrow(no_transfer)){
            adj_matrix[no_transfer$from[i], no_transfer$to[i]] = 0
            adj_matrix[no_transfer$to[i], no_transfer$from[i]] = 0
        }
    }

    return(adj_matrix)
}

# Helper function to convert vectors to matrices
vec_2_mat = function(x){
    n = trunc(0.5 * (1 + sqrt(1 + 4*length(x))))
    X = matrix(0, n, n)
    X[lower.tri(X) | upper.tri(X)] = x
    return(X)
}

# Construct a constraint matrix for the optimization
get_const_matrix = function(n, M = NULL){
    stopifnot(n>=2);

    mask = rep(1L, n*(n-1))
    if(!is.null(M)){
        mask = as.numeric(M[lower.tri(M) | upper.tri(M)] != 0);
    }

    block = rbind(rep(1L, n-1),
                  -diag(n-1));

    b_list = vector(mode = "list", length = n);
    b_list[[1]] = block;

    for(i in 2:n){
        block[(i-1):i, ] = block[i:(i-1), ];
        b_list[[i]] = block;
    }

    return(sweep(do.call(cbind, b_list), MARGIN=2, mask, "*"));
}

# Run linear program to get optimal transfer
solve_transfer <- function(move_list, adj_matrix){

    n = length(move_list)
    k = n*(n-1)
    C = get_const_matrix(n, adj_matrix)

    # solve lp
    fit_lp = lpSolve::lp(
        objective.in = rep(1, k),
        const.mat = C,
        const.dir = rep("==", n),
        const.rhs = move_list
    )

    return(vec_2_mat(fit_lp$solution))
}
