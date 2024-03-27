#' Selects Precincts and Blocks to Move to Equalize Population Between Districts
#'
#' This function, given a map and a plan, finds a set of precincts and blocks
#' to move such that it equalizes the population between the districts defined
#' in the plan. It returns the new assignment of blocks and precincts to a plan.
#'
#' @param map_prec a [redist_map] that is made up of precincts
#' @param map a [redist_map] that is made up of blocks
#' @param plan an integer vector containing the plan to be equalized with one
#' value per precinct
#' @param tolerance an integer value that indicates at what deviation the
#' algorithm should no longer attempt to find a new set of district swaps (this
#' is only used very narrowly, does not apply to the whole algorithm)
#'
#' @returns a dataframe that contains the GEOID of every block, the GEOID of the
#' precinct that the block belongs to, and the district that it is assigned to
#'
#' @examples
#' data(iowa)
#' iowa <- redist_map(iowa, existing_plan = cd_2010)
#' iowa_blk <- get_block_map("IA", 4)
#' #new_plan <- suppressWarnings(walnuts_blk(iowa, iowa_blk, iowa$cd_2010))
#'
#' @concept analyze
#' @md
#' @export
walnuts_blk <- function(map_prec, map_blk, plan, tolerance = 5) {

    map_prec <- validate_redist_map(map_prec)
    map_blk <- validate_redist_map(map_blk)

    if (!is.numeric(plan) && all(plan > 0) && length(plan) == nrow(map_prec))
        cli_abort("{.arg plan} must be a positive integer vector with one entry per precinct.")

    # Run precinct method
    map_prec$plan <- walnuts_prec(map_prec, plan)

    # Calculate transfer
    transfer <- optimal_transfer(map_prec, map_prec$plan)

    # Get block-level plan
    map_blk$plan <- walnuts_plan_to_block(map_prec, map_blk, map_prec$plan)

    # If no more changes to be made return block-level plan
    if (nrow(transfer) == 0 | sum(transfer$pop) %in% c(-1, 0, 1)){
        return(tibble::as_tibble(map_blk) %>% select(GEOID, BLOCKID, plan))
    }

    # Order transfer by largest population
    transfer <- transfer %>%
        arrange(desc(pop))

    line <- 1
    no_transfer <- data.frame()
    rerun <- FALSE

    while (line <= nrow(transfer)){

        # Set pair of districts to transfer population between
        districts <- c(transfer$from[line], transfer$to[line])

        # Get redist_map with just the pair of districts
        new_map <- map_blk %>%
            filter(plan %in% districts)

        # Get all of the precincts on the boundary
        on_boundary <- walnuts_find_boundary_prec(map_prec$adj, map_prec$plan,
                                                  districts[1], districts[2],
                                                  nrow(map_prec))
        map_prec$on_boundary <- on_boundary

        # Order precincts on boundary by size
        prec_on_boundary <- order_gpp_2dists(map_prec, map_blk,
                                             map_blk$plan,
                                             districts[1], districts[2])

        # Run the block method for all the precincts on the boundary
        best_result <- c()
        best_transfer <- 100000
        for (i in prec_on_boundary){

            # Run block method
            result <- walnuts_blk_2dists_mult(new_map, new_map$plan,
                                              transfer$pop[line],
                                              districts[1], districts[2], i)

            # Get deviation, update if better
            if (abs(result$transfer) < abs(best_transfer)){
                best_result <- result
                best_transfer <- result$transfer
            }

            if (best_transfer == 0){
                break
            }
        }

        result <- best_result

        new_map$new_plan <- result$plan

        # Get new plan
        map_blk$new_plan <- NULL
        map_blk <- map_blk %>%
            left_join(tibble::as_tibble(new_map)
                      %>% select(BLOCKID, new_plan), by = "BLOCKID") %>%
            mutate(plan = ifelse(is.na(new_plan), plan, new_plan))

        # Calculate new transfer
        new_transfer <- optimal_transfer(map_blk, map_blk$plan)

        if (nrow(new_transfer) == 0){
            return(tibble::as_tibble(map_blk) %>% select(GEOID, BLOCKID, plan))
        }

        new_transfer <- new_transfer %>%
            arrange(desc(pop))

        # Join transfer
        temp_transfer <- transfer %>%
            left_join(new_transfer, by = c("from", "to")) %>%
            mutate(pop_change = ifelse(pop.x == pop.y, FALSE, TRUE)) %>%
            mutate(pop_change = replace_na(pop_change, TRUE))

        if (sum(temp_transfer$pop_change[1:(line-1)]) != 0){
            line <- 1
        }
        else if (sum(temp_transfer$pop_change) == 0){
            line <- line + 1
        }

        # Set new transfer as transfer
        transfer <- new_transfer

        if(sum(transfer$pop) == 1){
            break
        }

        # Re-run optimal transfer without the offending districts (only try once)
        if (line > nrow(transfer) & max(transfer$pop) > tolerance & !rerun){
            no_transfer <- bind_rows(no_transfer, transfer)
            new_transfer <- optimal_transfer(map_blk, map_blk$plan,
                                             no_transfer = no_transfer)

            if (nrow(new_transfer) != 0){
                line <- 1
                transfer <- new_transfer %>% arrange(desc(pop))
            }
            rerun <- TRUE
        }
    }

    # Returns precinct GEOIDs, block GEOIDs, and plan
    return(tibble::as_tibble(map_blk) %>% select(GEOID, BLOCKID, plan))
}

# Walnuts functions ---------------------------------------------------------

# Walnuts at the precinct-level between 2 districts
walnuts_prec_2dists <- function(map, plan, transfer, dist_1, dist_2, shattered = NA){

    # Store original map
    original_map <- map

    # Convert redist_map to df
    map <- tibble::as_tibble(map)

    # Join plan to redist_map
    map$plan <- plan

    # Add indicator for precincts that will shatter the district
    map$will_shatter <- rep(FALSE, nrow(map))
    map$will_shatter[shattered] <- TRUE

    # Identify the precincts on the boundary of dist_1 to be moved to dist_2
    map$prec_to_swap <- walnuts_find_boundary_prec(map$adj, map$plan,
                                                   dist_1, dist_2, nrow(map))

    # Get list of precinct populations that can be swapped
    prec <- map %>%
        filter(prec_to_swap == TRUE & will_shatter == FALSE) %>%
        filter(pop < transfer) %>%
        arrange(desc(pop)) %>%
        filter(pop != 0)

    prec_pops <- prec %>% pull(pop)

    if (length(prec_pops) == 0){
        return(special_add_prec(map, plan, transfer, dist_1, dist_2, shattered))
    }

    # Subset sum
    target = sss_test(prec_pops, transfer)
    ss = subsetsum(prec_pops, target)$inds
    swapped_prec <- rep(FALSE, length(prec_pops))
    swapped_prec[ss] <- TRUE
    prec$swapped_prec <- swapped_prec

    # Get new plan
    map <- map %>%
        left_join(prec %>% select(GEOID, swapped_prec), by = "GEOID") %>%
        mutate(swapped_prec = replace_na(swapped_prec, FALSE)) %>%
        mutate(plan = ifelse(swapped_prec == TRUE, dist_2, plan))

    plan <- map %>% pull(plan)

    # Check contiguity
    shattered_prec <- walnuts_check_contiguity(map, plan,
                                               c(which(map$swapped_prec == TRUE)))

    if (!is.null(shattered_prec)){

        # Correct the shattered precincts
        map$plan[shattered_prec] <- dist_1
        map$will_shatter[shattered_prec] <- TRUE
        map$swapped_prec[shattered_prec] <- FALSE

        # Find new transfer amount
        transfer <- transfer - (map %>%
                                    filter(will_shatter == FALSE & swapped_prec == TRUE) %>%
                                    pull(pop) %>% sum())

        # Redo subset sum without the precincts that would shatter the district
        return(walnuts_prec_2dists(original_map, map$plan, transfer,
                                   dist_1, dist_2,
                                   shattered = c(shattered, shattered_prec)))
    }

    if(calc_diff(map, plan) == 0){
        return(list(imp = FALSE, plan = plan, transfer = 0))
    }

    # Find new transfer amount
    transfer <- transfer - (map %>%
                                filter(swapped_prec == TRUE) %>%
                                pull(pop) %>%
                                sum())

    # Return that an improvement can still be made to the plan
    return(list(imp = TRUE, plan = plan, transfer = transfer))
}

# Runs walnuts_prec_2dists multiple times until no more improvement
walnuts_prec_2dists_mult <- function(map, plan, transfer, dist_1, dist_2, run_limit = 10){

    result <- walnuts_prec_2dists(map, plan, transfer, dist_1, dist_2)
    alg_run <- 1
    d1 <- dist_1
    d2 <- dist_2

    # Keep running while there is still improvement, 10 runs max as default
    while (result$imp[1] == TRUE & alg_run < run_limit){

        # Swap if the district to transfer from swaps
        if (result$transfer < 0){
            temp <- d1
            d1 <- d2
            d2 <- temp
        }
        result <- walnuts_prec_2dists(map, result$plan, abs(result$transfer),
                                      d1, d2)

        alg_run <- alg_run + 1
    }

    return(result)
}

# Runs walnuts method for just precincts
walnuts_prec <- function(map, plan){

    # Join plan to redist_map
    map$plan <- plan

    # Calculate transfer
    transfer <- optimal_transfer(map, map$plan)

    # Return the plan if there is no more improvement that can be made
    if (nrow(transfer) == 0 | sum(transfer$pop) == 0){
        return(map$plan)
    }

    # Order transfer descending
    transfer <- transfer %>% arrange(desc(pop), desc(from), desc(to))

    line <- 1

    while (line <= nrow(transfer)){

        # Set pair of districts to transfer population between
        districts <- c(transfer$from[line], transfer$to[line])

        # Get redist_map with just the pair of districts
        new_map <- map %>%
            filter(plan %in% districts)

        # Run precinct method on the two districts
        new_map$new_plan <- walnuts_prec_2dists_mult(new_map, new_map$plan,
                                                     transfer$pop[line],
                                                     districts[1], districts[2],
                                                     run_limit = 5)$plan

        # Get new plan
        map$new_plan <- NULL
        map <- map %>%
            left_join(tibble::as_tibble(new_map) %>%
                          select(GEOID, new_plan), by = "GEOID") %>%
            mutate(plan = ifelse(is.na(new_plan), plan, new_plan))

        # Calculate new transfer
        new_transfer <- optimal_transfer(map, map$plan)

        if (nrow(new_transfer) == 0){
            return(map$plan)
        }

        # Join transfer
        temp_transfer <- transfer %>%
            left_join(new_transfer, by = c("from", "to")) %>%
            mutate(pop_change = ifelse(pop.x == pop.y, FALSE, TRUE)) %>%
            mutate(pop_change = replace_na(pop_change, TRUE))

        if (sum(temp_transfer$pop_change[1:(line-1)]) != 0){
            line <- 1
        }
        else if (sum(temp_transfer$pop_change) == 0){
            line <- line + 1
        }

        # Set new transfer as transfer
        transfer <- new_transfer %>% arrange(desc(pop), desc(from), desc(to))
    }

    return(map$plan)
}

# Walnuts at the block level for just one precinct
walnuts_blk_2dists <- function(map, plan, transfer, dist_1, dist_2, gpp, shattered = NA){

    # Store original map
    original_map <- map

    # Create dataframe
    map <- tibble::as_tibble(map)

    # Join plan to redist_map
    map$plan <- plan

    # Add indicator for blocks that will shatter the district
    map$will_shatter <- rep(FALSE, nrow(map))
    map$will_shatter[shattered] <- TRUE

    # Identify the blocks on the boundary of dist_1 to be moved to dist_2
    map$blks_to_swap <- walnuts_find_boundary_blk(map$adj, map$plan,
                                                  dist_1, dist_2,
                                                  nrow(map), map$GEOID, gpp)

    # Get block populations of largest precinct on the boundary
    blks <- map %>%
        filter(GEOID == gpp & blks_to_swap == TRUE) %>%
        filter(will_shatter == FALSE) %>%
        filter(pop != 0) %>%
        filter(pop < transfer) %>%
        arrange(desc(pop))

    blk_pops <- blks %>%
        pull(pop)

    # Case where no blocks that can be swapped are smaller than the transfer value
    if (length(blk_pops) == 0){
        return(special_add_blk(map, plan, transfer, dist_1, dist_2, gpp, shattered))
    }

    # Subset sum
    target = sss_test(blk_pops, transfer)
    ss = subsetsum(blk_pops, target)$inds
    swapped_blks <- rep(FALSE, length(blk_pops))
    swapped_blks[ss] <- TRUE
    blks$swapped_blks <- swapped_blks

    # Get new plan
    map <- map %>%
        left_join(blks %>% select(BLOCKID, swapped_blks), by = "BLOCKID") %>%
        mutate(swapped_blks = replace_na(swapped_blks, FALSE)) %>%
        mutate(plan = ifelse(swapped_blks == TRUE, dist_2, plan))

    plan <- map$plan

    # Check contiguity
    shattered_blks <- walnuts_check_contiguity(map, plan,
                                               c(which(map$swapped_blks == TRUE)))

    if (!is.null(shattered_blks)){

        # Correct the shattered precincts
        map$plan[shattered_blks] <- dist_1
        map$will_shatter[shattered_blks] <- TRUE
        map$swapped_blks[shattered_blks] <- FALSE

        # Find new transfer amount
        transfer <- transfer - (map %>%
                                    filter(will_shatter == FALSE & swapped_blks == TRUE) %>%
                                    pull(pop) %>%
                                    sum())

        # Redo subset sum without the precincts that would shatter the district
        return(walnuts_blk_2dists(original_map, map$plan, transfer,
                                  dist_1, dist_2,
                                  shattered = c(shattered, shattered_blks),
                                  gpp = gpp))
    }

    # Find new transfer amount
    transfer <- transfer - (map %>%
                                filter(swapped_blks == TRUE) %>%
                                pull(pop) %>%
                                sum())

    # Return that the plan no longer needs to be improved
    if(transfer == 0){
        return(list(imp = FALSE, plan = plan, transfer = transfer))
    }

    # Return that an improvement can still be made to the plan
    return(list(imp = TRUE, plan = plan, transfer = transfer))
}

# Runs walnuts_blk_2dists multiple times
walnuts_blk_2dists_mult <- function(map, plan, transfer, dist_1, dist_2, gpp, run_limit = 10){

    result <- walnuts_blk_2dists(map, plan, transfer, dist_1, dist_2, gpp)
    alg_run <- 1
    d1 <- dist_1
    d2 <- dist_2

    # Keep running while there is still improvement, 5 runs max as default
    while (!(result$transfer == 0) & (result$imp[1] == TRUE & alg_run < run_limit)){
        # Swap if the district to transfer from swaps
        if (result$transfer < 0){
            temp <- d1
            d1 <- d2
            d2 <- temp
        }
        result <- walnuts_blk_2dists(map, result$plan, abs(result$transfer),
                                     d1, d2, gpp)
        alg_run <- alg_run + 1
    }

    # Flip transfer value if districts swapped
    if (d1 != dist_1){
        result$transfer <- result$transfer*(-1)
    }

    return(result)
}


# Helper functions ---------------------------------------------------------

# Calculates difference in population for 2 districts
calc_diff <- function(map, plan){

    # Join plan to redist_map
    map$plan <- plan

    # Get population difference between the two districts
    pop_diff <- map %>%
        group_by(plan) %>%
        arrange(plan) %>%
        summarize(pop = sum(pop)) %>%
        pull(pop)

    return(pop_diff[1]-pop_diff[2])
}

# Optimal transfer for 2 districts - like calc_dff but halved and floored
optimal_transfer_2dists <- function(map, plan){
    # Get population difference between the two districts
    pop_diff <- calc_diff(map, plan)

    return(floor((pop_diff)/2))
}

# Find the greatest population precinct in the larger district of two that can be swapped
find_gpp_2dists <- function(map_prec, map_blk, plan, dist_1, dist_2){
    map_blk$to_swap <- walnuts_find_boundary_prec(map_blk$adj, plan, dist_1, dist_2, nrow(map_blk))

    scale <- function(x){(x-min(x))/(max(x)-min(x))}

    gpp <- tibble::as_tibble(map_blk) %>%
        filter(to_swap == TRUE) %>%
        group_by(GEOID) %>%
        summarize(count = n(), pop_on_boundary = sum(pop)) %>%
        left_join(tibble::as_tibble(map_prec) %>% select(GEOID, pop), by = "GEOID") %>%
        mutate_at(c('count', 'pop_on_boundary', 'pop'), scale) %>%
        mutate(score = count*pop_on_boundary*pop) %>%
        arrange(desc(score)) %>%
        head(1) %>%
        pull(GEOID)

    return(gpp)
}

# Scale values - used for ordering gpp
scale <- function(x){(x-min(x))/(max(x)-min(x))}

# Order districts by size in the larger district of two that can be swapped
order_gpp_2dists <- function(map_prec, map_blk, plan, dist_1, dist_2){
    map_blk$to_swap <- walnuts_find_boundary_prec(map_blk$adj, plan, dist_1, dist_2, nrow(map_blk))

    scale <- function(x){(x-min(x))/(max(x)-min(x))}

    gpp <- tibble::as_tibble(map_blk) %>%
        filter(to_swap == TRUE) %>%
        group_by(GEOID) %>%
        summarize(count = n(), pop_on_boundary = sum(pop)) %>%
        left_join(tibble::as_tibble(map_prec) %>% select(GEOID, pop), by = "GEOID") %>%
        mutate_at(c('count', 'pop_on_boundary', 'pop'), scale) %>%
        mutate(score = count*pop_on_boundary*pop) %>%
        arrange(desc(score)) %>%
        pull(GEOID)

    return(gpp)
}


# Creates a census block-level redist map for a state
get_block_map <- function(state_abbr, ndists){

    # Get block-assignment-file
    baf <- PL94171::pl_get_baf(state_abbr, 'VTD')$VTD %>%
        rename(vtd = DISTRICT) %>%
        select(-c(COUNTYFP))

    # Get block file from decennial census and join baf
    map_blk <- censable::build_dec('block', state_abbr) %>%
        rename(BLOCKID = GEOID) %>%
        left_join(baf, by = "BLOCKID") %>%
        mutate(GEOID = paste0(substr(BLOCKID, 1,5), vtd))

    # Get adjacency for census blocks
    adj <- geomander::adjacency(map_blk)

    if (geomander::ccm(adj) != 1){
        suggests <- suggest_component_connection(map_blk, adj)
        adj <- adj %>% add_edge(v1 = suggests$x, v2 = suggests$y)
    }

    # Turn into redist_map
    map_blk <- redist_map(map_blk, adj = adj, pop_tol = 0.01, ndists = ndists)

    return(map_blk)
}

# Creates block-level plan from a precinct-level plan
walnuts_plan_to_block <- function(map_prec, map_blk, plan){

    # Join plan to precinct map then to block-level map
    map_prec$plan <- plan
    plan_blk <- map_blk %>%
        left_join(st_drop_geometry(map_prec), by = "GEOID") %>%
        pull(plan)

    return(plan_blk)
}

# Check preservation of contiguity with swapped precincts
walnuts_check_contiguity <- function(map, plan, prec_to_swap){
    if (is.na(cct(map$adj, plan)[2])){
        return(c())
    }
    else {
        disconnected_prec <- check_contiguity(map$adj, plan) %>%
            cbind(index = 1:nrow(map)) %>%
            filter(component >= 2) %>%
            pull(index)

        prec_to_remove <- c()
        for (i in 1:length(disconnected_prec)){
            for (j in 1:length(prec_to_swap)){
                if (prec_to_swap[j] %in% (map[disconnected_prec[i],]$adj[[1]]+1)){
                    prec_to_remove <- c(prec_to_swap[j], prec_to_remove)
                }
            }
        }
        return(unique(prec_to_remove))
    }
}

# Add the the precinct closest to the population threshold on the boundary of dist_1 if it improves deviation
special_add_prec <- function(map, plan, transfer, dist_1, dist_2, shattered = NA){

    prec <- map %>%
        filter(prec_to_swap == TRUE & pop != 0 & pop < 2*transfer) %>%
        arrange(pop)

    # If there are no more precincts to swap, end
    if (nrow(prec) == 0){
        return(list(imp = FALSE, plan = plan, transfer = transfer))
    }

    index = 1

    while(index < nrow(prec)){
        # Add precinct closest to population threshold if it improves the deviation
        map$swapped_prec <- rep(FALSE, nrow(map))
        map$swapped_prec[which(map$GEOID == prec$GEOID[1])] <- TRUE
        plan[which(map$GEOID == prec$GEOID[1])] <- dist_2

        # Check contiguity
        if (is.null(walnuts_check_contiguity(map, plan, c(which(map$swapped_prec == TRUE))))){

            transfer <- transfer - (map %>% filter(swapped_prec == TRUE) %>% pull(pop) %>% sum())

            # Return that an improvement was made to the the plan
            return(list(imp = TRUE, plan = plan, transfer = transfer))
        }
        else {
            plan[which(map$GEOID == prec$GEOID[1])] <- dist_1
            index = index + 1
        }
        index <- index + nrow(prec)
    }

    # Return that no improvement was made to the plan
    return(list(imp = FALSE, plan = plan, transfer = transfer))
}

# Add the the block closest to the population threshold on the boundary of dist_1 if it improves deviation
special_add_blk <- function(map, plan, transfer, dist_1, dist_2, gpp, shattered = NA){

    blks <- map %>%
        filter(GEOID == gpp & blks_to_swap == TRUE & will_shatter == FALSE) %>%
        arrange(pop)

    # If there are no more blocks to swap, end
    if (nrow(blks) == 0){
        return(list(imp = FALSE, plan = plan, transfer = transfer))
    }

    # If all adjacent blocks have population 0, add them all
    if (sum(blks %>% pull(pop)) == 0){

        for (i in 1:length(blks$pop)){
            plan[which(map$BLOCKID == blks$BLOCKID[i])] <- dist_2
        }

        # Check contiguity
        contig <- walnuts_check_contiguity(map, plan, c(which(map$BLOCKID == blks$BLOCKID[i])))
        if (!is.null(contig)){
            for (i in 1:length(contig)){
                plan[contig[i]] <- dist_1
            }
        }

        # Return that an improvement was made to the the plan
        return(list(imp = TRUE, plan = plan, transfer = transfer))
    }

    # If there are no more blocks to swap, end
    if (nrow(blks %>% filter(pop != 0 & pop < 2*transfer)) == 0){

        # If there are blocks with population 0, swap
        if (0 %in% blks$pop){
            for (i in 1:length(blks$pop)){
                if (blks$pop[i] == 0){
                    plan[which(map$BLOCKID == blks$BLOCKID[i])] <- dist_2
                }
            }

            # Check contiguity
            contig <- walnuts_check_contiguity(map, plan, c(which(map$BLOCKID == blks$BLOCKID[i])))
            if (!is.null(contig)){
                for (i in 1:length(contig)){
                    plan[contig[i]] <- dist_1
                }
            }

            # Return that an improvement was made to the the plan
            return(list(imp = TRUE, plan = plan, transfer = transfer))
        }

        return(list(imp = FALSE, plan = plan, transfer = transfer))
    }

    blks <- blks %>%
        filter(pop != 0 & pop < 2*transfer)

    index = 1

    while(index < nrow(blks)){

        # Add precinct closest to population threshold if it improves the deviation
        map$swapped_blks <- rep(FALSE, nrow(map))
        map$swapped_blks[which(map$BLOCKID == blks$BLOCKID[1])] <- TRUE
        plan[which(map$BLOCKID == blks$BLOCKID[1])] <- dist_2

        # Check contiguity
        if (is.null(walnuts_check_contiguity(map, plan, c(which(map$swapped_blks == TRUE))))){

            transfer <- transfer - (map %>% filter(swapped_blks == TRUE) %>% pull(pop) %>% sum())

            # Return that an improvement was made to the the plan
            return(list(imp = TRUE, plan = plan, transfer = transfer))
        }
        else {
            map$swapped_blks[which(map$BLOCKID == blks$BLOCKID[1])] <- FALSE
            plan[which(map$BLOCKID == blks$BLOCKID[1])] <- dist_1
            index = index + 1
        }

        #index <- index + nrow(blks)
    }

    # Return that no improvement was made to the plan
    return(list(imp = FALSE, plan = plan, transfer = transfer))
}

# Get the max deviation from pop parity for the map and the plan
walnuts_plan_parity <- function(map, plan){
    map$plan <- plan
    pop <- map %>% group_by(plan) %>% summarize(pop = sum(pop)) %>% pull(pop)
    parity <- sum(pop)/length(pop)
    return(max(abs(pop - parity) / parity))
}

# Get the max deviation from pop parity for the map and the plan
walnuts_plan_deviation <- function(map, plan){
    ndists <- length(unique(plan))
    correct_pop <- round(sum(map$pop)/ndists)
    pops <- pop_tally(matrix(plan, ncol = 1), map$pop, length(unique(plan))) - correct_pop
    pops <- abs(pops)
    return(max(pops))
}

# Download data for a state and run n_tests for it based on the alarm_50_state simulated plans
test_parity <- function(state_abbrev, n_dists, n_tests = 100){

    map <- alarm_50state_map(state_abbrev)
    plans <- alarm_50state_plans(state_abbrev)
    plans_mat <- get_plans_matrix(plans)
    map_blk <- get_block_map(state_abbrev, n_dists)

    pb <- progress_bar$new(
        format = "  running [:bar] :percent eta: :eta",
        total = n_tests, clear = FALSE, width= 60)

    parity <- c()

    start_time <- Sys.time()

    for (i in 2:(n_tests+1)){
        pb$tick()
        parity <- c(parity, suppressWarnings(walnuts_plan_deviation(map_blk, walnuts_blk(map, map_blk, plans_mat[,i]))))
    }

    end_time <- Sys.time()

    print((end_time-start_time)/n_tests)
    return(parity)
}
