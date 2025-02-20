/********************************************************
 * Author: Cory McCartan
 * Institution: Harvard University
 * Date Created: 2021/02
 * Purpose: Merge-split MCMC redistricting sampler
 * (like Carter et al. 2019 but using the SMC proposal)
 ********************************************************/

#include "merge_split.h"

/*
 * Main entry point.
 *
 * USING MCMC
 * Sample `N` redistricting plans on map `g`, ensuring that the maximum
 * population deviation is between `lower` and `upper` (and ideally `target`)
 */
Rcpp::List ms_plans(int N, List l, const uvec init, const uvec &counties, const uvec &pop,
              int n_distr, double target, double lower, double upper, double rho,
              List constraints, List control, int k, int thin, int verbosity) {
    // re-seed MT
    seed_rng((int) Rcpp::sample(INT_MAX, 1)[0]);

    // unpack control params
    double thresh = (double) control["adapt_k_thresh"];
    bool do_mh = (bool) control["do_mh"];

    Graph g = list_to_graph(l);
    Multigraph cg = county_graph(g, counties);
    int V = g.size();
    int n_cty = max(counties);

    int rounded_up_N_over_thin = std::ceil(static_cast<double>(std::max(1, N-1))/thin);
    int n_out = rounded_up_N_over_thin + 2;
    umat districts(V, n_out, fill::zeros);
    districts.col(0) = init;
    districts.col(1) = init;

    Rcpp::IntegerVector mh_decisions(rounded_up_N_over_thin + 1);
    double mha;

    double tol = std::max(target - lower, upper - target) / target;

    if (verbosity >= 1) {
        Rcout.imbue(std::locale::classic());
        Rcout << "MARKOV CHAIN MONTE CARLO\n";
        Rcout << std::fixed << std::setprecision(0);
        Rcout << "Sampling " << N << " " << V << "-unit maps with " << n_distr
              << " districts and population between " << lower << " and " << upper << ".\n";
        if (cg.size() > 1)
            Rcout << "Sampling hierarchically with respect to the "
                  << cg.size() << " administrative units.\n";
    }

    // find k and multipliers
    if (k <= 0) {
        adapt_ms_parameters(g, n_distr, k, thresh, tol, init, counties, cg, pop, target);
    }
    if (verbosity >= 3)
        Rcout << "Using k = " << k << "\n";

    Graph dist_g = district_graph(g, init, n_distr);
    Graph new_dist_g;
    int distr_1, distr_2;
    select_pair(n_distr, dist_g, distr_1, distr_2);
    int n_accept = 0;
    int reject_ct;
    CharacterVector psi_names = CharacterVector::create(
        "pop_dev", "splits", "multisplits", "total_splits",
        "segregation", "grp_pow", "grp_hinge", "grp_inv_hinge",
        "compet", "status_quo", "incumbency",
        "polsby", "fry_hold", "log_st", "edges_removed",
        "qps", "custom"
    );
    NumericVector new_psi(psi_names.size());
    std::vector<int> distr_1_2;
    new_psi.names() = psi_names;
    RObject bar = cli_progress_bar(N, cli_config(false));
    int idx = 1;
    Tree ust = init_tree(V);
    std::vector<bool> visited(V);
    std::vector<bool> ignore(V);
    for (int i = 1; i <= N; i++) {
        // make the proposal
        double prop_lp = 0.0;
        mh_decisions(idx - 1) = 0;
        // copy old map to 'working' memory in `idx+1`
        districts.col(idx+1) = districts.col(idx);

        select_pair(n_distr, dist_g, distr_1, distr_2);
        prop_lp = split_map_ms(g, ust, counties, cg, districts.col(idx+1),
                               distr_1, distr_2, visited, ignore,
                               pop, lower, upper, target, k);

        if (!std::isfinite(prop_lp)) {
            districts.col(idx+1) = districts.col(idx);
            if (i % thin == 0) idx++;
            continue; // reject
        }

        // tau calculations
        if (rho != 1) {
            double log_st = 0;
            for (int j = 1; j <= n_cty; j++) {
                log_st += log_st_distr(g, districts, counties, idx, distr_1, j);
                log_st += log_st_distr(g, districts, counties, idx, distr_2, j);
                log_st -= log_st_distr(g, districts, counties, idx+1, distr_1, j);
                log_st -= log_st_distr(g, districts, counties, idx+1, distr_2, j);
            }
            log_st += log_st_contr(g, districts, counties, n_cty, idx, distr_1);
            log_st += log_st_contr(g, districts, counties, n_cty, idx, distr_2);
            log_st -= log_st_contr(g, districts, counties, n_cty, idx+1, distr_1);
            log_st -= log_st_contr(g, districts, counties, n_cty, idx+1, distr_2);

            prop_lp += (1 - rho) * log_st;
        }

        // add gibbs target
        // NOTE: different signs than above b/c of how Metropolis proposal has
        // transition ratio flipped relative to the target density ratio
        distr_1_2 = {distr_1, distr_2};

        prop_lp -= calc_gibbs_tgt(districts.col(idx+1), n_distr, V, distr_1_2, new_psi,
                                  pop, target, g, constraints);
        prop_lp += calc_gibbs_tgt(districts.col(idx), n_distr, V, distr_1_2, new_psi,
                                  pop, target, g, constraints);

        // adjust for prob of picking district pair
        new_dist_g = district_graph(g, districts.col(idx+1), n_distr); // update district graph
        prop_lp -= std::log(
            1.0/dist_g[distr_1 - 1].size() + 1.0/dist_g[distr_2 - 1].size()
        );
        prop_lp += std::log(
            1.0/new_dist_g[distr_1 - 1].size() + 1.0/new_dist_g[distr_2 - 1].size()
        );

        if (!do_mh || prop_lp >= 0 || std::log(r_unif()) <= prop_lp) { // ACCEPT
                n_accept++;
                districts.col(idx) = districts.col(idx+1); // copy over new map
                dist_g = new_dist_g;
                mh_decisions(idx - 1) = 1;
        } else { // reject
            districts.col(idx+1) = districts.col(idx);
        }

        if (i % thin == 0) idx++;

        if (verbosity >= 1 && CLI_SHOULD_TICK) {
            cli_progress_set(bar, i - 1);
            mha = (double) n_accept / (i - 1);
            cli_progress_set_format(bar, "{cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta} | MH Acceptance: %.2f", mha);
        }
        if (idx == n_out - 1) { // thin doesn't divide N and we are done early
            cli_progress_set(bar, N);
            break;
        }
        Rcpp::checkUserInterrupt();
    }
    cli_progress_done(bar);


    if (verbosity >= 1) {
        Rcout << "Acceptance rate: " << std::setprecision(2) << (100.0 * n_accept) / (N-1) << "%\n";
    }

    Rcpp::List out;
    out["plans"] = districts;
    out["est_k"] = k;
    out["mhdecisions"] = mh_decisions;

    return out;
}

/*
 * Split a map into two pieces with population lying between `lower` and `upper`
 */
double split_map_ms(const Graph &g, Tree &ust, const uvec &counties, Multigraph &cg,
                    subview_col<uword> districts, int distr_1, int distr_2,
                     std::vector<bool> &visited, std::vector<bool> &ignore,
                    const uvec &pop, double lower, double upper, double target,
                    int k) {
    int V = g.size();
    double orig_lb = log_boundary(g, districts, distr_1, distr_2);

    clear_tree(ust);
    double total_pop = 0;
    for (int i = 0; i < V; i++) {
        if (districts(i) == distr_1 || districts(i) == distr_2) {
            total_pop += pop(i);
            ignore[i] = false;
        } else {
            ignore[i] = true;
        }
    }

    int root;
    int result = sample_sub_ust(g, ust, V, root, visited, ignore,
                                pop, lower, upper, counties, cg);
    if (result != 0) return -log(0.0);

    // set `lower` as a way to return population of new district
    bool success = cut_districts_ms(ust, k, root, districts, distr_1, distr_2,
                                    pop, total_pop, lower, upper, target);

    if (!success) return -log(0.0); // reject sample

    return orig_lb - log_boundary(g, districts, distr_1, distr_2);
}


/*
 * Cut district into two pieces of roughly equal population
 */
// TESTED
bool cut_districts_ms(Tree &ust, int k, int root, subview_col<uword> &districts,
                      int distr_1, int distr_2, const uvec &pop, double total_pop,
                      double lower, double upper, double target) {
    int V = ust.size();
    // in case we pick a small-V district
    k = std::max(std::min(k, V-3), 1);
    // create list that points to parents & computes population below each vtx
    std::vector<int> pop_below(V, 0);
    std::vector<int> parent(V);
    parent[root] = -1;
    tree_pop(ust, root, pop, pop_below, parent);
    // compile a list of:
    std::vector<int> candidates; // candidate edges to cut,
    std::vector<double> deviances; // how far from target pop.
    std::vector<bool> is_ok; // whether they meet constraints
    int distr_root = districts(root);
    for (int i = 0; i < V; i++) {
        if (districts(i) != distr_root || i == root) continue;
        double below = pop_below.at(i);
        double dev1 = std::abs(below - target);
        double dev2 = std::abs(total_pop - below - target);
        candidates.push_back(i);
        deviances.push_back(std::max(dev1, dev2));
        is_ok.push_back(lower <= below && below <= upper &&
            lower <= total_pop - below && total_pop - below <= upper);
    }
    if ((int) candidates.size() < k) return false;

    int idx = r_int(k);
    idx = select_k(deviances, idx + 1);
    int cut_at = candidates[idx];
    // reject sample
    if (!is_ok[idx]) return false;

    // find index of node to cut at
    std::vector<int> *siblings = &ust[parent[cut_at]];
    int length = siblings->size();
    int j;
    for (j = 0; j < length; j++) {
        if ((*siblings)[j] == cut_at) break;
    }

    siblings->erase(siblings->begin()+j); // remove edge
    parent[cut_at] = -1;

    if (distr_root == distr_1) {
        assign_district(ust, districts, root, distr_1);
        assign_district(ust, districts, cut_at, distr_2);
    } else {
        assign_district(ust, districts, root, distr_2);
        assign_district(ust, districts, cut_at, distr_1);
    }

    return true;
}


/*
 * Choose k and multiplier for efficient, accurate sampling
 */
void adapt_ms_parameters(const Graph &g, int n_distr, int &k, double thresh,
                         double tol, const uvec &plan, const uvec &counties,
                         Multigraph &cg, const uvec &pop, double target) {
    // sample some spanning trees and compute deviances
    int V = g.size();
    Graph dist_g = district_graph(g, plan, n_distr);
    int k_max = std::min(20 + ((int) std::sqrt(V)), V - 1); // heuristic
    int N_adapt = (int) std::floor(4000.0 / sqrt((double) V));

    double lower = target * (1 - tol);
    double upper = target * (1 + tol);

    std::vector<std::vector<double>> devs;
    vec distr_ok(k_max+1, fill::zeros);
    int root;
    int max_ok = 0;
    std::vector<bool> ignore(V);
    std::vector<bool> visited(V);
    int distr_1, distr_2;
    int max_V = 0;
    Tree ust = init_tree(V);
    for (int i = 0; i < N_adapt; i++) {
        double joint_pop = 0;
        select_pair(n_distr, dist_g, distr_1, distr_2);
        int n_vtx = 0;
        for (int j = 0; j < V; j++) {
            if (plan(j) == distr_1 || plan(j) == distr_2) {
                joint_pop += pop(j);
                ignore[j] = false;
                n_vtx++;
            } else {
                ignore[j] = true;
            }
        }
        if (n_vtx > max_V) max_V = n_vtx;

        clear_tree(ust);
        int result = sample_sub_ust(g, ust, V, root, visited, ignore,
                                    pop, lower, upper, counties, cg);
        if (result != 0) {
            i--;
            continue;
        }

        devs.push_back(tree_dev(ust, root, pop, joint_pop, target));
        int n_ok = 0;
        for (int j = 0; j < V-1; j++) {
            if (ignore[j]) devs.at(i).at(j) = 2; // force not to work
            n_ok += devs.at(i).at(j) <= tol;
        }

        if (n_ok <= k_max)
            distr_ok(n_ok) += 1.0 / N_adapt;
        if (n_ok > max_ok && n_ok < k_max)
            max_ok = n_ok;
    }

    // For each k, compute pr(selected edge within top k),
    // among maps where valid edge was selected
    uvec idxs(N_adapt);
    for (k = 1; k <= k_max; k++) {
        idxs = as<uvec>(Rcpp::sample(k, N_adapt, true, R_NilValue, false));
        double sum_within = 0;
        int n_ok = 0;
        for (int i = 0; i < N_adapt; i++) {
            double dev = devs.at(i).at(idxs[i]);
            if (dev > tol) continue;
            else n_ok++;
            for (int j = 0; j < N_adapt; j++) {
                sum_within += ((double) (dev <= devs.at(j).at(k-1))) / N_adapt;
            }
        }
        if (sum_within / n_ok >= thresh) break;
    }

    if (k == k_max + 1) {
        Rcerr << "Warning: maximum hit; falling back to naive k estimator.\n";
        k = max_ok + 1;
    }

    k = std::min(k, max_V - 1);
}

/*
 * Select a pair of neighboring districts i, j
 */
void select_pair(int n_distr, const Graph &dist_g, int &i, int &j) {
    i = r_int(n_distr);
    std::vector<int> nbors = dist_g[i];
    j = nbors[r_int(nbors.size())] + 1;
    i++;
}
