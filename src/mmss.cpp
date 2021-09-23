
#include "merge_split.h"

Rcpp::List mmss_plans(int N, List l, const uvec init, const uvec &counties, const uvec &pop,
                    int n_distr, int n_merge, double target, double lower, double upper, double rho,
                    double beta_sq, const uvec &current, int n_current,
                    double beta_vra, double tgt_min, double tgt_other,
                    double pow_vra, const uvec &min_pop,
                    double beta_vra_hinge, const vec &tgts_min,
                    double beta_inc, const uvec &incumbents, double beta_splits,
                    double beta_fractures, double thresh, int k, int verbosity) {

    // re-seed MT
    generator.seed((int) Rcpp::sample(INT_MAX, 1)[0]);

    Graph g = list_to_graph(l);
    Multigraph cg = county_graph(g, counties);
    int V = g.size();
    int n_cty = max(counties);

    umat districts(V, N, fill::zeros);
    districts.col(0) = init;

    Rcpp::IntegerVector mh_decisions(N - 1);

    double tol = std::max(target - lower, upper - target) / target;

    if (verbosity >= 1) {
        Rcout << "MARKOV CHAIN MONTE CARLO\n";
        Rcout << "Sampling " << N-1 << " " << V << "-unit maps with " << n_distr
              << " districts and population between " << lower << " and " << upper << ".\n";
        if (cg.size() > 1)
            Rcout << "Sampling hierarchically with respect to the "
                  << cg.size() << " administrative units.\n";
    }

    // find k and multipliers
    if (k <= 0) {
        adapt_ms_parameters(g, n_distr, k, thresh, tol, init, counties, cg, pop, target);
    }
    if (verbosity >= 2)
        Rcout << "Using k = " << k << "\n";

    IntegerVector distrs(n_merge);
    select_merge(n_distr, g, init, n_merge, distrs);
    int refresh = std::max(N / 20, 1);
    int n_accept = 0;
    int reject_ct;
    for (int i = 1; i < N; i++) {
        districts.col(i) = districts.col(i - 1); // copy over old map

        // make the proposal
        double prop_lp = 0.0;
        reject_ct = 0;
        do {
            select_pair(n_distr, g, districts.col(i), distr_1, distr_2);
            prop_lp = split_map_ms(g, counties, cg, districts.col(i), distr_1,
                                   distr_2, pop, lower, upper, target, k);
            if (reject_ct % 200 == 0) Rcpp::checkUserInterrupt();
            reject_ct++;
        } while (!std::isfinite(prop_lp));

        // tau calculations
        if (rho != 1) {
            double log_st = 0;
            for (int j = 1; j <= n_cty; j++) {
                log_st += log_st_distr(g, districts, counties, i-1, distr_1, j);
                log_st += log_st_distr(g, districts, counties, i-1, distr_2, j);
                log_st -= log_st_distr(g, districts, counties, i, distr_1, j);
                log_st -= log_st_distr(g, districts, counties, i, distr_2, j);
            }
            log_st += log_st_contr(g, districts, counties, n_cty, i-1, distr_1);
            log_st += log_st_contr(g, districts, counties, n_cty, i-1, distr_2);
            log_st -= log_st_contr(g, districts, counties, n_cty, i, distr_1);
            log_st -= log_st_contr(g, districts, counties, n_cty, i, distr_2);

            prop_lp += (1 - rho) * log_st;
        }

        // add gibbs target
        // NOTE: different signs than above b/c of how Metropolis proposal has
        // transition ratio flipped relative to the target density ratio
        prop_lp -= calc_gibbs_tgt(districts.col(i), n_distr, V, distr_1, distr_2,
                                  pop, beta_sq, current, n_current, beta_vra,
                                  tgt_min, tgt_other, pow_vra, min_pop,
                                  beta_vra_hinge, tgts_min,
                                  beta_inc, incumbents,
                                  beta_splits, beta_fractures, counties, n_cty);
        prop_lp += calc_gibbs_tgt(districts.col(i-1), n_distr, V, distr_1, distr_2,
                                  pop, beta_sq, current, n_current, beta_vra,
                                  tgt_min, tgt_other, pow_vra, min_pop,
                                  beta_vra_hinge, tgts_min,
                                  beta_inc, incumbents,
                                  beta_splits, beta_fractures, counties, n_cty);

        double alpha = exp(prop_lp);
        if (alpha >= 1 || unif(generator) <= alpha) { // ACCEPT
            n_accept++;
            // map already stored in districts.col(i);
            mh_decisions(i - 1) = 1;
        } else { // REJECT
            districts.col(i) = districts.col(i - 1); // copy over old map
            mh_decisions(i - 1) = 0;
        }

        if (verbosity >= 2 && refresh > 0 && (i+1) % refresh == 0) {
            Rprintf("Iteration %'6d / %'d\n", i+1, N-1);
        }
        Rcpp::checkUserInterrupt();
    }

    if (verbosity >= 1) {
        Rprintf("Acceptance rate: %.1f%%.\n", (100.0 * n_accept) / (N-1));
    }

    Rcpp::List out;
    out["plans"] = districts;
    out["mhdecisions"] = mh_decisions;

    return out;
}


/*
 * Select a pair of neighboring districts i, j
 */
void select_merge(int n, const Graph &g, const uvec &plan, int &i, IntegerVector &j) {
    int V = g.size();
    i = 1 + rint(n);

    std::set<int> neighboring;
    for (int k = 0; k < V; k++) {
        if (plan(k) != i) continue;
        std::vector<int> nbors = g[k];
        int length = nbors.size();
        for (int l = 0; l < length; l++) {
            int nbor = nbors[l];
            if (plan(nbor) == i) continue;
            neighboring.insert(plan[nbor]);
        }
    }

    int n_nbor = neighboring.size();
    j = *std::next(neighboring.begin(), rint(n_nbor));

    return;
}
