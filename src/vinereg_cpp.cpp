#include "select_dvine.hpp"
#include <wrappers.hpp>
#include <RcppThread.h>

using namespace vinereg;

// [[Rcpp::export]]
Rcpp::List select_dvine_cpp(const Eigen::MatrixXd& data,
                            const std::vector<std::string>& var_types)
{
    auto fc = FitControlsVinecop({BicopFamily::gaussian});
    fc.set_selection_criterion("aic");
    DVineSelectStatus current_fit(data, var_types, fc);
    std::vector<std::vector<Bicop>> pcs;

    std::mutex m;
    RcppThread::ThreadPool pool;
    while (current_fit.get_selected_vars().size() < var_types.size() - 1) {
        auto old_fit = current_fit;
        auto fit_replace_if_better = [&](size_t var) {
            auto new_fit = old_fit;
            new_fit.extend_fit(var);
            std::lock_guard<std::mutex> lk(m);
            if (new_fit.get_crit() > current_fit.get_crit()) {
                current_fit = std::move(new_fit);
            }
        };

        auto rv = old_fit.get_remaining_vars();
        pool.map(fit_replace_if_better, rv);
        pool.wait();

        if (current_fit.get_selected_vars() == old_fit.get_selected_vars()) {
            break;
        }

        auto p_sel = current_fit.get_selected_vars().size();
        pcs.push_back(std::vector<Bicop>{ current_fit.get_pc(p_sel - 1) });
        for (size_t t = 0; t < p_sel - 1; t++) {
            pcs[t].push_back(current_fit.get_pc(t));
        }
    }
    pool.join();

    Rcpp::List vinecop_r;
    if (current_fit.get_selected_vars().size() > 0) {
        auto new_order = tools_stl::cat({ 0 }, current_fit.get_selected_vars());
        new_order = tools_stl::get_order(new_order);
        for (auto& o : new_order)
            o += 1;
        auto new_struct = DVineStructure(new_order);

        auto vine_structure = rvine_structure_wrap(new_struct);
        auto pair_copulas = pair_copulas_wrap(pcs, new_order.size(), true);
        double npars = Vinecop(pcs, new_struct).get_npars();
        double threshold = 0;
        double loglik = NAN;

        auto sel_vt = std::vector<std::string>();
        sel_vt.push_back(var_types[0]);
        for (auto v : current_fit.get_selected_vars())
            sel_vt.push_back(var_types[v]);

        vinecop_r = Rcpp::List::create(
            Rcpp::Named("pair_copulas") = pair_copulas,
            Rcpp::Named("structure")    = vine_structure,
            Rcpp::Named("var_types")    = sel_vt,
            Rcpp::Named("npars")        = npars,
            Rcpp::Named("nobs")         = data.rows(),
            Rcpp::Named("loglik")       = loglik,
            Rcpp::Named("threshold")    = threshold
        );
        vinecop_r.attr("class") = Rcpp::CharacterVector{"vinecop", "vinecop_dist"};
    }
    return Rcpp::List::create(
        Rcpp::Named("vine") = vinecop_r,
        Rcpp::Named("selected_vars") = current_fit.get_selected_vars()
    );
}

// [[Rcpp::export]]
std::vector<Eigen::VectorXd> cond_quantile_cpp(const Eigen::VectorXd& alpha,
                                               const Eigen::MatrixXd& u,
                                               const Rcpp::List& vinecop_r,
                                               size_t num_threads)
{
    tools_eigen::check_if_in_unit_cube(u);
    auto vinecop_cpp = vinecop_wrap(vinecop_r);
    auto vine_struct_ = vinecop_cpp.get_rvine_structure();
    auto d = vine_struct_.get_dim();
    auto var_types_ = vinecop_cpp.get_var_types();
    if ((static_cast<size_t>(u.cols()) != d) &&
        (static_cast<size_t>(u.cols()) != 2 * d))
        throw std::runtime_error("data dimension is incompatible with model.");

    auto trunc_lvl = vine_struct_.get_trunc_lvl();
    auto order = vine_struct_.get_order();
    auto inverse_order = tools_stl::invert_permutation(order);
    std::vector<Eigen::VectorXd> q(alpha.size());
    for (auto& qq : q)
        qq.resize(u.rows());

    auto do_batch = [&](const tools_batch::Batch& b) {
        TriangularArray<Eigen::VectorXd> hfunc2(d + 1, trunc_lvl + 1);
        TriangularArray<Eigen::VectorXd> hfunc1(d + 1, trunc_lvl + 1);
        TriangularArray<Eigen::VectorXd> hfunc2_sub(d + 1, trunc_lvl + 1);
        TriangularArray<Eigen::VectorXd> hfunc1_sub(d + 1, trunc_lvl + 1);

        // data have to be reordered to correspond to natural order
        for (size_t j = 0; j < d; ++j) {
            hfunc2(0, j) = u.col(order[j] - 1).segment(b.begin, b.size);
            hfunc1(0, j) = u.col(order[j] - 1).segment(b.begin, b.size);
            if (var_types_[order[j] - 1] == "d") {
                hfunc2_sub(0, j) = u.col(d + order[j] - 1).segment(b.begin, b.size);
                hfunc1_sub(0, j) = u.col(d + order[j] - 1).segment(b.begin, b.size);
            }
        }

        Eigen::MatrixXd u_e, u_e_sub;
        for (size_t tree = 0; tree < trunc_lvl; ++tree) {
            tools_interface::check_user_interrupt(u.rows() * d > 1e5);
            for (size_t edge = 1; edge < d - tree - 1; ++edge) {
                tools_interface::check_user_interrupt(d * u.rows() > 1e5);
                auto edge_copula = vinecop_cpp.get_pair_copula(tree, edge);
                auto var_types = edge_copula.get_var_types();
                size_t m = vine_struct_.min_array(tree, edge);

                u_e = Eigen::MatrixXd(b.size, 2);
                u_e.col(0) = hfunc2(tree, edge);
                u_e.col(1) = hfunc1(tree, m - 1);
                if ((var_types[0] == "d") | (var_types[1] == "d")) {
                    u_e.conservativeResize(b.size, 4);
                    u_e.col(2) =
                        (var_types[0] == "d") ? hfunc2_sub(tree, edge) : hfunc2(tree, edge);
                    u_e.col(3) = (var_types[1] == "d") ? hfunc1_sub(tree, m - 1)
                        : hfunc1(tree, m - 1);
                }

                if (vine_struct_.needed_hfunc1(tree, edge)) {
                    hfunc1(tree + 1, edge) = edge_copula.hfunc1(u_e);
                    if (var_types[1] == "d") {
                        u_e_sub = u_e;
                        u_e_sub.col(1) = u_e.col(3);
                        hfunc1_sub(tree + 1, edge) = edge_copula.hfunc1(u_e_sub);
                    }
                }
                hfunc2(tree + 1, edge) = edge_copula.hfunc2(u_e);
                if (var_types[0] == "d") {
                    u_e_sub = u_e;
                    u_e_sub.col(0) = u_e.col(2);
                    hfunc2_sub(tree + 1, edge) = edge_copula.hfunc2(u_e_sub);
                }
            }
        }

        for (size_t a = 0; a < static_cast<size_t>(alpha.size()); a++) {
            hfunc2(trunc_lvl, 0) = Eigen::VectorXd::Constant(b.size, alpha[a]);
            for (ptrdiff_t tree = trunc_lvl - 1; tree >= 0; --tree) {
                tools_interface::check_user_interrupt(d * u.rows() > 1e5);
                Bicop edge_copula =
                    vinecop_cpp.get_pair_copula(tree, 0).as_continuous();
                Eigen::MatrixXd U_e(b.size, 2);
                U_e.col(0) = hfunc2(tree + 1, 0);
                U_e.col(1) = hfunc1(tree, 1);
                hfunc2(tree, 0) = edge_copula.hinv2(U_e);
            }
            q[a].segment(b.begin, b.size) = hfunc2(0, 0);
        }
    };

    if (trunc_lvl > 0) {
        tools_thread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
        pool.map(do_batch, tools_batch::create_batches(u.rows(), num_threads));
        pool.join();
    }

    return q;
}

// [[Rcpp::export]]
Eigen::VectorXd cond_dist_cpp(const Eigen::MatrixXd& u,
                              const Rcpp::List& vinecop_r,
                              size_t num_threads)
{
    tools_eigen::check_if_in_unit_cube(u);
    auto vinecop_cpp = vinecop_wrap(vinecop_r);
    auto vine_struct_ = vinecop_cpp.get_rvine_structure();
    auto d = vine_struct_.get_dim();
    auto var_types_ = vinecop_cpp.get_var_types();
    if ((static_cast<size_t>(u.cols()) != d) &&
        (static_cast<size_t>(u.cols()) != 2 * d))
        throw std::runtime_error("data dimension is incompatible with model.");

    auto trunc_lvl = vine_struct_.get_trunc_lvl();
    auto order = vine_struct_.get_order();

    Eigen::VectorXd p(u.rows());
    auto do_batch = [&](const tools_batch::Batch& b) {
        Eigen::MatrixXd hfunc1, hfunc2, u_e, hfunc1_sub, hfunc2_sub, u_e_sub;
        hfunc1 = Eigen::MatrixXd::Zero(b.size, d);
        hfunc2 = Eigen::MatrixXd::Zero(b.size, d);
        hfunc1_sub = hfunc1;
        hfunc2_sub = hfunc2;

        // fill first row of hfunc2 matrix with evaluation points;
        // points have to be reordered to correspond to natural order
        for (size_t j = 0; j < d; ++j) {
            hfunc2.col(j) = u.block(b.begin, order[j] - 1, b.size, 1);
            if (var_types_[order[j] - 1] == "d") {
                hfunc2_sub.col(j) = u.block(b.begin, d + order[j] - 1, b.size, 1);
            }
        }

        for (size_t tree = 0; tree < trunc_lvl; ++tree) {
            tools_interface::check_user_interrupt(u.rows() * d > 1e5);
            for (size_t edge = 0; edge < d - tree - 1; ++edge) {
                tools_interface::check_user_interrupt(edge % 100 == 0);
                // extract evaluation point from hfunction matrices (have been
                // computed in previous tree level)
                Bicop edge_copula = vinecop_cpp.get_pair_copula(tree, edge);
                auto var_types = edge_copula.get_var_types();
                size_t m = vine_struct_.min_array(tree, edge);

                u_e = Eigen::MatrixXd(b.size, 2);
                u_e.col(0) = hfunc2.col(edge);
                if (m == vine_struct_.struct_array(tree, edge, true)) {
                    u_e.col(1) = hfunc2.col(m - 1);
                } else {
                    u_e.col(1) = hfunc1.col(m - 1);
                }

                if ((var_types[0] == "d") | (var_types[1] == "d")) {
                    u_e.conservativeResize(b.size, 4);
                    u_e.col(2) = hfunc2_sub.col(edge);
                    if (m == vine_struct_.struct_array(tree, edge, true)) {
                        u_e.col(3) = hfunc2_sub.col(m - 1);
                    } else {
                        u_e.col(3) = hfunc1_sub.col(m - 1);
                    }
                }

                if (vine_struct_.needed_hfunc1(tree, edge)) {
                    hfunc1.col(edge) = edge_copula.hfunc1(u_e);
                    if (var_types[1] == "d") {
                        u_e_sub = u_e;
                        u_e_sub.col(1) = u_e.col(3);
                        hfunc1_sub.col(edge) = edge_copula.hfunc1(u_e_sub);
                    }
                }
                if (vine_struct_.needed_hfunc2(tree, edge)) {
                    hfunc2.col(edge) = edge_copula.hfunc2(u_e);
                    if (var_types[0] == "d") {
                        u_e_sub = u_e;
                        u_e_sub.col(0) = u_e.col(2);
                        hfunc2_sub.col(edge) = edge_copula.hfunc2(u_e_sub);
                    }
                }
            }
        }
        p.segment(b.begin, b.size) = hfunc2.col(0);
    };

    if (trunc_lvl > 0) {
        tools_thread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
        pool.map(do_batch, tools_batch::create_batches(u.rows(), num_threads));
        pool.join();
    }

    return p;
}
