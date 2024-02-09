#ifndef BOOST_NO_AUTO_PTR
#define BOOST_NO_AUTO_PTR
#endif

#ifndef BOOST_MATH_PROMOTE_DOUBLE_POLICY
#define BOOST_MATH_PROMOTE_DOUBLE_POLICY false
#else
#undef BOOST_MATH_PROMOTE_DOUBLE_POLICY
#define BOOST_MATH_PROMOTE_DOUBLE_POLICY false
#endif

#ifndef BOOST_ALLOW_DEPRECATED_HEADERS
#define BOOST_ALLOW_DEPRECATED_HEADERS
#endif

#include "dvine_reg_selector.hpp"
#include <RcppThread.h>
#include <kde1d-wrappers.hpp>
#include <vinecopulib-wrappers.hpp>
#include <vinecopulib/bicop/fit_controls.hpp>

using namespace vinecopulib;

// [[Rcpp::export]]
std::vector<Rcpp::List>
fit_margins_cpp(const Eigen::MatrixXd& data,
                const Eigen::VectorXi& nlevels,
                const Eigen::VectorXd& mult,
                const Eigen::VectorXd& xmin,
                const Eigen::VectorXd& xmax,
                const Eigen::VectorXd& bw,
                const Eigen::VectorXi& deg,
                const Eigen::VectorXd& weights,
                size_t num_threads)
{
    size_t d = data.cols();
    std::vector<kde1d::Kde1d> fits_cpp(d);
    num_threads = (num_threads > 1) ? num_threads : 0;
    RcppThread::parallelFor(
      0,
      d,
      [&](const size_t& k) {
          fits_cpp[k] = kde1d::Kde1d(data.col(k),
                                     nlevels(k),
                                     bw(k),
                                     mult(k),
                                     xmin(k),
                                     xmax(k),
                                     deg(k),
                                     weights);
      },
      num_threads);

    // we can't do the following in parallel because it calls R API
    std::vector<Rcpp::List> fits_r(d);
    for (size_t k = 0; k < d; ++k) {
        fits_r[k] = kde1d::kde1d_wrap(fits_cpp[k]);
    }
    return fits_r;
}

// [[Rcpp::export]]
Rcpp::List
select_dvine_cpp(const Eigen::MatrixXd& data,
                 std::vector<std::string> family_set,
                 std::string par_method,
                 std::string nonpar_method,
                 double mult,
                 std::string selcrit,
                 const Eigen::VectorXd& weights,
                 double psi0,
                 bool preselect_families,
                 size_t cores,
                 const std::vector<std::string>& var_types)
{
    // set up the cpp fit controls from all the arguments ------
    std::vector<BicopFamily> fam_set(family_set.size());
    for (unsigned int fam = 0; fam < fam_set.size(); ++fam) {
        fam_set[fam] = to_cpp_family(family_set[fam]);
    }
    FitControlsBicop controls(fam_set,
                              par_method,
                              nonpar_method,
                              mult,
                              selcrit,
                              weights,
                              psi0,
                              preselect_families,
                              cores);

    // select the model -----------------------------------------
    vinereg::DVineRegSelector selector(data, var_types, controls);
    selector.select_model();
    auto selected_vars = selector.get_selected_vars();
    auto pcs = selector.get_pcs();

    // make results R-compatible -------------------------------
    Rcpp::List vinecop_r;
    // rank ensures that vars are 1, ..., p_sel
    auto order = tools_stl::cat(static_cast<size_t>(0), selected_vars);
    auto new_struct = DVineStructure(tools_stl::rank(order));

    auto sv = selected_vars;
    std::sort(sv.begin(), sv.end());
    auto vt = std::vector<std::string>();
    vt.push_back(var_types[0]);
    for (auto v : sv)
        vt.push_back(var_types[v]);

    vinecop_r = Rcpp::List::create(
      Rcpp::Named("pair_copulas") = pair_copulas_wrap(pcs, order.size(), true),
      Rcpp::Named("structure") = rvine_structure_wrap(new_struct),
      Rcpp::Named("var_types") = vt,
      Rcpp::Named("npars") = Vinecop(new_struct, pcs).get_npars(),
      Rcpp::Named("nobs") = data.rows(),
      Rcpp::Named("loglik") = NAN,
      Rcpp::Named("threshold") = 0);
    vinecop_r.attr("class") =
      Rcpp::CharacterVector{ "vinecop", "vinecop_dist" };
    for (auto& v : selected_vars) // R indexing starts at 1
        v++;
    return Rcpp::List::create(Rcpp::Named("vine") = vinecop_r,
                              Rcpp::Named("selected_vars") = selected_vars);
}

// [[Rcpp::export]]
std::vector<Eigen::VectorXd>
cond_quantile_cpp(const Eigen::VectorXd& alpha,
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
                hfunc2_sub(0, j) =
                  u.col(d + order[j] - 1).segment(b.begin, b.size);
                hfunc1_sub(0, j) =
                  u.col(d + order[j] - 1).segment(b.begin, b.size);
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
                if ((var_types[0] == "d") || (var_types[1] == "d")) {
                    u_e.conservativeResize(b.size, 4);
                    u_e.col(2) = (var_types[0] == "d") ? hfunc2_sub(tree, edge)
                                                       : hfunc2(tree, edge);
                    u_e.col(3) = (var_types[1] == "d") ? hfunc1_sub(tree, m - 1)
                                                       : hfunc1(tree, m - 1);
                }

                if (vine_struct_.needed_hfunc1(tree, edge)) {
                    hfunc1(tree + 1, edge) = edge_copula.hfunc1(u_e);
                    if (var_types[1] == "d") {
                        u_e_sub = u_e;
                        u_e_sub.col(1) = u_e.col(3);
                        hfunc1_sub(tree + 1, edge) =
                          edge_copula.hfunc1(u_e_sub);
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
                tools_interface::check_user_interrupt(d * u.rows() > 1e3);
                Bicop edge_copula = vinecop_cpp.get_pair_copula(tree, 0);
                Eigen::MatrixXd U_e(b.size, 2);
                U_e.col(0) = hfunc2(tree + 1, 0);
                U_e.col(1) = hfunc1(tree, 1);
                if (edge_copula.get_var_types()[1] == "d") {
                    U_e.conservativeResize(b.size, 4);
                    U_e.col(2) = U_e.col(0);
                    U_e.col(3) = hfunc1_sub(tree, 1);
                }
                hfunc2(tree, 0) = edge_copula.hinv2(U_e);
            }
            q[a].segment(b.begin, b.size) = hfunc2(0, 0);
        }
    };

    RcppThread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
    pool.map(do_batch, tools_batch::create_batches(u.rows(), num_threads));
    pool.join();

    return q;
}

// [[Rcpp::export]]
Eigen::VectorXd
cond_dist_cpp(const Eigen::MatrixXd& u,
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

    Eigen::VectorXd p(u.rows());
    auto do_batch = [&](const tools_batch::Batch& b) {
        Eigen::MatrixXd hfunc1, hfunc2, u_e, hfunc1_sub, hfunc2_sub, u_e_sub;
        hfunc1 = Eigen::MatrixXd::Zero(b.size, d);
        hfunc2 = Eigen::MatrixXd::Zero(b.size, d);
        hfunc1_sub = hfunc1;
        hfunc2_sub = hfunc2;

        // data have to be reordered to correspond to natural order
        for (size_t j = 0; j < d; ++j) {
            hfunc2.col(j) = u.block(b.begin, order[j] - 1, b.size, 1);
            if (var_types_[order[j] - 1] == "d") {
                hfunc2_sub.col(j) =
                  u.block(b.begin, d + order[j] - 1, b.size, 1);
            }
        }

        for (size_t tree = 0; tree < trunc_lvl; ++tree) {
            for (size_t edge = 0; edge < d - tree - 1; ++edge) {
                tools_interface::check_user_interrupt();
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

                if ((var_types[0] == "d") || (var_types[1] == "d")) {
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

                hfunc2.col(edge) = edge_copula.hfunc2(u_e);
                if (var_types[0] == "d" &&
                    vine_struct_.needed_hfunc2(tree, edge)) {
                    u_e_sub = u_e;
                    u_e_sub.col(0) = u_e.col(2);
                    hfunc2_sub.col(edge) = edge_copula.hfunc2(u_e_sub);
                }
            }
        }
        p.segment(b.begin, b.size) = hfunc2.col(0);
    };

    RcppThread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
    pool.map(do_batch, tools_batch::create_batches(u.rows(), num_threads));
    pool.join();

    return p;
}

// [[Rcpp::export]]
Eigen::VectorXd
cond_dens_cpp(const Eigen::MatrixXd& u,
              const Rcpp::List& vinecop_r,
              size_t num_threads)
{
    tools_eigen::check_if_in_unit_cube(u);
    auto vinecop_cpp = vinecop_wrap(vinecop_r);
    auto rvine_structure_ = vinecop_cpp.get_rvine_structure();
    auto d_ = rvine_structure_.get_dim();
    auto var_types_ = vinecop_cpp.get_var_types();
    auto pair_copulas_ = vinecop_cpp.get_all_pair_copulas();
    if ((static_cast<size_t>(u.cols()) != d_) &&
        (static_cast<size_t>(u.cols()) != 2 * d_))
        throw std::runtime_error("data dimension is incompatible with model.");

    int n_discrete = 0;
    for (auto t : var_types_) {
      n_discrete += (t == "d");
    }

    // info about the vine structure (reverse rows (!) for more natural
    // indexing)
    size_t trunc_lvl = rvine_structure_.get_trunc_lvl();
    auto order = rvine_structure_.get_order();
    auto disc_cols = tools_select::get_disc_cols(var_types_);

    // initial value must be 1.0 for multiplication
    Eigen::VectorXd pdf = Eigen::VectorXd::Constant(u.rows(), 1.0);

    auto do_batch = [&](const tools_batch::Batch& b) {
        // temporary storage objects (all data must be in (0, 1))
        Eigen::MatrixXd hfunc1, hfunc2, u_e, hfunc1_sub, hfunc2_sub, u_e_sub;
        hfunc1 = Eigen::MatrixXd::Zero(b.size, d_);
        hfunc2 = Eigen::MatrixXd::Zero(b.size, d_);
        if (n_discrete > 0) {
            hfunc1_sub = hfunc1;
            hfunc2_sub = hfunc2;
        }

        // fill first row of hfunc2 matrix with evaluation points;
        // points have to be reordered to correspond to natural order
        for (size_t j = 0; j < d_; ++j) {
            hfunc2.col(j) = u.block(b.begin, order[j] - 1, b.size, 1);
            if (var_types_[order[j] - 1] == "d") {
                hfunc2_sub.col(j) =
                  u.block(b.begin, d_ + disc_cols[order[j] - 1], b.size, 1);
            }
        }

        for (size_t tree = 0; tree < trunc_lvl; ++tree) {
            tools_interface::check_user_interrupt(
              static_cast<double>(u.rows()) * static_cast<double>(d_) > 1e5);
            for (size_t edge = 0; edge < d_ - tree - 1; ++edge) {
                tools_interface::check_user_interrupt(edge % 100 == 0);
                // extract evaluation point from hfunction matrices (have been
                // computed in previous tree level)
                Bicop* edge_copula = &pair_copulas_[tree][edge];
                auto var_types = edge_copula->get_var_types();
                size_t m = rvine_structure_.min_array(tree, edge);

                u_e = Eigen::MatrixXd(b.size, 2);
                u_e.col(0) = hfunc2.col(edge);
                if (m == rvine_structure_.struct_array(tree, edge, true)) {
                    u_e.col(1) = hfunc2.col(m - 1);
                } else {
                    u_e.col(1) = hfunc1.col(m - 1);
                }

                if ((var_types[0] == "d") || (var_types[1] == "d")) {
                    u_e.conservativeResize(b.size, 4);
                    u_e.col(2) = hfunc2_sub.col(edge);
                    if (m == rvine_structure_.struct_array(tree, edge, true)) {
                        u_e.col(3) = hfunc2_sub.col(m - 1);
                    } else {
                        u_e.col(3) = hfunc1_sub.col(m - 1);
                    }
                }

                if (edge == 0) {
                    pdf.segment(b.begin, b.size) =
                      pdf.segment(b.begin, b.size)
                        .cwiseProduct(edge_copula->pdf(u_e));
                }

                // h-functions are only evaluated if needed in next step
                if (rvine_structure_.needed_hfunc1(tree, edge)) {
                    hfunc1.col(edge) = edge_copula->hfunc1(u_e);
                    if (var_types[1] == "d") {
                        u_e_sub = u_e;
                        u_e_sub.col(1) = u_e.col(3);
                        hfunc1_sub.col(edge) = edge_copula->hfunc1(u_e_sub);
                    }
                }
                if (rvine_structure_.needed_hfunc2(tree, edge)) {
                    hfunc2.col(edge) = edge_copula->hfunc2(u_e);
                    if (var_types[0] == "d") {
                        u_e_sub = u_e;
                        u_e_sub.col(0) = u_e.col(2);
                        hfunc2_sub.col(edge) = edge_copula->hfunc2(u_e_sub);
                    }
                }
            }
        }
    };

    if (trunc_lvl > 0) {
        tools_thread::ThreadPool pool((num_threads == 1) ? 0 : num_threads);
        pool.map(do_batch, tools_batch::create_batches(u.rows(), num_threads));
        pool.join();
    }

    return pdf;
}

// [[Rcpp::export]]
double
cond_loglik_cpp(const Eigen::MatrixXd& u,
                const Rcpp::List& vinecop_r,
                size_t num_threads)
{
  return cond_dens_cpp(u, vinecop_r, num_threads).array().log().sum();
}

// [[Rcpp::export()]]
Rcpp::List
with_parameters_cop_cpp(const Rcpp::List& vinecop_r,
                        const Eigen::VectorXd parameters)
{
  auto vc = vinecop_wrap(vinecop_r);
  auto d = vc.get_dim();
  auto pcs = vc.get_all_pair_copulas();

  size_t i = 0;
  for (size_t t = 0; t < d - 1; t++) {
    for (size_t e = 0; e < std::min(d, d - t - 1); e++) {
      if (pcs[t][e].get_family() == BicopFamily::indep)
        continue;
      auto lb = pcs[t][e].get_parameters_lower_bounds();
      auto ub = pcs[t][e].get_parameters_upper_bounds();
      auto new_pars =
        parameters.segment(i, lb.size()).cwiseMax(lb).cwiseMin(ub);
      pcs[t][e].set_parameters(new_pars);
      i += lb.size();
    }
  }

  vc.set_all_pair_copulas(pcs);
  return vinecop_wrap(vc, false);
}
