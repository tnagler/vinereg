#pragma once
#include <RcppThread.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/bicop/class.hpp>

namespace vinecopulib {
  class Bicop;

  namespace tools_stl {
    template<typename T>
    std::vector<T> rank(const std::vector<T>& x)
    {
      std::vector<size_t> r(x.size());
      auto order = tools_stl::get_order(x);
      for (auto i : order)
        r[order[i]] = static_cast<double>(i + 1);
      return r;
    }
  }
}

namespace vinereg {
using namespace vinecopulib;

struct DVineFitTemporaries
{
  std::vector<Eigen::VectorXd> hfunc1;
  std::vector<Eigen::VectorXd> hfunc2;
  std::vector<Eigen::VectorXd> hfunc1_sub;
  std::vector<Eigen::VectorXd> hfunc2_sub;
  std::vector<Bicop> pcs;
  std::vector<size_t> remaining_vars;
  std::vector<size_t> selected_vars;
  double crit;
};

class DVineRegSelector
{
public:
  DVineRegSelector(const Eigen::MatrixXd& data,
                   const std::vector<std::string>& var_types,
                   const FitControlsBicop& controls);

  void select_model();
  std::vector<size_t> get_selected_vars() const { return fit_.selected_vars; }
  std::vector<std::vector<Bicop>> get_pcs() const { return pcs_; }

private:
  void extend_fit(DVineFitTemporaries& fit, size_t var) const;
  void initialize_var(DVineFitTemporaries& fit, size_t var) const;
  std::vector<std::string> get_edge_types(const DVineFitTemporaries& fit,
                                          size_t t) const;
  Eigen::MatrixXd get_edge_data(const DVineFitTemporaries& fit, size_t t) const;
  void fit_pair_copula(DVineFitTemporaries& fit,
                       size_t t,
                       const Eigen::MatrixXd& u_e) const;
  void update_hfunc1(DVineFitTemporaries& fit,
                     size_t t,
                     const Eigen::MatrixXd& u_e) const;
  void update_hfunc2(DVineFitTemporaries& fit,
                     size_t t,
                     const Eigen::MatrixXd& u_e) const;
  void update_hfuncs(DVineFitTemporaries& fit,
                     size_t t,
                     const Eigen::MatrixXd& u_e) const;
  void update_selcrit(DVineFitTemporaries& fit) const;
  void update_vars(DVineFitTemporaries& fit, size_t var) const;
  void update_status(DVineFitTemporaries& fit, size_t var) const;

  size_t p_;
  Eigen::MatrixXd data_;
  std::vector<std::string> var_types_;
  FitControlsBicop controls_;
  DVineFitTemporaries fit_;
  std::vector<std::vector<Bicop>> pcs_;
};

inline DVineRegSelector::DVineRegSelector(
    const Eigen::MatrixXd& data,
    const std::vector<std::string>& var_types,
    const FitControlsBicop& controls)
  : p_(var_types.size() - 1)
  , data_(data)
  , var_types_(var_types)
  , controls_(controls)
{
  fit_.hfunc1.resize(p_);
  fit_.hfunc2.resize(p_);
  fit_.hfunc1_sub.resize(p_);
  fit_.hfunc2_sub.resize(p_);
  fit_.pcs.resize(p_);
  fit_.remaining_vars = tools_stl::seq_int(1, p_);
  fit_.selected_vars.reserve(p_);
  fit_.crit = 0.0;

  fit_.hfunc2[0] = data_.col(0);
  if (var_types_[0] == "d") {
    fit_.hfunc2_sub[0] = data_.col(p_ + 1);
  }
}

inline void DVineRegSelector::select_model()
{
  std::vector<std::vector<Bicop>> pcs;
  std::mutex m; // required to synchronize write/reads to the selector
  auto num_threads = controls_.get_num_threads();
  RcppThread::ThreadPool pool(num_threads > 1 ? num_threads : 0);

  while (fit_.selected_vars.size() < p_) {
    auto old_fit = fit_;  // fix current model (fit_ will be modified below)
    auto fit_replace_if_better = [&](size_t var) {
      DVineFitTemporaries new_fit = old_fit;
      this->extend_fit(new_fit, var);
      std::lock_guard<std::mutex> lk(m);  // synchronize
      if (new_fit.crit > fit_.crit)
        fit_ = std::move(new_fit);
    };
    pool.map(fit_replace_if_better, old_fit.remaining_vars);
    pool.wait();

    if (fit_.selected_vars == old_fit.selected_vars)
      break;  // could not improve the selection criterion

    // model improved; store pair copulas for new variable
    auto p_sel = fit_.selected_vars.size();
    pcs_.push_back(std::vector<Bicop>{ fit_.pcs[p_sel - 1] });
    for (size_t t = 0; t < p_sel - 1; t++)
      pcs_[t].push_back(fit_.pcs[t]);
  }
  pool.join();
}

inline void DVineRegSelector::extend_fit(DVineFitTemporaries& fit,
                                         size_t var) const
{
  this->initialize_var(fit, var);

  for (size_t t = 0; t < fit.selected_vars.size() + 1; t++) {
    auto u_e = this->get_edge_data(fit, t);
    this->fit_pair_copula(fit, t, u_e);
    this->update_hfuncs(fit, t, u_e);
  }
  this->update_status(fit, var);
}

inline void DVineRegSelector::initialize_var(DVineFitTemporaries& fit,
                                             size_t var) const
{
  fit.hfunc1[0] = data_.col(var);
  fit.hfunc1_sub[0] =
    (var_types_[var] == "d") ? data_.col(p_ + 1 + var) : Eigen::VectorXd();
}

// obtain variable types for the new edge in tree t
inline std::vector<std::string> DVineRegSelector::get_edge_types(
    const DVineFitTemporaries& fit, size_t t) const
{
  // the variable type can be inferred from the existence of _sub data
  std::vector<std::string> var_types(2);
  var_types[0] = fit.hfunc2_sub[t].size() ? "d" : "c";
  var_types[1] = fit.hfunc1_sub[t].size() ? "d" : "c";
  return var_types;
}

// obtain data for the new edge in tree t
inline Eigen::MatrixXd DVineRegSelector::get_edge_data(
    const DVineFitTemporaries& fit, size_t t) const
{
  // (hfunc2 has been computed in previous fit, hfunc1 in previous tree)
  Eigen::MatrixXd u_e(data_.rows(), 2);
  u_e.col(0) = fit.hfunc2[t];
  u_e.col(1) = fit.hfunc1[t];

  if (fit.hfunc2_sub[t].size() | fit.hfunc1_sub[t].size()) {
    u_e.conservativeResize(u_e.rows(), 4);
    // use dummys for _sub data if variable is not discrete
    u_e.col(2) =
      fit.hfunc2_sub[t].size() ? fit.hfunc2_sub[t] : fit.hfunc2[t];
    u_e.col(3) =
      fit.hfunc1_sub[t].size() ? fit.hfunc1_sub[t] : fit.hfunc1[t];
  }

  return u_e;
}

inline void DVineRegSelector::fit_pair_copula(DVineFitTemporaries& fit,
                                              size_t t,
                                              const Eigen::MatrixXd& u_e) const
{
  auto var_types = this->get_edge_types(fit, t);
  fit.pcs[t].set_var_types(var_types);
  fit.pcs[t].select(u_e, controls_);
}

inline void DVineRegSelector::update_hfunc1(DVineFitTemporaries& fit,
                                            size_t t,
                                            const Eigen::MatrixXd& u_e) const
{
  if (p_ == t + 1) // selection is complete
    return;
  fit.hfunc1[t + 1] = fit.pcs[t].hfunc1(u_e);
  if (fit.hfunc1_sub[t].size()) { // second variable is discrete
    auto u_e_sub = u_e;
    u_e_sub.col(1) = u_e.col(3);
    fit.hfunc1_sub[t + 1] = fit.pcs[t].hfunc1(u_e_sub);
  } else {
    fit.hfunc1_sub[t + 1] = Eigen::VectorXd();
  }
}

inline void DVineRegSelector::update_hfunc2(DVineFitTemporaries& fit,
                                            size_t t,
                                            const Eigen::MatrixXd& u_e) const
{
  // use unneeded space to store hfunc2 (required when we add the next var);
  // will be shifted one up at the end
  fit.hfunc2[t] = fit.pcs[t].hfunc2(u_e);
  if (u_e.cols() > 2) {
    if (fit.hfunc2_sub[t].size()) { // first variable is discrete
      auto u_e_sub = u_e;
      u_e_sub.col(0) = u_e.col(2);
      fit.hfunc2_sub[t] = fit.pcs[t].hfunc2(u_e_sub);
    } else {
      fit.hfunc2_sub[t] = Eigen::VectorXd();
    }
  }

  if (t == fit.selected_vars.size()) { // all trees have been fit
    // shift hfunc2 entries into correct tree level
    std::rotate(
      fit.hfunc2.begin(), fit.hfunc2.end() - 1, fit.hfunc2.end());
    std::rotate(
      fit.hfunc2_sub.begin(), fit.hfunc2_sub.end() - 1, fit.hfunc2_sub.end());

    // fill first tree with actual observations
    fit.hfunc2[0] = fit.hfunc1[0];
    if (fit.hfunc1_sub[0].size()) {
      fit.hfunc2_sub[0] = fit.hfunc1_sub[0];
    } else {
      fit.hfunc2_sub[0] = Eigen::VectorXd();
    }
  }
}

void DVineRegSelector::update_hfuncs(DVineFitTemporaries& fit,
                                     size_t t,
                                     const Eigen::MatrixXd& u_e) const
{
  this->update_hfunc1(fit, t, u_e);
  this->update_hfunc2(fit, t, u_e);
}


// update value of the criterion used for variable and family selection
inline void DVineRegSelector::update_selcrit(DVineFitTemporaries& fit) const
{
  if (controls_.get_selection_criterion() == "loglik")
    fit.crit += fit.pcs[fit.selected_vars.size()].get_loglik();
  if (controls_.get_selection_criterion() == "aic")
    fit.crit -= fit.pcs[fit.selected_vars.size()].get_aic();
  if (controls_.get_selection_criterion() == "bic")
    fit.crit -= fit.pcs[fit.selected_vars.size()].get_bic();
}

// remove var from remaining variables; add to selected variables
inline void DVineRegSelector::update_vars(DVineFitTemporaries& fit, size_t var)
  const
{
  fit.remaining_vars.erase(
    std::remove(fit.remaining_vars.begin(), fit.remaining_vars.end(), var));
  fit.selected_vars.push_back(var);
}

inline void DVineRegSelector::update_status(DVineFitTemporaries& fit,size_t var)
  const
{
  this->update_selcrit(fit);
  this->update_vars(fit, var);
}

}
