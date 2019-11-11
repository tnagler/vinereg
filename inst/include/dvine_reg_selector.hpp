#pragma once
#include <RcppEigen.h>
#include <algorithm>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/vinecop/class.hpp>

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

  void extend_fit(size_t var);

  double get_crit() { return fit_.crit; }
  std::vector<size_t> get_remaining_vars() { return fit_.remaining_vars; }
  std::vector<size_t> get_selected_vars() { return fit_.selected_vars; }
  Bicop get_pc(size_t t) { return fit_.pcs[t]; }

private:
  void initialize_var(size_t var);
  std::vector<std::string> get_edge_types(size_t t);
  Eigen::MatrixXd get_edge_data(size_t t);
  void update_hfunc1(size_t t, const Eigen::MatrixXd& u_e);
  void update_hfunc2(size_t t, const Eigen::MatrixXd& u_e);
  void update_selcrit();
  void update_vars(size_t var);

  size_t p_;
  Eigen::MatrixXd data_;
  std::vector<std::string> var_types_;
  FitControlsBicop controls_;

  DVineFitTemporaries fit_;
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
  fit_.hfunc2[0] = data_.col(0);
  fit_.crit = 0.0;
  if (var_types_[0] == "d") {
    fit_.hfunc2_sub[0] = data_.col(p_ + 1);
  }
}

inline void DVineRegSelector::extend_fit(size_t var)
{
  this->initialize_var(var);

  for (size_t t = 0; t < fit_.selected_vars.size() + 1; t++) {
    auto u_e = this->get_edge_data(t);

    fit_.pcs[t].set_var_types(this->get_edge_types(t));
    fit_.pcs[t].select(u_e, controls_);

    this->update_hfunc1(t, u_e);
    this->update_hfunc2(t, u_e);
  }

  this->update_selcrit();
  this->update_vars(var);
}

inline void DVineRegSelector::initialize_var(size_t var)
{
  fit_.hfunc1[0] = data_.col(var);
  if (var_types_[var] == "d") {
    fit_.hfunc1_sub[0] = data_.col(p_ + 1 + var);
  }
}

// obtain variable types for the new edge in tree t
inline std::vector<std::string>
  DVineRegSelector::get_edge_types(size_t t)
  {
    // the variable type can be inferred from the existence of _sub data
    std::vector<std::string> var_types{ "c", "c" };
    if (fit_.hfunc2_sub[t].size())
      var_types[0] = "d";
    if (fit_.hfunc1_sub[t].size())
      var_types[1] = "d";

    return var_types;
  }

// obtain data for the new edge in tree t
inline Eigen::MatrixXd DVineRegSelector::get_edge_data(size_t t)
{
  // (hfunc2 has been computed in previous fit, hfunc1 in previous tree)
  Eigen::MatrixXd u_e(data_.rows(), 2);
  u_e.col(0) = fit_.hfunc2[t];
  u_e.col(1) = fit_.hfunc1[t];

  if (fit_.hfunc2_sub[t].size() | fit_.hfunc1_sub[t].size()) {
    u_e.conservativeResize(u_e.rows(), 4);
    // use dummys for _sub data if variable is not discrete
    u_e.col(2) =
      fit_.hfunc2_sub[t].size() ? fit_.hfunc2_sub[t] : fit_.hfunc2[t];
    u_e.col(3) =
      fit_.hfunc1_sub[t].size() ? fit_.hfunc1_sub[t] : fit_.hfunc1[t];
  }

  return u_e;
}

inline void DVineRegSelector::update_hfunc1(size_t t, const Eigen::MatrixXd& u_e)
{
  if (p_ == t + 1) // selection is complete
    return;
  fit_.hfunc1[t + 1] = fit_.pcs[t].hfunc1(u_e);
  if (u_e(0, 1) != u_e(0, 3)) { // second variable is discrete
    auto u_e_sub = u_e;
    u_e_sub.col(1) = u_e.col(3);
    fit_.hfunc1_sub[t + 1] = fit_.pcs[t].hfunc1(u_e_sub);
  } else {
    fit_.hfunc1_sub[t + 1] = Eigen::VectorXd();
  }
}

inline void DVineRegSelector::update_hfunc2(size_t t, const Eigen::MatrixXd& u_e)
{
  // use unneeded space to store hfunc2 (required when we add the next var);
  // will be shifted one up at the end
  fit_.hfunc2[t] = fit_.pcs[t].hfunc2(u_e);
  if (u_e.cols() > 2) {
    if (u_e(0, 0) != u_e(0, 2)) { // first variable is discrete
      auto u_e_sub = u_e;
      u_e_sub.col(0) = u_e.col(2);
      fit_.hfunc2_sub[t] = fit_.pcs[t].hfunc2(u_e_sub);
    } else {
      fit_.hfunc2_sub[t] = Eigen::VectorXd();
    }
  }

  if (t == fit_.selected_vars.size()) { // all trees have been fit
    // shift hfunc2 entries into correct tree level
    std::rotate(fit_.hfunc2.begin(), fit_.hfunc2.end() - 1, fit_.hfunc2.end());

    // fill first tree with actual observations
    fit_.hfunc2[0] = fit_.hfunc1[0];
    if (fit_.hfunc1_sub[0].size() > 0) {
      fit_.hfunc2_sub[0] = fit_.hfunc1_sub[0];
    }
  }
}

// update value of the criterion used for variable and family selection
inline void DVineRegSelector::update_selcrit()
{
  if (controls_.get_selection_criterion() == "loglik")
    fit_.crit += fit_.pcs[fit_.selected_vars.size()].get_loglik();
  if (controls_.get_selection_criterion() == "aic")
    fit_.crit -= fit_.pcs[fit_.selected_vars.size()].get_aic();
  if (controls_.get_selection_criterion() == "bic")
    fit_.crit -= fit_.pcs[fit_.selected_vars.size()].get_bic();
}

// remove var from remaining variables; add to selected variables
inline void DVineRegSelector::update_vars(size_t var)
{
  fit_.remaining_vars.erase(
    std::remove(fit_.remaining_vars.begin(), fit_.remaining_vars.end(), var));
  fit_.selected_vars.push_back(var);
}

}
