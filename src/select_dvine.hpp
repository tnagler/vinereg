#pragma once
#include <algorithm>
#include <vinecopulib/misc/tools_stl.hpp>
#include <vinecopulib/vinecop/class.hpp>
#include <RcppEigen.h>

namespace vinereg {
using namespace vinecopulib;
class DVineSelectStatus
{

public:
    DVineSelectStatus(const Eigen::MatrixXd& data,
                      const std::vector<std::string>& var_types,
                      const FitControlsBicop& controls)
        : p_(var_types.size() - 1)
        , data_(data)
        , var_types_(var_types)
        , controls_(controls)
    {
        hfunc1_.resize(p_);
        hfunc2_.resize(p_);
        hfunc1_sub_.resize(p_);
        hfunc2_sub_.resize(p_);
        pcs_.resize(p_);
        remaining_vars_ = tools_stl::seq_int(1, p_);
        selected_vars_.reserve(p_);

        hfunc2_[0] = data_.col(0);
        if (var_types_[0] == "d") {
            hfunc2_sub_[0] = data_.col(p_ + 1);
        }
    }

    void extend_fit(size_t var)
    {
        Eigen::MatrixXd u_e(data_.rows(), 2);
        hfunc1_[0] = data_.col(var);
        if (var_types_[var] == "d") {
            hfunc1_sub_[0] = data_.col(p_ + 1 + var);
        }

        for (size_t t = 0; t < selected_vars_.size() + 1; t++) {
            std::vector<std::string> var_types{ "c", "c" };
            if (hfunc2_sub_[t].size())
                var_types[0] = "d";
            if (hfunc1_sub_[t].size())
                var_types[1] = "d";

            u_e.col(0) = hfunc2_[t];
            u_e.col(1) = hfunc1_[t];
            if ((var_types[0] == "d") | (var_types[1] == "d")) {
                u_e.conservativeResize(u_e.rows(), 4);
                u_e.col(2) = (var_types[0] == "d") ? hfunc2_sub_[t] : u_e.col(0);
                u_e.col(3) = (var_types[1] == "d") ? hfunc1_sub_[t] : u_e.col(1);
            }
            pcs_[t].set_var_types(var_types);
            pcs_[t].select(u_e, controls_);

            if (p_ == t + 1)
                break;

            hfunc1_[t + 1] = pcs_[t].hfunc1(u_e);
            if (var_types[1] == "d") {
                auto u_e_sub = u_e;
                u_e_sub.col(1) = u_e.col(3);
                hfunc1_sub_[t + 1] = pcs_[t].hfunc1(u_e_sub);
            } else {
                hfunc1_sub_[t + 1] = Eigen::VectorXd();
            }
            // use unneeded space to store hfunc1 for addition of next variable;
            // will be shifted one up after the fit.
            hfunc2_[t] = pcs_[t].hfunc2(u_e);
            if (var_types[0] == "d") {
                auto u_e_sub = u_e;
                u_e_sub.col(0) = u_e.col(2);
                hfunc2_sub_[t] = pcs_[t].hfunc2(u_e_sub);
            } else {
                hfunc2_sub_[t] = Eigen::VectorXd();
            }
        }

        // shift hfunc2 entries and fill first tree with observations
        std::rotate(hfunc2_.begin(), hfunc2_.end() - 1, hfunc2_.end());
        hfunc2_[0] = data_.col(var);
        if (var_types_[var] == "d") {
            hfunc2_sub_[0] = data_.col(p_ + 1 + var);
        }

        if (controls_.get_selection_criterion() == "loglik")
            crit_ += pcs_[selected_vars_.size()].get_loglik();
        if (controls_.get_selection_criterion() == "aic")
            crit_ -= pcs_[selected_vars_.size()].get_aic();
        if (controls_.get_selection_criterion() == "bic")
            crit_ -= pcs_[selected_vars_.size()].get_bic();

        remaining_vars_.erase(
            std::remove(remaining_vars_.begin(), remaining_vars_.end(), var));
        selected_vars_.push_back(var);
    }

    double get_crit() { return crit_; }

    std::vector<size_t> get_remaining_vars() { return remaining_vars_; }

    std::vector<size_t> get_selected_vars() { return selected_vars_; }

    Bicop get_pc(size_t t) { return pcs_[t]; }

private:
    size_t p_;
    std::vector<Eigen::VectorXd> hfunc1_;
    std::vector<Eigen::VectorXd> hfunc2_;
    std::vector<Eigen::VectorXd> hfunc1_sub_;
    std::vector<Eigen::VectorXd> hfunc2_sub_;
    std::vector<Bicop> pcs_;
    std::vector<size_t> remaining_vars_;
    std::vector<size_t> selected_vars_;
    Eigen::MatrixXd data_;
    std::vector<std::string> var_types_;
    FitControlsBicop controls_;
    double crit_{ 0.0 };
};

}
