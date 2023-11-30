/**
 * @file inset.h
 * @brief Implemente the Insetting Formation (IF) algorithm
 */

#ifndef _INSET_H
#define _INSET_H

#include <exception>
#include <cstdio>
#include <mosek.h>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <vector>
#include <casadi/casadi.hpp>
#include <inset_basis/data_utils.h>

enum class Insetting_Formation_Algorithm
{
    NO_ALGORITHM_USED = 0,
    ELLIPSE_BASED_INSETTING_FORMATION,
    MVB_BASED_INSETTING_FORMATION,
    HEIGHT_FREE_CONSTRAINED_INSETTING_FORMATION
};

namespace inset {
/// formation class
class formation
{
private:
    Eigen::VectorXd quadrotor_1_;
    Eigen::VectorXd quadrotor_2_;
    Eigen::VectorXd payload_;
    Eigen::VectorXd quadrotor_1_bound_;
    Eigen::VectorXd quadrotor_2_bound_;
    Eigen::VectorXd payload_bound_;
    bool change_lable_;
    bool change_bound_lable_;
    Insetting_Formation_Algorithm insetting_formation_type_;

public:
    formation(int dim = 2);
    formation(Eigen::VectorXd quadrotor_1, Eigen::VectorXd quadrotor_2, Eigen::VectorXd payload);
    ~formation() {}
    /// get quadrotor_1 position
    const Eigen::VectorXd &get_quadrotor_1() const;
    /// get quadrotor_2 position
    const Eigen::VectorXd &get_quadrotor_2() const;
    /// get payload position
    const Eigen::VectorXd &get_payload() const;
    /// get quadrotor_1 Virtual Boundary Point(VBP)
    const Eigen::VectorXd &get_quadrotor_1_bound() const;
    /// get quadrotor_2 Virtual Boundary Point(VBP)
    const Eigen::VectorXd &get_quadrotor_2_bound() const;
    /// get payload Virtual Boundary Point(VBP)
    const Eigen::VectorXd &get_payload_bound() const;
    /// get Inset success lable, true represent success, false represent fail
    const bool &get_change_lable() const;
    /// get construct Virtual Boundary Point(VBP) success lable, true represent success, false represent fail
    const bool &get_change_bound_lable() const;
    /// get insetting formation(IF) algorithm type
    const Insetting_Formation_Algorithm &get_insetting_formation_type() const;

    /// set quadrotor_1 position
    void set_quadrotor_1(const Eigen::VectorXd &quadrotor_1);
    /// set quadrotor_2 position
    void set_quadrotor_2(const Eigen::VectorXd &quadrotor_2);
    /// set payload position
    void set_payload(const Eigen::VectorXd &payload);
    /// set quadrotor_1 Virtual Boundary Point(VBP)
    void set_quadrotor_1_bound(const Eigen::VectorXd &quadrotor_1_bound);
    /// set quadrotor_2 Virtual Boundary Point(VBP)
    void set_quadrotor_2_bound(const Eigen::VectorXd &quadrotor_2_bound);
    /// set payload Virtual Boundary Point(VBP)
    void set_payload_bound(const Eigen::VectorXd &payload_bound);
    /// set Inset success lable, true represent success, false represent fail
    void set_change_lable(const bool &change_lable);
    /// set construct Virtual Boundary Point(VBP) success lable, true represent success, false represent fail
    void set_change_bound_lable(const bool &change_bound_lable);
    /// set Insetting Formation(IF) algorithm type
    void set_insetting_formation_type(const Insetting_Formation_Algorithm &insetting_formation_type);

    /// Maximum Virtual Boundary(MVB)-based Insetting Formation(IF)
    void construct_formation_high_equ(const double length, const double low_rate);
    /// Maximum Virtual Boundary(MVB)-based Insetting Formation(IF)
    void construct_formation_high_equ(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, const Eigen::Vector4d &bound,
                                      const Eigen::Vector2d &limit_rate, const double length);
    /// Ellipse-based IF schematic
    void construct_formation_ellipse(const Eigen::MatrixXd &C, const Eigen::VectorXd &D, const Eigen::Vector4d &bound,
                                     const Eigen::Vector2d &limit_rate, const double length);
    /// construct Maximum Virtual Boundary(MVB)
    void construct_formation_bound(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, const Eigen::Vector4d &bound,
                                   const Eigen::Vector2d &limit_rate, const double length,
                                   MSKenv_t *existing_env = NULL);
    /// Height-free constrained-based Insetting Formation(IF)
    void construct_formation_high_diff(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, const Eigen::Vector4d &bound,
                                       const Eigen::Vector2d &limit_rate, const double length,
                                       const Eigen::Vector2d &initial_point);
    /// Height-free constrained-based Insetting Formation(IF), Consider dynamics
    void construct_formation_high_diff(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, const Eigen::Vector4d &bound,
                                       const double low_rate, const double length, const Eigen::Vector2d &initial_point,
                                       const double T_max, const double G_payload);
    /// Insetting Formation(IF)
    void construct_formation(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, const Eigen::MatrixXd &C,
                             const Eigen::VectorXd &D, const Eigen::Vector4d &bound, const Eigen::Vector2d &limit_rate,
                             const double length, const Eigen::Vector2d &initial_point);
    /// Insetting Formation(IF), Consider dynamics
    void construct_formation(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, const Eigen::MatrixXd &C,
                             const Eigen::VectorXd &D, const Eigen::Vector4d &bound, const Eigen::Vector2d &limit_rate,
                             const double length, const Eigen::Vector2d &initial_point, const double T_max,
                             const double G_payload);
};

///  based Mosek INSET Error class
class INSETMosekError : public std::exception
{
private:
    std::string message;

public:
    explicit INSETMosekError(MSKrescodee res)
    {
        /* In case of an error print error code and description. */
        char symname[MSK_MAX_STR_LEN];
        char desc[MSK_MAX_STR_LEN];
        MSK_getcodedesc(res, symname, desc);
        message = std::string(symname) + ": " + std::string(desc);
    }

    const char *what() const throw() { return message.c_str(); }

    ~INSETMosekError() throw() {}
};
/// Inset Infeasible Error class
class InsetInfeasibleError : public std::exception
{
    const char *what() const throw() { return "Inset formation problem is infeasible"; }
};

/// Determine if the solution is correctly completed
void check_res(MSKrescodee res);
/// Setting the geometric boundaries of the quadrotor and the payload
Eigen::Vector4d set_bound(const double quadrotor_geometric_bound_x, const double quadrotor_geometric_bound_y,
                          const double payload_geometric_bound_x, const double payload_geometric_bound_y);
/// Set the distance limit between UAVs, T_max and G_payload in N
Eigen::Vector2d set_limit_dis_rate(const double limit_low_dis_rate, const double T_max, const double G_payload);
} // namespace inset

#endif
