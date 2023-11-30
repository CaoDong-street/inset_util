#include "inset_util/inset.h"

using namespace Eigen;
using namespace casadi;

namespace inset {

Eigen::Vector4d set_bound(const double quadrotor_geometric_bound_x, const double quadrotor_geometric_bound_y,
                          const double payload_geometric_bound_x, const double payload_geometric_bound_y)
{
    Eigen::Vector4d geometric_bound(quadrotor_geometric_bound_x, quadrotor_geometric_bound_y, payload_geometric_bound_x,
                                    payload_geometric_bound_y);
    return geometric_bound;
}

Eigen::Vector2d set_limit_dis_rate(const double limit_low_dis_rate, const double T_max, const double G_payload)
{
    double limit_up_dis_rate = 2 * sqrt(1 - pow(G_payload / (2.0 * T_max), 2));
    Eigen::Vector2d limit_dis(limit_low_dis_rate, limit_up_dis_rate);
    return limit_dis;
}

#define ADD_VAR(x)                     \
    std::vector<int> ndx_##x(num_##x); \
    for (int i = 0; i < num_##x; i++) ndx_##x[i] = nvar++;

void check_res(MSKrescodee res)
{
    if (res != MSK_RES_OK)
    {
        throw INSETMosekError(res);
    }
}

#ifndef NDEBUG
/* This function prints log output from MOSEK to the terminal. */
static void MSKAPI printstr(void *handle, MSKCONST char str[]) { printf("%s", str); } /* printstr */
#endif

formation::formation(int dim)
    : quadrotor_1_(Eigen::VectorXd(dim)),
      quadrotor_2_(Eigen::VectorXd(dim)),
      payload_(Eigen::VectorXd(dim)),
      quadrotor_1_bound_(Eigen::VectorXd(dim)),
      quadrotor_2_bound_(Eigen::VectorXd(dim)),
      payload_bound_(Eigen::VectorXd(dim)),
      change_lable_(false),
      change_bound_lable_(false),
      insetting_formation_type_(Insetting_Formation_Algorithm::NO_ALGORITHM_USED)
{
}
formation::formation(Eigen::VectorXd quadrotor_1, Eigen::VectorXd quadrotor_2, Eigen::VectorXd payload)
    : quadrotor_1_(quadrotor_1),
      quadrotor_2_(quadrotor_2),
      payload_(payload),
      quadrotor_1_bound_(Eigen::VectorXd(quadrotor_1.size())),
      quadrotor_2_bound_(Eigen::VectorXd(quadrotor_2.size())),
      payload_bound_(Eigen::VectorXd(payload_.size())),
      change_lable_(false),
      change_bound_lable_(false),
      insetting_formation_type_(Insetting_Formation_Algorithm::NO_ALGORITHM_USED)
{
}

const Eigen::VectorXd &formation::get_quadrotor_1() const { return quadrotor_1_; }
const Eigen::VectorXd &formation::get_quadrotor_2() const { return quadrotor_2_; }
const Eigen::VectorXd &formation::get_payload() const { return payload_; }
const Eigen::VectorXd &formation::get_quadrotor_1_bound() const { return quadrotor_1_bound_; }
const Eigen::VectorXd &formation::get_quadrotor_2_bound() const { return quadrotor_2_bound_; }
const Eigen::VectorXd &formation::get_payload_bound() const { return payload_bound_; }
const bool &formation::get_change_lable() const { return change_lable_; }
const bool &formation::get_change_bound_lable() const { return change_bound_lable_; }
const Insetting_Formation_Algorithm &formation::get_insetting_formation_type() const
{
    return insetting_formation_type_;
}

void formation::set_quadrotor_1(const Eigen::VectorXd &quadrotor_1) { quadrotor_1_ = quadrotor_1; }
void formation::set_quadrotor_2(const Eigen::VectorXd &quadrotor_2) { quadrotor_2_ = quadrotor_2; }
void formation::set_payload(const Eigen::VectorXd &payload) { payload_ = payload; }
void formation::set_quadrotor_1_bound(const Eigen::VectorXd &quadrotor_1_bound)
{
    quadrotor_1_bound_ = quadrotor_1_bound;
}
void formation::set_quadrotor_2_bound(const Eigen::VectorXd &quadrotor_2_bound)
{
    quadrotor_2_bound_ = quadrotor_2_bound;
}
void formation::set_payload_bound(const Eigen::VectorXd &payload_bound) { payload_bound_ = payload_bound; }
void formation::set_change_lable(const bool &change_lable) { change_lable_ = change_lable; }
void formation::set_change_bound_lable(const bool &change_bound_lable) { change_bound_lable_ = change_bound_lable; }
void formation::set_insetting_formation_type(const Insetting_Formation_Algorithm &insetting_formation_type)
{
    insetting_formation_type_ = insetting_formation_type;
}

void formation::construct_formation_high_equ(const double length, const double low_rate)
{
    if (change_bound_lable_ == false)
    {
        change_lable_ = false;
        insetting_formation_type_ = Insetting_Formation_Algorithm::NO_ALGORITHM_USED;
#ifndef NDEBUG
        std::cout << "Maximum Virtual Boundary(MVB)-based insetting formation fail : "
                     "construct Maximum Virtual Boundary(MVB) fail "
                  << std::endl;
#endif
        return;
    }
    const double xq2_xq1_dis_origin = (quadrotor_1_ - quadrotor_2_).norm();
    const double xq2_xq1_dis_bound = (quadrotor_1_bound_ - quadrotor_2_bound_).norm();
    const Eigen::VectorXd quadrotor1_2_mid = (quadrotor_1_bound_ + quadrotor_2_bound_) / 2;
    const double bound_y = (quadrotor1_2_mid - payload_bound_).norm();
    double up_dis = xq2_xq1_dis_bound;
    double low_dis;
    double length_x;
    double length_y;
    if (length <= bound_y)
    {
        low_dis = low_rate * length;
    }
    else if (low_rate * length >= 2 * sqrt(pow(length, 2) - pow(bound_y, 2)))
    {
        low_dis = low_rate * length;
        ;
    }
    else
    {
        low_dis = 2 * sqrt(pow(length, 2) - pow(bound_y, 2));
    }

    if (bound_y < sqrt(pow(length, 2) - pow(xq2_xq1_dis_bound / 2.0, 2)))
    {
        change_lable_ = false;
        insetting_formation_type_ = Insetting_Formation_Algorithm::NO_ALGORITHM_USED;
#ifndef NDEBUG
        std::cout << "Maximum Virtual Boundary(MVB)-based insetting formation fail : "
                     "The height difference between the Virtual Boundary Point(VBP) of the"
                     "payload and the Virtual Boundary Point(VBP) of the quadrotor is too small "
                  << std::endl;
#endif
    }
    else
    {
        if (xq2_xq1_dis_origin <= up_dis && xq2_xq1_dis_origin >= low_dis)
        {
            length_x = xq2_xq1_dis_origin / 2;
            length_y = sqrt(pow(length, 2) - pow(length_x, 2));
        }
        else if (xq2_xq1_dis_origin > up_dis)
        {
            length_x = up_dis / 2;
            length_y = sqrt(pow(length, 2) - pow(length_x, 2));
        }
        else
        {
            length_x = low_dis / 2;
            length_y = sqrt(pow(length, 2) - pow(length_x, 2));
        }
        quadrotor_1_(0) = quadrotor1_2_mid(0) + length_x;
        quadrotor_2_(0) = quadrotor1_2_mid(0) - length_x;
        payload_(0) = quadrotor1_2_mid(0);
        quadrotor_1_(1) = quadrotor1_2_mid(1);
        quadrotor_2_(1) = quadrotor1_2_mid(1);
        payload_(1) = quadrotor1_2_mid(1) - length_y;
        change_lable_ = true;
        insetting_formation_type_ = Insetting_Formation_Algorithm::MVB_BASED_INSETTING_FORMATION;
    }
}
void formation::construct_formation_high_equ(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
                                             const Eigen::Vector4d &bound, const Eigen::Vector2d &limit_rate,
                                             const double length)
{
    formation::construct_formation_bound(A, b, bound, limit_rate, length);
    formation::construct_formation_high_equ(length, limit_rate(0));
}

void formation::construct_formation_ellipse(const Eigen::MatrixXd &C, const Eigen::VectorXd &D,
                                            const Eigen::Vector4d &bound, const Eigen::Vector2d &limit_rate,
                                            const double length)
{
    if (payload_(0) <= quadrotor_1_(0) && payload_(0) >= quadrotor_2_(0) && payload_(1) < quadrotor_1_(1) &&
        payload_(1) < quadrotor_2_(1) && quadrotor_1_(0) - quadrotor_2_(0) >= limit_rate(0) * length &&
        quadrotor_1_(0) - quadrotor_2_(0) <= limit_rate(1) * length)
    {
        Eigen::MatrixXd center_x_numerator(3, 3);
        Eigen::MatrixXd center_x_denominator(3, 3);
        Eigen::MatrixXd center_y_numerator(3, 3);
        Eigen::MatrixXd center_y_denominator(3, 3);
        Eigen::VectorXd xq1_center(2);
        Eigen::VectorXd xq2_center(2);
        Eigen::VectorXd xp_center(2);
        center_x_numerator << 1, pow(quadrotor_1_(0), 2) + pow(quadrotor_1_(1), 2), quadrotor_1_(1), 1,
            pow(quadrotor_2_(0), 2) + pow(quadrotor_2_(1), 2), quadrotor_2_(1), 1,
            pow(payload_(0), 2) + pow(payload_(1), 2), payload_(1);
        center_x_denominator << 1, quadrotor_1_(0), quadrotor_1_(1), 1, quadrotor_2_(0), quadrotor_2_(1), 1,
            payload_(0), payload_(1);
        center_y_numerator << 1, quadrotor_1_(0), pow(quadrotor_1_(0), 2) + pow(quadrotor_1_(1), 2), 1, quadrotor_2_(0),
            pow(quadrotor_2_(0), 2) + pow(quadrotor_2_(1), 2), 1, payload_(0),
            pow(payload_(0), 2) + pow(payload_(1), 2);
        center_y_denominator << 1, quadrotor_1_(0), quadrotor_1_(1), 1, quadrotor_2_(0), quadrotor_2_(1), 1,
            payload_(0), payload_(1);
        double center_x = 0.5 * center_x_numerator.determinant() / center_x_denominator.determinant();
        double center_y = 0.5 * center_y_numerator.determinant() / center_y_denominator.determinant();
        xq1_center << quadrotor_1_(0) - center_x, quadrotor_1_(1) - center_y;
        xq2_center << quadrotor_2_(0) - center_x, quadrotor_2_(1) - center_y;
        xp_center << payload_(0) - center_x, payload_(1) - center_y;
        double radius = xq1_center.norm();
        Eigen::EigenSolver<Eigen::MatrixXd> Cs(C);
        Eigen::MatrixXcd Cvals = Cs.eigenvalues();
        Eigen::MatrixXd CvalsReal = Cvals.real();
        MatrixXf::Index CsMin;
        MatrixXf::Index boundsmax;
        double short_axis = CvalsReal.rowwise().sum().minCoeff(&CsMin);
        double bound_max = bound.maxCoeff(&boundsmax);
        if (short_axis > radius + bound_max)
        {
            quadrotor_1_ = D + xq1_center;
            quadrotor_2_ = D + xq2_center;
            payload_ = D + xp_center;
            change_lable_ = true;
            insetting_formation_type_ = Insetting_Formation_Algorithm::ELLIPSE_BASED_INSETTING_FORMATION;
        }
        else
        {
            change_lable_ = false;
            insetting_formation_type_ = Insetting_Formation_Algorithm::NO_ALGORITHM_USED;
#ifndef NDEBUG
            std::cout << "Ellipse-based insetting foramtion fail : The formation cannot be inset into the ellipse"
                      << std::endl;
#endif
        }
    }
    else
    {
        change_lable_ = false;
        insetting_formation_type_ = Insetting_Formation_Algorithm::NO_ALGORITHM_USED;
#ifndef NDEBUG
        std::cout << "Ellipse-based insetting foramtion fail : "
                     "The original formation relationship does not meet the dynamic requirements"
                  << std::endl;
#endif
    }
}
void formation::construct_formation_bound(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
                                          const Eigen::Vector4d &bound, const Eigen::Vector2d &limit_rate,
                                          const double length, MSKenv_t *existing_env)
{
    MSKenv_t *env;
    if (existing_env)
    {
        env = existing_env;
    }
    else
    {
        env = (MSKenv_t *)malloc(sizeof(MSKenv_t));
        check_res(MSK_makeenv(env, NULL));
    }
    MSKtask_t task = NULL;
    const double quadrotor_bound_x = bound(0);
    const double quadrotor_bound_y = bound(1);
    const double payload_bound_x = bound(2);
    const double payload_bound_y = bound(3);
    const double low_rate = limit_rate(0);
    const double up_rate = limit_rate(1);
    const int m = A.rows();
    const int n = A.cols();
    const int p_sum = 3;

    const int num_p = p_sum * n;
    const int num_b = p_sum * n * 2;
    const int num_s = (p_sum - 1) * n;
    const int num_g = p_sum - 1;

    int nvar = 0;
    ADD_VAR(p)
    ADD_VAR(b)
    ADD_VAR(s)
    ADD_VAR(g)

    const int ncon =
        p_sum * n * 2 + m * p_sum * 2 * 2 + p_sum - 1 + p_sum - 2 + p_sum - 1 + (p_sum - 1) * n + p_sum - 1;

    check_res(MSK_maketask(*env, ncon, 0, &task));

#ifndef NDEBUG
    /* Directs the log task stream to the 'printstr' function. */
    check_res(MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr));
#endif

    check_res(MSK_appendcons(task, ncon));

    check_res(MSK_appendvars(task, nvar));
    for (unsigned int i = 0; i < ndx_g.size(); i++)
    {
        check_res(MSK_putcj(task, ndx_g[i], 1.0));
    }

    for (int i = 0; i < nvar; i++)
    {
        check_res(MSK_putvarbound(task, i, MSK_BK_FR, -MSK_INFINITY, MSK_INFINITY));
    }
    int con_ndx = 0;

    for (int p_num = 0; p_num < p_sum - 1; p_num++)
    {
        MSKint32t subi[] = {ndx_p[n * p_num], ndx_b[2 * n * p_num]};
        MSKrealt vali[] = {1, -1};
        check_res(MSK_putarow(task, con_ndx, 2, subi, vali));
        check_res(MSK_putconbound(task, con_ndx, MSK_BK_FX, quadrotor_bound_x, quadrotor_bound_x));
        con_ndx = con_ndx + 1;
    }
    for (int p_num = 0; p_num < p_sum - 1; p_num++)
    {
        MSKint32t subi[] = {ndx_p[n * p_num], ndx_b[2 * n * p_num + 1]};
        MSKrealt vali[] = {-1, 1};
        check_res(MSK_putarow(task, con_ndx, 2, subi, vali));
        check_res(MSK_putconbound(task, con_ndx, MSK_BK_FX, quadrotor_bound_x, quadrotor_bound_x));
        con_ndx = con_ndx + 1;
    }
    for (int p_num = 0; p_num < p_sum - 1; p_num++)
    {
        MSKint32t subi[] = {ndx_p[n * p_num + 1], ndx_b[2 * n * p_num + 2]};
        MSKrealt vali[] = {1, -1};
        check_res(MSK_putarow(task, con_ndx, 2, subi, vali));
        check_res(MSK_putconbound(task, con_ndx, MSK_BK_FX, quadrotor_bound_y, quadrotor_bound_y));
        con_ndx = con_ndx + 1;
    }
    for (int p_num = 0; p_num < p_sum - 1; p_num++)
    {
        MSKint32t subi[] = {ndx_p[n * p_num + 1], ndx_b[2 * n * p_num + 3]};
        MSKrealt vali[] = {-1, 1};
        check_res(MSK_putarow(task, con_ndx, 2, subi, vali));
        check_res(MSK_putconbound(task, con_ndx, MSK_BK_FX, quadrotor_bound_y, quadrotor_bound_y));
        con_ndx = con_ndx + 1;
    }

    {
        MSKint32t subi[] = {ndx_p[n * (p_sum - 1)], ndx_b[2 * n * (p_sum - 1)]};
        MSKrealt vali[] = {1, -1};
        check_res(MSK_putarow(task, con_ndx, 2, subi, vali));
        check_res(MSK_putconbound(task, con_ndx, MSK_BK_FX, payload_bound_x, payload_bound_x));
        con_ndx = con_ndx + 1;
    }

    {
        MSKint32t subi[] = {ndx_p[n * (p_sum - 1)], ndx_b[2 * n * (p_sum - 1) + 1]};
        MSKrealt vali[] = {-1, 1};
        check_res(MSK_putarow(task, con_ndx, 2, subi, vali));
        check_res(MSK_putconbound(task, con_ndx, MSK_BK_FX, payload_bound_x, payload_bound_x));
        con_ndx = con_ndx + 1;
    }

    {
        MSKint32t subi[] = {ndx_p[n * (p_sum - 1) + 1], ndx_b[2 * n * (p_sum - 1) + 2]};
        MSKrealt vali[] = {1, -1};
        check_res(MSK_putarow(task, con_ndx, 2, subi, vali));
        check_res(MSK_putconbound(task, con_ndx, MSK_BK_FX, payload_bound_y, payload_bound_y));
        con_ndx = con_ndx + 1;
    }

    {
        MSKint32t subi[] = {ndx_p[n * (p_sum - 1) + 1], ndx_b[2 * n * (p_sum - 1) + 3]};
        MSKrealt vali[] = {-1, 1};
        check_res(MSK_putarow(task, con_ndx, 2, subi, vali));
        check_res(MSK_putconbound(task, con_ndx, MSK_BK_FX, payload_bound_y, payload_bound_y));
        con_ndx = con_ndx + 1;
    }

    for (int p_num = 0; p_num < p_sum; p_num++)
    {
        for (int i = 0; i < m; i++)
        {
            {
                MSKint32t subi[] = {ndx_b[2 * n * p_num], ndx_b[2 * n * p_num + 2]};
                MSKrealt vali[] = {A(i, n - 2), A(i, n - 1)};
                check_res(MSK_putarow(task, con_ndx, 2, subi, vali));
                check_res(MSK_putconbound(task, con_ndx, MSK_BK_UP, -MSK_INFINITY, b(i)));
                con_ndx = con_ndx + 1;
            }
            {
                MSKint32t subi[] = {ndx_b[2 * n * p_num], ndx_b[2 * n * p_num + 3]};
                MSKrealt vali[] = {A(i, n - 2), A(i, n - 1)};
                check_res(MSK_putarow(task, con_ndx, 2, subi, vali));
                check_res(MSK_putconbound(task, con_ndx, MSK_BK_UP, -MSK_INFINITY, b(i)));
                con_ndx = con_ndx + 1;
            }
        }
    }

    for (int p_num = 0; p_num < p_sum; p_num++)
    {
        for (int i = 0; i < m; i++)
        {
            {
                MSKint32t subi[] = {ndx_b[2 * n * p_num + 1], ndx_b[2 * n * p_num + 2]};
                MSKrealt vali[] = {A(i, n - 2), A(i, n - 1)};
                check_res(MSK_putarow(task, con_ndx, 2, subi, vali));
                check_res(MSK_putconbound(task, con_ndx, MSK_BK_UP, -MSK_INFINITY, b(i)));
                con_ndx = con_ndx + 1;
            }
            {
                MSKint32t subi[] = {ndx_b[2 * n * p_num + 1], ndx_b[2 * n * p_num + 3]};
                MSKrealt vali[] = {A(i, n - 2), A(i, n - 1)};
                check_res(MSK_putarow(task, con_ndx, 2, subi, vali));
                check_res(MSK_putconbound(task, con_ndx, MSK_BK_UP, -MSK_INFINITY, b(i)));
                con_ndx = con_ndx + 1;
            }
        }
    }

    {
        MSKint32t subi[] = {ndx_p[n * 1 - 1], ndx_p[n * p_sum - 1]};
        MSKrealt vali[] = {1, -1};
        check_res(MSK_putarow(task, con_ndx, 2, subi, vali));
        check_res(MSK_putconbound(task, con_ndx, MSK_BK_LO, 0, MSK_INFINITY));
        con_ndx = con_ndx + 1;
    }

    {
        MSKint32t subi[] = {ndx_p[n * 1 - 1], ndx_p[n * 2 - 1]};
        MSKrealt vali[] = {1, -1};
        check_res(MSK_putarow(task, con_ndx, 2, subi, vali));
        check_res(MSK_putconbound(task, con_ndx, MSK_BK_FX, 0, 0));
        con_ndx = con_ndx + 1;
    }

    {
        MSKint32t subi[] = {ndx_p[n * (p_sum - 1)], ndx_p[(n - 2) * 2]};
        MSKrealt vali[] = {-1, 1};
        check_res(MSK_putarow(task, con_ndx, 2, subi, vali));
        check_res(MSK_putconbound(task, con_ndx, MSK_BK_LO, 0, MSK_INFINITY));
        con_ndx = con_ndx + 1;
    }
    {
        MSKint32t subi[] = {ndx_p[n * (p_sum - 1)], ndx_p[(n - 1) * 2]};
        MSKrealt vali[] = {1, -1};
        check_res(MSK_putarow(task, con_ndx, 2, subi, vali));
        check_res(MSK_putconbound(task, con_ndx, MSK_BK_LO, 0, MSK_INFINITY));
        con_ndx = con_ndx + 1;
    }

    {
        MSKint32t subi[] = {ndx_p[n * (p_sum - 3)], ndx_p[n * (p_sum - 2)]};
        MSKrealt vali[] = {1, -1};
        check_res(MSK_putarow(task, con_ndx, 2, subi, vali));
        check_res(MSK_putconbound(task, con_ndx, MSK_BK_RA, low_rate * length, up_rate * length));
        con_ndx = con_ndx + 1;
    }

    {
        MSKint32t subi[] = {ndx_s[0], ndx_p[0], ndx_p[4]};
        MSKrealt vali[] = {1, -1, 1};
        check_res(MSK_putarow(task, con_ndx, 3, subi, vali));
        check_res(MSK_putconbound(task, con_ndx, MSK_BK_FX, 0, 0));
        con_ndx = con_ndx + 1;
    }
    {
        MSKint32t subi[] = {ndx_s[1], ndx_p[1], ndx_p[5]};
        MSKrealt vali[] = {1, -1, 1};
        check_res(MSK_putarow(task, con_ndx, 3, subi, vali));
        check_res(MSK_putconbound(task, con_ndx, MSK_BK_FX, 0, 0));
        con_ndx = con_ndx + 1;
    }
    {
        MSKint32t subi[] = {ndx_s[2], ndx_p[4], ndx_p[2]};
        MSKrealt vali[] = {1, -1, 1};
        check_res(MSK_putarow(task, con_ndx, 3, subi, vali));
        check_res(MSK_putconbound(task, con_ndx, MSK_BK_FX, 0, 0));
        con_ndx = con_ndx + 1;
    }
    {
        MSKint32t subi[] = {ndx_s[3], ndx_p[3], ndx_p[5]};
        MSKrealt vali[] = {1, -1, 1};
        check_res(MSK_putarow(task, con_ndx, 3, subi, vali));
        check_res(MSK_putconbound(task, con_ndx, MSK_BK_FX, 0, 0));
        con_ndx = con_ndx + 1;
    }

    {
        MSKint32t subi[] = {ndx_g[0], ndx_g[1]};
        MSKrealt vali[] = {-1, 1};
        check_res(MSK_putarow(task, con_ndx, 2, subi, vali));
        check_res(MSK_putconbound(task, con_ndx, MSK_BK_FX, 0, 0));
        con_ndx = con_ndx + 1;
    }
    {
        MSKint32t subi[] = {ndx_g[1]};
        MSKrealt vali[] = {1};
        check_res(MSK_putarow(task, con_ndx, 1, subi, vali));
        check_res(MSK_putconbound(task, con_ndx, MSK_BK_LO, length, MSK_INFINITY));
        con_ndx = con_ndx + 1;
    }
    assert(con_ndx == ncon);

    MSKint32t csub[3];
    for (int p_num = 0; p_num < p_sum - 1; p_num++)
    {
        csub[0] = ndx_s[p_num * 2];
        csub[1] = ndx_s[p_num * 2 + 1];
        csub[2] = ndx_g[p_num];
        check_res(MSK_appendcone(task, MSK_CT_RQUAD, 0.0, 3, csub));
    }

    check_res(MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MAXIMIZE));

    double *xx;
    Eigen::VectorXd xq1_b(2);
    Eigen::VectorXd xq2_b(2);
    Eigen::VectorXd xp_b(2);
    MSKrescodee trmcode;
    MSKrescodee res = MSK_optimizetrm(task, &trmcode);
    MSK_solutionsummary(task, MSK_STREAM_MSG);

    check_res(res);

    MSKsolstae solsta;
    MSK_getsolsta(task, MSK_SOL_ITR, &solsta);

    switch (solsta)
    {
        case MSK_SOL_STA_OPTIMAL:
        case MSK_SOL_STA_NEAR_OPTIMAL:
            xx = (double *)MSK_calloctask(task, nvar, sizeof(MSKrealt));

            MSK_getxx(task, MSK_SOL_ITR, xx);
            xq1_b(0) = xx[0];
            xq1_b(1) = xx[1];
            xq2_b(0) = xx[2];
            xq2_b(1) = xx[3];
            xp_b(0) = xx[4];
            xp_b(1) = xx[5];

            formation::set_quadrotor_1_bound(xq1_b);
            formation::set_quadrotor_2_bound(xq2_b);
            formation::set_payload_bound(xp_b);
            change_bound_lable_ = true;

            MSK_freetask(task, xx);
            break;
        case MSK_SOL_STA_DUAL_INFEAS_CER:
        case MSK_SOL_STA_PRIM_INFEAS_CER:
        case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
        case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
#ifndef NDEBUG
            std::cout << "construct Maximum Virtual Boundary(MVB) fail : "
                         "Primal or dual infeasibility certificate found."
                      << std::endl;
#endif
            change_bound_lable_ = false;
            break;
        case MSK_SOL_STA_UNKNOWN:
#ifndef NDEBUG
            std::cout << "construct Maximum Virtual Boundary(MVB) fail : "
                         "The status of the solution could not be determined."
                      << std::endl;
#endif
            change_bound_lable_ = false;
            break;
        default:
#ifndef NDEBUG
            printf("construct Maximum Virtual Boundary(MVB) fail : other solution status: %d\n", solsta);
#endif
            change_bound_lable_ = false;
            break;
    }
    MSKrealt obj_val;
    MSK_getprimalobj(task, MSK_SOL_ITR, &obj_val);

    MSK_deletetask(&task);
    if (!existing_env)
    {
        MSK_deleteenv(env);
        free(env);
    }
}
void formation::construct_formation_high_diff(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
                                              const Eigen::Vector4d &bound, const double low_rate, const double length,
                                              const Eigen::Vector2d &initial_point, const double T_max,
                                              const double G_payload)
{
    const double tolerance_value = 1e-6;
    const double quadrotor_bound_x = bound(0);
    const double quadrotor_bound_y = bound(1);
    const double payload_bound_x = bound(2);
    const double payload_bound_y = bound(3);
    const int m = A.rows();
    int ncon;
    bool flag = true;
    std::vector<double> res_g{};
    Eigen::VectorXd q1_p, p_q2;
    q1_p = quadrotor_1_ - payload_;
    p_q2 = payload_ - quadrotor_2_;
    std::vector<double> vq1_p, vp_q2;
    vq1_p.resize(q1_p.size());
    VectorXd::Map(vq1_p.data(), vq1_p.size()) = q1_p;
    vp_q2.resize(p_q2.size());
    VectorXd::Map(vp_q2.data(), vp_q2.size()) = p_q2;
    SX x = SX::sym("x", 8);
    SX xq1_p = vertcat(x(0) - x(4), x(1) - x(5));
    SX xp_q2 = vertcat(x(4) - x(2), x(5) - x(3));
    SX f = -(dot(vq1_p, xq1_p) + dot(vp_q2, xp_q2));
    SX g;
    std::vector<double> lbg{};
    std::vector<double> ubg{};

    for (int i = 0; i < m; i++)
    {
        g = vertcat(g, A(i, 0) * x(0) + A(i, 1) * x(1));
        lbg.push_back(-inf);
        ubg.push_back(b(i) - A(i, 0) * quadrotor_bound_x - A(i, 1) * quadrotor_bound_y);
        g = vertcat(g, A(i, 0) * x(0) + A(i, 1) * x(1));
        lbg.push_back(-inf);
        ubg.push_back(b(i) - A(i, 0) * quadrotor_bound_x + A(i, 1) * quadrotor_bound_y);
        g = vertcat(g, A(i, 0) * x(0) + A(i, 1) * x(1));
        lbg.push_back(-inf);
        ubg.push_back(b(i) + A(i, 0) * quadrotor_bound_x + A(i, 1) * quadrotor_bound_y);
        g = vertcat(g, A(i, 0) * x(0) + A(i, 1) * x(1));
        lbg.push_back(-inf);
        ubg.push_back(b(i) + A(i, 0) * quadrotor_bound_x - A(i, 1) * quadrotor_bound_y);
    }

    for (int i = 0; i < m; i++)
    {
        g = vertcat(g, A(i, 0) * x(2) + A(i, 1) * x(3));
        lbg.push_back(-inf);
        ubg.push_back(b(i) - A(i, 0) * quadrotor_bound_x - A(i, 1) * quadrotor_bound_y);
        g = vertcat(g, A(i, 0) * x(2) + A(i, 1) * x(3));
        lbg.push_back(-inf);
        ubg.push_back(b(i) - A(i, 0) * quadrotor_bound_x + A(i, 1) * quadrotor_bound_y);
        g = vertcat(g, A(i, 0) * x(2) + A(i, 1) * x(3));
        lbg.push_back(-inf);
        ubg.push_back(b(i) + A(i, 0) * quadrotor_bound_x + A(i, 1) * quadrotor_bound_y);
        g = vertcat(g, A(i, 0) * x(2) + A(i, 1) * x(3));
        lbg.push_back(-inf);
        ubg.push_back(b(i) + A(i, 0) * quadrotor_bound_x - A(i, 1) * quadrotor_bound_y);
    }

    for (int i = 0; i < m; i++)
    {
        g = vertcat(g, A(i, 0) * x(4) + A(i, 1) * x(5));
        lbg.push_back(-inf);
        ubg.push_back(b(i) - A(i, 0) * payload_bound_x - A(i, 1) * payload_bound_y);
        g = vertcat(g, A(i, 0) * x(4) + A(i, 1) * x(5));
        lbg.push_back(-inf);
        ubg.push_back(b(i) - A(i, 0) * payload_bound_x + A(i, 1) * payload_bound_y);
        g = vertcat(g, A(i, 0) * x(4) + A(i, 1) * x(5));
        lbg.push_back(-inf);
        ubg.push_back(b(i) + A(i, 0) * payload_bound_x + A(i, 1) * payload_bound_y);
        g = vertcat(g, A(i, 0) * x(4) + A(i, 1) * x(5));
        lbg.push_back(-inf);
        ubg.push_back(b(i) + A(i, 0) * payload_bound_x - A(i, 1) * payload_bound_y);
    }

    g = vertcat(g, pow((x(2) - x(4)), 2) + pow((x(3) - x(5)), 2) - pow(length, 2));
    lbg.push_back(0);
    ubg.push_back(0);
    g = vertcat(g, pow((x(0) - x(4)), 2) + pow((x(1) - x(5)), 2) - pow(length, 2));
    lbg.push_back(0);
    ubg.push_back(0);

    g = vertcat(g, pow((x(0) - x(2)), 2) + pow((x(1) - x(3)), 2) - pow(low_rate * length, 2));
    lbg.push_back(0);
    ubg.push_back(inf);

    g = vertcat(g, x(0) - x(4));
    lbg.push_back(0);
    ubg.push_back(inf);
    g = vertcat(g, x(2) - x(4));
    lbg.push_back(-inf);
    ubg.push_back(0);

    g = vertcat(g, x(5) - x(1));
    lbg.push_back(-inf);
    ubg.push_back(0);
    g = vertcat(g, x(5) - x(3));
    lbg.push_back(-inf);
    ubg.push_back(0);

    g = vertcat(g, x(6) * (x(0) - x(4)) / length - x(7) * (x(4) - x(2)) / length);
    lbg.push_back(0);
    ubg.push_back(0);
    g = vertcat(g, x(6) * (x(1) - x(5)) / length + x(7) * (x(3) - x(5)) / length - G_payload);
    lbg.push_back(0);
    ubg.push_back(0);

    ncon = ubg.size();

    std::vector<double> x0{};
    for (int i = 0; i < 6; i++)
    {
        x0.push_back(initial_point(i % 2));
    }
    x0.push_back(T_max);
    x0.push_back(T_max);
    std::vector<double> lbx = {-inf, -inf, -inf, -inf, -inf, -inf, 0, 0};
    std::vector<double> ubx = {inf, inf, inf, inf, inf, inf, T_max, T_max};

    SXDict nlp = {{"x", x}, {"f", f}, {"g", g}};

    Dict opts;
    opts["verbose_init"] = false;
    opts["verbose"] = false;
    opts["print_time"] = false;
    opts["ipopt.print_level"] = 0;

#ifndef NDEBUG
    opts["ipopt.print_level"] = 5;
#endif

    Function solver = nlpsol("solver", "ipopt", nlp, opts);
    std::map<std::string, DM> arg, res;

    arg["lbx"] = lbx;
    arg["ubx"] = ubx;
    arg["lbg"] = lbg;
    arg["ubg"] = ubg;
    arg["x0"] = x0;
    res = solver(arg);

    for (int i = 0; i < ncon; i++)
    {
        res_g.push_back((double)res["g"](i));
    }

    for (int i = 0; i < ncon; i++)
    {
        if (res_g[i] < lbg[i] - tolerance_value || res_g[i] > ubg[i] + tolerance_value)
        {
            flag = false;
            break;
        }
    }

    if (flag)
    {
        Eigen::VectorXd xq1(2), xq2(2), xq3(2);
        xq1(0) = (double)res["x"](0);
        xq1(1) = (double)res["x"](1);
        xq2(0) = (double)res["x"](2);
        xq2(1) = (double)res["x"](3);
        xq3(0) = (double)res["x"](4);
        xq3(1) = (double)res["x"](5);
        formation::set_quadrotor_1(xq1);
        formation::set_quadrotor_2(xq2);
        formation::set_payload(xq3);
        formation::set_change_lable(flag);
        insetting_formation_type_ = Insetting_Formation_Algorithm::HEIGHT_FREE_CONSTRAINED_INSETTING_FORMATION;
    }
    else
    {
        formation::set_change_lable(flag);
        insetting_formation_type_ = Insetting_Formation_Algorithm::NO_ALGORITHM_USED;
#ifndef NDEBUG
        std::cout << "Height-free constrained-based insetting formation fail : "
                     "The obstacle avoidance environment is too tough"
                  << std::endl;
#endif
    }
}
void formation::construct_formation_high_diff(const Eigen::MatrixXd &A, const Eigen::VectorXd &b,
                                              const Eigen::Vector4d &bound, const Eigen::Vector2d &limit_rate,
                                              const double length, const Eigen::Vector2d &initial_point)
{
    const double tolerance_value = 1e-6;
    const double quadrotor_bound_x = bound(0);
    const double quadrotor_bound_y = bound(1);
    const double payload_bound_x = bound(2);
    const double payload_bound_y = bound(3);
    const double low_rate = limit_rate(0);
    const double up_rate = limit_rate(1);
    const int m = A.rows();
    int ncon;
    bool flag = true;
    std::vector<double> res_g{};

    Eigen::VectorXd q1_p, p_q2;
    q1_p = quadrotor_1_ - payload_;
    p_q2 = payload_ - quadrotor_2_;

    std::vector<double> vq1_p, vp_q2;
    vq1_p.resize(q1_p.size());
    VectorXd::Map(vq1_p.data(), vq1_p.size()) = q1_p;
    vp_q2.resize(p_q2.size());
    VectorXd::Map(vp_q2.data(), vp_q2.size()) = p_q2;

    SX x = SX::sym("x", 6);
    SX xq1_p = vertcat(x(0) - x(4), x(1) - x(5));
    SX xp_q2 = vertcat(x(4) - x(2), x(5) - x(3));

    SX f = -(dot(vq1_p, xq1_p) + dot(vp_q2, xp_q2));

    SX g;

    std::vector<double> lbg{};
    std::vector<double> ubg{};

    for (int i = 0; i < m; i++)
    {
        g = vertcat(g, A(i, 0) * x(0) + A(i, 1) * x(1));
        lbg.push_back(-inf);
        ubg.push_back(b(i) - A(i, 0) * quadrotor_bound_x - A(i, 1) * quadrotor_bound_y);
        g = vertcat(g, A(i, 0) * x(0) + A(i, 1) * x(1));
        lbg.push_back(-inf);
        ubg.push_back(b(i) - A(i, 0) * quadrotor_bound_x + A(i, 1) * quadrotor_bound_y);
        g = vertcat(g, A(i, 0) * x(0) + A(i, 1) * x(1));
        lbg.push_back(-inf);
        ubg.push_back(b(i) + A(i, 0) * quadrotor_bound_x + A(i, 1) * quadrotor_bound_y);
        g = vertcat(g, A(i, 0) * x(0) + A(i, 1) * x(1));
        lbg.push_back(-inf);
        ubg.push_back(b(i) + A(i, 0) * quadrotor_bound_x - A(i, 1) * quadrotor_bound_y);
    }
    for (int i = 0; i < m; i++)
    {
        g = vertcat(g, A(i, 0) * x(2) + A(i, 1) * x(3));
        lbg.push_back(-inf);
        ubg.push_back(b(i) - A(i, 0) * quadrotor_bound_x - A(i, 1) * quadrotor_bound_y);
        g = vertcat(g, A(i, 0) * x(2) + A(i, 1) * x(3));
        lbg.push_back(-inf);
        ubg.push_back(b(i) - A(i, 0) * quadrotor_bound_x + A(i, 1) * quadrotor_bound_y);
        g = vertcat(g, A(i, 0) * x(2) + A(i, 1) * x(3));
        lbg.push_back(-inf);
        ubg.push_back(b(i) + A(i, 0) * quadrotor_bound_x + A(i, 1) * quadrotor_bound_y);
        g = vertcat(g, A(i, 0) * x(2) + A(i, 1) * x(3));
        lbg.push_back(-inf);
        ubg.push_back(b(i) + A(i, 0) * quadrotor_bound_x - A(i, 1) * quadrotor_bound_y);
    }
    for (int i = 0; i < m; i++)
    {
        g = vertcat(g, A(i, 0) * x(4) + A(i, 1) * x(5));
        lbg.push_back(-inf);
        ubg.push_back(b(i) - A(i, 0) * payload_bound_x - A(i, 1) * payload_bound_y);
        g = vertcat(g, A(i, 0) * x(4) + A(i, 1) * x(5));
        lbg.push_back(-inf);
        ubg.push_back(b(i) - A(i, 0) * payload_bound_x + A(i, 1) * payload_bound_y);
        g = vertcat(g, A(i, 0) * x(4) + A(i, 1) * x(5));
        lbg.push_back(-inf);
        ubg.push_back(b(i) + A(i, 0) * payload_bound_x + A(i, 1) * payload_bound_y);
        g = vertcat(g, A(i, 0) * x(4) + A(i, 1) * x(5));
        lbg.push_back(-inf);
        ubg.push_back(b(i) + A(i, 0) * payload_bound_x - A(i, 1) * payload_bound_y);
    }

    g = vertcat(g, pow((x(2) - x(4)), 2) + pow((x(3) - x(5)), 2) - pow(length, 2));
    lbg.push_back(0);
    ubg.push_back(0);
    g = vertcat(g, pow((x(0) - x(4)), 2) + pow((x(1) - x(5)), 2) - pow(length, 2));
    lbg.push_back(0);
    ubg.push_back(0);

    g = vertcat(g, pow((x(0) - x(2)), 2) + pow((x(1) - x(3)), 2) - pow(low_rate * length, 2));
    lbg.push_back(0);
    ubg.push_back(inf);
    g = vertcat(g, pow((x(0) - x(2)), 2) + pow((x(1) - x(3)), 2) - pow(up_rate * length, 2));
    lbg.push_back(-inf);
    ubg.push_back(0);

    g = vertcat(g, x(0) - x(4));
    lbg.push_back(0);
    ubg.push_back(inf);
    g = vertcat(g, x(2) - x(4));
    lbg.push_back(-inf);
    ubg.push_back(0);

    g = vertcat(g, x(5) - x(1));
    lbg.push_back(-inf);
    ubg.push_back(0);
    g = vertcat(g, x(5) - x(3));
    lbg.push_back(-inf);
    ubg.push_back(0);

    ncon = ubg.size();

    std::vector<double> x0{};
    for (int i = 0; i < 6; i++)
    {
        x0.push_back(initial_point(i % 2));
    }
    std::vector<double> lbx = {-inf, -inf, -inf, -inf, -inf, -inf};
    std::vector<double> ubx = {inf, inf, inf, inf, inf, inf};

    SXDict nlp = {{"x", x}, {"f", f}, {"g", g}};

    Dict opts;
    opts["verbose_init"] = false;
    opts["verbose"] = false;
    opts["print_time"] = false;
    opts["ipopt.print_level"] = 0;

#ifndef NDEBUG
    opts["ipopt.print_level"] = 5;
#endif

    Function solver = nlpsol("solver", "ipopt", nlp, opts);
    std::map<std::string, DM> arg, res;

    arg["lbx"] = lbx;
    arg["ubx"] = ubx;
    arg["lbg"] = lbg;
    arg["ubg"] = ubg;
    arg["x0"] = x0;
    res = solver(arg);

    for (int i = 0; i < ncon; i++)
    {
        res_g.push_back((double)res["g"](i));
    }

    for (int i = 0; i < ncon; i++)
    {
        if (res_g[i] < lbg[i] - tolerance_value || res_g[i] > ubg[i] + tolerance_value)
        {
            flag = false;
            break;
        }
    }

    if (flag)
    {
        Eigen::VectorXd xq1(2), xq2(2), xq3(2);
        xq1(0) = (double)res["x"](0);
        xq1(1) = (double)res["x"](1);
        xq2(0) = (double)res["x"](2);
        xq2(1) = (double)res["x"](3);
        xq3(0) = (double)res["x"](4);
        xq3(1) = (double)res["x"](5);
        formation::set_quadrotor_1(xq1);
        formation::set_quadrotor_2(xq2);
        formation::set_payload(xq3);
        formation::set_change_lable(flag);
        insetting_formation_type_ = Insetting_Formation_Algorithm::HEIGHT_FREE_CONSTRAINED_INSETTING_FORMATION;
    }
    else
    {
        formation::set_change_lable(flag);
        insetting_formation_type_ = Insetting_Formation_Algorithm::NO_ALGORITHM_USED;
#ifndef NDEBUG
        std::cout << "Height-free constrained-based insetting formation fail : "
                     "The obstacle avoidance environment is too tough"
                  << std::endl;
#endif
    }
}

void formation::construct_formation(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, const Eigen::MatrixXd &C,
                                    const Eigen::VectorXd &D, const Eigen::Vector4d &bound,
                                    const Eigen::Vector2d &limit_rate, const double length,
                                    const Eigen::Vector2d &initial_point)
{
    const double tolerance_value = 1e-6;
    formation::construct_formation_ellipse(C, D, bound, limit_rate, length);
    if (!change_lable_)
    {
        if (quadrotor_1_(1) - quadrotor_2_(1) < tolerance_value && quadrotor_1_(1) - quadrotor_2_(1) > -tolerance_value)
        {
            formation::construct_formation_high_equ(A, b, bound, limit_rate, length);
            if (!change_lable_)
            {
                formation::construct_formation_high_diff(A, b, bound, limit_rate, length, initial_point);
            }
        }
        else
        {
            formation::construct_formation_high_diff(A, b, bound, limit_rate, length, initial_point);
        }
    }
}

void formation::construct_formation(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, const Eigen::MatrixXd &C,
                                    const Eigen::VectorXd &D, const Eigen::Vector4d &bound,
                                    const Eigen::Vector2d &limit_rate, const double length,
                                    const Eigen::Vector2d &initial_point, const double T_max, const double G_payload)
{
    const double tolerance_value = 1e-6;
    formation::construct_formation_ellipse(C, D, bound, limit_rate, length);
    if (!change_lable_)
    {
        if (quadrotor_1_(1) - quadrotor_2_(1) < tolerance_value && quadrotor_1_(1) - quadrotor_2_(1) > -tolerance_value)
        {
            formation::construct_formation_high_equ(A, b, bound, limit_rate, length);
            if (!change_lable_)
            {
                formation::construct_formation_high_diff(A, b, bound, limit_rate(0), length, initial_point, T_max,
                                                         G_payload);
            }
        }
        else
        {
            formation::construct_formation_high_diff(A, b, bound, limit_rate(0), length, initial_point, T_max,
                                                     G_payload);
        }
    }
}
} // namespace inset