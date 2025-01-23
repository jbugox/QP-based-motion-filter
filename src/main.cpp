#include <iostream>
#include "QP_based_motion_filter.hpp"


int main() {

    QP_based_motion_filter filter;

    // From here, assume that you can find robot dynamic parameters:
    // M, C, g, tau, qddot, qdot, J ...
    // Also, all matrices must be arranged in following order.
    // e.g. [ 1 3 2 ]
    //      [ 0 1 5 ]
    //      [ 2 0 1 ] => [1 3 2 0 1 5 2 0 1]


    // This is just arbitrary values
    // DOF of manipulator is 6, 4 sets of constraints
    qpOASES::real_t tau[QP_based_motion_filter::DOF] = {-1, 0, 1.5, 0, -3.2, 2};
    qpOASES::real_t M[QP_based_motion_filter::DOF*QP_based_motion_filter::DOF] = {};

    // By applying the constraints, the equation will be organized in the following form
    // case of joint constraints: lbA_in < inv(M)*tau < ubA_in
    // case of task constraints:  lbA_in < J*inv(M)*tau < ubA_in
    qpOASES::real_t Minv[QP_based_motion_filter::DOF] = {};
    qpOASES::real_t JMinv[QP_based_motion_filter::DOF] = {};
    qpOASES::real_t lbA_in[QP_based_motion_filter::DOF*QP_based_motion_filter::CONSTRAINTS_NUM] = {};
    qpOASES::real_t ubA_in[QP_based_motion_filter::DOF*QP_based_motion_filter::CONSTRAINTS_NUM] = {};

    constexpr bool constraints_choices[4] = {false, false, false, false};
    // false: joint constraints, true: task constraints
    filter.solve(tau, Minv, JMinv, lbA_in, ubA_in, constraints_choices,
        25,
        true,
        true);



    std::cout << "DONE" << std::endl;
    return 0;
}
