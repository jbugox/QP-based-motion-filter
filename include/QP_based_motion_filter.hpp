#ifndef QP_based_motion_filter_hpp
#define QP_based_motion_filter_hpp

#include <qpOASES.hpp>


USING_NAMESPACE_QPOASES


class QP_based_motion_filter
{

public:

	enum
	{
		 DOF = 6, // degree of freedom of the manipulator
		 CONSTRAINTS_NUM = 4,  // constraints set num
		 k_slack = 6, // slack variable num

	};

private:

	real_t const omega_j = 1e6; // slack variable size

	real_t const filter_eps =1.192092896e-07F; // epsilons

	// arrays
	bool * apply_epsilon_flags;
	real_t* delta_i;
	real_t* tau_n;
	real_t* H;
	real_t* g;
	real_t* lb;
	real_t* ub;
	real_t* A;
	real_t* lbA;
	real_t* ubA;
	real_t* xOpt;
	real_t* yOpt;


	// qpOASES solving classes
	SQProblem example;
	Options options;

	// Epsilon apply
	void applyEpsilon(qpOASES::real_t lbA_in[DOF * CONSTRAINTS_NUM], qpOASES::real_t ubA_in[DOF * CONSTRAINTS_NUM]) const;


	void get_solution(const bool get_dual_sol = false) const;

	// printing
	static void print_tau(qpOASES::real_t tau[DOF], const bool & is_tau_c);
	static void print_Minv(qpOASES::real_t Minv[DOF*DOF]);

public:

	QP_based_motion_filter();
	virtual ~QP_based_motion_filter();

	qpOASES::real_t getEPS();

	// lb, ub setting (optional)
	void setlb(qpOASES::real_t lb_input[DOF+k_slack]) const;
	void setub(qpOASES::real_t ub_input[DOF+k_slack]) const;

	// printing
	void print_A() const;
	void print_xOpt() const;
	void print_yOpt() const;
	void print_tau() const;
	void print_rho() const;


	// solving
	void solve(qpOASES::real_t tau[DOF], const qpOASES::real_t Minv[DOF*DOF], const qpOASES::real_t JMinv[DOF*DOF], qpOASES::real_t lbA_in[DOF*CONSTRAINTS_NUM],  qpOASES::real_t ubA_in[DOF*CONSTRAINTS_NUM], const bool constraints_choices[CONSTRAINTS_NUM], qpOASES::int_t nWSR_input, const bool is_apply_epsilon, const bool is_print_tau);


};

#endif





