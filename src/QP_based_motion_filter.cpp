
#include "QP_based_motion_filter.hpp"


QP_based_motion_filter::QP_based_motion_filter() : example(DOF + k_slack, DOF*CONSTRAINTS_NUM)
{
	// applying epsilons
	apply_epsilon_flags = new bool[CONSTRAINTS_NUM]{true, true, true, true};

	delta_i = new real_t[CONSTRAINTS_NUM]{ 0.0, 0.0, 1.0, 0.0}; // apply rou setting: CLF
	// 0: not apply rou
	// 1: apply rou

	// input tau
	tau_n = new real_t[DOF]{};


	// QP problem:
	// (tau_c)^2 - 2*tau_n*tau_c
	// = 0.5*tau_c^2 - tau_n*tau_c
	// = 0.5*tau_c^2 - tau_n*tau_c + 0.5*W*lou^2
    //  = 0.5*transpose(x)*H*x + transpose(g)*x;
	H = new real_t[(DOF + k_slack) * (DOF + k_slack)]
	{   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, omega_j, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, omega_j, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, omega_j, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, omega_j, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, omega_j, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, omega_j };

	// g = [-transpose(tau_n), 0]
	g = new real_t[DOF + k_slack]{ -tau_n[0], -tau_n[1],-tau_n[2],-tau_n[3],-tau_n[4],-tau_n[5], 0, 0, 0, 0, 0, 0 };

	// constraints bound, initial value: inf
	lb = new real_t[DOF + k_slack]{ -INFTY, -INFTY, -INFTY, -INFTY ,-INFTY ,-INFTY ,-INFTY ,-INFTY ,-INFTY ,-INFTY ,-INFTY ,-INFTY };
	ub = new real_t[DOF + k_slack]{ INFTY,  INFTY,  INFTY,  INFTY,  INFTY,  INFTY,  INFTY,  INFTY,  INFTY,  INFTY,  INFTY,  INFTY };


	A = new real_t[DOF * CONSTRAINTS_NUM * (DOF + k_slack)]
	{
	/*1  2  3  4  5  6  1  2  3  4  5  6*/

		/*constraints 1 */
		/* 1 */      0, 0, 0, 0, 0, 0, delta_i[0],          0,          0,          0,          0,          0,
		/* 2 */      0, 0, 0, 0, 0, 0,          0, delta_i[0],          0,          0,          0,          0,
		/* 3 */      0, 0, 0, 0, 0, 0,          0,          0, delta_i[0],          0,          0,          0,
		/* 4 */      0, 0, 0, 0, 0, 0,          0,          0,          0, delta_i[0],          0,          0,
		/* 5 */      0, 0, 0, 0, 0, 0,          0,          0,          0,          0, delta_i[0],          0,
		/* 6 */      0, 0, 0, 0, 0, 0,          0,          0,          0,          0,          0, delta_i[0],

		/* constraints 2 */
		/* 7 */      0, 0, 0, 0, 0, 0, delta_i[1],          0,          0,          0,          0,          0,
		/* 8 */      0, 0, 0, 0, 0, 0,          0, delta_i[1],          0,          0,          0,          0,
		/* 9 */      0, 0, 0, 0, 0, 0,          0,          0, delta_i[1],          0,          0,          0,
		/* 10*/      0, 0, 0, 0, 0, 0,          0,          0,          0, delta_i[1],          0,          0,
		/* 11*/      0, 0, 0, 0, 0, 0,          0,          0,          0,          0, delta_i[1],          0,
		/* 12*/      0, 0, 0, 0, 0, 0,          0,          0,          0,          0,          0, delta_i[1],

		/* constraints 3 */
		/* 13*/      0, 0, 0, 0, 0, 0, delta_i[2],          0,          0,          0,          0,          0,
		/* 14*/      0, 0, 0, 0, 0, 0,          0, delta_i[2],          0,          0,          0,          0,
		/* 15*/      0, 0, 0, 0, 0, 0,          0,          0, delta_i[2],          0,          0,          0,
		/* 16*/      0, 0, 0, 0, 0, 0,          0,          0,          0, delta_i[2],          0,          0,
		/* 17*/      0, 0, 0, 0, 0, 0,          0,          0,          0,          0, delta_i[2],          0,
		/* 18*/      0, 0, 0, 0, 0, 0,          0,          0,          0,          0,          0, delta_i[2],

		/* constraints 4 */
		/* 19*/      0, 0, 0, 0, 0, 0, delta_i[3],          0,          0,          0,          0,          0,
		/* 20*/      0, 0, 0, 0, 0, 0,          0, delta_i[3],          0,          0,          0,          0,
		/* 21*/      0, 0, 0, 0, 0, 0,          0,          0, delta_i[3],          0,          0,          0,
		/* 22*/      0, 0, 0, 0, 0, 0,          0,          0,          0, delta_i[3],          0,          0,
		/* 23*/      0, 0, 0, 0, 0, 0,          0,          0,          0,          0, delta_i[3],          0,
		/* 24*/      0, 0, 0, 0, 0, 0,          0,          0,          0,          0,          0, delta_i[3],

	};


	lbA = new real_t[DOF * CONSTRAINTS_NUM]{};
	ubA = new real_t[DOF * CONSTRAINTS_NUM]{};


	xOpt = new real_t[DOF + k_slack]{};
	yOpt = new real_t[(DOF + k_slack) + CONSTRAINTS_NUM * DOF]{};


	// Function to disable output
	options.printLevel = qpOASES::PrintLevel::PL_NONE;
	
	example.setOptions(options);

	int_t nwsr = 15;
	example.init(H, g, A, lb, ub, lbA, ubA, nwsr);


}

QP_based_motion_filter::~QP_based_motion_filter()
{
	delete[] apply_epsilon_flags;
	delete[] delta_i;
	delete[] tau_n;
	delete[] H;
	delete[] g;
	delete[] lb;
	delete[] ub;
	delete[] A;
	delete[] lbA;
	delete[] ubA;
	delete[] xOpt;
	delete[] yOpt;

}


qpOASES::real_t QP_based_motion_filter::getEPS()
{
	return filter_eps;
}

void QP_based_motion_filter::setlb(qpOASES::real_t lb_input[]) const {
	for (int i = 0; i < DOF + k_slack; i++) lb[i] = lb_input[i];

}
void QP_based_motion_filter::setub(qpOASES::real_t ub_input[]) const {
	for (int i = 0; i < DOF + k_slack; i++) ub[i] = ub_input[i];
}

void QP_based_motion_filter::applyEpsilon(qpOASES::real_t lbA_in[DOF * CONSTRAINTS_NUM], qpOASES::real_t ubA_in[DOF * CONSTRAINTS_NUM]) const {
	for (int j = 0; j < CONSTRAINTS_NUM; j++)
	{
		if (apply_epsilon_flags[j])
		{
			for (int i = 0; i < DOF; i++)
			{
				lbA_in[i+DOF*j] += filter_eps;
				ubA_in[i+DOF*j] -= filter_eps;
			}
		}
	}

}


void QP_based_motion_filter::get_solution(const bool get_dual_sol) const {
	example.getPrimalSolution(xOpt);
	if (get_dual_sol) example.getDualSolution(yOpt);

}

void QP_based_motion_filter::print_tau(qpOASES::real_t tau[DOF], const bool & is_tau_c)
{
	if (is_tau_c) printf("\n tau_c: ");
	else          printf("\n tau_n: ");
	for (int i = 0; i < DOF; i++) printf("%.16lf ", tau[i]);
	if (is_tau_c) printf("\n");
}
void QP_based_motion_filter::print_tau() const
{
	printf("\nOptimized Tau: ");
	for (int i = 0; i < DOF; i++) {

		printf("%.16lf, ", xOpt[i]);
	}

}

void QP_based_motion_filter::print_Minv(qpOASES::real_t Minv[DOF*DOF])
{
	printf("\n Minv(input): [\n");
	for (int i = 0; i < DOF; i++)
	{
		for (int j = 0; j < DOF; j++)
		{
			printf("%.2lf, ", Minv[i*DOF+j]);
		}
		printf("\n");

	}
	printf("]\n\n");

}

void QP_based_motion_filter::print_A() const
{
	printf("\n A: [ \n");
	for (int i = 0; i < DOF*CONSTRAINTS_NUM; i++)
	{
		for (int j = 0; j < (DOF + k_slack); j++)
		{
			printf("%.2lf, ", A[i*(DOF + k_slack) + j]);
		}
		printf("\n");

	}

	printf("]\n\n");

}

void QP_based_motion_filter::print_rho() const
{
	printf("\nSlack Var Rho: [ ");
	for (int i = 0; i < k_slack; i++) {

		printf("%.16lf, ", xOpt[DOF+i]);
	}
	printf("]\n\n");
}


void QP_based_motion_filter::print_xOpt() const
{
	printf("\nxOpt: [ ");
	for (int i = 0; i < DOF + k_slack; i++) {
		if (i % DOF == 0) printf("|\n");

		printf("%e, ", xOpt[i]);
	}
	printf("]\n\n");
}

void QP_based_motion_filter::print_yOpt() const
{

	printf("\nyOpt: [ ");
	for (int i = 0; i < (DOF + k_slack) + CONSTRAINTS_NUM * (DOF); i++) {
		if (i % DOF == 0) printf("|\n");

		printf("%e, ", yOpt[i]);
	}
	printf("]\n\n");
}

void QP_based_motion_filter::solve(qpOASES::real_t tau[DOF], const qpOASES::real_t Minv[DOF*DOF], const qpOASES::real_t JMinv[DOF*DOF], qpOASES::real_t lbA_in[DOF*CONSTRAINTS_NUM],  qpOASES::real_t ubA_in[DOF*CONSTRAINTS_NUM], const bool constraints_choices[CONSTRAINTS_NUM], qpOASES::int_t nWSR_input, const bool is_apply_epsilon, const bool is_print_tau)
{
	/*
	 *	tau : input & output
	 *	Minv: inv(M), M is Mass Matrix (M(q))
	 *	JMinv : J*inv(M), J is Jacobian Matrix. For taskspace constraints.
	 *	lbA_in, ubA_in : lbA_in < inv(M)*tau < ubA_in
	 *  constraints_choices: Jointspace(false) or Taskspace(true).
	 *  e.g. joint, joint, joint, task: {false, false, false, true}
	 * nWSR_input : Maximum calculation count for QP sloving.
	 * is_apply_epsilon : true -> a < x < b
	 *                    false -> a <= x <= b
	 * is_print_tau : print tau before, and tau after.
	 */


	// tau_n setting
	for (int i = 0; i < DOF; i++) g[i] = -tau[i];

	// epsilon applying
	if (is_apply_epsilon) applyEpsilon(lbA_in, ubA_in);

	// A setting
	for (int i = 0; i < DOF; i++)
	{
		for (int j = 0; j < DOF; j++)
		{
			for (int k = 0; k < CONSTRAINTS_NUM; k++)
			{
				if (constraints_choices[k] == true) A[(k * DOF + i)*(DOF + k_slack) + j] = JMinv[i * DOF + j]; // q constraints
				else A[(k * DOF + i)*(DOF + k_slack) + j] = Minv[i * DOF + j]; // p constraints
			}

		}
	}

	// tau_n print
	if (is_print_tau) print_tau(tau, false);

	// QP setting & solving
	example.hotstart(H, g, A, lb, ub, lbA_in, ubA_in, nWSR_input);

	// solution get
	get_solution(false);

	//output
	for (int i = 0; i < DOF; i++) tau[i] = xOpt[i];

	// tau_c print
	if (is_print_tau) print_tau(tau, true);



}

