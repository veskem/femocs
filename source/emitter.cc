

#include "../include/currents_and_heating.h"

/*
int main() {
	Laplace laplace_problem;
	laplace_problem.run();

	return 0;
}
*/


int main() {
	try {
		using namespace dealii;
		using namespace CurrentsAndHeating2d;

		deallog.depth_console(1);

		CurrentsAndHeating problem;
		problem.run();

	} catch (std::exception &exc) {
		std::cerr << std::endl << std::endl
				<< "----------------------------------------------------"
				<< std::endl;
		std::cerr << "Exception on processing: " << std::endl
				<< exc.what() << std::endl
				<< "Aborting!" << std::endl
				<< "----------------------------------------------------"
				<< std::endl;
		return 1;
	} catch (...) {
		std::cerr << std::endl << std::endl
				<< "----------------------------------------------------"
				<< std::endl;
		std::cerr << "Unknown exception!" << std::endl
				<< "Aborting!" << std::endl
				<< "----------------------------------------------------"
				<< std::endl;
		return 1;
	}

	return 0;
}
