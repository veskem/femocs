
#include <iostream>
#include <cstdlib>


#include "utility.h"
#include "currents_and_heating.h"
#include "physical_quantities.h"


int main() {

	PhysicalQuantities physical_quantities;

	if (!physical_quantities.load_emission_data("../data/gtf_grid_1000x1000.dat")) return EXIT_FAILURE;
	if (!physical_quantities.load_resistivity_data("../data/cu_res_mod.dat")) return EXIT_FAILURE;

	//physical_quantities.output_to_files();

	try {
		using namespace dealii;

		deallog.depth_console(1);

		CurrentsAndHeating problem(physical_quantities);
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
		return EXIT_FAILURE;
	} catch (...) {
		std::cerr << std::endl << std::endl
				<< "----------------------------------------------------"
				<< std::endl;
		std::cerr << "Unknown exception!" << std::endl
				<< "Aborting!" << std::endl
				<< "----------------------------------------------------"
				<< std::endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;

}
