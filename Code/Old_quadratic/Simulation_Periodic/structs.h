#ifndef STRUCTS_H_
#define STRUCTS_H_

/*-------------------------------------------------------------------------------------------------------*/
// define a struct for retrieving parameters from a given file
struct parameters {
	// file to be written to
	char 	outfile[1024];

	// system properties
	double 	Gamma_A;
	double 	m;
	double 	shear_A;
	double 	D_rat;

	// simulation and computer system properties
	double 	timestep;
	int		max_timesteps;
	int		write_step;
};


/*-------------------------------------------------------------------------------------------------------*/
// define a struct for handing parameters to the simulation
struct sim_struct {

	// variable parameters set by the user
	double	Gamma_A;
	double 	m;
	double 	shear;
	double 	D_rat;

	// variables relevant to the simulation
	double 	timestep;
	int		max_timesteps;
	int		write_step;
	int		sim_number;

	// file to write to
	char	outfile[1024];
};


#endif /* STRUCTS_H_ */
