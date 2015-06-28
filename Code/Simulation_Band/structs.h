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

	double 	gamma_shear;
	double 	D_rat;

	// simulation and computer system properties
	double 	timestep;
	int		tau;
	int		write_step;
	int		sim_number;
};

/*-------------------------------------------------------------------------------------------------------*/
// define a struct for passing all relevant system parameters to file
struct attributes {
	// system parameters
	int 	Num;
	double 	Gamma_A;
	double	m;
	double 	X;
	double 	gamma_shear;
	double 	D_rat;

	// box parameters
	double 	L_x;
	double 	L_y_attr;

	// simulation parameters
	int 	max_writeouts;
	int 	tau;
	int 	write_step;
	double 	timestep;
};


#endif /* STRUCTS_H_ */
