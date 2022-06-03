# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <math.h>
# include <omp.h>

int main ( int argc, char *argv[] );
void compute ( int np, int nd, double pos[], double vel[],
  double mass, double f[], double *pot, double *kin );
double dist ( int nd, double r1[], double r2[], double dr[] );
void initialize ( int np, int nd, double box[], int *seed, double pos[],
  double vel[], double acc[] );
void r8mat_uniform_ab ( int m, int n, double a, double b, int *seed, double r[] );
void timestamp ( void );
void update ( int np, int nd, double pos[], double vel[], double f[],
  double acc[], double mass, double dt );

/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:
    MAIN is the main program for MD_OPENMP.
  Discussion:
    MD implements a simple molecular dynamics simulation.
    The program uses Open MP directives to allow parallel computation.
    The velocity Verlet time integration scheme is used.
    The particles interact with a central pair potential.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    30 July 2009
  Author:
    Original FORTRAN77 version by Bill Magro.
    C version by John Burkardt.
  Parameters:
    None
*/
{
  double *acc;
  double *box;
  double dt;
  double e0;
  double *force;
  int nbrThreads;
  int i;
  int id;
  double kinetic;
  double mass = 1.0;
  int nd;
  int np;
  double *pos;
  double potential;
  int proc_num;
  int seed = 123456789;
  int step;
  int step_num;
  int step_print;
  int step_print_index;
  int step_print_num;
  double *vel;
  double wtime;

  timestamp ( );

  printf ( "\nC/OpenMP\n\n" );
  printf ( "A molecular dynamics program.\n\n" );
  
  proc_num = omp_get_num_procs ( );

	/*
	 Get the spatial dimension.
	*/
	 if ( 1 < argc )
	 {
	 nd = atoi ( argv[1] );
	 }
	 else
	 {
	 printf ( "\n" );
	 printf ( " Enter ND, the spatial dimension : " );
	 scanf ( "%d", &nd );
	 printf ( "\n" );
	 }
	//
	// Get the number of particles.
	//
	 if ( 2 < argc )
	 {
	 np = atoi ( argv[2] );
	 }
	 else
	 {
	 printf ( "\n" );
	 printf ( " Enter NP, the number of particles  : " );
	 scanf ( "%d", &np );
	 printf ( "\n" );
	 }
	//
	// Get the number of time steps.
	//
	 if ( 3 < argc )
	 {
	 step_num = atoi ( argv[3] );
	 }
	 else
	 {
	 printf ( "\n" );
	 printf ( " Enter step_num, the number of time steps : " );
	 scanf ( "%d", &step_num );
	 printf ( "\n" );
	 }
	//
	// Get the time steps.
	//
	 if ( 4 < argc )
	 {
	 dt = atof ( argv[4] );
	 }
	 else
	 {
	 printf ( "\n" );
	 printf ( " Enter DT, the size of the time step (0.1, for instance) : " );
	 scanf ( "%lf", &dt );
	 printf ( "\n" );
	 }
	 //Get number of threads
	 if ( 5 < argc )
	 {
	 dt = atof ( argv[5] );
	 }
	 else
	 {
	 printf ( "\n" );
	 printf ( " Enter nbrThreads, the number of threads to use : ");
	 scanf ( "%d", &nbrThreads );
	 printf ( "\n" );
	 }
	 
	/*
	 Report.
	*/
	 printf ( "\n" );
	 printf ( " ND, the spatial dimension, is %d\n", nd );
	 printf ( " NP, the number of particles in the simulation, is %d\n", np );
	 printf ( " STEP_NUM, the number of time steps, is %d\n", step_num );
	 printf ( " DT, the size of each time step, is %f\n", dt );
	/*
	 Allocate memory.
	*/


  acc = ( double * ) malloc ( nd * np * sizeof ( double ) );
  box = ( double * ) malloc ( nd * sizeof ( double ) );
  force = ( double * ) malloc ( nd * np * sizeof ( double ) );
  pos = ( double * ) malloc ( nd * np * sizeof ( double ) );
  vel = ( double * ) malloc ( nd * np * sizeof ( double ) );
/*
  Set the dimensions of the box.
*/
  for ( i = 0; i < nd; i++ )
  {
    box[i] = 10.0;
  }

  printf ( "\n" );
  printf ( "  Initializing positions, velocities, and accelerations.\n" );
/*
  Set initial positions, velocities, and accelerations.
*/
  initialize ( np, nd, box, &seed, pos, vel, acc );
/*
  Compute the forces and energies.
*/
  printf ( "\n" );
  printf ( "  Computing initial forces and energies.\n" );

  compute ( np, nd, pos, vel, mass, force, &potential, &kinetic );

  e0 = potential + kinetic;
/*
  This is the main time stepping loop:
    Compute forces and energies,
    Update positions, velocities, accelerations.
*/
  printf ( "\n" );
  printf ( "  At each step, we report the potential and kinetic energies.\n" );
  printf ( "  The sum of these energies should be a constant.\n" );
  printf ( "  As an accuracy check, we also print the relative error\n" );
  printf ( "  in the total energy.\n" );
  printf ( "\n" );
  printf ( "      Step      Potential       Kinetic        (P+K-E0)/E0\n" );
  printf ( "                Energy P        Energy K       Relative Energy Error\n" );
  printf ( "\n" );

  step_print = 0;
  step_print_index = 0;
  step_print_num = 10;

  step = 0;
  printf ( "  %8d  %14f  %14f  %14e\n",
    step, potential, kinetic, ( potential + kinetic - e0 ) / e0 );
  step_print_index = step_print_index + 1;
  step_print = ( step_print_index * step_num ) / step_print_num;
	
omp_set_num_threads(nbrThreads);
  wtime = omp_get_wtime ( );

  for ( step = 1; step <= step_num; step++ )
  {
    compute ( np, nd, pos, vel, mass, force, &potential, &kinetic );

    if ( step == step_print )
    {
      printf ( "  %8d  %14f  %14f  %14e\n", step, potential, kinetic,
       ( potential + kinetic - e0 ) / e0 );
      step_print_index = step_print_index + 1;
      step_print = ( step_print_index * step_num ) / step_print_num;
    }
    update ( np, nd, pos, vel, force, acc, mass, dt );
  }
  wtime = omp_get_wtime ( ) - wtime;

  printf ( "\n" );
  printf ( "  Elapsed time for main computation:\n" );
  printf ( "  %f seconds.\n", wtime );

  free ( acc );
  free ( box );
  free ( force );
  free ( pos );
  free ( vel );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "  Normal end of execution.\n" );

  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void compute ( int np, int nd, double pos[], double vel[],
  double mass, double f[], double *pot, double *kin )

/******************************************************************************/
/*
  Purpose:
    COMPUTE computes the forces and energies.
  Discussion:
    The computation of forces and energies is fully parallel.
    The potential function V(X) is a harmonic well which smoothly
    saturates to a maximum value at PI/2:
      v(x) = ( sin ( min ( x, PI2 ) ) )**2
    The derivative of the potential is:
      dv(x) = 2.0 * sin ( min ( x, PI2 ) ) * cos ( min ( x, PI2 ) )
            = sin ( 2.0 * min ( x, PI2 ) )
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 November 2007
  Author:
    Original FORTRAN77 version by Bill Magro.
    C version by John Burkardt.
  Parameters:
    Input, int NP, the number of particles.
    Input, int ND, the number of spatial dimensions.
    Input, double POS[ND*NP], the position of each particle.
    Input, double VEL[ND*NP], the velocity of each particle.
    Input, double MASS, the mass of each particle.
    Output, double F[ND*NP], the forces.
    Output, double *POT, the total potential energy.
    Output, double *KIN, the total kinetic energy.
*/
{
  double d;
  double d2;
  int i;
  int j;
  int k;
  double ke;
  double pe;
  double PI2 = 3.141592653589793 / 2.0;
  double rij[3];

  pe = 0.0;
  ke = 0.0;

# pragma omp parallel \
  shared ( f, nd, np, pos, vel ) \
  private ( i, j, k, rij, d, d2 )


  for ( k = 0; k < np; k++ )
  {
/*
  Compute the potential energy and forces.
*/
# pragma omp for reduction ( + : ke )
    for ( i = 0; i < nd; i++ )
    {
      f[i+k*nd] = 0.0;
    }

    for ( j = 0; j < np; j++ )
    {
      if ( k != j )
      {
        d = dist ( nd, pos+k*nd, pos+j*nd, rij );
/*
  Attribute half of the potential energy to particle J.
*/
        if ( d < PI2 )
        {
          d2 = d;
        }
        else
        {
          d2 = PI2;
        }

        pe = pe + 0.5 * pow ( sin ( d2 ), 2 );

        for ( i = 0; i < nd; i++ )
        {
          f[i+k*nd] = f[i+k*nd] - rij[i] * sin ( 2.0 * d2 ) / d;
        }
      }
    }
/*
  Compute the kinetic energy.
*/
    for ( i = 0; i < nd; i++ )
    {
      ke = ke + vel[i+k*nd] * vel[i+k*nd];
    }
  }

  ke = ke * 0.5 * mass;

  *pot = pe;
  *kin = ke;

  return;
}
/******************************************************************************/

double dist ( int nd, double r1[], double r2[], double dr[] )

/******************************************************************************/
/*
  Purpose:
    DIST computes the displacement (and its norm) between two particles.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 November 2007
  Author:
    Original FORTRAN77 version by Bill Magro.
    C version by John Burkardt.
  Parameters:
    Input, int ND, the number of spatial dimensions.
    Input, double R1[ND], R2[ND], the positions of the particles.
    Output, double DR[ND], the displacement vector.
    Output, double D, the Euclidean norm of the displacement.
*/
{
  double d;
  int i;

  d = 0.0;
  for ( i = 0; i < nd; i++ )
  {
    dr[i] = r1[i] - r2[i];
    d = d + dr[i] * dr[i];
  }
  d = sqrt ( d );

  return d;
}
/******************************************************************************/

void initialize ( int np, int nd, double box[], int *seed, double pos[],
  double vel[], double acc[] )

/******************************************************************************/
/*
  Purpose:
    INITIALIZE initializes the positions, velocities, and accelerations.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    21 November 2007
  Author:
    Original FORTRAN77 version by Bill Magro.
    C version by John Burkardt.
  Parameters:
    Input, int NP, the number of particles.
    Input, int ND, the number of spatial dimensions.
    Input, double BOX[ND], specifies the maximum position
    of particles in each dimension.
    Input, int *SEED, a seed for the random number generator.
    Output, double POS[ND*NP], the position of each particle.
    Output, double VEL[ND*NP], the velocity of each particle.
    Output, double ACC[ND*NP], the acceleration of each particle.
*/
{
	 int i;
	 int j;
	/*
	 Set positions.
	*/
	 r8mat_uniform_ab ( nd, np, 0.0, 10.0, seed, pos );
	/*
	 Set velocities.
	*/
	 for ( j = 0; j < np; j++ )
	 {
	 for ( i = 0; i < nd; i++ )
	 {
	 vel[i+j*nd] = 0.0;
	 }
	 }
	/*
	 Set accelerations.
	*/
	 for ( j = 0; j < np; j++ )
	 {
	 for ( i = 0; i < nd; i++ )
	 {
	 acc[i+j*nd] = 0.0;
	 }
	 }
	
	 return;
	}
/******************************************************************************/

void r8mat_uniform_ab ( int m, int n, double a, double b, int *seed, double r[] )
	
	/******************************************************************************/
	/*
	 Purpose:
	
	 R8MAT_UNIFORM_AB returns a scaled pseudorandom R8MAT.
	
	 Discussion:
	
	 This routine implements the recursion
	
	 seed = 16807 * seed mod ( 2^31 - 1 )
	 unif = seed / ( 2^31 - 1 )
	
	 The integer arithmetic never requires more than 32 bits,
	 including a sign bit.
	
	 Parameters:
	
	 Input, int M, N, the number of rows and columns.
	
	 Input, double A, B, the limits of the pseudorandom values.
	
	 Input/output, int *SEED, the "seed" value. Normally, this
	 value should not be 0. On output, SEED has
	 been updated.
	
	 Output, double R[M*N], a matrix of pseudorandom values.
	*/
	{
	 int i;
	 const int i4_huge = 2147483647;
	 int j;
	 int k;
	
	 if ( *seed == 0 )
	 {
	 fprintf ( stderr, "\n" );
	 fprintf ( stderr, "R8MAT_UNIFORM_AB - Fatal error!\n" );
	 fprintf ( stderr, " Input value of SEED = 0.\n" );
	 exit ( 1 );
	 }
	
	 for ( j = 0; j < n; j++ )
	 {
	 for ( i = 0; i < m; i++ )
	 {
	 k = *seed / 127773;
	
	 *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;
	
	 if ( *seed < 0 )
	 {
	 *seed = *seed + i4_huge;
	 }
	 r[i+j*m] = a + ( b - a ) * ( double ) ( *seed ) * 4.656612875E-10;
	 }
	 }
	
	 return;
	}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:
    TIMESTAMP prints the current YMDHMS date as a time stamp.
  Example:
    31 May 2001 09:45:54 AM
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 September 2003
  Author:
    John Burkardt
  Parameters:
    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
/******************************************************************************/

void update ( int np, int nd, double pos[], double vel[], double f[],
  double acc[], double mass, double dt )

/******************************************************************************/
/*
  Purpose:
    UPDATE updates positions, velocities and accelerations.
  Discussion:
    The time integration is fully parallel.
    A velocity Verlet algorithm is used for the updating.
    x(t+dt) = x(t) + v(t) * dt + 0.5 * a(t) * dt * dt
    v(t+dt) = v(t) + 0.5 * ( a(t) + a(t+dt) ) * dt
    a(t+dt) = f(t) / m
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    17 April 2009
  Author:
    Original FORTRAN77 version by Bill Magro.
    C version by John Burkardt.
  Parameters:
    Input, int NP, the number of particles.
    Input, int ND, the number of spatial dimensions.
    Input/output, double POS[ND*NP], the position of each particle.
    Input/output, double VEL[ND*NP], the velocity of each particle.
    Input, double F[ND*NP], the force on each particle.
    Input/output, double ACC[ND*NP], the acceleration of each particle.
    Input, double MASS, the mass of each particle.
    Input, double DT, the time step.
*/
{
  int i;
  int j;
  double rmass;

  rmass = 1.0 / mass;

  for ( j = 0; j < np; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      pos[i+j*nd] = pos[i+j*nd] + vel[i+j*nd] * dt + 0.5 * acc[i+j*nd] * dt * dt;
      vel[i+j*nd] = vel[i+j*nd] + 0.5 * dt * ( f[i+j*nd] * rmass + acc[i+j*nd] );
      acc[i+j*nd] = f[i+j*nd] * rmass;
    }
  }

  return;
}
