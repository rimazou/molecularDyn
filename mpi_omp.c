#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h> 
#include <mpi.h>
#include <omp.h>

int  myrank, numproc;

int main ( int argc, char *argv[] );
void compute ( int np, int nd, double pos[], double vel[], double mass,
double f[], double *pot, double *kin, int myrank, int nbproc );
double cpu_time ( );
double dist ( int nd, double r1[], double r2[], double dr[] );
void initialize ( int np, int nd, double pos[], double vel[], double acc[] );
void r8mat_uniform_ab ( int m, int n, double a, double b, int *seed, double r[] );
void timestamp ( );
void update ( int np, int nd,double loc_pos[], double loc_vel[], double loc_acc[], double pos[], double vel[], double f[],
  double acc[], double mass, double dt, int myrank, int nbproc );
/******************************************************************************/

int main ( int argc, char *argv[] )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MD.

  Discussion:

    MD implements a simple molecular dynamics simulation.

    The velocity Verlet time integration scheme is used.

    The particles interact with a central pair potential.

    This program is based on a FORTRAN90 program by Bill Magro.

  Usage:

    md nd np step_num dt

    where:

    * nd is the spatial dimension (2 or 3);
    * np is the number of particles (500, for instance);
    * step_num is the number of time steps (500, for instance).
    * dt is the time step (0.1 for instance)


*/
{
  double *acc, *loc_acc;
  double ctime;
  double dt;
  double e0;
  double *force, *loc_force;
  double kinetic, loc_kinetic;
  double mass = 1.0;
  int nd;
  int np;
  double *pos, *loc_pos;
  double potential, loc_potential;
  int step;
  int step_num;
  int step_print;
  int step_print_index;
  int step_print_num;
  double *vel, *loc_vel;

/*
  MPI initialization and getting number of processes and the rank
*/
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);

 if (myrank==0){

    timestamp ( );
    printf ( "\n" );
    printf ( "MD\n" );
    printf ( "  C version\n" );
    printf ( "  A molecular dynamics program.\n" );
 }

    if ( 4 < argc )
    {
      // Get the spatial dimension.
      nd = atoi ( argv[1] );

      //  Get the number of particles.
      np = atoi ( argv[2] );
    
      //  Get the number of time steps.
      step_num = atoi ( argv[3] );

      //  Get the time steps.
      dt = atof ( argv[4] );

    }
    else
    {
      printf ( "\n" );
      printf ( "ERROR, the program needs 4 parameters \n" );
      return 0 ;
    }
/*
  Report.
*/
if (myrank==0)
{
  printf ( "\n" );
  printf ( "  ND, the spatial dimension, is %d\n", nd );
  printf ( "  NP, the number of particles in the simulation, is %d\n", np );
  printf ( "  STEP_NUM, the number of time steps, is %d\n", step_num );
  printf ( "  DT, the size of each time step, is %f\n", dt );

}

  
/*
  Allocate memory.
*/
  acc = ( double * ) malloc ( nd * np * sizeof ( double ) );
  force = ( double * ) malloc ( nd * np * sizeof ( double ) );
  pos = ( double * ) malloc ( nd * np * sizeof ( double ) );
  vel = ( double * ) malloc ( nd * np * sizeof ( double ) );
  loc_acc = ( double * ) malloc ( nd * np * sizeof ( double ) );
  loc_force = ( double * ) malloc ( nd * np * sizeof ( double ) );
  loc_pos = ( double * ) malloc ( nd * np * sizeof ( double ) );
  loc_vel = ( double * ) malloc ( nd * np * sizeof ( double ) );


/*
  This is the main time stepping loop:
    Compute forces and energies,
    Update positions, velocities, accelerations.
*/
if (myrank==0)
{
  printf ( "\n" );
  printf ( "  At each step, we report the potential and kinetic energies.\n" );
  printf ( "  The sum of these energies should be a constant.\n" );
  printf ( "  As an accuracy check, we also print the relative error\n" );
  printf ( "  in the total energy.\n" );
  printf ( "\n" );
  printf ( "      Step      Potential       Kinetic        (P+K-E0)/E0\n" );
  printf ( "                Energy P        Energy K       Relative Energy Error\n" );
  printf ( "\n" );

 
}
/*
initialisation de local var
*/
  step_print = 0;
  step_print_index = 0;
  step_print_num = 10;

 for (int o = 0; o < np; o++ )
  {
    for ( int l = 0; l < nd; l++ )
    {
     
    loc_acc[l+ o*nd] = 0.0 ;
    loc_force[l+ o*nd] = 0.0;
    force[l+ o*nd] = 0.0;
    loc_pos[l+o*nd] = 0.0;
    loc_vel[l+o*nd] = 0.0;

    }
  }
  ctime = MPI_Wtime(); 
  
  for ( step = 0; step <= step_num; step++ )
  {
    if ( step == 0 )
    {
         
      initialize ( np, nd, pos, vel, acc );
      
    }
    else
    {
      update ( np, nd,loc_pos, loc_vel, loc_acc, pos, vel, force, acc, mass, dt, myrank, numproc );
      
      MPI_Allreduce(loc_pos,pos, nd*np, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(loc_vel,vel,nd*np , MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(loc_acc,acc, nd*np, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    }


    compute ( np, nd, pos, vel, mass, loc_force, &loc_potential, &loc_kinetic, myrank, numproc);
    MPI_Allreduce(loc_force,force, nd*np, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_potential,&potential,1 , MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_kinetic,&kinetic, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 

    if ( step == 0 )
    {
      e0 = potential + kinetic;
    }
    if ( step == step_print && myrank==0)
    {
  
      printf ( "  %8d  %14f  %14f  %14e\n", step, potential, kinetic,
       ( potential + kinetic - e0 ) / e0 );
      step_print_index = step_print_index + 1;
      step_print = ( step_print_index * step_num ) / step_print_num;
    }

  }
  
/*
  Report timing.
*/
  ctime = MPI_Wtime( ) - ctime;
  printf ( "\n" );
  printf ( "  Elapsed cpu time: %f seconds for processus %d.\n", ctime, myrank);
/*
  Free memory.
*/
  free ( acc );
  free ( force );
  free ( pos );
  free ( vel );
  free ( loc_acc );
  free ( loc_force );
  free ( loc_pos );
  free ( loc_vel );

  MPI_Finalize();
/*
  Terminate.
*/

  printf ( "\n" );
  printf ( "MD\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );
  return 0;
}
/******************************************************************************/

void compute ( int np, int nd, double pos[], double vel[], double mass,
  double f[], double *pot, double *kin, int myrank, int nbproc )

/******************************************************************************/
/*
  Purpose:

    COMPUTE computes the forces and energies.

  Discussion:

    The computation of forces and energies is fully parallel.

    The potential function V(X) is a harmonic well which smoothly
    saturates to a maximum value at PI/2:

      v(x) = ( sin ( min ( x, PI/2 ) ) )^2

    The derivative of the potential is:

      dv(x) = 2.0 * sin ( min ( x, PI/2 ) ) * cos ( min ( x, PI/2 ) )
            = sin ( 2.0 * min ( x, PI/2 ) )



  Parameters:

    Input, int NP, the number of particles.

    Input, int ND, the number of spatial dimensions.

    Input, double POS[ND*NP], the positions.

    Input, double VEL[ND*NP], the velocities.

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
  
  int size = np / nbproc;
  int start = size * myrank;
  int end = size * (myrank+1) ;

 # pragma omp parallel for private(j, k,d, i, d2, rij) shared(np,nd, pos, f, vel) reduction ( + : pe ) reduction (+ : ke) schedule(static)
  for ( k = start ; k < end; k++ )
  {
/*
  Compute the potential energy and forces.
*/

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
//# pragma omp parallel for  private(i) shared(nd, k, vel) reduction ( + : ke )
    for (int i = 0; i < nd; i++ )
    {
      ke = ke + vel[i+k*nd] * vel[i+k*nd];
    }

  }

  ke = ke * 0.5 * mass;

  *pot = pe;
  *kin = ke;


  return;
}
/*******************************************************************************/

double cpu_time ( )

/*******************************************************************************/
/*
  Purpose:

    CPU_TIME reports the total CPU time for a program.


  Parameters:

    Output, double CPU_TIME, the current total elapsed CPU time in second.
*/
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
/******************************************************************************/

double dist ( int nd, double r1[], double r2[], double dr[] )

/******************************************************************************/
/*

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

void initialize ( int np, int nd, double pos[], double vel[], double acc[] )

/******************************************************************************/
/*
  Purpose:

    INITIALIZE initializes the positions, velocities, and accelerations.


  Parameters:

    Input, int NP, the number of particles.

    Input, int ND, the number of spatial dimensions.

    Output, double POS[ND*NP], the positions.

    Output, double VEL[ND*NP], the velocities.

    Output, double ACC[ND*NP], the accelerations.
*/
{
  int i;
  int j;
  int seed;
/*
  Set positions.
*/
  seed = 123456789;
  r8mat_uniform_ab ( nd, np, 0.0, 10.0, &seed, pos );
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

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has
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
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
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

void timestamp ( )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
/******************************************************************************/

void update ( int np, int nd,double loc_pos[], double loc_vel[], double loc_acc[], double pos[], double vel[], double f[],
  double acc[], double mass, double dt, int myrank, int nbproc )

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


  Parameters:

    Input, int NP, the number of particles.

    Input, int ND, the number of spatial dimensions.

    Input/output, double POS[ND*NP], the positions.

    Input/output, double VEL[ND*NP], the velocities.

    Input, double F[ND*NP], the forces.

    Input/output, double ACC[ND*NP], the accelerations.

    Input, double MASS, the mass.

    Input, double DT, the time step.
*/
{
  int i;
  int j;
  double rmass;

  rmass = 1.0 / mass;
 int size = np / nbproc;
  int start = size * myrank;
  int end = size * (myrank+1) ;

 # pragma omp parallel for private(j,i) shared(np,nd, pos,f, vel, rmass, acc,dt,loc_pos,loc_vel,loc_acc)  schedule(static)
  for ( j = start; j < end; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      loc_pos[i+j*nd] = pos[i+j*nd] + vel[i+j*nd] * dt + 0.5 * acc[i+j*nd] * dt * dt;
      loc_vel[i+j*nd] = vel[i+j*nd] + 0.5 * dt * ( f[i+j*nd] * rmass + acc[i+j*nd] );
      loc_acc[i+j*nd] = f[i+j*nd] * rmass;
    }
  }
  return;
}