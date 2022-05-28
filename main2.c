#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <pthread.h>

#define NUM_THREADS 4
#define n 100
pthread_t tid[NUM_THREADS];

double *acc;
double *force;
double *pos;
double *vel;

double e0;
double kinetic;
double potential;

//double r[200];
double a= 0.0;
double b=10.0;
// r8mat_uniform_ab function arguments
typedef struct thread_args
{
	int my_first;
	int my_last;
        double seed;
        double r[200];
} thread_args;
typedef struct thread_args thread_args;

// update function arguments
typedef struct update_struct
{
  int np_start;
  int np_end;
  int nd;
  double mass;
  double dt;
} UpdateStruct;
// compute function arguments
typedef struct compute_struct
{
  int np;
  int np_start;
  int np_end;
  int nd;
  double mass;
  double *pot;
  double *kin;
} ComputeStruct;
// update function return type
typedef struct pot_kin
{
  double pot;
  double kin;
} PotKin;

int main(int argc, char *argv[]);

void *compute(void *data);

double cpu_time();

double dist(int nd, double r1[], double r2[], double dr[]);

void initialize(int np, int nd, double pos[], double vel[], double acc[]);

void *r8mat_uniform_ab(void *input);

void timestamp();

void *update(void *data);

/******************************************************************************/

int main(int argc, char *argv[])

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

  double ctime;
  double dt;
  double mass = 1.0;
  int nd;
  int np;
  int step;
  int step_num;
  int step_print;
  int step_print_index;
  int step_print_num;

  timestamp();
  printf("\n");
  printf("MD\n");
  printf("  C version\n");
  printf("  A molecular dynamics program.\n");
  /*
    Get the spatial dimension.
  */
  if (1 < argc)
  {
    nd = atoi(argv[1]);
  }
  else
  {
    printf("\n");
    printf("  Enter ND, the spatial dimension (2 or 3).\n");
    scanf("%d", &nd);
  }
  //
  //  Get the number of particles.
  //
  if (2 < argc)
  {
    np = atoi(argv[2]);
  }
  else
  {
    printf("\n");
    printf("  Enter NP, the number of particles (500, for instance).\n");
    scanf("%d", &np);
  }
  //
  //  Get the number of time steps.
  //
  if (3 < argc)
  {
    step_num = atoi(argv[3]);
  }
  else
  {
    printf("\n");
    printf("  Enter ND, the number of time steps (500 or 1000, for instance).\n");
    scanf("%d", &step_num);
  }
  //
  //  Get the time steps.
  //
  if (4 < argc)
  {
    dt = atof(argv[4]);
  }
  else
  {
    printf("\n");
    printf("  Enter DT, the size of the time step (0.1, for instance).\n");
    scanf("%lf", &dt);
  }
  /*
    Report.
  */
  printf("\n");
  printf("  ND, the spatial dimension, is %d\n", nd);
  printf("  NP, the number of particles in the simulation, is %d\n", np);
  printf("  STEP_NUM, the number of time steps, is %d\n", step_num);
  printf("  DT, the size of each time step, is %f\n", dt);
  /*
    Allocate memory.
  */
  acc = (double *)malloc(nd * np * sizeof(double));
  force = (double *)malloc(nd * np * sizeof(double));
  pos = (double *)malloc(nd * np * sizeof(double));
  vel = (double *)malloc(nd * np * sizeof(double));
  /*
    This is the main time stepping loop:
      Compute forces and energies,
      Update positions, velocities, accelerations.
  */
  printf("\n");
  printf("  At each step, we report the potential and kinetic energies.\n");
  printf("  The sum of these energies should be a constant.\n");
  printf("  As an accuracy check, we also print the relative error\n");
  printf("  in the total energy.\n");
  printf("\n");
  printf("      Step      Potential       \n");
  printf("                Energy P        \n");
  printf("\n");

  step_print = 0;
  step_print_index = 0;
  step_print_num = 10;

  ctime = cpu_time();

  for (step = 0; step <= step_num; step++)
  {
    if (step == 0)
    {
      initialize(np, nd, pos, vel, acc);
    }
    else
    {

      // Start update function
      for (int i = 0; i < NUM_THREADS; i++)
      {
        UpdateStruct *p = malloc(sizeof(UpdateStruct));
        p->dt = dt;
        p->mass = mass;
        p->nd = nd;
        int taux;
        if (np % NUM_THREADS == 0)
        {
          taux = (np / NUM_THREADS);
        }
        else
        {
          taux = (np / NUM_THREADS) + np % NUM_THREADS;
        }
        p->np_start = taux * i;
        p->np_end = (i + 1) * taux;

        if (p->np_end > np)
          p->np_end = np;

        int rc = pthread_create(&tid[i], NULL, update,
                                p);
        if (rc)
        {
          printf("ERROR; return code from pthread_create() is %d\n", rc);
          exit(-1);
        }
      }

      for (int t = 0; t < NUM_THREADS; t++)
      {
        pthread_join(tid[t], NULL);
      }
      // End update function
    }

    // compute function
    for (int i = 0; i < NUM_THREADS; i++)
    {
      ComputeStruct *c = malloc(sizeof(ComputeStruct));
      int taux = 0;
      if (np % NUM_THREADS == 0)
      {
        taux = (np / NUM_THREADS);
      }
      else
      {
        taux = (np / NUM_THREADS) + np % NUM_THREADS;
      }
      c->np_start = taux * i;
      c->np_end = (i + 1) * taux;
      if (c->np_end > np)
        c->np_end = np;

      c->nd = nd;
      c->mass = mass;
      c->pot = &potential;
      c->kin = &kinetic;
      c->np = np;

      int rc = pthread_create(&tid[i], NULL, compute,
                              c);
      if (rc)
      {
        printf("ERROR; return code from pthread_create() is %d\n", rc);
        exit(-1);
      }
    }
    PotKin *potkin[NUM_THREADS];
    for (int t = 0; t < NUM_THREADS; t++)
    {
      pthread_join(tid[t], (void **)&potkin[t]);
    }
    potential = 0.0;
    kinetic = 0.0;
    for (int k = 0; k < NUM_THREADS; k++)
    {
      potential += (double)(potkin[k])->pot;
      kinetic += (double)(potkin[k])->kin;
    }
    // end compute

    if (step == 0)
    {
      e0 = potential + kinetic;
    
    }

    if (step == step_print)
    {
      //printf("e0= %f, potential= %f,kinetic=%f ", e0, potential, kinetic);
      //printf("  %8d  %14f  %14f  %14e\n", step, potential, kinetic, (potential + kinetic - e0) / e0);
      printf("  %8d  %14f  \n", step, potential);
      step_print_index = step_print_index + 1;
      step_print = (step_print_index * step_num) / step_print_num;
    }
  }
  /*
    Report timing.
  */
  ctime = cpu_time() - ctime;
  printf("\n");
  printf("  Elapsed cpu time: %f seconds.\n", ctime);
  /*
    Free memory.  */
  free(acc);
  free(force);
  free(pos);
  free(vel);
  /*
    Terminate.
  */
  printf("\n");
  printf("MD\n");
  printf("  Normal end of execution.\n");
  printf("\n");
  timestamp();

  return 0;
}

/******************************************************************************/

void *compute(void *data)
// int np, int nd, double pos[], double vel[], double mass, double f[], double *pot, double *kin
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

  // get the arguments
  ComputeStruct *args = (ComputeStruct *)data;
  int np_start = args->np_start;
  int np_end = args->np_end;
  double mass = args->mass;
  double *pot = args->pot;
  double *kin = args->kin;
  int nd = args->nd;
  int np = args->np;

  pe = 0.0;
  ke = 0.0;

  for (k = np_start; k < np_end; k++)
  {
    /*
      Compute the potential energy and forces.
    */
    for (i = 0; i < nd; i++)
    {
      force[i + k * nd] = 0.0;
    }

    for (j = 0; j < np; j++)
    {
      if (k != j)
      {
        d = dist(nd, pos + k * nd, pos + j * nd, rij);
        /*
          Attribute half of the potential energy to particle J.
        */
        if (d < PI2)
        {
          d2 = d;
        }
        else
        {
          d2 = PI2;
        }

        pe = pe + 0.5 * pow(sin(d2), 2);

        for (i = 0; i < nd; i++)
        {
          force[i + k * nd] = force[i + k * nd] - rij[i] * sin(2.0 * d2) / d;
        }
      }
    }
    /*
      Compute the kinetic energy.
    */
    for (i = 0; i < nd; i++)
    {
      ke = ke + vel[i + k * nd] * vel[i + k * nd];
      
    }
  }

  ke = ke * 0.5 * mass;

  //    pthread_mutex_lock(&mutex_compute);
  //    *pot = pe;
  //    *kin = ke;
  //    pthread_mutex_unlock(&mutex_compute);

  PotKin *potKin = malloc(sizeof(PotKin));
  potKin->pot = pe;
  potKin->kin = ke;
  pthread_exit(potKin);
}

/*******************************************************************************/

double cpu_time()

/*******************************************************************************/
/*
  Purpose:

    CPU_TIME reports the total CPU time for a program.


  Parameters:

    Output, double CPU_TIME, the current total elapsed CPU time in second.
*/
{
  double value;

  value = (double)clock() / (double)CLOCKS_PER_SEC;

  return value;
}

/******************************************************************************/

double dist(int nd, double r1[], double r2[], double dr[])

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
  for (i = 0; i < nd; i++)
  {
    dr[i] = r1[i] - r2[i];
    d = d + dr[i] * dr[i];
  }
  d = sqrt(d);

  return d;
}

/******************************************************************************/

void initialize(int np, int nd, double pos[], double vel[], double acc[])

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
  //seed = 123456789;
 pthread_t threads[NUM_THREADS];
 int rs;
  
 
   	

 	for (i = 0; i < NUM_THREADS; i++)
	{
		thread_args *args = malloc(sizeof(thread_args));
		args->my_first = i * (n / NUM_THREADS);
		args->my_last = (i + 1) * (n / NUM_THREADS) - 1;
        	args->seed = 123456789;
                args ->r[200] = pos[200];
		rs = pthread_create(&threads[i], NULL, r8mat_uniform_ab, (void*)&args);
       
		if (rs)
		{
			printf("ERROR; return code from pthread_create() is %d\n", rs);
			exit(-1);
		}
     
	}
	for (i = 0; i < NUM_THREADS; i++)
	{
		pthread_join(threads[i], NULL);
       
	}
 // r8mat_uniform_ab(nd, np, 0.0, 10.0, &seed, pos);
  /*
    Set velocities.
  */
  for (j = 0; j < np; j++)
  {
    for (i = 0; i < nd; i++)
    {
      vel[i + j * nd] = 0.0;
      
    }
  }
  /*
    Set accelerations.
  */
  for (j = 0; j < np; j++)
  {
    for (i = 0; i < nd; i++)
    {
      acc[i + j * nd] = 0.0;
    }
  }

  return;
}

/******************************************************************************/

void *r8mat_uniform_ab(void *input)
{
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

 thread_args *args = (thread_args *)input;
	int my_first = args->my_first;
	int my_last = args->my_last;
    double seed = args->seed;
    double r[200];
   r[200]= args->r[200];
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;

  if ( seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_UNIFORM_AB - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }  

  for ( j = 0; j < n; j++)
  {   
    for (i = my_first; i < my_last; ++i)
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;
       
      if ( seed < 0 )
      {
        seed = seed + i4_huge;
      }
    
      r[i+j*my_last] = a + ( b - a ) * ( double ) ( seed ) * 4.656612875E-10;
    }
  }
}

/******************************************************************************/

void timestamp()

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
#define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time(NULL);
  tm = localtime(&now);

  strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

  printf("%s\n", time_buffer);

  return;
#undef TIME_SIZE
}

/******************************************************************************/

void *update(void *data)
// int np, int nd, double mass, double dt
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

  // function arguments
  UpdateStruct *args = (UpdateStruct *)data;

  double mass = args->mass;
  rmass = 1.0 / mass;

  int np_start = args->np_start;
  int np_end = args->np_end;

  int nd = args->nd;
  double dt = args->dt;

  for (j = np_start; j < np_end; j++)
  {
    for (i = 0; i < nd; i++)
    {
      pos[i + j * nd] = pos[i + j * nd] + vel[i + j * nd] * dt + 0.5 * acc[i + j * nd] * dt * dt;
      vel[i + j * nd] = vel[i + j * nd] + 0.5 * dt * (force[i + j * nd] * rmass + acc[i + j * nd]);
      acc[i + j * nd] = force[i + j * nd] * rmass;
    }
  }

  pthread_exit(NULL);
}
