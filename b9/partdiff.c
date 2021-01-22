/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                 TU München - Institut für Informatik                   **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff.c                                                  **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauß-Seidel and    **/
/**            Jacobi method.                                              **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* ************************************************************************ */
/* Include standard header file.                                            */
/* ************************************************************************ */
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>

#include <mpi.h>

#include "partdiff.h"

struct calculation_arguments
{
    uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
    uint64_t  num_matrices;   /* number of matrices                             */
    double    h;              /* length of a space between two lines            */
    double    ***Matrix;      /* index matrix used for addressing M             */
    double    *M;             /* two matrices with real values                  */
    uint64_t      numberOfRows;   /* Number of rows to calculate for each rank  */
    uint64_t       start;         /* two matrices with real values                  */
    uint64_t       end;           /* two matrices with real values                  */
};

struct calculation_results
{
    uint64_t  m;
    uint64_t  stat_iteration; /* number of current iteration                    */
    double    stat_precision; /* actual precision of all slaves in iteration    */
};

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time; /* time when program started                      */
struct timeval comp_time;  /* time when calculation completed                */

int nprocs, rank;


/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static
void
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
{
    arguments->N = (options->interlines * 8) + 9 - 1;
    arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
    arguments->h = 1.0 / arguments->N;

        int size = (arguments->N+1) / nprocs;
        int rest = (arguments->N+1) % nprocs;

        if(rank < rest)
        {
            arguments->start = (rank * (size + 1));
            arguments->end = arguments->start + size;
        }
        else
        {
            arguments->start = (rank * size) + 1;
            arguments->end = arguments->start + size - 1;
        }
        if(rank == 0)
        {
            arguments->start = 1;
        }
        if(rank == nprocs -1)
        {
            arguments->end = arguments->N-1;
        }
        arguments->numberOfRows = arguments->end - arguments->start + 1;
    results->m = 0;
    results->stat_iteration = 0;
    results->stat_precision = 0;
}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static
void
freeMatrices (struct calculation_arguments* arguments)
{
    uint64_t i;

    for (i = 0; i < arguments->num_matrices; i++)
    {
        free(arguments->Matrix[i]);
    }

    free(arguments->Matrix);
    free(arguments->M);
}

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static
void*
allocateMemory (size_t size)
{
    void *p;

    if ((p = malloc(size)) == NULL)
    {
        printf("Speicherprobleme! (%" PRIu64 " Bytes angefordert)\n", size);
        exit(1);
    }

    return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static
void
allocateMatrices (struct calculation_arguments* arguments, struct options const* options)
{
    uint64_t i, j;
    uint64_t numberOfMatrixRows = arguments->numberOfRows+2;
    uint64_t const N = arguments->N;

    arguments->M = allocateMemory(arguments->num_matrices * (numberOfMatrixRows) * (N+1) * sizeof(double));
    arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

    for (i = 0; i < arguments->num_matrices; i++)
    {
        arguments->Matrix[i] = allocateMemory((numberOfMatrixRows) * sizeof(double*));

        for (j = 0; j < numberOfMatrixRows; j++)
        {
            arguments->Matrix[i][j] = arguments->M + (i * (numberOfMatrixRows) * (N + 1)) + (j * (N + 1));
        }
    }
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options)
{
    uint64_t g, i, j; /* local variables for loops */

    uint64_t const N = arguments->N;
    double const h = arguments->h;
    double*** Matrix = arguments->Matrix;
    int globalRow;

    uint64_t numberOfMatrixRows = arguments->numberOfRows+2;

    /* initialize matrix/matrices with zeros */
    for (g = 0; g < arguments->num_matrices; g++)
    {
        for (i = 0; i < numberOfMatrixRows; i++)
        {
            for (j = 0; j <= N; j++)
            {
                Matrix[g][i][j] = 0.0;
            }
        }

        if(options->inf_func == FUNC_F0)
        {
            for (i = 1; i < numberOfMatrixRows; i++)
            {
                globalRow = (i + arguments->start - 1);
                Matrix[g][i][0] = 1.0 - (h * globalRow);
                Matrix[g][i][N] = h * globalRow;
            }

            if (rank == 0)
            {
                for (i = 0; i <= N; i++)
                {
                    Matrix[g][0][i] = 1.0 - (h * i);
                }

                Matrix[g][0][N] = 0.0;
            }

            if (rank == nprocs - 1)
            {
                for (i = 1; i <= N; i++)
                {
                    Matrix[g][numberOfMatrixRows - 1][i] = h * i;
                }
                Matrix[g][numberOfMatrixRows - 1][0] = 0.0;
            }
        }
    }
}


/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculateGS (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
    int i, j;           /* local variables for loops */
    int m1, m2;         /* used as indices for old and new matrices */
    double star;        /* four times center value minus 4 neigh.b values */
    double residuum;    /* residuum of current iteration */
    double maxResiduum; /* maximum residuum value of a slave in iteration */

    int flag = 0; // flag for stoping using iprobe with precision
    MPI_Status status; // status for iprobe
    int sentStopAlready = 0; // last rank send stop message only once

    //MPI_Request request; // request for isend and i probe (abbruch prec)

    int const N = arguments->N;
    double const h = arguments->h;

    double pih = 0.0;
    double fpisin = 0.0;

    int term_iteration = options->term_iteration;

    /* initialize m1 and m2 depending on algorithm */
    if (options->method == METH_JACOBI)
    {
        m1 = 0;
        m2 = 1;
    }
    else
    {
        m1 = 0;
        m2 = 0;
    }

    if (options->inf_func == FUNC_FPISIN)
    {
        pih = PI * h;
        fpisin = 0.25 * TWO_PI_SQUARE * h * h;
    }

    while (term_iteration > 0)
    {


        double** Matrix_Out = arguments->Matrix[m1];
        double** Matrix_In  = arguments->Matrix[m2];

        maxResiduum = 0;

        if (rank != 0)
        {

          // recive first row from rank beneath
          MPI_Recv(Matrix_In[0], arguments->N+1 , MPI_DOUBLE, rank - 1, rank ,MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          // receive maxResiduum from rank beneath
          MPI_Recv(&maxResiduum, 1 , MPI_DOUBLE, rank - 1, rank + 1 ,MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        }

        /* over all rows */
        for (i = 1; i < arguments->numberOfRows + 1; i++)
        {
          // erste echte row wurde gerade berechnet
          // also erstmal nach oben senden
          // rank 0 muss nichts nach oben senden
          if (i == 2 && rank != 0)
          {
            // send first real line of just calculated matrix out
            // to rank - 1
            // tag = sender (rank)
            MPI_Send(Matrix_In[1], N + 1, MPI_DOUBLE, rank - 1, rank, MPI_COMM_WORLD);
          }

            double fpisin_i = 0.0;

            if (options->inf_func == FUNC_FPISIN)
            {
                fpisin_i = fpisin * sin(pih * (double)i);
            }

            /* over all columns */
            for (j = 1; j < N; j++)
            {
                star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

                if (options->inf_func == FUNC_FPISIN)
                {
                    star += fpisin_i * sin(pih * (double)j);
                }

                if (options->termination == TERM_PREC || term_iteration == 1)
                {
                    residuum = Matrix_In[i][j] - star;
                    residuum = (residuum < 0) ? -residuum : residuum;
                    maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
                }
                Matrix_Out[i][j] = star;
            }
            // hat gearde letzte row berechnet
            // und ist nicht der letzte rank (da der nichts runter sendet)
            if (i == arguments->numberOfRows && rank != nprocs - 1)
            {
              // sende letzte echte row  einen rank runter
              // sent to: rank + 1 // tag: rank + 1
              MPI_Send(Matrix_In[arguments->numberOfRows], arguments->N+1 , MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD);

              // sende maxResiduum and rank + 1
              // tag ist rank + 2
              MPI_Send(&maxResiduum,1,MPI_DOUBLE,rank + 1, rank + 2,MPI_COMM_WORLD);
          
              // warte auf erste berechnete row in dem darunter berechnetem rank
              MPI_Recv(Matrix_In[arguments->numberOfRows + 1], arguments->N+1 , MPI_DOUBLE, rank + 1,rank + 1,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        results->stat_iteration++;
        results->stat_precision = maxResiduum;

        /* exchange m1 and m2 */
        i = m1;
        m1 = m2;
        m2 = i;


        if (options->termination == TERM_ITER || sentStopAlready)
        {
            term_iteration--;
        }
        
        /* check for stopping calculation depending on termination method */
        if (options->termination == TERM_PREC)
        {   // wenn letzter rank dann check abbruch
            if (!sentStopAlready)
            {
            if (rank == nprocs - 1) {
              if (maxResiduum < options->term_precision)
              {
                int send = 1;
                // send a 1 to all ranks
                for (int allranks = 0; allranks < nprocs - 1; allranks++) {
                  MPI_Send(&send, 1, MPI_INT, allranks , nprocs, MPI_COMM_WORLD);
                  term_iteration = rank;
                  sentStopAlready = 1;
                }
              }
            }
            else
            {
              //flag =  0 = not done / 1 = done

              // check ob letzter rank abbrechen will
              // source = nprocs - 1
              // tag = nprocs
              MPI_Iprobe(nprocs - 1,nprocs,MPI_COMM_WORLD,&flag,&status);
              if (flag)
              {
                term_iteration = rank;
                sentStopAlready = 1;
              }
            }
            }

        }

    }

    if (nprocs > 1)
    {
      if (rank == 0)
      {
      // empfange maxResiduum vom letzten rank
      MPI_Recv(&maxResiduum, 1 , MPI_DOUBLE, nprocs - 1, 0 ,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      results->stat_precision = maxResiduum;
      }
      else if (rank == nprocs - 1)
      {
      //sende maxResiduum an rank 0, da dieser die stats ausgibt
      MPI_Send(&maxResiduum,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
      }
    }
    results->m = m2;
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculateJacobi (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
    int i, j;           /* local variables for loops */
    int m1, m2;         /* used as indices for old and new matrices */
    double star;        /* four times center value minus 4 neigh.b values */
    double residuum;    /* residuum of current iteration */
    double maxResiduum; /* maximum residuum value of a slave in iteration */
    double maxResiduumGlobal;

    int const N = arguments->N;
    double const h = arguments->h;

    double pih = 0.0;
    double fpisin = 0.0;

    int term_iteration = options->term_iteration;

    m1 = 0;
    m2 = 1;


    if (options->inf_func == FUNC_FPISIN)
    {
        pih = PI * h;
        fpisin = 0.25 * TWO_PI_SQUARE * h * h;
    }

    while (term_iteration > 0)
    {
        double** Matrix_Out = arguments->Matrix[m1];
        double** Matrix_In  = arguments->Matrix[m2];

        maxResiduum = 0;


        if(rank > 0)
        {
            MPI_Sendrecv(Matrix_In[1], arguments->N+1, MPI_DOUBLE, rank - 1, rank, Matrix_In[0], arguments->N + 1, MPI_DOUBLE, rank - 1, rank - 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if(rank < nprocs - 1)
        {
            MPI_Sendrecv(Matrix_In[arguments->numberOfRows], arguments->N+1, MPI_DOUBLE, rank + 1, rank, Matrix_In[arguments->numberOfRows + 1], arguments->N + 1, MPI_DOUBLE, rank + 1, rank + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        /* over all rows */
        for (i = 1; i < arguments->numberOfRows + 1; i++)
        {
            double fpisin_i = 0.0;

            if (options->inf_func == FUNC_FPISIN)
            {
                // i + start is the overall row number
                fpisin_i = fpisin * sin(pih * ((double)i + arguments->start - 1));
            }


            /* over all columns */
            for (j = 1; j < N; j++)
            {
                star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

                if (options->inf_func == FUNC_FPISIN)
                {
                    star += fpisin_i * sin(pih * (double)j);
                }

                if (options->termination == TERM_PREC || term_iteration == 1)
                {
                    residuum = Matrix_In[i][j] - star;
                    residuum = (residuum < 0) ? -residuum : residuum;
                    maxResiduum = (residuum < maxResiduum) ? maxResiduum : residuum;
                }

                Matrix_Out[i][j] = star;
            }
        }

        MPI_Allreduce(&maxResiduum, &maxResiduumGlobal, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        results->stat_iteration++;
        results->stat_precision = maxResiduumGlobal;
        maxResiduum = maxResiduumGlobal;
        /* exchange m1 and m2 */
        i = m1;
        m1 = m2;
        m2 = i;

        /* check for stopping calculation depending on termination method */
        if (options->termination == TERM_PREC)
        {
            if (maxResiduumGlobal < options->term_precision)
            {
                term_iteration = 0;
            }
        }
        else if (options->termination == TERM_ITER)
        {
            term_iteration--;
        }
    }

    results->m = m2;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static
void
displayStatistics (struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options)
{
    int N = arguments->N;
    double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

    printf("Berechnungszeit:    %f s \n", time);
    printf("Speicherbedarf:     %f MiB\n", (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
    printf("Berechnungsmethode: ");

    if (options->method == METH_GAUSS_SEIDEL)
    {
        printf("Gauß-Seidel");
    }
    else if (options->method == METH_JACOBI)
    {
        printf("Jacobi");
    }

    printf("\n");
    printf("Interlines:         %" PRIu64 "\n",options->interlines);
    printf("Stoerfunktion:      ");

    if (options->inf_func == FUNC_F0)
    {
        printf("f(x,y) = 0");
    }
    else if (options->inf_func == FUNC_FPISIN)
    {
        printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
    }

    printf("\n");
    printf("Terminierung:       ");

    if (options->termination == TERM_PREC)
    {
        printf("Hinreichende Genaugkeit");
    }
    else if (options->termination == TERM_ITER)
    {
        printf("Anzahl der Iterationen");
    }

    printf("\n");
    printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
    printf("Norm des Fehlers:   %e\n", results->stat_precision);
    printf("\n");
}

/****************************************************************************/
/** Beschreibung der Funktion displayMatrix:                               **/
/**                                                                        **/
/** Die Funktion displayMatrix gibt eine Matrix                            **/
/** in einer "ubersichtlichen Art und Weise auf die Standardausgabe aus.   **/
/**                                                                        **/
/** Die "Ubersichtlichkeit wird erreicht, indem nur ein Teil der Matrix    **/
/** ausgegeben wird. Aus der Matrix werden die Randzeilen/-spalten sowie   **/
/** sieben Zwischenzeilen ausgegeben.                                      **/
/****************************************************************************/
static
void
DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to)
{
  int const elements = 8 * options->interlines + 9;

  int x, y;
  double** Matrix = arguments->Matrix[results->m];
  MPI_Status status;

  /* first line belongs to rank 0 */
  if (rank == 0)
    from--;

  /* last line belongs to rank size - 1 */
  if (rank + 1 == size)
    to++;

  if (rank == 0)
    printf("Matrix:\n");

  for (y = 0; y < 9; y++)
  {
    int line = y * (options->interlines + 1);

    if (rank == 0)
    {
      /* check whether this line belongs to rank 0 */
      if (line < from || line > to)
      {
        /* use the tag to receive the lines in the correct order
         * the line is stored in Matrix[0], because we do not need it anymore */
        MPI_Recv(Matrix[0], elements, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y, MPI_COMM_WORLD, &status);
      }
    }
    else
    {
      if (line >= from && line <= to)
      {
        /* if the line belongs to this process, send it to rank 0
         * (line - from + 1) is used to calculate the correct local address */
        MPI_Send(Matrix[line - from + 1], elements, MPI_DOUBLE, 0, 42 + y, MPI_COMM_WORLD);
      }
    }

    if (rank == 0)
    {
      for (x = 0; x < 9; x++)
      {
        int col = x * (options->interlines + 1);

        if (line >= from && line <= to)
        {
          /* this line belongs to rank 0 */
          printf("%7.4f", Matrix[line][col]);
        }
        else
        {
          /* this line belongs to another rank and was received above */
          printf("%7.4f", Matrix[0][col]);
        }
      }

      printf("\n");
    }
  }

  fflush(stdout);
}



/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{
    struct options options;
    struct calculation_arguments arguments;
    struct calculation_results results;

    int valid_input = 0;

	MPI_Init(&argc, &argv);
    // get nr of processes
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    // get rank of this process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // prozess 0 macht die eingabe
    if(rank==0)
    {
    	askParams(&options, argc, argv);
	}
    
    /* 
    HIER FEHLT NOCH DER CHECK OB ES GEKLAPPT HAT, ABER askparams HAT 
    KEINEN RETURN, DAHER WEISS ICH NOCH NICHT WIE DAS ZU LOESEN IST OHNE 
    QUASI DIE GESAMTE INPUTVALIDATION AUS askparams.c HIER HIN ZU KOPIEREN 
    */
    // abbrechen, falls die eingabe nicht klappt
    MPI_Bcast(&valid_input, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(valid_input != 0)
    {
    	MPI_Finalize();
    	exit(1);
    }

    // die eigentliche eingabe an die anderen prozesse senden
    MPI_Bcast(&options.number, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&options.method, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&options.interlines, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&options.inf_func, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&options.termination, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&options.term_iteration, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&options.term_precision, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    
    initVariables(&arguments, &results, &options);
    allocateMatrices(&arguments, &options);
    initMatrices(&arguments, &options);

    gettimeofday(&start_time, NULL);
    if (options.method == METH_JACOBI)
    {
      calculateJacobi(&arguments, &results, &options);
    }
    else
    {
      calculateGS(&arguments, &results, &options);
    }
    gettimeofday(&comp_time, NULL);

    if(rank == 0)
    {
      displayStatistics(&arguments, &results, &options);
    }
    DisplayMatrix (&arguments,&results, &options,rank, nprocs, (int)arguments.start, (int)arguments.end);
    freeMatrices(&arguments);
    MPI_Finalize();
    return 0;
}
