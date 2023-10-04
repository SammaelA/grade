#include "optimization.h"
#include "common_utils/blk.h"
#include "third_party/CMA_ES/cmaes_interface.h"

namespace opt
{
  double fitfun(double const *x, int dim); 

    /* the objective (fitness) function to be minimized */
    double fitfun(double const *x, int N) { /* function "cigtab" */
      int i; 
      double sum = 1e4*x[0]*x[0] + 1e-4*x[1]*x[1];
      for(i = 2; i < N; ++i)  
        sum += x[i]*x[i]; 
      return sum;
  }

  void CMA_ES::optimize(opt_func_with_grad_vector &F, const std::vector<float> &min_X, const std::vector<float> &max_X, Block &settings,
                        init_params_func &get_init_params)
  {
    cmaes_t evo; /* an CMA-ES type struct or "object" */
    double *arFunvals, *const*pop, *xfinal;
    int i; 

    /* Initialize everything into the struct evo, 0 means default */
    arFunvals = cmaes_init(&evo, 0, NULL, NULL, 0, 0, "src/third_party/CMA_ES/cmaes_initials.par");
    printf("%s\n", cmaes_SayHello(&evo));
    cmaes_ReadSignals(&evo, "src/third_party/CMA_ES/cmaes_signals.par");  /* write header and initial values */

    /* Iterate until stop criterion holds */
    while(!cmaes_TestForTermination(&evo))
      { 
        /* generate lambda new search points, sample population */
        pop = cmaes_SamplePopulation(&evo); /* do not change content of pop */

        /* Here we may resample each solution point pop[i] until it
      becomes feasible. function is_feasible(...) needs to be
      user-defined.
      Assumptions: the feasible domain is convex, the optimum is
      not on (or very close to) the domain boundary, initialX is
      feasible and initialStandardDeviations are sufficiently small
      to prevent quasi-infinite looping. */
        /* for (i = 0; i < cmaes_Get(&evo, "popsize"); ++i)
            while (!is_feasible(pop[i]))
              cmaes_ReSampleSingle(&evo, i);
        */

        /* evaluate the new search points using fitfun */
        for (i = 0; i < cmaes_Get(&evo, "lambda"); ++i) {
          arFunvals[i] = fitfun(pop[i], (int) cmaes_Get(&evo, "dim"));
        }

        /* update the search distribution used for cmaes_SamplePopulation() */
        cmaes_UpdateDistribution(&evo, arFunvals);  

        /* read instructions for printing output or changing termination conditions */ 
        cmaes_ReadSignals(&evo, "src/third_party/CMA_ES/cmaes_signals.par");
        fflush(stdout); /* useful in MinGW */
      }
    printf("Stop:\n%s\n",  cmaes_TestForTermination(&evo)); /* print termination reason */
    cmaes_WriteToFile(&evo, "all", "src/third_party/CMA_ES/allcmaes.dat");         /* write final results */

    /* get best estimator for the optimum, xmean */
    xfinal = cmaes_GetNew(&evo, "xmean"); /* "xbestever" might be used as well */
    cmaes_exit(&evo); /* release memory */ 

    /* do something with final solution and finally release memory */
    free(xfinal); 
  }
}