#include "biom.hpp"
#include "skbio_alt.hpp"
#include "api.hpp"
#include <unistd.h>
#include <vector>
#include "test_helper.hpp"

void test_center_mat() {
    SUITE_START("test center mat");

    // unweighted unifrac of test.biom
    double matrix[] = { 
      0.0000000000, 0.2000000000, 0.5714285714, 0.6000000000, 0.5000000000, 0.2000000000,
      0.2000000000, 0.0000000000, 0.4285714286 ,0.6666666667, 0.6000000000, 0.3333333333,
      0.5714285714, 0.4285714286, 0.0000000000, 0.7142857143, 0.8571428571, 0.4285714286,
      0.6000000000, 0.6666666667, 0.7142857143, 0.0000000000, 0.3333333333, 0.4000000000,
      0.5000000000, 0.6000000000, 0.8571428571, 0.3333333333, 0.0000000000, 0.6000000000,
      0.2000000000, 0.3333333333, 0.4285714286, 0.4000000000, 0.6000000000, 0.0000000000};

    const uint32_t n_samples = 6;

    double exp[] = { 0.05343726,  0.04366213, -0.0329743 , -0.07912698, -0.00495654, 0.01995843,
                     0.04366213,  0.073887  ,  0.04867914, -0.11112434, -0.04973167,-0.00537226,
                    -0.0329743 ,  0.04867914,  0.20714475, -0.07737528, -0.17044974, 0.02497543,
                    -0.07912698, -0.11112434, -0.07737528,  0.14830877,  0.11192366, 0.00739418,
                    -0.00495654, -0.04973167, -0.17044974,  0.11192366,  0.18664966,-0.07343537,
                     0.01995843, -0.00537226,  0.02497543,  0.00739418, -0.07343537, 0.02647959 };

    {
      double *centered = (double *) malloc(6*6*sizeof(double)); 

      su::mat_to_centered(matrix, n_samples, centered);

      for(int i = 0; i < (6*6); i++) {
        //printf("%i %f %f\n",i,float(centered[i]),float(exp[i]));
        ASSERT(fabs(centered[i] - exp[i]) < 0.000001);
      }

      free(centered);
    }

    float *matrix_fp32 = (float *) malloc(6*6*sizeof(float));
    for(int i = 0; i < (6*6); i++) matrix_fp32[i] = matrix[i];


    {
      float *centered_fp32 = (float *) malloc(6*6*sizeof(float));

      su::mat_to_centered(matrix_fp32, n_samples, centered_fp32);

      for(int i = 0; i < (6*6); i++) {
        //printf("%i %f %f\n",i,float(centered_fp32[i]),float(exp[i]));
        ASSERT(fabs(centered_fp32[i] - exp[i]) < 0.000001);
      }

      free(centered_fp32);
    }

    free(matrix_fp32);


    SUITE_END();
}

void test_pcoa() {
    SUITE_START("test pcoa");

    // unweighted unifrac of crawford.biom
    double matrix[] = { 
      0.         , 0.71836067 , 0.71317361 , 0.69746044 , 0.62587207 , 0.72826674
    , 0.72065895 , 0.72640581 , 0.73606053, 
      0.71836067 , 0.         , 0.70302967 , 0.73407301 , 0.6548042  , 0.71547381
    , 0.78397813 , 0.72318399 , 0.76138933,
      0.71317361 , 0.70302967 , 0.         , 0.61041275 , 0.62331299 , 0.71848305
    , 0.70416337 , 0.75258475 , 0.79249029,
      0.69746044 , 0.73407301 , 0.61041275 , 0.         , 0.6439278  , 0.70052733
    , 0.69832716 , 0.77818938 , 0.72959894,
      0.62587207 , 0.6548042  , 0.62331299 , 0.6439278  , 0.         , 0.75782689
    , 0.71005144 , 0.75065046 , 0.78944369,
      0.72826674 , 0.71547381 , 0.71848305 , 0.70052733 , 0.75782689 , 0.
    , 0.63593642 , 0.71283615 , 0.58314638,
      0.72065895 , 0.78397813 , 0.70416337 , 0.69832716 , 0.71005144 , 0.63593642
    , 0.         , 0.69200762 , 0.68972056,
      0.72640581 , 0.72318399 , 0.75258475 , 0.77818938 , 0.75065046 , 0.71283615
    , 0.69200762 , 0.         , 0.71514083,
      0.73606053 , 0.76138933 , 0.79249029 , 0.72959894 , 0.78944369 , 0.58314638
    , 0.68972056 , 0.71514083 , 0. };


    const uint32_t n_samples = 9;

    // Test centering

    double exp2[] = { 
        0.22225336 , -0.025481   , -0.03491711 , -0.02614685  , 0.01899659 , -0.05103589
     , -0.03971812 , -0.02699392 , -0.03695706,
       -0.025481   ,  0.24282669 , -0.01744751 , -0.04206624  , 0.0107569  , -0.03151438
     , -0.07706765 , -0.0143721  , -0.0456347,
       -0.03491711 , -0.01744751 ,  0.21652901 ,  0.02791465  , 0.0177328  , -0.04682078
     , -0.03082866 , -0.04921529 , -0.08294711,
       -0.02614685 , -0.04206624 ,  0.02791465 ,  0.21190401  , 0.00235833 , -0.03639361
     , -0.02904855 , -0.07112525 , -0.03739649,
        0.01899659 ,  0.0107569  ,  0.0177328  ,  0.00235833  , 0.20745566 , -0.08039931
     , -0.03952883 , -0.05229812 , -0.08507403,
       -0.05103589 , -0.03151438 , -0.04682078 , -0.03639361  , -0.08039931 , 0.20604732
     ,  0.00964595 , -0.02533192 ,  0.05580262,
       -0.03971812 , -0.07706765 , -0.03082866 , -0.02904855  , -0.03952883 , 0.00964595
     ,  0.21765972 , -0.00489531 , -0.00621856,
       -0.02699392 , -0.0143721  , -0.04921529 , -0.07112525  , -0.05229812 , -0.02533192
     , -0.00489531 ,  0.2514242  , -0.00719229,
       -0.03695706 , -0.0456347  , -0.08294711 , -0.03739649  , -0.08507403 , 0.05580262
     , -0.00621856 , -0.00719229 , 0.24561762 };

    double *centered = (double *) malloc(9*9*sizeof(double));

    su::mat_to_centered(matrix, n_samples, centered);

    for(int i = 0; i < (9*9); i++) {
      //printf("%i %f %f\n",i,float(centered[i]),float(exp[i]));
      ASSERT(fabs(centered[i] - exp2[i]) < 0.00001);
    }

    // Test eigens

    double exp3a[] = {0.45752162, 0.3260088 , 0.2791141 , 0.26296948, 0.20924533};
    double exp3b[] = {
       -0.17316152,  0.17579996,  0.23301609, -0.74519625, -0.05194624,
       -0.19959264,  0.53235665, -0.53370018,  0.2173474 ,  0.26736004,
       -0.35794942, -0.27956624,  0.01114096,  0.40488848, -0.13121464,
       -0.2296467 , -0.47494333, -0.12571292,  0.02313551, -0.46916459,
       -0.44501584,  0.05597451,  0.07717711, -0.15881922,  0.24594442,
        0.40335552, -0.16290597, -0.30327343,  0.03778646,  0.21664806,
        0.23769142, -0.29034629,  0.46813757,  0.13858945,  0.58624834,
        0.2407584 ,  0.51300752,  0.48211607,  0.34422672, -0.42416046,
        0.52356078, -0.0693768 , -0.30890127, -0.26195855, -0.23971493};

    {
      double *eigenvalues;
      double *eigenvectors;
      su::find_eigens_fast(n_samples, 5, centered, eigenvalues, eigenvectors);

      for(int i = 0; i < 5; i++) {
        //printf("%i %f %f\n",i,float(eigenvalues[i]),float(exp3a[i]));
        ASSERT(fabs(eigenvalues[i] - exp3a[i]) < 0.00001);
      }

      // signs may flip, that's normal
      for(int i = 0; i < (5*9); i++) {
        //printf("%i %f %f %f\n",i,float(eigenvectors[i]),float(exp3b[i]),float(fabs(eigenvectors[i]) - fabs(exp3b[i])));
        ASSERT( fabs(fabs(eigenvectors[i]) - fabs(exp3b[i])) < 0.00001);
      }
    
      free(eigenvectors);
      free(eigenvalues);
    }
    free(centered);


    // Test PCoA (incudes the above calls

    double *exp4a = exp3a; // same eigenvals;
    double exp4b[] = {
       -0.11712705,  0.10037682, -0.12310531, -0.38214073, -0.02376195,
       -0.13500515,  0.30396064,  0.28196047,  0.11145694,  0.12229942,
       -0.24211822, -0.15962444, -0.00588591,  0.20762904, -0.06002196,
       -0.15533382, -0.27117925,  0.06641571,  0.01186402, -0.21461156,
       -0.30101024,  0.03195987, -0.04077363, -0.08144337,  0.1125032 ,
        0.27283106, -0.09301471,  0.16022314,  0.0193771 ,  0.09910206,
        0.16077529, -0.16577955, -0.24732293,  0.07106943,  0.26816958,
        0.16284981,  0.29291283, -0.25470794,  0.17652136, -0.19402517,
        0.35413832, -0.0396122 ,  0.1631964 , -0.13433378, -0.10965362};
    double exp4c[] = {0.22630343, 0.16125338, 0.13805791, 0.13007231, 0.10349879};

    {
      double *eigenvalues;
      double *samples;
      double *proportion_explained;

      su::pcoa(matrix, n_samples, 5, eigenvalues, samples, proportion_explained); 

      for(int i = 0; i < 5; i++) {
        //printf("%i %f %f\n",i,float(eigenvalues[i]),float(exp4a[i]));
        ASSERT(fabs(eigenvalues[i] - exp4a[i]) < 0.00001);
      }

      // signs may flip, that's normal
      for(int i = 0; i < (5*9); i++) {
        //printf("%i %f %f %f\n",i,float(samples[i]),float(exp4b[i]),float(fabs(samples[i]) - fabs(exp4b[i])));
        ASSERT( fabs(fabs(samples[i]) - fabs(exp4b[i])) < 0.00001);
      }

      for(int i = 0; i < 5; i++) {
        //printf("%i %f %f\n",i,float(proportion_explained[i]),float(exp4c[i]));
        ASSERT(fabs(proportion_explained[i] - exp4c[i]) < 0.00001);
      } 

      free(eigenvalues);
      free(samples);
      free(proportion_explained);
    }

    // Test PCoA mixed mode
    {
      float *eigenvalues_fp32;
      float *samples_fp32;
      float *proportion_explained_fp32;

      su::pcoa(matrix, n_samples, 5, eigenvalues_fp32, samples_fp32, proportion_explained_fp32);

      for(int i = 0; i < 5; i++) {
        //printf("%i %f %f\n",i,float(eigenvalues_fp32[i]),float(exp4a[i]));
       ASSERT(fabs(eigenvalues_fp32[i] - exp4a[i]) < 0.00001);
      }

      // signs may flip, that's normal
      for(int i = 0; i < (5*9); i++) {
        //printf("%i %f %f %f\n",i,float(samples_fp32[i]),float(exp4b[i]),float(fabs(samples_fp32[i]) - fabs(exp4b[i])));
        ASSERT( fabs(fabs(samples_fp32[i]) - fabs(exp4b[i])) < 0.00001);
      }

      for(int i = 0; i < 5; i++) {
        //printf("%i %f %f\n",i,float(proportion_explained_fp32[i]),float(exp4c[i]));
        ASSERT(fabs(proportion_explained_fp32[i] - exp4c[i]) < 0.00001);
      }

      free(eigenvalues_fp32);
      free(samples_fp32);
      free(proportion_explained_fp32);
    }

    // test in-place
    {
      double *eigenvalues;
      double *samples;
      double *proportion_explained;

      su::pcoa_inplace(matrix, n_samples, 5, eigenvalues, samples, proportion_explained);
      // Note: matrix content has been destroyed

      for(int i = 0; i < 5; i++) {
        //printf("%i %f %f\n",i,float(eigenvalues[i]),float(exp4a[i]));
        ASSERT(fabs(eigenvalues[i] - exp4a[i]) < 0.00001);
      }

      // signs may flip, that's normal
      for(int i = 0; i < (5*9); i++) {
        //printf("%i %f %f %f\n",i,float(samples[i]),float(exp4b[i]),float(fabs(samples[i]) - fabs(exp4b[i])));
        ASSERT( fabs(fabs(samples[i]) - fabs(exp4b[i])) < 0.00001);
      }

      for(int i = 0; i < 5; i++) {
        //printf("%i %f %f\n",i,float(proportion_explained[i]),float(exp4c[i]));
        ASSERT(fabs(proportion_explained[i] - exp4c[i]) < 0.00001);
      }

      free(eigenvalues);
      free(samples);
      free(proportion_explained);
    }

    SUITE_END();
}

void test_pcoa_big() {
    SUITE_START("test pcoa big");

    //too big to inline, use support file
    FILE *fptr = fopen("test_ska_pcoa_big.dat", "r");
    ASSERT(fptr != NULL)

    unsigned int n_samples = 0;
    ASSERT(fscanf(fptr,"# unifrac %u\n",&n_samples)==1);
    ASSERT(n_samples == 57);


    // first 57 rows/cols of unweighted unifrac of EMP
    double matrix[57*57];
    for (unsigned int i=0; i<(57*57); i++)
     ASSERT(fscanf(fptr,"%lf\n",&(matrix[i]))==1);

    unsigned int n_dims = 0;
    ASSERT(fscanf(fptr,"# pcoa %u\n",&n_dims)==1);
    ASSERT(n_dims == 7);

    double exp1[7];
    for (unsigned int i=0; i<(7); i++)
     ASSERT(fscanf(fptr,"%lf\n",&(exp1[i]))==1);

    double exp2[7*57];
    for (unsigned int i=0; i<(7*57); i++)
     ASSERT(fscanf(fptr,"%lf\n",&(exp2[i]))==1);

    double exp3[7];
    for (unsigned int i=0; i<(7); i++)
     ASSERT(fscanf(fptr,"%lf\n",&(exp3[i]))==1);
    
    fclose(fptr);

    {
      double *eigenvalues;
      double *samples;
      double *proportion_explained;

      su::pcoa(matrix, n_samples, n_dims, eigenvalues, samples, proportion_explained);

      // last three eignes are very close to each other and could come back in reverse order

      for(unsigned int i = 0; i < n_dims ; i++) {
        //printf("%i %f %f %f\n",i,float(eigenvalues[i]),float(exp1[i]),float(fabs(eigenvalues[i] - exp1[i])));
        const double max_err = (i<4) ? 0.01 : 0.1; // the values are approximate, based on a random number in the algo
        ASSERT(fabs(eigenvalues[i] - exp1[i]) < max_err);
      }

      // signs may flip, that's normal
      for(unsigned int i = 0; i < (n_samples*n_dims); i++) {
        if ((i%n_dims)<4) {
          //printf("%i %f %f %f\n",i,float(samples[i]),float(exp2[i]),float(fabs(samples[i]) - fabs(exp2[i])));
          ASSERT( fabs(fabs(samples[i]) - fabs(exp2[i])) < 0.05)
        } else {
          // any of the 3 will do
          unsigned int ibase = (i/n_dims)*n_dims;
          //printf("%i %f %f %f %f\n",i,float(samples[i]),float(exp2[ibase+4]),float(exp2[ibase+5]),float(exp2[ibase+6]));
          ASSERT( (fabs(fabs(samples[i]) - fabs(exp2[ibase+4])) < 0.15) || (fabs(fabs(samples[i]) - fabs(exp2[ibase+5])) < 0.15) || (fabs(fabs(samples[i]) - fabs(exp2[ibase+6])) < 0.15) );
        }
      }

      for(unsigned int i = 0; i < n_dims; i++) {
        //printf("%i %f %f\n",i,float(proportion_explained[i]),float(exp3[i]),float(fabs(proportion_explained[i] - exp3[i])));
        const double max_err = (i<4) ? 0.001 : 0.01;
        ASSERT(fabs(proportion_explained[i] - exp3[i]) < max_err);
      }

      free(proportion_explained);
      free(samples);
      free(eigenvalues);
    }

    {
      float *eigenvalues;
      float *samples;
      float *proportion_explained;

      su::pcoa(matrix, n_samples, n_dims, eigenvalues, samples, proportion_explained);

      // last three eignes are very close to each other and could come back in reverse order

      for(unsigned int i = 0; i < n_dims ; i++) {
        //printf("%i %f %f %f\n",i,float(eigenvalues[i]),float(exp1[i]),float(fabs(eigenvalues[i] - exp1[i])));
        const float max_err = (i<4) ? 0.01 : 0.1; // the values are approximate, based on a random number in the algo
        ASSERT(fabs(eigenvalues[i] - exp1[i]) < max_err);
      }

      // signs may flip, that's normal
      for(unsigned int i = 0; i < (n_samples*n_dims); i++) {
        if ((i%n_dims)<4) {
          //printf("%i %f %f %f\n",i,float(samples[i]),float(exp2[i]),float(fabs(samples[i]) - fabs(exp2[i])));
          ASSERT( fabs(fabs(samples[i]) - fabs(exp2[i])) < 0.1)
        } else {
          // any of the 3 will do
          unsigned int ibase = (i/n_dims)*n_dims;
          //printf("%i %f %f %f %f\n",i,float(samples[i]),float(exp2[ibase+4]),float(exp2[ibase+5]),float(exp2[ibase+6]));
          ASSERT( (fabs(fabs(samples[i]) - fabs(exp2[ibase+4])) < 0.15) || (fabs(fabs(samples[i]) - fabs(exp2[ibase+5])) < 0.15) || (fabs(fabs(samples[i]) - fabs(exp2[ibase+6])) < 0.15) );
        }
      }

      for(unsigned int i = 0; i < n_dims; i++) {
        //printf("%i %f %f\n",i,float(proportion_explained[i]),float(exp3[i]),float(fabs(proportion_explained[i] - exp3[i])));
        const float max_err = (i<4) ? 0.001 : 0.01;
        ASSERT(fabs(proportion_explained[i] - exp3[i]) < max_err);
      }

      free(proportion_explained);
      free(samples);
      free(eigenvalues);
    }

    {
      float matrix_fp32[57*57];
      for (unsigned int i=0; i<(57*57); i++)
        matrix_fp32[i] = matrix[i];

      {
        float *eigenvalues;
        float *samples;
        float *proportion_explained;

        su::pcoa(matrix_fp32, n_samples, n_dims, eigenvalues, samples, proportion_explained);

        // last three eignes are very close to each other and could come back in reverse order

        for(unsigned int i = 0; i < n_dims ; i++) {
          //printf("%i %f %f %f\n",i,float(eigenvalues[i]),float(exp1[i]),float(fabs(eigenvalues[i] - exp1[i])));
          const float max_err = (i<4) ? 0.01 : 0.1; // the values are approximate, based on a random number in the algo
          ASSERT(fabs(eigenvalues[i] - exp1[i]) < max_err);
        }

        // signs may flip, that's normal
        for(unsigned int i = 0; i < (n_samples*n_dims); i++) {
          if ((i%n_dims)<4) {
            //printf("%i %f %f %f\n",i,float(samples[i]),float(exp2[i]),float(fabs(samples[i]) - fabs(exp2[i])));
            ASSERT( fabs(fabs(samples[i]) - fabs(exp2[i])) < 0.1)
          } else {
            // any of the 3 will do
            unsigned int ibase = (i/n_dims)*n_dims;
            //printf("%i %f %f %f %f\n",i,float(samples[i]),float(exp2[ibase+4]),float(exp2[ibase+5]),float(exp2[ibase+6]));
            ASSERT( (fabs(fabs(samples[i]) - fabs(exp2[ibase+4])) < 0.15) || (fabs(fabs(samples[i]) - fabs(exp2[ibase+5])) < 0.15) || (fabs(fabs(samples[i]) - fabs(exp2[ibase+6])) < 0.15) );
          }
        }

        for(unsigned int i = 0; i < n_dims; i++) {
          //printf("%i %f %f\n",i,float(proportion_explained[i]),float(exp3[i]),float(fabs(proportion_explained[i] - exp3[i])));
          const float max_err = (i<4) ? 0.001 : 0.01;
          ASSERT(fabs(proportion_explained[i] - exp3[i]) < max_err);
        }

        free(proportion_explained);
        free(samples);
        free(eigenvalues);
      }

      {
        float *eigenvalues;
        float *samples;
        float *proportion_explained;

        su::pcoa_inplace(matrix_fp32, n_samples, n_dims, eigenvalues, samples, proportion_explained);
        // Note: content of matrix_fp32 has been destroyed 

        // last three eignes are very close to each other and could come back in reverse order

        for(unsigned int i = 0; i < n_dims ; i++) {
          //printf("%i %f %f %f\n",i,float(eigenvalues[i]),float(exp1[i]),float(fabs(eigenvalues[i] - exp1[i])));
          const float max_err = (i<4) ? 0.01 : 0.1; // the values are approximate, based on a random number in the algo
          ASSERT(fabs(eigenvalues[i] - exp1[i]) < max_err);
        }

        // signs may flip, that's normal
        for(unsigned int i = 0; i < (n_samples*n_dims); i++) {
          if ((i%n_dims)<4) {
            //printf("%i %f %f %f\n",i,float(samples[i]),float(exp2[i]),float(fabs(samples[i]) - fabs(exp2[i])));
            ASSERT( fabs(fabs(samples[i]) - fabs(exp2[i])) < 0.1)
          } else {
            // any of the 3 will do
            unsigned int ibase = (i/n_dims)*n_dims;
            //printf("%i %f %f %f %f\n",i,float(samples[i]),float(exp2[ibase+4]),float(exp2[ibase+5]),float(exp2[ibase+6]));
            ASSERT( (fabs(fabs(samples[i]) - fabs(exp2[ibase+4])) < 0.15) || (fabs(fabs(samples[i]) - fabs(exp2[ibase+5])) < 0.15) || (fabs(fabs(samples[i]) - fabs(exp2[ibase+6])) < 0.15) );
          }
        }

        for(unsigned int i = 0; i < n_dims; i++) {
          //printf("%i %f %f\n",i,float(proportion_explained[i]),float(exp3[i]),float(fabs(proportion_explained[i] - exp3[i])));
          const float max_err = (i<4) ? 0.001 : 0.01;
          ASSERT(fabs(proportion_explained[i] - exp3[i]) < max_err);
        }

        free(proportion_explained);
        free(samples);
        free(eigenvalues);
      }

    }



    SUITE_END();
}

/* PCoA draws a random matrix for the randomized SVD. With an explicit seed that
 * draw comes from a generator local to the call, so the answer is reproducible
 * without touching -- or being disturbed by -- the process-global generator.
 * That is what makes several PCoAs in flight at once possible. A negative seed
 * keeps the legacy behaviour of drawing from the global generator that
 * su::set_random_seed() sets.
 */
static double max_abs_diff(const std::vector<double> &a, const std::vector<double> &b) {
    if (a.size() != b.size()) return 1.0/0.0;
    double m = 0.0;
    for (size_t i = 0; i < a.size(); i++) m = std::max(m, fabs(a[i] - b[i]));
    return m;
}

static std::vector<double> run_seeded_pcoa(const double *mat, uint32_t n_samples,
                                           uint32_t n_dims, int seed) {
    double *eigenvalues = NULL, *samples = NULL, *proportion_explained = NULL;
    su::pcoa(mat, n_samples, n_dims, eigenvalues, samples, proportion_explained, seed);

    std::vector<double> out;
    out.insert(out.end(), eigenvalues, eigenvalues + n_dims);
    out.insert(out.end(), samples, samples + (size_t(n_dims) * n_samples));
    out.insert(out.end(), proportion_explained, proportion_explained + n_dims);
    free(eigenvalues);
    free(samples);
    free(proportion_explained);
    return out;
}

static double run_seeded_permanova(const double *mat, uint32_t n_samples,
                                   const uint32_t *grouping, unsigned int n_perm, int seed) {
    double fstat = 0., pvalue = 0.;
    su::permanova(mat, n_samples, grouping, n_perm, fstat, pvalue, seed);
    return pvalue;
}

void test_pcoa_seeded() {
    SUITE_START("test pcoa seeded");

    // unweighted unifrac of crawford.biom, as in test_pcoa
    const double matrix[] = {
      0.         , 0.71836067 , 0.71317361 , 0.69746044 , 0.62587207 , 0.72826674
    , 0.72065895 , 0.72640581 , 0.73606053,
      0.71836067 , 0.         , 0.70302967 , 0.73407301 , 0.6548042  , 0.71547381
    , 0.78397813 , 0.72318399 , 0.76138933,
      0.71317361 , 0.70302967 , 0.         , 0.61041275 , 0.62331299 , 0.71848305
    , 0.70416337 , 0.75258475 , 0.79249029,
      0.69746044 , 0.73407301 , 0.61041275 , 0.         , 0.6439278  , 0.70052733
    , 0.69832716 , 0.77818938 , 0.72959894,
      0.62587207 , 0.6548042  , 0.62331299 , 0.6439278  , 0.         , 0.75782689
    , 0.71005144 , 0.75065046 , 0.78944369,
      0.72826674 , 0.71547381 , 0.71848305 , 0.70052733 , 0.75782689 , 0.
    , 0.63593642 , 0.71283615 , 0.58314638,
      0.72065895 , 0.78397813 , 0.70416337 , 0.69832716 , 0.71005144 , 0.63593642
    , 0.         , 0.69200762 , 0.68972056,
      0.72640581 , 0.72318399 , 0.75258475 , 0.77818938 , 0.75065046 , 0.71283615
    , 0.69200762 , 0.         , 0.71514083,
      0.73606053 , 0.76138933 , 0.79249029 , 0.72959894 , 0.78944369 , 0.58314638
    , 0.68972056 , 0.71514083 , 0. };

    const uint32_t n_samples = 9;
    const uint32_t n_dims    = 5;

    /* Tolerance, not bit-exact -- see "Ordination reproduces to a tolerance" in
     * README.md. SAME is still far tighter than the signal it discriminates:
     * changing the seed moves the answer by ~0.5 here, which SEED_DIFF pins.
     */
    const double SAME      = 1e-12;
    const double SEED_DIFF = 1e-6;

    // an explicit seed reproduces, with no seeding call in between
    std::vector<double> a = run_seeded_pcoa(matrix, n_samples, n_dims, 7);
    std::vector<double> b = run_seeded_pcoa(matrix, n_samples, n_dims, 7);
    ASSERT(max_abs_diff(a, b) < SAME);

    // ... and is not perturbed by the global generator moving underneath it
    su::set_random_seed(999);
    std::vector<double> c = run_seeded_pcoa(matrix, n_samples, n_dims, 7);
    ASSERT(max_abs_diff(a, c) < SAME);

    // a different seed gives a different draw, so the seed is really consumed
    std::vector<double> d = run_seeded_pcoa(matrix, n_samples, n_dims, 8);
    ASSERT(max_abs_diff(a, d) > SEED_DIFF);

    // a negative seed is the legacy path: same global seed, same answer
    su::set_random_seed(42);
    std::vector<double> g1 = run_seeded_pcoa(matrix, n_samples, n_dims, -1);
    su::set_random_seed(42);
    std::vector<double> g2 = run_seeded_pcoa(matrix, n_samples, n_dims, -1);
    ASSERT(max_abs_diff(g1, g2) < SAME);

    // ... and the legacy path really does follow the global generator: a
    // different global seed gives a different answer, which is precisely why it
    // cannot be called concurrently
    su::set_random_seed(4242);
    std::vector<double> g3 = run_seeded_pcoa(matrix, n_samples, n_dims, -1);
    ASSERT(max_abs_diff(g1, g3) > SEED_DIFF);

    SUITE_END();
}

void test_permanova_ties() {
    SUITE_START("test permanova ties");

    // Same as skbio tests
    const double matrix_fp64[] = { 
      0., 1., 1., 4.,
      1., 0., 3., 2.,
      1., 3., 0., 3.,
      4., 2., 3., 0.};
    const float matrix_fp32[] = { 
      0., 1., 1., 4.,
      1., 0., 3., 2.,
      1., 3., 0., 3.,
      4., 2., 3., 0.};

    const uint32_t grouping_equal[] = { 0, 0, 1, 1};
    const uint32_t grouping_equal_swapped[] = { 1, 1, 0, 0};

    const uint32_t n_samples = 4;


    // value from skbio
    const float exp_stat = 2.0;
    // using a different random than skbio, so result different
    const float exp_pvalue = 0.66;

    double stat_fp64, pvalue_fp64;
    float stat_fp32, pvalue_fp32;

    su::permanova(matrix_fp64, n_samples,
                  grouping_equal, 999,
                  stat_fp64, pvalue_fp64);
    ASSERT(fabs(stat_fp64 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp64 - exp_pvalue) < 0.05);

    su::permanova(matrix_fp32, n_samples,
                  grouping_equal, 999,
                  stat_fp32, pvalue_fp32);
    ASSERT(fabs(stat_fp64 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp64 - exp_pvalue) < 0.05);

    su::permanova(matrix_fp32, n_samples,
                  grouping_equal_swapped, 999,
                  stat_fp32, pvalue_fp32);
    ASSERT(fabs(stat_fp32 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp32 - exp_pvalue) < 0.05);

    su::permanova(matrix_fp64, n_samples,
                  grouping_equal_swapped, 999,
                  stat_fp64, pvalue_fp64);
    ASSERT(fabs(stat_fp64 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp64 - exp_pvalue) < 0.05);

    SUITE_END();
}

void test_permanova_noties() {
    SUITE_START("test permanova noties");

    // Same as skbio tests
    const double matrix_fp64[] = { 
      0., 1., 5., 4.,
      5., 0., 3., 2.,
      1., 3., 0., 3.,
      4., 2., 3., 0.};
    const float matrix_fp32[] = { 
      0., 1., 5., 4.,
      1., 0., 3., 2.,
      5., 3., 0., 3.,
      4., 2., 3., 0.};

    const uint32_t grouping_equal[] = { 0, 0, 1, 1};
    const uint32_t grouping_equal_swapped[] = { 1, 1, 0, 0};

    const uint32_t n_samples = 4;


    // value from skbio
    const float exp_stat = 4.4;
    // using a different random than skbio, so result different
    const float exp_pvalue = 0.34;

    double stat_fp64, pvalue_fp64;
    float stat_fp32, pvalue_fp32;

    su::permanova(matrix_fp64, n_samples,
                  grouping_equal, 999,
                  stat_fp64, pvalue_fp64);
    ASSERT(fabs(stat_fp64 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp64 - exp_pvalue) < 0.05);

    su::permanova(matrix_fp32, n_samples,
                  grouping_equal, 999,
                  stat_fp32, pvalue_fp32);
    ASSERT(fabs(stat_fp32 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp32 - exp_pvalue) < 0.05);

    su::permanova(matrix_fp32, n_samples,
                  grouping_equal_swapped, 999,
                  stat_fp32, pvalue_fp32);
    ASSERT(fabs(stat_fp32 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp32 - exp_pvalue) < 0.05);

    su::permanova(matrix_fp64, n_samples,
                  grouping_equal_swapped, 999,
                  stat_fp64, pvalue_fp64);
    ASSERT(fabs(stat_fp64 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp64 - exp_pvalue) < 0.05);

    SUITE_END();
}

void test_permanova_unequal() {
    SUITE_START("test permanova unequal");

    // Same as skbio tests
    const double matrix_fp64[] = { 
      0.0,    1.0,   0.1,   0.5678, 1.0,   1.0,
      1.0,    0.0,   0.002, 0.42,   0.998, 0.0,
      0.1,    0.002, 0.0,   1.0,    0.123, 1.0,
      0.5678, 0.42,  1.0,   0.0,    0.123, 0.43,
      1.0,    0.998, 0.123, 0.123,  0.0,   0.5,
      1.0,    0.0,   1.0,   0.43,   0.5,   0.0 };
    const float matrix_fp32[] = { 
      0.0,    1.0,   0.1,   0.5678, 1.0,   1.0,
      1.0,    0.0,   0.002, 0.42,   0.998, 0.0,
      0.1,    0.002, 0.0,   1.0,    0.123, 1.0,
      0.5678, 0.42,  1.0,   0.0,    0.123, 0.43,
      1.0,    0.998, 0.123, 0.123,  0.0,   0.5,
      1.0,    0.0,   1.0,   0.43,   0.5,   0.0 };

    const uint32_t grouping_equal[] = { 0, 1, 2, 1, 0 , 0};
    const uint32_t grouping_equal_swapped[] = { 1, 2, 0, 2, 1, 1};

    const uint32_t n_samples = 6;


    // value from skbio
    const float exp_stat = 0.578848;
    // using a different random than skbio, so result different
    const float exp_pvalue = 0.65;

    double stat_fp64, pvalue_fp64;
    float stat_fp32, pvalue_fp32;

    su::permanova(matrix_fp64, n_samples,
                  grouping_equal, 999,
                  stat_fp64, pvalue_fp64);
    ASSERT(fabs(stat_fp64 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp64 - exp_pvalue) < 0.05);

    su::permanova(matrix_fp32, n_samples,
                  grouping_equal, 999,
                  stat_fp32, pvalue_fp32);
    ASSERT(fabs(stat_fp32 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp32 - exp_pvalue) < 0.05);

    su::permanova(matrix_fp32, n_samples,
                  grouping_equal, 999,
                  stat_fp32, pvalue_fp32);
    ASSERT(fabs(stat_fp32 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp32 - exp_pvalue) < 0.05);

    su::permanova(matrix_fp64, n_samples,
                  grouping_equal_swapped, 999,
                  stat_fp64, pvalue_fp64);
    ASSERT(fabs(stat_fp64 - exp_stat) < 0.00001);
    ASSERT(fabs(pvalue_fp64 - exp_pvalue) < 0.05);

    SUITE_END();
}

/* Same seeding contract as test_pcoa_seeded, for the permutation draw. Only the
 * p-value is at stake: the F statistic is computed from the data alone, so it
 * does not move with the seed.
 *
 * The tolerance is looser here, and deliberately so: a p-value is a rank within
 * the permutation distribution, over n_perm+1 values counting the unpermuted
 * one, so ULP drift in the observed F either leaves it alone or steps it by
 * 1/(n_perm+1).
 */
void test_permanova_seeded() {
    SUITE_START("test permanova seeded");

    // Same as skbio tests, as in test_permanova_noties
    const double matrix[] = {
      0., 1., 5., 4.,
      5., 0., 3., 2.,
      1., 3., 0., 3.,
      4., 2., 3., 0.};

    const uint32_t grouping[] = { 0, 0, 1, 1};
    const uint32_t n_samples  = 4;
    const unsigned int n_perm = 999;

    const double SAME = 1e-2;

    // an explicit seed reproduces, with no seeding call in between
    const double a = run_seeded_permanova(matrix, n_samples, grouping, n_perm, 7);
    const double b = run_seeded_permanova(matrix, n_samples, grouping, n_perm, 7);
    ASSERT(fabs(a - b) < SAME);

    // ... and is not perturbed by the global generator moving underneath it
    su::set_random_seed(999);
    const double c = run_seeded_permanova(matrix, n_samples, grouping, n_perm, 7);
    ASSERT(fabs(a - c) < SAME);

    /* The seed is really consumed. Checked across a spread of seeds rather than
     * against one alternative: a p-value is a rank out of n_perm+1, so any two
     * seeds can legitimately land on the same one. Requiring that *some* seed
     * disagrees keeps the assertion honest without making it a coin flip.
     */
    {
        unsigned int differing = 0;
        for (int s = 8; s <= 12; s++)
            if (fabs(run_seeded_permanova(matrix, n_samples, grouping, n_perm, s) - a) > SAME)
                differing++;
        ASSERT(differing > 0);
    }

    // a negative seed is the legacy path: same global seed, same answer
    su::set_random_seed(42);
    const double g1 = run_seeded_permanova(matrix, n_samples, grouping, n_perm, -1);
    su::set_random_seed(42);
    const double g2 = run_seeded_permanova(matrix, n_samples, grouping, n_perm, -1);
    ASSERT(fabs(g1 - g2) < SAME);

    SUITE_END();
}

void test_subsample_replacement() {
    SUITE_START("test subsample with replacement");

    su::biom org_table("test.biom");
    uint32_t exp_org_n_samples = 6;
    uint32_t exp_org_n_obs = 5;

    std::string sids_org[] = {"Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"};
    std::vector<std::string> exp_org_sids = _string_array_to_vector(sids_org, exp_org_n_samples);

    std::string oids_org[] = {"GG_OTU_1", "GG_OTU_2","GG_OTU_3", "GG_OTU_4", "GG_OTU_5"};
    std::vector<std::string> exp_org_oids = _string_array_to_vector(oids_org, exp_org_n_obs);

    su::skbio_biom_subsampled table(org_table,true,5);
    uint32_t exp_n_samples = 2;
    uint32_t exp_n_obs = 3;

    std::string sids[] = {"Sample1", "Sample4"};
    std::vector<std::string> exp_sids = _string_array_to_vector(sids, exp_n_samples);

    std::string oids[] = {"GG_OTU_2","GG_OTU_3", "GG_OTU_4"};
    std::vector<std::string> exp_oids = _string_array_to_vector(oids, exp_n_obs);

    // verify it did not destroy the original table
    ASSERTINTEQ(org_table.n_samples , exp_org_n_samples);
    ASSERTINTEQ(org_table.n_obs , exp_org_n_obs);
    ASSERT(org_table.get_sample_ids() == exp_org_sids);
    ASSERT(org_table.get_obs_ids() == exp_org_oids);

    // now check basic rarefaction properties
    ASSERTINTEQ(table.n_samples , exp_n_samples);
    ASSERTINTLE(table.n_obs , exp_n_obs);
    ASSERT(table.get_sample_ids() == exp_sids);

    // all columns should add to n==5
    {
      double exp_data[2] = {5.0, 5.0};
      std::vector<double> exp_vec = _double_array_to_vector(exp_data, 2);

      double data_sum[2] = {0.0, 0.0};
      for (auto obs_id : table.get_obs_ids()) {
         double line[2];
         table.get_obs_data(obs_id, line);
         for (int j=0; j<2; j++) data_sum[j] += line[j];
      }
      std::vector<double> sum_vec = _double_array_to_vector(data_sum, 2);
      ASSERT(sum_vec == exp_vec);
    }

    su::skbio_biom_subsampled table_empty(org_table,true,8);
    uint32_t exp_empty_n_samples = 0;
    uint32_t exp_empty_n_obs = 0;

    ASSERT(table_empty.n_samples == exp_empty_n_samples);
    ASSERT(table_empty.n_obs == exp_empty_n_obs);

    SUITE_END();
}

void test_subsample_replacement_limit() {
    SUITE_START("test subsample with replacement limit");

    // test case when n == sum(vector) - 1

    const char *t_obs_ids[] = {"OTU0","OTU1","OTU2","OTU3","OTU4",
                               "OTU5","OTU6","OTU7","OTU8","OTU9"};
    const char *t_samp_ids[] = {"S1"};
    uint32_t t_index[] = {0,0,0,0,0,0,0,0,0,0};
    uint32_t t_idxptr[] = {0,1,2,3,4,5,6,7,8,9,10}; 
    double t_data[] = {2., 1., 2., 1., 7., 6., 3., 3., 5., 5.};
    su::biom_inmem org_table(t_obs_ids,t_samp_ids,t_index,t_idxptr,t_data,10,1);

    const uint32_t exp_org_n_samples = 1;
    const uint32_t exp_org_n_obs = 10;

    // verify we did the right thing
    ASSERTINTEQ(org_table.n_samples , exp_org_n_samples);
    ASSERTINTEQ(org_table.n_obs , exp_org_n_obs);

    // column should add to n==35
    {
      double exp_data[1] = {35.0};
      std::vector<double> exp_vec = _double_array_to_vector(exp_data, 1);

      double data_sum[1] = {0.0};
      for (auto obs_id : org_table.get_obs_ids()) {
         double line[1];
         org_table.get_obs_data(obs_id, line);
         for (int j=0; j<1; j++) data_sum[j] += line[j];
      }
      std::vector<double> sum_vec = _double_array_to_vector(data_sum, 1);
      ASSERT(sum_vec == exp_vec);
    }

    su::skbio_biom_subsampled table(org_table,true,34);
    uint32_t exp_n_samples = 1;
    uint32_t exp_n_obs = 10;

    // now check basic rarefaction properties
    ASSERTINTEQ(table.n_samples , exp_n_samples);
    ASSERTINTLE(table.n_obs , exp_n_obs);

    // column should add to n==34
    {
      double exp_data[1] = {34.0};
      std::vector<double> exp_vec = _double_array_to_vector(exp_data, 1);

      double data_sum[1] = {0.0};
      for (auto obs_id : table.get_obs_ids()) {
         double line[1];
         table.get_obs_data(obs_id, line);
         for (int j=0; j<1; j++) data_sum[j] += line[j];
      }
      std::vector<double> sum_vec = _double_array_to_vector(data_sum, 1);
      ASSERT(sum_vec == exp_vec);
    }

    // As in scikit-bio and biom-format
    // when with replacement
    // the number of permutations should be more than n
    {
      std::unordered_set<uint64_t> perm_set;
      // will pick column 1
      for (int i=0; i<1000; i++) {
         su::skbio_biom_subsampled table2(org_table,true,34);

         uint64_t val = 0;
         for (auto obs_id : table2.get_obs_ids()) {
            double line[1];
            table2.get_obs_data(obs_id, line);
            val = val*8 + uint64_t(line[0]);
         }
         perm_set.insert(val);
      }
      ASSERT(uint32_t(perm_set.size()) > uint32_t(10));
    }

    SUITE_END();
}

void test_subsample_woreplacement() {
    SUITE_START("test subsample without replacement");

    su::biom org_table("test.biom");
    uint32_t exp_org_n_samples = 6;
    uint32_t exp_org_n_obs = 5;

    std::string sids_org[] = {"Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"};
    std::vector<std::string> exp_org_sids = _string_array_to_vector(sids_org, exp_org_n_samples);

    std::string oids_org[] = {"GG_OTU_1", "GG_OTU_2","GG_OTU_3", "GG_OTU_4", "GG_OTU_5"};
    std::vector<std::string> exp_org_oids = _string_array_to_vector(oids_org, exp_org_n_obs);

    su::skbio_biom_subsampled table(org_table,false,5);
    uint32_t exp_n_samples = 2;
    uint32_t exp_n_obs = 3;

    std::string sids[] = {"Sample1", "Sample4"};
    std::vector<std::string> exp_sids = _string_array_to_vector(sids, exp_n_samples);

    std::string oids[] = {"GG_OTU_2","GG_OTU_3", "GG_OTU_4"};
    std::vector<std::string> exp_oids = _string_array_to_vector(oids, exp_n_obs);

    // verify it did not destroy the original table
    ASSERTINTEQ(org_table.n_samples , exp_org_n_samples);
    ASSERTINTEQ(org_table.n_obs , exp_org_n_obs);
    ASSERT(org_table.get_sample_ids() == exp_org_sids);
    ASSERT(org_table.get_obs_ids() == exp_org_oids);

    // now check basic rarefaction properties
    ASSERTINTEQ(table.n_samples , exp_n_samples);
    ASSERTINTLE(table.n_obs , exp_n_obs);
    ASSERT(table.get_sample_ids() == exp_sids);

    // all columns should add to n==5
    {
      double exp_data[2] = {5.0, 5.0};
      std::vector<double> exp_vec = _double_array_to_vector(exp_data, 2);

      double data_sum[2] = {0.0, 0.0};
      for (auto obs_id : table.get_obs_ids()) {
         double line[2];
         table.get_obs_data(obs_id, line);
         for (int j=0; j<2; j++) data_sum[j] += line[j];
      }
      std::vector<double> sum_vec = _double_array_to_vector(data_sum, 2);
      ASSERT(sum_vec == exp_vec);
    }

    su::skbio_biom_subsampled table_empty(org_table,false,8);
    uint32_t exp_empty_n_samples = 0;
    uint32_t exp_empty_n_obs = 0;

    ASSERT(table_empty.n_samples == exp_empty_n_samples);
    ASSERT(table_empty.n_obs == exp_empty_n_obs);

    SUITE_END();
}

void test_subsample_woreplacement_limit() {
    SUITE_START("test subsample without replacement limit");

    // test case when n == sum(vector) - 1

    const char *t_obs_ids[] = {"OTU0","OTU1","OTU2","OTU3","OTU4",
                               "OTU5","OTU6","OTU7","OTU8","OTU9"};
    const char *t_samp_ids[] = {"S1"};
    uint32_t t_index[] = {0,0,0,0,0,0,0,0,0,0};
    uint32_t t_idxptr[] = {0,1,2,3,4,5,6,7,8,9,10}; 
    double t_data[] = {2., 1., 2., 1., 7., 6., 3., 3., 5., 5.};
    su::biom_inmem org_table(t_obs_ids,t_samp_ids,t_index,t_idxptr,t_data,10,1);

    const uint32_t exp_org_n_samples = 1;
    const uint32_t exp_org_n_obs = 10;

    // verify we did the right thing
    ASSERTINTEQ(org_table.n_samples , exp_org_n_samples);
    ASSERTINTEQ(org_table.n_obs , exp_org_n_obs);

    // column should add to n==35
    {
      double exp_data[1] = {35.0};
      std::vector<double> exp_vec = _double_array_to_vector(exp_data, 1);

      double data_sum[1] = {0.0};
      for (auto obs_id : org_table.get_obs_ids()) {
         double line[1];
         org_table.get_obs_data(obs_id, line);
         for (int j=0; j<1; j++) data_sum[j] += line[j];
      }
      std::vector<double> sum_vec = _double_array_to_vector(data_sum, 1);
      ASSERT(sum_vec == exp_vec);
    }

    su::skbio_biom_subsampled table(org_table,false,34);
    uint32_t exp_n_samples = 1;
    uint32_t exp_n_obs = 10;

    // now check basic rarefaction properties
    ASSERTINTEQ(table.n_samples , exp_n_samples);
    ASSERTINTLE(table.n_obs , exp_n_obs);

    // column should add to n==34
    {
      double exp_data[1] = {34.0};
      std::vector<double> exp_vec = _double_array_to_vector(exp_data, 1);

      double data_sum[1] = {0.0};
      for (auto obs_id : table.get_obs_ids()) {
         double line[1];
         table.get_obs_data(obs_id, line);
         for (int j=0; j<1; j++) data_sum[j] += line[j];
      }
      std::vector<double> sum_vec = _double_array_to_vector(data_sum, 1);
      ASSERT(sum_vec == exp_vec);
    }

    // As in scikit-bio and biom-format
    // when without replacement
    // the number of permutations should be exactly n
    {
      std::unordered_set<uint64_t> perm_set;
      // will pick column 1
      for (int i=0; i<1000; i++) {
         su::skbio_biom_subsampled table2(org_table,false,34);

         uint64_t val = 0;
         for (auto obs_id : table2.get_obs_ids()) {
            double line[1];
            table2.get_obs_data(obs_id, line);
            val = val*8 + uint64_t(line[0]);
         }
         perm_set.insert(val);
      }
      ASSERT(uint32_t(perm_set.size()) == uint32_t(10));
    }

    SUITE_END();
}

int main(int argc, char** argv) {
    test_center_mat();
    test_pcoa();
    // The following test is unstable, disable for now
    // test_pcoa_big();
    test_pcoa_seeded();
    test_permanova_ties();
    test_permanova_noties();
    test_permanova_unequal();
    test_permanova_seeded();
    test_subsample_replacement();
    test_subsample_replacement_limit();
    test_subsample_woreplacement();
    test_subsample_woreplacement_limit();

    printf("\n");
    printf(" %i / %i suites failed\n", suites_failed, suites_run);
    printf(" %i / %i suites empty\n", suites_empty, suites_run);
    printf(" %i / %i tests failed\n", tests_failed, tests_run);
  
    printf("\n THE END.\n");
    
    return tests_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
