#include <iostream>
#include <omp.h>
#include <cmath>
#include <iomanip>

class SimpleTaskRunner {
    private:
        // number of iterations the task should be run
        int iterCount;

        // number of cores used to run the task
        int coreCount;

        // number of caches in the system used to run the task
        int cacheCount;

        vector<vector<set<int>>> simplifiedCacheHierarchy;

        // keep track of the iteration cutoff points for static and dynamic tasks
        vector<int> iterations;
    
    public:
        SimpleTaskRunner(int numIters, int numCores, int numCaches, const vector<vector<set<int>>>& cacheHierarchy)
            : iterCount(numIters), coreCount(numCores), cacheCount(numCaches), simplifiedCacheHierarchy(cacheHierarchy) {
            
            for (int i = numCores; i <= numCaches ; i++) {
                iterations.push_back(iterCount * i / numCaches);
            }
        }

        void runTaskHybrid() {
            // Set number of workers using core count
            omp_set_num_threads(coreCount);

            // start time measurement of running the task
            double start_time = omp_get_wtime();

            #pragma omp parallel
            {
                int thread_id = omp_get_thread_num();

                // Static scheduling for each of the cores
                #pragma omp for schedule(static) nowait
                for (int i = 0; i < iterations[0]; ++i) {
                    std::cout << "Thread " << thread_id << " is processing static iteration " << i << std::endl;
                }

                // Iterate through each set of each level for dynamic scheduling
                int cur_iter = 1;
                for (int level = 1; level < simplifiedCacheHierarchy.size(); ++level) {
                    int tmp_cur_iter = cur_iter;
                    for (const auto& cacheSet : simplifiedCacheHierarchy[level]) {
                        // check if thread id is in the set
                        if (cacheSet.find(thread_id) != cacheSet.end()) {
                            #pragma omp for schedule(dynamic) nowait
                            for (int i = iterations[tmp_cur_iter - 1]; i < iterations[tmp_cur_iter]; ++i) {
                                std::cout << "Thread " << thread_id << " is processing dynamic iteration " << i << std::endl;
                            }
                            break;
                        } else {
                            tmp_cur_iter += 1;
                        }
                    }
                    cur_iter += simplifiedCacheHierarchy[level].size();
                }
            }

            // end time measurement of running the task
            double end_time = omp_get_wtime();
            double total_time = end_time - start_time;

            std::cout << "Total execution time: " << total_time << " seconds" << std::endl;
        }

        void runTaskStatic() {
            // Set number of workers using core count
            omp_set_num_threads(coreCount);

            // start time measurement of running the task
            double start_time = omp_get_wtime();

            #pragma omp parallel
            {
                int thread_id = omp_get_thread_num();

                // Static scheduling for all of the cores
                #pragma omp for schedule(static) nowait
                for (int i = 0; i < iterations.back(); ++i) {
                    std::cout << "Thread " << thread_id << " is processing static iteration " << i << std::endl;
                }
            }

            // end time measurement of running the task
            double end_time = omp_get_wtime();
            double total_time = end_time - start_time;

            std::cout << "Total execution time: " << total_time << " seconds" << std::endl;
        }

        void runTaskDynamic() {
            // Set number of workers using core count
            omp_set_num_threads(coreCount);

            // start time measurement of running the task
            double start_time = omp_get_wtime();

            #pragma omp parallel
            {
                int thread_id = omp_get_thread_num();

                // Static scheduling for all of the cores
                #pragma omp for schedule(dynamic) nowait
                for (int i = 0; i < iterations.back(); ++i) {
                    std::cout << "Thread " << thread_id << " is processing dynamic iteration " << i << std::endl;
                }
            }

            // end time measurement of running the task
            double end_time = omp_get_wtime();
            double total_time = end_time - start_time;

            std::cout << "Total execution time: " << total_time << " seconds" << std::endl;
        }
};

class LINPACKTaskRunner {
    private:
        // number of iterations the task should be run
        int iterCount;

        // number of cores used to run the task
        int coreCount;

        // number of caches in the system used to run the task
        int cacheCount;

        vector<vector<set<int>>> simplifiedCacheHierarchy;

        // variables used for LINPACK solver
        float *a;
        float *b;
        float *x;
        int lda;
        int *ipvt;
        int info;
        int n;
        float err;
        int job;
    
    public:
        LINPACKTaskRunner(int numIters, int numCores, int numCaches, const vector<vector<set<int>>>& cacheHierarchy)
            : iterCount(numIters), coreCount(numCores), cacheCount(numCaches), simplifiedCacheHierarchy(cacheHierarchy) {
            
            n = iterCount;
            lda = n;
            a = new float[lda * n];
            b = new float[n];
            x = new float[n];

            matgen ( lda, n, a, x, b );

            ipvt = new int[n];
        }

        void runTaskHybrid() {
            // Set number of workers using core count
            omp_set_num_threads(coreCount);

            // start time measurement of running the task
            double start_time = omp_get_wtime();

            int k,kp1,l,nm1;
            float t;

            info = 0;
            nm1 = n - 1;
            for ( k = 0; k < nm1; k++ )
            {
                kp1 = k + 1;
                l = isamax ( n-k, a+k+k*lda, 1 ) + k - 1;
                ipvt[k] = l + 1;

                if ( a[l+k*lda] == 0.0 )
                {
                info = k + 1;
                break;
                }

                if ( l != k )
                {
                t          = a[l+k*lda];
                a[l+k*lda] = a[k+k*lda];
                a[k+k*lda] = t;
                }
                t = -1.0 / a[k+k*lda]; 
                sscal ( n-k-1, t, a+kp1+k*lda, 1 );
            //
            //  Interchange the pivot row and the K-th row.
            //
                if ( l != k )
                {
                sswap ( n-k-1, a+l+kp1*lda, lda, a+k+kp1*lda, lda );
                }
            //
            //  Add multiples of the K-th row to rows K+1 through N.
            //
                int i,j;
                // keep track of the iteration cutoff points for static and dynamic tasks
                vector<int> iterations;

                for (int i = coreCount; i <= cacheCount ; i++) {
                    iterations.push_back(n-k-1 * i / cacheCount);
                }

                // # pragma omp parallel \
                // shared ( a+k+kp1*lda, n-k-1, n-k-1, a+kp1+k*lda, a+kp1+kp1*lda) \
                // private ( i, j ) 

                // int thread_id = omp_get_thread_num();
                
                // // Static scheduling for each of the cores
                // #pragma omp for schedule(static) nowait
                // for (int j = 0; j < iterations[0]; ++j) {
                //     for ( int i = 0; i < n-k-1; i++ )
                //     {
                //     (a+kp1+kp1*lda)[i+j*n] += (a+k+kp1*lda)[j*n] * (a+kp1+k*lda)[i];
                //     }
                // }

                // // Iterate through each set of each level for dynamic scheduling
                // int cur_iter = 1;
                // for (int level = 1; level < simplifiedCacheHierarchy.size(); ++level) {
                //     int tmp_cur_iter = cur_iter;
                //     for (const auto& cacheSet : simplifiedCacheHierarchy[level]) {
                //         // check if thread id is in the set
                //         if (cacheSet.find(thread_id) != cacheSet.end()) {
                //             #pragma omp for schedule(dynamic) nowait
                //             for (int j = iterations[tmp_cur_iter - 1]; j < iterations[tmp_cur_iter]; ++j) {
                //                 for ( int i = 0; i < n-k-1; i++ )
                //                 {
                //                 (a+kp1+kp1*lda)[i+j*n] += (a+k+kp1*lda)[j*n] * (a+kp1+k*lda)[i];
                //                 }
                //             }
                //             break;
                //         } else {
                //             tmp_cur_iter += 1;
                //         }
                //     }
                //     cur_iter += simplifiedCacheHierarchy[level].size();
                // }

                #pragma omp parallel shared(a, k, kp1, lda, n, simplifiedCacheHierarchy, iterations)
                {
                    int thread_id = omp_get_thread_num();
                    
                    // Static scheduling for each of the cores
                    #pragma omp for schedule(static) nowait
                    for (int j = 0; j < iterations[0]; ++j) {
                        for (int i = 0; i < n - k - 1; i++) {
                            (a+kp1+kp1*lda)[i+j*n] += (a+k+kp1*lda)[j*n] * (a+kp1+k*lda)[i];
                        }
                    }

                    // Iterate through each set of each level for dynamic scheduling
                    int cur_iter = 1;
                    for (int level = 1; level < simplifiedCacheHierarchy.size(); ++level) {
                        int tmp_cur_iter = cur_iter;
                        for (const auto& cacheSet : simplifiedCacheHierarchy[level]) {
                            // check if thread id is in the set
                            if (cacheSet.find(thread_id) != cacheSet.end()) {
                                #pragma omp for schedule(dynamic) nowait
                                for (int j = iterations[tmp_cur_iter - 1]; j < iterations[tmp_cur_iter]; ++j) {
                                    for (int i = 0; i < n - k - 1; i++) {
                                        (a+kp1+kp1*lda)[i+j*n] += (a+k+kp1*lda)[j*n] * (a+kp1+k*lda)[i];
                                    }
                                }
                                break;
                            } else {
                                tmp_cur_iter += 1;
                            }
                        }
                        cur_iter += simplifiedCacheHierarchy[level].size();
                    }
                }
            }

            ipvt[n-1] = n;

            if ( a[n-1+(n-1)*lda] == 0.0 )
            {
                info = n;
            }

            // end time measurement of running the task
            double end_time = omp_get_wtime();
            double total_time = end_time - start_time;

            std::cout << "Total execution time: " << total_time << " seconds" << std::endl;

            if ( info != 0 )
            {
                cout << "\n";
                cout << "TEST02 - Fatal error!\n";
                cout << "  MSGEFA reports the matrix is singular.\n";
                exit ( 1 );
            }
            //
            //  Solve the linear system.
            //
            job = 0;
            sgesl ( a, lda, n, ipvt, b, job );

            err = 0.0;
            for ( int i = 0; i < n; i++ )
            {
                err = err + fabs ( x[i] - b[i] );
            }

            cout << "N: " << n << endl;
            cout << "Error: " << err << endl;
            
            delete [] a;
            delete [] b;
            delete [] ipvt;
            delete [] x;

            return;
        }

        // MATGEN generates a "random" matrix for testing.
        void matgen ( int lda, int n, float a[], float x[], float b[] ) {
            int i;
            int j;
            int seed;
            float value;

            seed = 1325;
            //
            //  Set the matrix A.
            //
            for ( j = 0; j < n; j++ )
            {
                for ( i = 0; i < n; i++ )
                {
                seed = ( 3125 * seed ) % 65536;
                value = ( ( float ) seed - 32768.0 ) / 16384.0;
                a[i+j*lda] = value;
                }
            }
            //
            //  Set x.
            //
            for ( i = 0; i < n; i++ )
            {
                x[i] = ( float ) ( i + 1 ) / ( ( float ) n );
            }
            //
            //  Set b = A * x.
            //
            for ( i = 0; i < n; i++ ) 
            {
                b[i] = 0.0;
                for ( j = 0; j < n; j++ )
                {
                b[i] = b[i] + a[i+j*lda] * x[j];
                }
            }
            return;
        }

        //    ISAMAX finds the index of the vector element of maximum absolute value.
        int isamax ( int n, float x[], int incx ) {
            float xmax;
            int i;
            int ix;
            int value;

            value = 0;

            if ( n < 1 || incx <= 0 )
            {
                return value;
            }

            value = 1;

            if ( n == 1 )
            {
                return value;
            }

            if ( incx == 1 )
            {
                xmax = fabs ( x[0] );

                for ( i = 1; i < n; i++ )
                {
                if ( xmax < fabs ( x[i] ) )
                {
                    value = i + 1;
                    xmax = fabs ( x[i] );
                }
                }
            }
            else
            {
                ix = 0;
                xmax = fabs ( x[0] );
                ix = ix + incx;

                for ( i = 1; i < n; i++ )
                {
                if ( xmax < fabs ( x[ix] ) )
                {
                    value = i + 1;
                    xmax = fabs ( x[ix] );
                }
                ix = ix + incx;
                }
            }

            return value;
        }

        //    SSCAL scales a float vector by a constant.
        void sscal ( int n, float sa, float x[], int incx ) {
            int i;
            int ix;
            int m;

            if ( n <= 0 )
            {
            }
            else if ( incx == 1 )
            {
                m = n % 5;

                for ( i = 0; i < m; i++ )
                {
                x[i] = sa * x[i];
                }

                for ( i = m; i < n; i = i + 5 )
                {
                x[i]   = sa * x[i];
                x[i+1] = sa * x[i+1];
                x[i+2] = sa * x[i+2];
                x[i+3] = sa * x[i+3];
                x[i+4] = sa * x[i+4];
                }
            }
            else
            {
                if ( 0 <= incx )
                {
                ix = 0;
                }
                else
                {
                ix = ( - n + 1 ) * incx;
                }

                for ( i = 0; i < n; i++ )
                {
                x[ix] = sa * x[ix];
                ix = ix + incx;
                }

            }

            return;
        }

        //    SSWAP interchanges two float vectors.
        void sswap ( int n, float x[], int incx, float y[], int incy ) {
            int i;
            int ix;
            int iy;
            int m;
            float temp;

            if ( n <= 0 )
            {
            }
            else if ( incx == 1 && incy == 1 )
            {
                m = n % 3;

                for ( i = 0; i < m; i++ )
                {
                temp = x[i];
                x[i] = y[i];
                y[i] = temp;
                }

                for ( i = m; i < n; i = i + 3 )
                {
                temp = x[i];
                x[i] = y[i];
                y[i] = temp;

                temp = x[i+1];
                x[i+1] = y[i+1];
                y[i+1] = temp;

                temp = x[i+2];
                x[i+2] = y[i+2];
                y[i+2] = temp;
                }
            }
            else
            {
                if ( 0 <= incx )
                {
                ix = 0;
                }
                else
                {
                ix = ( - n + 1 ) * incx;
                }

                if ( 0 <= incy )
                {
                iy = 0;
                }
                else
                {
                iy = ( - n + 1 ) * incy;
                }

                for ( i = 0; i < n; i++ )
                {
                temp = x[ix];
                x[ix] = y[iy];
                y[iy] = temp;
                ix = ix + incx;
                iy = iy + incy;
                }
            }
            return;
        }

        //    SGESL solves a real general linear system A * X = B.
        void sgesl ( float a[], int lda, int n, int ipvt[], float b[], int job ) {
            int k;
            int l;
            float t;
            //
            //  Solve A * X = B.
            //

            if ( job == 0 )
            {
                for ( k = 1; k <= n-1; k++ )
                {
                l = ipvt[k-1];
                t = b[l-1];

                if ( l != k )
                {
                    b[l-1] = b[k-1];
                    b[k-1] = t;
                }
                saxpy ( n-k, t, a+k+(k-1)*lda, 1, b+k, 1 );
                }

                for ( k = n; 1 <= k; k-- )
                {
                b[k-1] = b[k-1] / a[k-1+(k-1)*lda];
                t = -b[k-1];
                saxpy ( k-1, t, a+0+(k-1)*lda, 1, b, 1 );
                }
            }
            //
            //  Solve A' * X = B.
            //
            else
            {
                for ( k = 1; k <= n; k++ )
                {
                t = sdot ( k-1, a+0+(k-1)*lda, 1, b, 1 );
                b[k-1] = ( b[k-1] - t ) / a[k-1+(k-1)*lda];
                }

                for ( k = n-1; 1 <= k; k-- )
                {
                b[k-1] = b[k-1] + sdot ( n-k, a+k+(k-1)*lda, 1, b+k, 1 );
                l = ipvt[k-1];

                if ( l != k )
                {
                    t = b[l-1];
                    b[l-1] = b[k-1];
                    b[k-1] = t;
                }
                }
            }
            return;
        }

        //    SAXPY computes float constant times a vector plus a vector.
        void saxpy ( int n, float a, float x[], int incx, float y[], int incy ) {
            int i;
            int ix;
            int iy;
            int m;

            if ( n <= 0 )
            {
                return;
            }

            if ( a == 0.0 )
            {
                return;
            }
            //
            //  Code for unequal increments or equal increments
            //  not equal to 1.
            //
            if ( incx != 1 || incy != 1 )
            {
                if ( 0 <= incx )
                {
                ix = 0;
                }
                else
                {
                ix = ( - n + 1 ) * incx;
                }

                if ( 0 <= incy )
                {
                iy = 0;
                }
                else
                {
                iy = ( - n + 1 ) * incy;
                }

                for ( i = 0; i < n; i++ )
                {
                y[iy] = y[iy] + a * x[ix];
                ix = ix + incx;
                iy = iy + incy;
                }
            }
            //
            //  Code for both increments equal to 1.
            //
            else
            {
                m = n % 4;

                for ( i = 0; i < m; i++ )
                {
                y[i] = y[i] + a * x[i];
                }

                for ( i = m; i < n; i = i + 4 )
                {
                y[i  ] = y[i  ] + a * x[i  ];
                y[i+1] = y[i+1] + a * x[i+1];
                y[i+2] = y[i+2] + a * x[i+2];
                y[i+3] = y[i+3] + a * x[i+3];
                }
            }

            return;
        }

        //    SDOT forms the dot product of two vectors.
        float sdot ( int n, float x[], int incx, float y[], int incy ) {
            int i;
            int ix;
            int iy;
            int m;
            float temp;

            temp = 0.0;

            if ( n <= 0 )
            {
                return temp;
            }
            //
            //  Code for unequal increments or equal increments
            //  not equal to 1.
            //
            if ( incx != 1 || incy != 1 )
            {
                if ( 0 <= incx )
                {
                ix = 0;
                }
                else
                {
                ix = ( - n + 1 ) * incx;
                }

                if ( 0 <= incy )
                {
                iy = 0;
                }
                else
                {
                iy = ( - n + 1 ) * incy;
                }

                for ( i = 0; i < n; i++ )
                {
                temp = temp + x[ix] * y[iy];
                ix = ix + incx;
                iy = iy + incy;
                }
            }
            //
            //  Code for both increments equal to 1.
            //
            else
            {
                m = n % 5;

                for ( i = 0; i < m; i++ )
                {
                temp = temp + x[i] * y[i];
                }

                for ( i = m; i < n; i = i + 5 )
                {
                temp = temp + x[i  ] * y[i  ] 
                            + x[i+1] * y[i+1] 
                            + x[i+2] * y[i+2] 
                            + x[i+3] * y[i+3] 
                            + x[i+4] * y[i+4];
                }
            }

            return temp;
        }
};