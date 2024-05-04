#include <iostream>
#include <omp.h>

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