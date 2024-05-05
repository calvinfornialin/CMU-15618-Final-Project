#include <iostream>
#include <fstream>
#include <numeric>
#include "compute_system.h"
#include "task_runner.h"

int main() {

    // initialize the underlying hardware system including 
    // cores number and cache topology
    ComputeSystem system(3, 8);

    system.addCache(1, 32 * 1024, 64, {0});     // core 0's private L1 cache 
    system.addCache(1, 32 * 1024, 64, {1});     // core 1's private L1 cache 
    system.addCache(1, 32 * 1024, 64, {2});     // core 2's private L1 cache 
    system.addCache(1, 32 * 1024, 64, {3});     // core 3's private L1 cache 
    system.addCache(1, 32 * 1024, 64, {4});     // core 4's private L1 cache 
    system.addCache(1, 32 * 1024, 64, {5});     // core 5's private L1 cache 
    system.addCache(1, 32 * 1024, 64, {6});     // core 6's private L1 cache 
    system.addCache(1, 32 * 1024, 64, {7});     // core 7's private L1 cache 
    system.addCache(2, 256 * 1024, 64, {0, 1}); // L2 cache shared by cores 0 and 1
    system.addCache(2, 256 * 1024, 64, {2, 3}); // L2 cache shared by cores 2 and 3
    system.addCache(2, 256 * 1024, 64, {4, 5}); // L2 cache shared by cores 4 and 5
    system.addCache(2, 256 * 1024, 64, {6, 7}); // L2 cache shared by cores 6 and 7
    system.addCache(3, 2048 * 1024, 64, {0, 1, 2, 3, 4, 5, 6, 7}); // L3 cache shared by all cores

    /*
    int cacheLevels, coreCount;
    int cacheLevel, cacheSize, blockSize, sharingCoreID;

    ifstream inputFile("config1.txt");
    if (!inputFile) {
        cerr << "Error opening file!" << endl;
        return 1;
    }

    inputFile >> cacheLevels >> coreCount;
    ComputeSystem system(cacheLevels, coreCount);

    while (inputFile >> cacheLevel >> cacheSize >> blockSize) {
        vector<int> sharingIDs;
        while (inputFile >> sharingCoreID) {
            sharingIDs.push_back(sharingCoreID);
        }
        system.addCache(cacheLevel, cacheSize * 1024, blockSize, sharingIDs);
    }
    */

    system.display();
    vector<vector<set<int>>> simplifiedCacheHierarchy = system.getSimplifiedCacheHierarchy();


    // Iterate over each level of the hierarchy
    for (int level = 0; level < simplifiedCacheHierarchy.size(); ++level) {
        cout << "Level " << level + 1 << ":" << endl;
        
        // Iterate over each set in the current level
        for (const auto& cacheSet : simplifiedCacheHierarchy[level]) {
            cout << "Set: { ";
            
            // Print the elements of the set
            for (const auto& id : cacheSet) {
                cout << id << " ";
            }
            
            cout << "}" << endl;
        }
    }

    int coreCount = system.getCoreCount();
    int cacheCount = system.getCacheCount();
    int numIters = 1000;

    cout << "coreCount: " << coreCount << endl;
    cout << "cacheCount: " << cacheCount << endl;

    // SimpleTaskRunner simpletaskrunner(numIters, coreCount, cacheCount, simplifiedCacheHierarchy);
    // simpletaskrunner.runTaskRefinedHybrid();
    // simpletaskrunner.runTaskHybrid();
    // simpletaskrunner.runTaskStatic();
    // simpletaskrunner.runTaskDynamic();


    // LINPACKTaskRunner linpacktaskrunner(numIters, coreCount, cacheCount, simplifiedCacheHierarchy);
    // linpacktaskrunner.runTaskRefinedHybrid();
    // linpacktaskrunner.runTaskHybrid();
    // linpacktaskrunner.runTaskStatic();
    // linpacktaskrunner.runTaskDynamic();

    vector<double> execution_times;

    // Run task 100 times
    for (int i = 0; i < 100; ++i) {
        LINPACKTaskRunner linpacktaskrunner(numIters, coreCount, cacheCount, simplifiedCacheHierarchy);
        double execution_time = linpacktaskrunner.runTaskRefinedHybrid();
        // double execution_time = linpacktaskrunner.runTaskHybrid();
        // double execution_time = linpacktaskrunner.runTaskStatic();
        // double execution_time = linpacktaskrunner.runTaskDynamic();
        execution_times.push_back(execution_time);
    }

    // Calculate average execution time
    double average_execution_time = accumulate(execution_times.begin(), execution_times.end(), 0.0) / execution_times.size();

    cout << "Average execution time: " << average_execution_time << " seconds" << endl;

    return 0;
}
