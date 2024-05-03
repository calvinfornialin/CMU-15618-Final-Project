#include <iostream>
#include "compute_system.h"

int main() {
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
    
    system.display();

    return 0;
}
