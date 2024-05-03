#include <iostream>
#include <vector>

using namespace std;
 
class Cache {
    private:
        // level of cache
        int level;

        // cache size in bytes
        int size;           

        // associativity of the cache
        // int associativity;  

        // block size in bytes
        int blockSize;      

        // number of cores sharing this cache
        int numSharingCore;

        // core ID that share this cache
        vector<int> sharingCoreID;

    public:
        // constructor
        Cache(int l, int s, int b, const vector<int>& sharingIDs) 
            : level(l), size(s), blockSize(b), sharingCoreID(sharingIDs) {
            numSharingCore = sharingCoreID.size();
        }

        // display cache info
        void display() {
            cout << "Cache details:" << endl;
            cout << "Cache level: L" << level << endl;
            cout << "Cache size: " << size << " bytes" << endl;
            cout << "Block size: " << blockSize << " bytes" << endl;
            cout << "Number of cores sharing: " << numSharingCore << endl;
            cout << "Core ID(s) sharing this cache: ";
            for (int id : sharingCoreID) {
                cout << id << " ";
            }
            cout << endl;
        }
};


class ComputeSystem {

    private:
        // number of cache levels
        int cacheLevels; 
        
        // total number of cores
        int coreCount;

        // number of cores sharing each cache level
        vector<int> coresSharingCache; 

        vector<Cache> cacheHierarchy;

    public:
        // constructor
        ComputeSystem(int l, int c) : cacheLevels(l), coreCount(c) {}

        // add a cache to the system
        void addCache(int level, int size, int blockSize, const vector<int>& sharingIDs) {
            Cache newCache(level, size, blockSize, sharingIDs);
            cacheHierarchy.push_back(newCache);
        }

        // Display compute system info
        void display() {
            cout << "Compute system details:" << endl;
            cout << "There are " << coreCount << " cores in the system" << endl;
            cout << "There are " << cacheLevels << " cache levels in the hierarchy" << endl;
            cout << "Each level cache info:" << endl;
            for (auto& cache : cacheHierarchy) {
                cache.display();
                cout << endl;
            }
            cout << endl;
        }
        
};
 

