#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>

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

        bool privateCache;

        // core ID that share this cache
        vector<int> sharingCoreID;

    public:
        // constructor
        Cache(int l, int s, int b, const vector<int>& sharingIDs) 
            : level(l), size(s), blockSize(b), sharingCoreID(sharingIDs) {
            numSharingCore = sharingCoreID.size();
            privateCache = (sharingCoreID.size() == 1)? true : false;
        }


        // return the level of the cache
        int getCacheLevel() const{
            return level;
        }


        // return the privacy of the cache
        bool getCachePrivacy() const{
            return privateCache;
        }


        vector<int> getsharingCoreID() const {
            return sharingCoreID;
        }


        // display cache info
        void display() {
            cout << "Cache details:" << endl;
            cout << "Cache level: L" << level << endl;
            cout << "Cache size: " << size << " bytes" << endl;
            cout << "Block size: " << blockSize << " bytes" << endl;
            cout << "Number of cores sharing: " << numSharingCore << endl;
            if (privateCache) {
                cout << "It's a private cache to: ";
                for (int id : sharingCoreID) {
                    cout << id << " ";
                }
                cout << endl;
            }
            else {
                cout << "It's shared by Core ID(s): ";
                for (int id : sharingCoreID) {
                    cout << id << " ";
                }
                cout << endl;
            }
        }
};


class ComputeSystem {

    private:
        // number of cache levels
        int cacheLevels; 
        
        // total number of cores
        int coreCount;

        // total number of caches
        int cacheCount;

        // number of cores sharing each cache level
        vector<int> coresSharingCache; 

        vector<Cache> cacheHierarchy;

        // a simplified version for easier access of cahce hierarchy
        vector<vector<set<int>>> simplifiedCacheHierarchy;

    public:
        // constructor
        ComputeSystem(int l, int c) : cacheLevels(l), coreCount(c) {
            simplifiedCacheHierarchy.resize(cacheLevels);
            cacheCount = 0;
        }


        // add a cache to the system
        void addCache(int level, int size, int blockSize, const vector<int>& sharingIDs) {
            Cache newCache(level, size, blockSize, sharingIDs);
            cacheHierarchy.push_back(newCache);
            cacheCount += 1;
            
            // also update the simplifed cache hierarchy
            set<int> sharingIDsSet(sharingIDs.begin(), sharingIDs.end());
            simplifiedCacheHierarchy[level - 1].push_back(sharingIDsSet);
        }


        // return number of shared caches at all levels
        /*int numSharedCaches(const vector<Cache>& cacheHierarchy) {
            int numSharedCaches = 0;
            for (auto& cache : cacheHierarchy) {
                // case it's a shared cache
                if (!cache.getCachePrivacy()) {
                    // case it's L1 cache
                    if (cache.getCacheLevel() == 1) {
                        numSharedCaches++;
                    } 
                    else {
                        // check if there's other higher level cache
                        // that's also shared and are shared by the same
                        // set of cores

                    }
                }
            }
        }*/
        int numSharedCaches(const vector<Cache>& cacheHierarchy) {
            map<set<int>, int> sharedCacheCount;  // Maps sets of core IDs to the count of shared caches.

            for (const auto& cache : cacheHierarchy) {
                set<int> coreSet(cache.getsharingCoreID().begin(), cache.getsharingCoreID().end());
                
                // only consider shared caches
                if (coreSet.size() > 1) {  
                    // increments the count for each unique core set
                    sharedCacheCount[coreSet]++;
                }
            }

            // each unique core set that has more than one cache sharing it counts as one
            int uniqueSharedCaches = 0;
            for (const auto& entry : sharedCacheCount) {
                // only count if there's at least one shared cache
                if (entry.second > 0) {  
                    uniqueSharedCaches++;
                }
            }

            return uniqueSharedCaches;
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

        vector<vector<set<int>>> getSimplifiedCacheHierarchy() const {
            return simplifiedCacheHierarchy;
        }

        int getCacheCount() const {
            return cacheCount;
        }

        int getCoreCount() const {
            return coreCount;
        }
};
 

