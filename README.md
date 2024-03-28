# CMU-15618-Final-Project
## Topic: Design and Implementation Of A Cache Hierarchy-Aware Task Scheduling On Multicore Architectures ##
URL: https://github.com/calvinfornialin/15618-Final-Project

Summary

This project aims to design and implement a cache hierarchy-aware task scheduling system for multicore architectures. The goal is to optimize task execution on parallel systems by leveraging knowledge of the cache hierarchy to reduce cache misses and improve cache utilization and overall system performance.

Background

Many modern applications, particularly those in scientific computing, data analytics, and machine learning, are compute-intensive and can benefit significantly from parallel execution. As learned in lectures, optimal performance on multicore systems requires careful consideration of how tasks are scheduled among multiple cores or processors. Several crucial factors that we should consider are data locality, workload balance, and synchronization/communication time. However, the static or dynamic scheduling approaches introduced in lectures have neglected the cache hierarchy factor. So, in this project, we will propose a novel scheduling approach that takes not only data locality and workload balance but also cache hierarchy into account, and focus on a compute-intensive component that involves frequent data accesses, such as matrix multiplication, where parallelism can be exploited through task partitioning and scheduling based on cache topology and behavior.

Challenges
Correctly reflecting what we learned about task scheduling from lecture in the actual implementation requires a firm understanding of crucial metrics and trade-offs between different scheduling approaches.
Although tasks such as matrix multiplication may seem easy to parallelize at first sight, the challenge arises when dealing with large matrices. How to decrease the latencies due to cache misses based on the cache hierarchy while reducing the dependencies between calculations and maintaining workload balance will be a huge challenge.
Understanding each CADSS module and how they are interconnected and interact with each other.
Dependencies: Some tasks might depend on the completion of others, which requires careful scheduling to avoid deadlocks and minimize synchronization/communication times.
System Constraints: Mapping tasks effectively to a multicore architecture, considering the data locality and workload balance, is already challenging; considering cache hierarchy and core interconnect bandwidth is even more difficult.
Resources

Machine used: We plan to develop the scheduling approach on local machines and then run tests and benchmarks on the GHC and PSC machines.
Starting framework: CADSS: We'll start implementing the computing system from scratch, building upon the modules exposed by CADSS.
Reference paper: Nader Khammassi, Jean-Christophe Le Lann, "Design and Implementation of a Cache Hierarchy-Aware Task Scheduling for Parallel Loops on Multicore Architectures".
Goals and Deliverables

Plan to Achieve: Implement a cache hierarchy-aware task scheduler that considers cache hierarchy for task execution. Demonstrate improved cache hit rates and reduced execution times for selected compute-intensive applications compared to several standard schedulers.
Hope to Achieve: If ahead of schedule, we’d like to explore factors beside cache hierarchy that might also affect task scheduling performance metrics and include that in our scheduler.
Demo

An interactive demonstration showcasing the performance improvements in task execution times and cache utilization compared to non-optimized scheduling. Speedup graphs and cache hit/miss statistics will be presented.

Analysis Goals

Understand how task granularity influences task performance metrics.
Understand how different cache hierarchy levels influence task scheduling performance metrics and to identify the optimal task granularity for maximizing cache utilization.
Understand how CPU core number influences task scheduling performance metrics and to identify the optimal task granularity for maximizing cache utilization.
Platform Choice

We will be doing most of our development and testing locally and on the GHC machines. For the scalability study, we will be using PSC machine.

Schedule

Week 1: Study CADSS’s basic modules code and how the simulation platform functions.
Week 2: Design and revise the algorithm of cache hierarchy-aware scheduler.
Week 3: Implementation of the scheduler.
Week 4: Perform analysis and gather data among several schedulers.
Week 5: Work on report and extending scheduler to be dynamically adjustable.
