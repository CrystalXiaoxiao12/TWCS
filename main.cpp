#include <iostream>
#include <vector>
#include <algorithm>
#include <queue>
#include <chrono>
#include <csignal>
#include <unistd.h>
#include <pthread.h>
#include "decomp.h"
#include "myGraph.h"
#include "tools/debug.h"
#include "file.h"
#include "myAlgorithm.h"
#include "common.h"
#include <thread>
#include <condition_variable>
#include <atomic>
#include <chrono>
#include <mutex>

myGraph G;

std::atomic<bool> terminateFlag{false}; // 定义全局的原子变量
std::atomic<bool> done{false};
std::condition_variable cv;
std::mutex cv_m;

void run(myAlgorithm &algorithm_obj, std::string algorithm_name)
{

    auto start = std::chrono::high_resolution_clock::now();
    algorithm_obj.run(algorithm_name);
   
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    if (terminateFlag.load())
    {
        cout << -1 <<" "<<-1<<endl;
    }
    else
    {
        std::cout << algorithm_obj.getNewEdgesSize() << " ";
        std::cout << elapsed.count() << std::endl;
    }

    std::lock_guard<std::mutex> lk(cv_m);
    done = true;
    cv.notify_one();
    algorithm_obj.verify(algorithm_name);
}
void monitorFunction(int timeoutSeconds)
{
    std::unique_lock<std::mutex> lk(cv_m);
    if (cv.wait_for(lk, std::chrono::seconds(timeoutSeconds), []
                    { return done.load(); }))
    {
        if (done.load())
        {
        }
    }
    else
    {

        terminateFlag.store(true); // 设置终止标志
        // 保证线程已经识别到终止标志，并设置完成标志。
        cv.wait(lk, []
                { return done.load(); });
    }
}
int main(int argc, char *argv[])
{
    DEBUG_ASSERT(argc == 6, "invalid input");
    std::string dataset_file_name = argv[1];
    int k = std::atoi(argv[2]);
    std::string k_truss_community_file_name = argv[3];
    int why_not_vertex = std::atoi(argv[4]);
    std::string algorithm_name = argv[5];

    file_ReadGraph(G, dataset_file_name);

    decomp decomp_obj(G);
    decomp_obj.trussDecomp();
    // G.printGraph();
    // G.printTau();
    // cout<<k_truss_community_file_name<<" "<<why_not_vertex<<endl;
    USetEdge k_truss_community;
    file_ReadKTrussCommunity(k_truss_community, k_truss_community_file_name);

    myAlgorithm algorithm_obj(G, k, k_truss_community, why_not_vertex,dataset_file_name);
    std::thread t([&algorithm_obj, algorithm_name]()
                  { run(algorithm_obj, algorithm_name); });

    std::thread monitor(monitorFunction, 7200);

    t.join();
    monitor.join();

    return 0;
}

