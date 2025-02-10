#pragma once
#include <bitset>
#include "decomp.h"
#include <numeric>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <queue>
#include "common.h"
#include "myGraph.h"
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include "tools/debug.h"
#include <set>
#include <climits>
#include <thread>
#include <condition_variable>
#include <atomic>
#include <chrono>
#include <mutex>
using namespace std;

extern std::atomic<bool> terminateFlag;
extern std::atomic<bool> done;
extern std::condition_variable cv;
extern std::mutex cv_m;
extern myGraph G;
class myAlgorithm
{
private:
    vector<UMapIntInt> &graph, &tau;
    int n, m;
    int k;
    USetEdge &k_truss_community;
    unordered_set<int> k_truss_vertices;
    int w;
    vector<Edge> newEdges;
    string dataSetFileName;

public:
    myAlgorithm(myGraph &G, int k, USetEdge &k_truss_community, int why_not_vertex, string dataSetFileName) : graph(G.graph), tau(G.tau), n(G.n), m(G.m), k(k), k_truss_community(k_truss_community), w(why_not_vertex)
    {
        for (auto e : k_truss_community)
        {
            k_truss_vertices.insert(e.u);
            k_truss_vertices.insert(e.v);
        }
        int i = 0, j = 0;
        for (int tt = dataSetFileName.size() - 1; tt >= 0; tt--)
        {
            if (dataSetFileName[tt] == '.')
                j = tt;
            if (dataSetFileName[tt] == '/')
            {
                i = tt + 1;
                break;
            }
        }
        this->dataSetFileName = dataSetFileName.substr(i, j - i);

        // cout<<"dataSetFileName: "<<dataSetFileName<<endl;
    };
    void run(string algorithm_name);

    void localPlus();


    void mDebug();

   
    void globalPlus();  
   
    vector<Edge> sort_community_local(USetEdge &k_truss_community);
    USetEdge calculateDependentEdgesNaive(Edge e, vector<TRI> &mRecoverTau);

    void gtm();
    bool checkWhetherTriangleConnected(vector<int> &outofFollowers);

    USetEdge calculateNewComers(USetEdge &activators, vector<TRI> &mRecoverTau);
    USetEdge calculateDependentEdges(Edge e, USetEdge &newKtrussEdges, vector<TRI> &mRecoverTau);
    USetEdge calculateDependentEdges(Edge e, USetEdge &newComers);
    void allTrussUpdateWithDeletion(Edge e0);
    void allTrussUpdateWithInsertion(Edge e0);
    void allTrussUpdateWithInsertion(Edge e0, UMapEdgeInt &mRecoverTau);
    void allTrussUpdateWithInsertion(Edge e0, UMapEdgeInt &mRecoverTau, USetEdge &newKtrussEdges, USetEdge &newK_1trussEdges);

    vector<Edge> kTrussUpdateWithInsertion(Edge e0);
    void verify(string algorithm_name);
    void printEdges(string message, vector<Edge> &edges);
    void printTau();
    void calbound(Edge e0, int &k1, int &k2);
    vector<int> getNeighbors(int u, int v);
    vector<Edge> queryProcess(int k);
    int getNewEdgesSize();
    vector<Edge> selectCandidates();
    bool checkWhetherTimeOut();
    void printNewCommunity(string &algorithm_name);
    size_t getSetSizeInBytes(const set<set<Edge>> &s);
};