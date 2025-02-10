#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "common.h"
#include "myGraph.h"
using namespace std;
class decomp
{
    int n,m;
    vector<UMapIntInt> &adj , &tau , &graph ;
    vector<int> deg,bin, mapto;
    vector<Edge> binEdge;
    vector<UMapIntInt> pos;
    vector<vector<int>> A;
    

public:
    decomp(myGraph &G):n(G.n), m(G.m), adj(G.adj), tau(G.tau), graph(G.graph),deg(G.deg){};
    decomp(myGraph &G, bool flag);
    void trussDecomp();
    void updateSupport(int u, int v, int delta);
    void reorder();
    void countTriangles();
    void binSort();
    bool compVertex(int i, int j);
    void orderPair(int &u, int &v);
    void intersect(const vector<int> &a, const vector<int> &b, vector<int> &c);
    void printClass(int u, int v, int t);
    void removeEdge(int u, int v);
    void updateEdge(int u, int v, int minsup);
};