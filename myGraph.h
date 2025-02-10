#pragma once
#include <iostream>
#include <vector>
#include <map>
#include "common.h"
using namespace std;

class myGraph
{
    int n, m;                      // the number of vertices and edges
    vector<UMapIntInt> graph, tau; // the support graph and tau
    vector<int> deg;               // degree of each vertex
    vector<UMapIntInt> adj;        // adjacent list of each vertex
    friend class decomp;
    friend void file_ReadGraph(myGraph &G, string dataset_file_name);
    friend class myAlgorithm;

public:
    myGraph(){};
    void printTau();
    void printGraph();
    void getTau(vector<UMapIntInt> &t);
    void getGraph(vector<UMapIntInt> g);
    int getN();
    int getM();
  
};

