#include "myGraph.h"
void myGraph::printTau()
{
    cout << "tau:" << endl;
    for (int i = 0; i < n; i++)
    {
        for (auto it = tau[i].begin(); it != tau[i].end(); it++)
        {
            cout << i << " " << it->first << " " << it->second << endl;
        }
    }
}
void myGraph::printGraph()
{
    cout << "graph:" << endl;
    for (int i = 0; i < n; i++)
    {
        for (auto it = graph[i].begin(); it != graph[i].end(); it++)
        {
            cout << i << " " << it->first << " " << it->second << endl;
        }
    }
}

void myGraph::getTau(vector<UMapIntInt> &t)
{
    t = this->tau;
}
void myGraph::getGraph(vector<UMapIntInt>g)
{
    g=this->graph;
}

int myGraph::getN()
{
    return n;
}
int myGraph::getM()
{
    return m;
}
