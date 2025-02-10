#include "file.h"
/*
input:
    dataset_file_name: the name of the dataset file
output:
    G: the graph read from the dataset file
*/
void file_ReadGraph(myGraph &G, string dataset_file_name)
{
    ifstream fin(dataset_file_name);

    int vMax = 0;
    int u, v;

    // read the edge list and count the number of vertices
    while (fin >> u >> v)
    {
        if (u == v)
            continue;
        G.m++;
        vMax = max(vMax, max(u, v));
    }
    G.n = vMax + 1;

    // read the edge list again and build the graph
    fin.clear();
    fin.seekg(0, ios::beg);

    // initialize the graph
    G.adj.resize(G.n);
    G.deg.resize(G.n, 0);

    while (fin >> u >> v)
    {
        if (u == v)
            continue;
        G.adj[u][v] = 0;
        G.adj[v][u] = 0;
        G.deg[u]++;
        G.deg[v]++;
    }
    fin.close();
}
/*
input:
    k_truss_community_file_name: the name of the k-truss community file
    k_truss_community file is of format:
        u v
        u v
       ...
output:
    k_truss_community: the k-truss community read from the file
*/
void file_ReadKTrussCommunity(USetEdge &k_truss_community, string k_truss_community_file_name)
{
    ifstream fin(k_truss_community_file_name);
    int u, v;
    while (fin >> u >> v)
    {
        k_truss_community.insert(makeEdge(u, v));
    }
    fin.close();
}


