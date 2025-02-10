#include "decomp.h"

void decomp::trussDecomp()
{

    tau.resize(n);
    for (int i = 0; i < n; ++i)
        tau[i].clear();
    reorder();
    countTriangles();
    binSort();
    for (int s = 0; s < m; ++s)
    {
        int u = binEdge[s].u;
        int v = binEdge[s].v;
        orderPair(u, v);
        int supuv = adj[u][v];
        printClass(u, v, supuv + 2);
        int nfound = 0;
        for (auto it = adj[u].begin(); it != adj[u].end(); ++it)
        { // 将（u，v）边从G中删除
            if (nfound == supuv)
                break;
            int w = it->first;
            if (w == v)
                continue;
            if (adj[v].find(w) != adj[v].end())
            { // w 是u、v共同的邻居
                ++nfound;
                updateEdge(u, w, supuv);
                updateEdge(v, w, supuv);
            }
        }
        removeEdge(u, v);
    }
}
void decomp::updateSupport(int u, int v, int delta)
{
    adj[u][v] += delta;
    adj[v][u] += delta;
}
void decomp::reorder()
{
    mapto.resize(n);
    for (int i = 0; i < n; ++i)
        mapto[i] = i;
    sort(mapto.begin(), mapto.end(), [this](int i, int j)
         { return deg[i] < deg[j] || (deg[i] == deg[j] && i < j); });
}
void decomp::countTriangles()
{
    A.resize(n);
    for (int i = 0; i < n; ++i)
        A[i].clear();
    int nDeltas = 0;
    for (int vi = n - 1; vi >= 0; --vi)
    {
        int v = mapto[vi];
        for (auto it = adj[v].begin(); it != adj[v].end(); ++it)
        {
            int u = it->first; // v的邻接点u
            if (!compVertex(u, v))
                continue; // 必须满足u<v
            vector<int> common;
            intersect(A[u], A[v], common); // 保证每个三角形只会被计算一遍，而不是三遍。由度数第2大和第三大的两点构成的边计算而来
            for (unsigned i = 0; i < common.size(); ++i)
            {
                int w = mapto[common[i]];
                ++nDeltas;
                updateSupport(u, v, 1);
                updateSupport(v, w, 1);
                updateSupport(w, u, 1);
            }
            A[u].push_back(vi);
        }
    }
    graph = adj;
}
void decomp::binSort()
{

    bin.clear();
    bin.resize(n, 0); // bin[i]:桶i的中元素的多少
    int nBins = 0;
    int mp = 0;
    for (int u = 0; u < n; ++u)
    {
        auto tadj = adj[u];
        for (auto it = tadj.begin(); it != tadj.end(); ++it)
        {
            int v = it->first;
            if (!compVertex(u, v))
                continue;
            int sup = it->second;
            if (sup == 0)
            {
                printClass(u, v, 2);
                removeEdge(u, v);
                continue;
            }
            ++mp;
            ++bin[sup];
            nBins = max(sup, nBins);
        }
    }
    m = mp;
    ++nBins;
    int count = 0;
    for (int i = 0; i < nBins; ++i)
    { // bin[i]:第i个桶中第一个元素的位置。
        int binsize = bin[i];
        bin[i] = count;
        count += binsize;
    }
    pos.clear();
    pos.resize(n);
    for (int i = 0; i < n; ++i)
        pos[i].clear();
    binEdge.resize(m);
    for (int u = 0; u < n; ++u)
        for (auto it = adj[u].begin(); it != adj[u].end(); ++it)
        {
            int v = it->first;
            if (!compVertex(u, v))
                continue;
            int sup = it->second;
            Edge e = {u, v};
            int &b = bin[sup]; //
            binEdge[b] = e;    // binEdge：按照边的sup从小到大排列，元素为边
            pos[u][v] = b++;   // 与binEdge 相反，记录边的位置，元素为int
        }
    for (int i = nBins; i > 0; --i)
        bin[i] = bin[i - 1]; // bin[i]:第i个桶里的元素位置为：[bin[i]-bin[i+1])
    bin[0] = 0;
}
bool decomp::compVertex(int i, int j)
{
    return deg[i] < deg[j] || (deg[i] == deg[j] && i < j);
}
/*
功能：使得边(u,v)的deg (u)<deg(v)
输入：(u,v)
输出：空
*/
void decomp::orderPair(int &u, int &v)
{
    if (!compVertex(u, v))
        swap(u, v);
}
void decomp::intersect(const vector<int> &a, const vector<int> &b, vector<int> &c)
{ // a,b vector 里面的元素都是逆序排列
    c.clear();
    unsigned j = 0;
    for (unsigned i = 0; i < a.size(); ++i)
    {
        while (j < b.size() && b[j] > a[i])
            ++j;
        if (j >= b.size())
            break;
        if (b[j] == a[i])
            c.push_back(a[i]);
    }
}
void decomp::printClass(int u, int v, int t)
{
    tau[u][v] = t;
    tau[v][u] = t;
    // fout << "(" << u << "," << v << "):" << tau << endl;
}
void decomp::removeEdge(int u, int v)
{
    adj[u].erase(v);
    adj[v].erase(u);
}
void decomp::updateEdge(int u, int v, int minsup)
{
    orderPair(u, v);
    int sup = adj[u][v];
    if (sup <= minsup)
        return;
    int p = pos[u][v];
    int posbin = bin[sup];
    Edge se = binEdge[posbin];
    Edge e = {u, v};
    if (p != posbin)
    {
        pos[u][v] = posbin;
        pos[se.u][se.v] = p;
        binEdge[p] = se;
        binEdge[posbin] = e;
    }
    ++bin[sup];
    updateSupport(u, v, -1);
}
// 2th truss decomposition: adj and deg have to be recomputed
decomp::decomp(myGraph &G, bool flag):tau(G.tau),graph(G.graph),adj(G.adj)
{
    this->n = G.n;
    this->m=G.m;
    
    this->tau.resize(n);
    this->adj.resize(n);
    this->deg.resize(n, 0);

    for(int i=0;i<n;i++)
    {
        for(auto it=G.graph[i].begin();it!=G.graph[i].end();it++)
        {
            int j=it->first;
            this->adj[i][j]=0;
            this->adj[j][i]=0;
            this->deg[i]++;
            this->deg[j]++;
        }
    }
    graph.resize(n);
}
