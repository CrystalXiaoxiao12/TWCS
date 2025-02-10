
#include "myAlgorithm.h"
/*
input: algorithm_name
output: none
description: This function is used to run the specified algorithm.
*/

void myAlgorithm::run(string algorithm_name)
{
   if (algorithm_name == "GlobalPlus")
    {
        globalPlus();
    }
    else if (algorithm_name == "LocalPlus")
    {
        localPlus();
    }
    else if (algorithm_name == "GTM")
    {
        gtm();
    }
}
/*
This function corresponds to the optimization O1 of local.
 */
vector<Edge> myAlgorithm::sort_community_local(USetEdge &k_truss_community)
{
    // Create a vector from the given edge set
    int k_truss_size = k_truss_community.size();
    vector<pair<int, Edge>> community_vector;

    // Calculate support for each edge in the community
    for (auto e : k_truss_community)
    {
        int u = e.u, v = e.v;
        vector<int> wList1 = getNeighbors(w, u);
        vector<int> wList2 = getNeighbors(w, v);
        int supk_2_1 = 0;
        for (int z : wList1)
        {
            if (tau[z][u] >= k - 2 && tau[z][w] >= k - 2)
            {
                supk_2_1++;
            }
        }
        int supk_2_2 = 0;
        for (int z : wList2)
        {
            if (tau[z][v] >= k - 2 && tau[z][w] >= k - 2)
            {
                supk_2_2++;
            }
        }
        int support = max(0, k - 2 - supk_2_1) + max(0, k - 2 - supk_2_2);
        community_vector.push_back(pair<int, Edge>(support, e));
    }

    // Sort edges based on their support values in increasing order
    sort(community_vector.begin(), community_vector.end());
    // cout<<"sorted"<<endl;
    // Clear previous community and insert sorted edges
    vector<Edge> res;
    for (auto it : community_vector)
    {
        Edge e = it.second;
        // cout<<e.u<<" "<<e.v<<" "<<it.first<<endl;
        // res.insert(e);
        res.push_back(e);
    }
    // for(auto e:res)
    // {
    //     cout<<e.u<<" "<<e.v<<endl;
    // }
    return res;
}

USetEdge myAlgorithm::calculateDependentEdgesNaive(Edge e0, vector<TRI> &mRecoverTau)
{
    USetEdge dependentEdges;
    dependentEdges.insert(e0);
    int u0 = e0.u, v0 = e0.v;
    // cout << "delete u,v: " << u0 << " " << v0 << endl;
    // printTau();
    int kMax = tau[u0][v0];
    mRecoverTau.push_back(TRI(u0, v0, tau[u0][v0]));
    // 改变e0
    graph[u0].erase(v0);
    graph[v0].erase(u0);
    tau[u0].erase(v0);
    tau[v0].erase(u0);

    queue<Edge> q; // 记录可能truss值减少1的边
    // 将与e0形成三角形的边，truss值<=tau[eo][v0]的边加入queue
    USetEdge inL; // 所有边只入一次队列
    vector<int> neighbors = getNeighbors(u0, v0);

    for (auto w : neighbors)
    {
        int minn = min(tau[w][u0], tau[w][v0]);
        if (minn == k)
        {
            if (tau[w][u0] == minn)
            {
                q.push(makeEdge(w, u0));
                inL.insert(makeEdge(w, u0));
            }
            if (tau[w][v0] == minn)
            {
                q.push(makeEdge(w, v0));
                inL.insert(makeEdge(w, v0));
            }
        }
    }

    UMapEdgeInt s; // support
    int tk = k;
    while (!q.empty())
    {
        int u = q.front().u;
        int v = q.front().v;
        q.pop();
        s[makeEdge(u, v)] = 0;
        vector<int> neighbors = getNeighbors(u, v);
        for (auto w : neighbors)
        {
            if (tau[w][v] >= tk && tau[w][u] >= tk)
            {
                s[makeEdge(u, v)]++;
            }
        }
        queue<Edge> supQueue;
        if (s[makeEdge(u, v)] < tk - 2)
        {
            supQueue.push(makeEdge(u, v));
        }

        while (!supQueue.empty())
        {
            int x = supQueue.front().u;
            int y = supQueue.front().v;

            dependentEdges.insert(makeEdge(x, y));
            mRecoverTau.push_back(TRI{x, y, k});
            supQueue.pop();

            tau[x][y]--;
            tau[y][x]--;
            s[makeEdge(x, y)] = -1;
            vector<int> neighbors = getNeighbors(x, y);
            for (auto z : neighbors)
            {
                if (tau[x][z] < tk || tau[y][z] < tk)
                    continue;
                if (tau[x][z] == tk)
                {
                    if (s.find(makeEdge(x, z)) != s.end() && s[makeEdge(x, z)] != -1)
                    {
                        s[makeEdge(x, z)]--;
                        if (s[makeEdge(x, z)] == tk - 3)
                        {
                            supQueue.push(makeEdge(x, z));
                        }
                    }
                    if (inL.find(makeEdge(x, z)) == inL.end())
                    {
                        inL.insert(makeEdge(x, z));
                        q.push(makeEdge(x, z));
                    }
                }
                if (tau[y][z] == tk)
                {
                    if (s.find(makeEdge(y, z)) != s.end() && s[makeEdge(y, z)] != -1)
                    {
                        s[makeEdge(y, z)]--;
                        if (s[makeEdge(y, z)] == tk - 3)
                        {
                            supQueue.push(makeEdge(y, z));
                        }
                    }
                    if (inL.find(makeEdge(y, z)) == inL.end())
                    {
                        inL.insert(makeEdge(y, z));
                        q.push(makeEdge(y, z));
                    }
                }
            }
        }
    }
    
    return dependentEdges;
}

USetEdge myAlgorithm::calculateNewComers(USetEdge &activators, vector<TRI> &mRecoverTau)
{
    // Timer timer("calculateNewComers");
    USetEdge L;

    UMapEdgeInt mIsInRecoverTau;
    // w原本存在graph中的邻居边：support>=k-2:所在的三角形除w邻边外的另外一边的truss值>=k-1（因为该种边truss值最多+1）
    for (auto iter : graph[w])
    {
        int u = iter.first;
        int sup = 0;
        for (auto iter2 : graph[u])
        {
            int v = iter2.first;
            if (graph[w].find(v) != graph[w].end() || activators.find(makeEdge(w, v)) != activators.end())
            {
                if (tau[u][v] >= k - 1)
                    sup++;
            }
        }
        if (sup >= k - 2)
        {
            L.insert(makeEdge(w, u));
            mRecoverTau.push_back(TRI(w, u, tau[w][u]));
            mIsInRecoverTau[makeEdge(w, u)] = 1;
            tau[w][u] = k - 1;
            tau[u][w] = k - 1;
        }
    }
    for (auto e : activators)
    {
        graph[e.u][e.v] = 1;
        graph[e.v][e.u] = 1;
        tau[e.u][e.v] = k;
        tau[e.v][e.u] = k;
    }
    // 边界边(u,v)，u，v都是w的邻居:truss值=k-1
    for (auto iter : graph[w])
    {
        int u = iter.first;
        for (auto iter2 : graph[w])
        {
            int v = iter2.first;
            if (u != v && graph[u].find(v) != graph[u].end() && tau[u][v] == k - 1)
            {
                L.insert(makeEdge(u, v));
            }
        }
    }
   
    queue<Edge> Q;
    UMapEdgeInt s;
    USetEdge inQueue;        // truss值待更新的边
    UMapEdgeInt mMovedToBin; // 记录是否移动到bin
    // cout<<"tk: "<<tk<<endl;
    for (auto e : L)
    {
        Q.push(e);
        inQueue.insert(e);
    }
    int tk = k - 1;
    while (!Q.empty())
    {

        int x = Q.front().u;
        int y = Q.front().v;
        Q.pop();

        vector<int> wList = getNeighbors(x, y);
        vector<Edge> tempEdges;

        for (int i = 0; i < wList.size(); i++)
        {
            int z = wList[i];
            if ((tau[z][x] < tk) || (tau[z][y] < tk))
                continue;
            if (tau[z][x] == tk && mMovedToBin[makeEdge(z, x)] == 2)
                continue;
            if (tau[z][y] == tk && mMovedToBin[makeEdge(z, y)] == 2)
                continue;
            s[{x, y}]++;
            if ((tau[z][x] == tk) && inQueue.find(makeEdge(z, x)) == inQueue.end())
            {
                tempEdges.push_back(makeEdge(z, x));
            }
            if ((tau[z][y] == tk) && inQueue.find(makeEdge(z, y)) == inQueue.end())
            {
                tempEdges.push_back(makeEdge(z, y));
            }
        }
        mMovedToBin[makeEdge(x, y)] = 1;
       
        if (s[{x, y}] > tk - 2)
        {

            for (auto edge : tempEdges)
            {
                if (inQueue.find(edge) == inQueue.end())
                {
                    Q.push(edge);
                    inQueue.insert(edge);
                }
            }
        }
        else
        {
            queue<Edge> Q1; // mToUpdatedEges里面s值<=tk-2的边
            Q1.push(makeEdge(x, y));
            // 入队又出队
            while (!Q1.empty())
            {
                int x1 = Q1.front().u, y1 = Q1.front().v;
                mMovedToBin[makeEdge(x1, y1)] = 2;
               
                Q1.pop();
                vector<int> wList = getNeighbors(x1, y1);
                for (int i = 0; i < wList.size(); i++)
                {
                    int z1 = wList[i];
                    if (tau[x1][z1] < tk || tau[y1][z1] < tk)
                        continue;
                    if (tau[x1][z1] == tk && mMovedToBin[makeEdge(x1, z1)] == 2)
                        continue;
                    if (tau[y1][z1] == tk && mMovedToBin[makeEdge(y1, z1)] == 2)
                        continue;
                    if (tau[x1][z1] == tk)
                    {
                        if (mMovedToBin[makeEdge(x1, z1)] == 1)
                        {
                            s[makeEdge(x1, z1)] = s[makeEdge(x1, z1)] - 1;
                            if (s[makeEdge(x1, z1)] == tk - 2)
                            {
                                Q1.push(makeEdge(x1, z1));
                            }
                        }
                    }
                    if (tau[y1][z1] == tk)
                    {
                        if (mMovedToBin[makeEdge(y1, z1)] == 1)
                        {
                            s[makeEdge(y1, z1)] = s[makeEdge(y1, z1)] - 1;
                            if (s[makeEdge(y1, z1)] == tk - 2)
                            {
                                Q1.push(makeEdge(y1, z1));
                            }
                        }
                    }
                }
            }
        }
    }
    USetEdge improvedEdges;
    for (auto edge : inQueue)
    {
        if (mMovedToBin[edge] == 1)
        {
            int u = edge.u;
            int v = edge.v;
            improvedEdges.insert(makeEdge(u, v));
            if (mIsInRecoverTau.find(edge) == mIsInRecoverTau.end())
            {
                mRecoverTau.push_back(TRI{edge.u, edge.v, tau[edge.u][edge.v]});
            }
            tau[u][v] += 1;
            tau[v][u] += 1;
            // cout << "truss increase: " << u << " " << v << " " << tau[u][v] << endl;
        }
    }
    for (auto it : mRecoverTau)
    {
        int u = get<0>(it);
        int v = get<1>(it);
        int truss = get<2>(it);
       
    }
    improvedEdges.insert(activators.begin(), activators.end());
    // end = clock();
    // if (double(end - start) / CLOCKS_PER_SEC * 1000 > 100)
    // cout << "边插入时间: " << double(end - start) / CLOCKS_PER_SEC << endl;

    return improvedEdges;
}

void myAlgorithm::globalPlus()
{

    // 将ver中所有点映射到唯一元素
    unordered_map<int, int> mp;
    int _ = 0;
    for (auto x : k_truss_vertices)
    {
        mp[x] = _++;
    }
   // 将set<edge>输入转换成唯一映射的string
    auto edgeToString = [&](USetEdge S) -> string
    {
        string res(k_truss_vertices.size(), '0');
        for (auto e : S)
        {
            int u = e.u == w ? e.v : e.u;
            res[mp[u]] = '1';
        }
        return res;
    };
    auto edgeToBitset = [&](USetEdge &S) -> DynamicBitset
    {
        DynamicBitset res(k_truss_vertices.size());
        for (auto e : S)
        {
            int u = e.u == w ? e.v : e.u;
            res.set(mp[u], 1);
        }
        return res;
    };
   
    // 两个set取交集
    auto setIntersection = [&](USetEdge a, USetEdge b) -> USetEdge
    {
        USetEdge res;
        for (auto e : a)
        {
            if (b.find(e) != b.end())
            {
                res.insert(e);
            }
        }
        return res;
    };

    RadixTree *trie = new RadixTree();
    auto whetherSkip = [&](USetEdge &S) -> bool
    {
        
        DynamicBitset state = edgeToBitset(S);
        return trie->search(state);
    };
    newEdges.clear();

    // // 按照连接边数由小到大的社区点
    for (auto it : k_truss_vertices)
    {

        if (checkWhetherTimeOut())
            return;
        int x = it;
        bool f = false; // 中途发现已被访问
       

        // 模仿x插入的边
        USetEdge S;
        for (auto iter : graph[x])
        {
            int y = iter.first;
            if (k_truss_vertices.find(y) != k_truss_vertices.end() && graph[w].find(y) == graph[w].end() && tau[x][y] >= k)
            {
                S.insert(makeEdge(w, y));
            }
        }
        if (whetherSkip(S))
        {
            
            continue;
        }
        else
        {
            DynamicBitset state = edgeToBitset(S);
            trie->insert(state);
        }

        // 将S复制给mRecoverGraph
        vector<TRI> mRecoverTau;
        vector<Edge> mRecoverGraph(S.begin(), S.end());
        USetEdge newComers = calculateNewComers(S, mRecoverTau);
        

        USetEdge tempNewEdges;
        UMapEdgeBool mRemain; // 标记是否存进tempNewEdges,该条边无需再处理
        int len = S.size();

        // 该轮最少要加的边数
        int minSize = 0;
        int lastlen = 0;
        // 直到所有插入边被删除或者加入到tempNewEdges
        while (len)
        {

            if (len == lastlen)
            {
                DEBUG_ASSERT(false, "infinite loop");
            }
            lastlen = len;
            // 每一轮删除followers最大的边
            int maxFollowers = 0;
            Edge maxEdge;
            UMapEdgeBool mCaled; // 计算该轮最大followers时，已经计算过的边
       
            for (auto e : S)
            {

                if (mRemain.find(e) != mRemain.end() || mCaled.find(e) != mCaled.end())
                    continue;

                int followers = 0;
                vector<TRI> mRecoverTau;
                // 计算e的followers
                USetEdge followersEdges = calculateDependentEdges(e, newComers, mRecoverTau); // 包括e

                // 取其中属于S也就是实际插入边的交集
                followersEdges = setIntersection(followersEdges, S);

                // 如果followersEdges.size()==S.size()-1，还需要检查一下最后一条边和k truss 社区是否是边连通的，不是的话加到followers里面
                if (followersEdges.size() < S.size())
                {
                    // Timer timer("checkWhetherTriangleConnected");
                    // S中不是followers的那条边
                    vector<int> outofFollowers;
                    for (auto e1 : S)
                    {
                        if (followersEdges.find(e1) == followersEdges.end())
                        {
                            int u = e1.u == w ? e1.v : e1.u;
                            outofFollowers.push_back(u);
                        }
                    }
                    bool ff = checkWhetherTriangleConnected(outofFollowers);
                    if (!ff)
                    {
                        followersEdges = S;
                    }
                }
    

                //  恢复graph和tau
                graph[e.u][e.v] = 1;
                graph[e.v][e.u] = 1;
                for (auto tri : mRecoverTau)
                {
                    int u = get<0>(tri), v = get<1>(tri), mTau = get<2>(tri);
                    tau[u][v] = mTau;
                    tau[v][u] = mTau;
                }
                if (followersEdges.size() != S.size() && followersEdges.size() >= maxFollowers)
                {
                    maxFollowers = followersEdges.size();
                    maxEdge = e;
                }
              
                // 等于当前全集则不删除，放入tempNewEdges
                if (followersEdges.size() == S.size())
                {
                    
                    USetEdge group;
                    group.insert(e);
                    for (auto edge1 : group)
                    {
                        // 一条边只会被加入到tempNewEdges一次
                        if (tempNewEdges.find(edge1) == tempNewEdges.end())
                        {
                            len--;
                            minSize++;
                            mRemain[edge1] = true;
                            tempNewEdges.insert(edge1);
                        }
                    }

                    // 下界超过newEdges.size()则终止
                    if (newEdges.size() > 0 && minSize >= newEdges.size())
                    {
                        f = true;
                       
                        break;
                    }
                }
                else
                {
                    // 标记所有followers为计算过
                    for (auto edge1 : followersEdges)
                    {
                        mCaled[edge1] = true;
                    }
                }
            }
            if (f)
                break;
            
            if (maxFollowers > 0)
            {
               
                // 删除最大followers的边（更新graph和tau）

                // 更新S和newKtrussedges

                USetEdge maxfollowersEdges = calculateDependentEdges(maxEdge, newComers);
                for (auto e : maxfollowersEdges)
                {
                    // 有一些是已经存在的边，不需要在S中删除
                    if (S.find(e) != S.end())
                    {
                        S.erase(e);
                        len--;
                    }
                    newComers.erase(e);
                }

                if (whetherSkip(S))
                {
                    
                    f = true;
                    break;
                }
                else
                {
                    DynamicBitset state = edgeToBitset(S);
                    trie->insert(state);
                }
            }
        }
       
        if (!f && (newEdges.size() == 0 || tempNewEdges.size() < newEdges.size()))
        {
            
            newEdges = vector<Edge>(tempNewEdges.begin(), tempNewEdges.end());
            
        }
        // 恢复graph和tau
        for (auto tri : mRecoverTau)
        {
            int u = get<0>(tri), v = get<1>(tri), mTau = get<2>(tri);
            tau[u][v] = mTau;
            tau[v][u] = mTau;
        }
        for (auto e : mRecoverGraph)
        {
            if (graph[e.u].find(e.v) != graph[e.u].end())
            {
                graph[e.u].erase(e.v);
                graph[e.v].erase(e.u);
            }
            if (tau[e.u].find(e.v) != tau[e.u].end())
            {
                tau[e.u].erase(e.v);
                tau[e.v].erase(e.u);
            }
        }

        // 已达到最小，无需再模仿
        if (newEdges.size() == 1)
        {
           
            break;
        }
    }
}

void myAlgorithm::localPlus()
{
   
    newEdges.clear();

    set<set<Edge>> visitedEdgeSets;
    auto visitEdgeSet = [&](set<Edge> &edgeSet) -> bool
    {
        if (visitedEdgeSets.find(edgeSet) != visitedEdgeSets.end())
        {
            return true;
        }
        else
        {
            visitedEdgeSets.insert(edgeSet);
            return false;
        }
    };

   
    vector<Edge> sortedCommunity = sort_community_local(k_truss_community);

    for (auto mSeedEdge : sortedCommunity)
    {
        if (checkWhetherTimeOut())
            return;
        bool mIsVisited = false;
        bool mTerminateFlag = false;
        // 插入种子三角形
        int mSeedu = mSeedEdge.u, mSeedv = mSeedEdge.v;
       
        set<Edge> mAllInsertedEdges; //
        set<Edge> mNotInKTrussEdges; // 这些边还需要其他边的加入才能加入k truss 社区

        if (graph[w].find(mSeedu) == graph[w].end())
            mAllInsertedEdges.insert(makeEdge(w, mSeedu));
        if (graph[w].find(mSeedv) == graph[w].end())
            mAllInsertedEdges.insert(makeEdge(w, mSeedv));

        if (mAllInsertedEdges.size() > 0 && visitEdgeSet(mAllInsertedEdges))
            continue;
        UMapEdgeInt mRecoverTau;
        for (auto e : mAllInsertedEdges)
            allTrussUpdateWithInsertion(e, mRecoverTau);

        if (tau[w][mSeedu] < k)
            mNotInKTrussEdges.insert(makeEdge(w, mSeedu));
        if (tau[w][mSeedv] < k)
            mNotInKTrussEdges.insert(makeEdge(w, mSeedv));
        // 计算p2:e至少还要多少个三角形才能加入k truss 社区：k-2-已经存在的三角形数（另外两边必须满足truss值>=k）
        auto calP2 = [&](Edge e) -> int
        {
            int p2 = 0;
            int u = e.u, v = e.v;
            for (auto iter : graph[u])
            {
                int x = iter.first;
                if (x != v && graph[v].find(x) != graph[v].end())
                {
                    if (tau[u][x] >= k - 1 && tau[v][x] >= k - 1)
                        p2++;
                }
            }
            return max(0, k - 2 - p2);
        };
        //  计算 p1:侯选边能让mNotInKTrussEdges中多少条边构成三角形，三角形内除侯选边，mNotInKTrussEdges的那条边之外，其他边的truss值>=k
        auto calP1 = [&](Edge e) -> int
        {
            int p1 = 0;
            int u = e.u, v = e.v;
            for (auto iter : graph[u])
            {
                int x = iter.first;
                // x为u,v的共同邻居
                if (x != v && graph[v].find(x) != graph[v].end())
                {
                    if (mNotInKTrussEdges.find(makeEdge(u, x)) != mNotInKTrussEdges.end() && graph[x].find(v) != graph[x].end() && tau[x][v] >= k)
                        p1++;
                    if (mNotInKTrussEdges.find(makeEdge(v, x)) != mNotInKTrussEdges.end() && graph[x].find(u) != graph[x].end() && tau[x][u] >= k)
                        p1++;
                }
            }
            return p1;
        };
        //  筛选侯选边：e所在的三角形中，除该边外，另一边tau>=k，剩下的那条边作为侯选边。并计算goodness, 将其加入侯选边goodness队列
        auto getCandidates = [&](Edge e, USetEdge &mCandidateEdges, priority_queue<pair<int, Edge>> &mCandidateEdgeQueue, UMapEdgeInt &mP1ofCandidateEdges, UMapEdgeInt &mP2ofCandidateEdges)
        {
            int u = e.u, v = e.v;
            // u的邻居
            for (auto iter : graph[u])
            {
                int x = iter.first;
                if (x != v && graph[v].find(x) == graph[v].end() && tau[u][x] >= k)
                {
                    Edge mEdge = makeEdge(v, x);
                    if (mCandidateEdges.find(mEdge) == mCandidateEdges.end())
                    {
                        mCandidateEdges.insert(mEdge);
                        // 计算所有候选边的p1-p2
                        int mP1 = calP1(mEdge), mP2 = calP2(mEdge);
                        mP1ofCandidateEdges[mEdge] = mP1;
                        mP2ofCandidateEdges[mEdge] = mP2;
                        
                        mCandidateEdgeQueue.push({mP1 - mP2, mEdge});
                    }
                }
            }
            // v的邻居
            for (auto iter : graph[v])
            {
                int x = iter.first;

                if (x != u && graph[u].find(x) == graph[u].end() && tau[v][x] >= k)
                {
                    Edge mEdge = makeEdge(u, x);
                    if (mCandidateEdges.find(mEdge) == mCandidateEdges.end())
                    {
                        mCandidateEdges.insert(mEdge);
                        // 计算所有候选边的p1-p2
                        int mP1 = calP1(mEdge), mP2 = calP2(mEdge);
                        mP1ofCandidateEdges[mEdge] = mP1;
                        mP2ofCandidateEdges[mEdge] = mP2;
                        
                        mCandidateEdgeQueue.push({mP1 - mP2, mEdge});
                    }
                }
            }
        };
        auto getCandidates_v2 = [&](Edge e, USetEdge &mCandidateEdges, priority_queue<pair<int, Edge>> &mCandidateEdgeQueue, UMapEdgeInt &mP1ofCandidateEdges, UMapEdgeInt &mP2ofCandidateEdges, UMapEdgeInt &toChangedEdges)
        {
            int u = e.u, v = e.v;
            // u的邻居
            for (auto iter : graph[u])
            {
                int x = iter.first;
                if (x != v && graph[v].find(x) == graph[v].end() && tau[u][x] >= k)
                {
                    Edge mEdge = makeEdge(v, x);
                    if (mCandidateEdges.find(mEdge) == mCandidateEdges.end())
                    {
                        mCandidateEdges.insert(mEdge);
                        // 计算所有候选边的p1-p2
                        int mP1 = calP1(mEdge), mP2 = calP2(mEdge);
                        mP1ofCandidateEdges[mEdge] = mP1;
                        mP2ofCandidateEdges[mEdge] = mP2;
                       
                        mCandidateEdgeQueue.push({mP1 - mP2, mEdge});
                    }
                    else
                    {
                        mP1ofCandidateEdges[mEdge]++;
                        toChangedEdges[mEdge]++;
                    }
                }
            }
            // v的邻居
            for (auto iter : graph[v])
            {
                int x = iter.first;

                if (x != u && graph[u].find(x) == graph[u].end() && tau[v][x] >= k)
                {
                    Edge mEdge = makeEdge(u, x);
                    if (mCandidateEdges.find(mEdge) == mCandidateEdges.end())
                    {
                        mCandidateEdges.insert(mEdge);
                        // 计算所有候选边的p1-p2
                        int mP1 = calP1(mEdge), mP2 = calP2(mEdge);
                        mP1ofCandidateEdges[mEdge] = mP1;
                        mP2ofCandidateEdges[mEdge] = mP2;
                       
                        mCandidateEdgeQueue.push({mP1 - mP2, mEdge});
                    }
                    else
                    {
                        mP1ofCandidateEdges[mEdge]++;
                        toChangedEdges[mEdge]++;
                    }
                }
            }
        };

        auto updateCandidatesP2_v1 = [&](Edge mInsertEdge, Edge e, USetEdge &mCandidateEdges, UMapEdgeInt &mP2ofCandidateEdges, USetEdge &mNewK_1trussEdges, UMapEdgeInt &toChangedEdges)
        {
            
            UMapEdgeBool mVisited; // 记录已经访问过的边,避免重复-1
            mVisited[e] = true;    // e已经更新过其侯选边
            int u = e.u, v = e.v;
            // u的邻居
            for (auto iter : graph[u])
            {
                int x = iter.first;
                if (x != v && graph[v].find(x) == graph[v].end())
                {
                    if (tau[u][x] < k - 1)
                        continue;
                    if (mNewK_1trussEdges.find(makeEdge(u, x)) != mNewK_1trussEdges.end() && mVisited.find(makeEdge(u, x)) != mVisited.end())
                        continue;
                    mVisited[makeEdge(u, x)] = true;
                    Edge mEdge = makeEdge(v, x);
                    if (mCandidateEdges.find(mEdge) != mCandidateEdges.end())
                    {
                       
                        if (mP2ofCandidateEdges[mEdge] > 0)
                        {
                            mP2ofCandidateEdges[mEdge]--;
                            toChangedEdges[mEdge]++;
                        }
                    }
                }
            }
            // v的邻居
            for (auto iter : graph[v])
            {
                int x = iter.first;
                if (x != u && graph[u].find(x) == graph[u].end())
                {
                    if (tau[v][x] < k - 1)
                        continue;
                    if (mNewK_1trussEdges.find(makeEdge(v, x)) != mNewK_1trussEdges.end() && mVisited.find(makeEdge(v, x)) != mVisited.end())
                        continue;
                    mVisited[makeEdge(v, x)] = true;
                    Edge mEdge = makeEdge(u, x);
                    if (mCandidateEdges.find(mEdge) != mCandidateEdges.end())
                    {
                       
                        if (mP2ofCandidateEdges[mEdge] > 0)
                        {
                            mP2ofCandidateEdges[mEdge]--;
                            toChangedEdges[mEdge]++;
                        }
                    }
                }
            }
        };
        // 当e退出truss值不为k的插入边集合时，将e原来的侯选边的p1值减1，如果p1-值为0，则从候选边集合中删除该边
        auto updateCandidatesP1_v1 = [&](Edge e, USetEdge &mCandidateEdges, UMapEdgeInt &mP1ofCandidateEdges, UMapEdgeInt &mP2ofCandidateEdges, USetEdge &mNewKtrussEdges, UMapEdgeInt &toChangedEdges)
        {
           
            int u = e.u, v = e.v;
            // u的邻居
            for (auto iter : graph[u])
            {
                int x = iter.first;
                // 另一边原本的trus值>=k
                if (x != v && graph[v].find(x) == graph[v].end() && tau[u][x] >= k && mNewKtrussEdges.find(makeEdge(u, x)) == mNewKtrussEdges.end())
                {
                    Edge mEdge = makeEdge(v, x);
                    if (mCandidateEdges.find(mEdge) != mCandidateEdges.end())
                    {
                        mP1ofCandidateEdges[mEdge]--;
                        toChangedEdges[mEdge]--;
                        if (mP1ofCandidateEdges[mEdge] == 0)
                        {
                            mCandidateEdges.erase(mEdge);
                            mP1ofCandidateEdges.erase(mEdge);
                            mP2ofCandidateEdges.erase(mEdge);
                        }
                    }
                }
            }
            // v的邻居
            for (auto iter : graph[v])
            {
                int x = iter.first;

                if (x != u && graph[u].find(x) == graph[u].end() && tau[v][x] >= k && mNewKtrussEdges.find(makeEdge(v, x)) == mNewKtrussEdges.end())
                {
                    Edge mEdge = makeEdge(u, x);
                    if (mCandidateEdges.find(mEdge) != mCandidateEdges.end())
                    {
                       
                        mP1ofCandidateEdges[mEdge]--;
                        toChangedEdges[mEdge]--;
                        if (mP1ofCandidateEdges[mEdge] == 0)
                        {
                            mCandidateEdges.erase(mEdge);
                            mP1ofCandidateEdges.erase(mEdge);
                            mP2ofCandidateEdges.erase(mEdge);
                        }
                    }
                }
            }
        };
        auto updateCandidatesP1_v2 = [&](Edge e, USetEdge &mCandidateEdges, UMapEdgeInt &mP1ofCandidateEdges, set<Edge> &mNotInKTrussEdges, UMapEdgeInt &toChangedEdges)
        {
            
            int u = e.u, v = e.v;
            // u的邻居
            for (auto iter : graph[u])
            {
                int x = iter.first;
                // 另一边还在notInKtrussedges中
                if (x != v && graph[v].find(x) == graph[v].end() && mNotInKTrussEdges.find(makeEdge(u, x)) != mNotInKTrussEdges.end())
                {
                    Edge mEdge = makeEdge(v, x);
                    if (graph[u].find(x) == graph[u].end())
                    {
                        
                        mP1ofCandidateEdges[mEdge]++;
                        toChangedEdges[mEdge]++;
                    }
                }
            }
            // v的邻居
            for (auto iter : graph[v])
            {
                int x = iter.first;

                if (x != u && graph[u].find(x) == graph[u].end() && mNotInKTrussEdges.find(makeEdge(v, x)) != mNotInKTrussEdges.end())
                {
                    Edge mEdge = makeEdge(u, x);
                    if (graph[u].find(x) == graph[u].end())
                    {
                      
                        mP1ofCandidateEdges[mEdge]++;
                        toChangedEdges[mEdge]++;
                    }
                }
            }
        };
        auto updateCandidates = [&](Edge mInsertEdge, set<Edge> &mNotInKTrussEdges, USetEdge &mCandidateEdges, priority_queue<pair<int, Edge>> &mCandidateEdgeQueue, UMapEdgeInt &mP1ofCandidateEdges, UMapEdgeInt &mP2ofCandidateEdges, USetEdge &mNewKtrussEdges, USetEdge &mNewK_1TrussEdges)
        {
            
            UMapEdgeInt toChangedEdges; // 记录需要更新的边的p1-p2,以便之后放到queue中
            if (tau[mInsertEdge.u][mInsertEdge.v] >= k)
            {
                // 使得mNotInKTrussEdges的插入边truss值达到k，从而退出mNotInKTrussEdges
                UMapEdgeBool mErased;
                for (auto e : mNotInKTrussEdges)
                {
                    if (tau[e.u][e.v] >= k)
                    {
                        mErased[e] = true;
                    }
                }

                for (auto e : mErased)
                {
                    mNotInKTrussEdges.erase(e.first);
                    updateCandidatesP1_v1(e.first, mCandidateEdges, mP1ofCandidateEdges, mP2ofCandidateEdges, mNewKtrussEdges, toChangedEdges);
                }
                // 使得侯选边的p+值增加
                for (auto e : mNewKtrussEdges)
                {
                    updateCandidatesP1_v2(e, mCandidateEdges, mP1ofCandidateEdges, mNotInKTrussEdges, toChangedEdges);
                }
                updateCandidatesP2_v1(mInsertEdge, mInsertEdge, mCandidateEdges, mP2ofCandidateEdges, mNewK_1TrussEdges, toChangedEdges);
            }
            else
            {
                mNotInKTrussEdges.insert(mInsertEdge);
                // 将该边的侯选边加入候选边集合，并计算goodness
                getCandidates_v2(mInsertEdge, mCandidateEdges, mCandidateEdgeQueue, mP1ofCandidateEdges, mP2ofCandidateEdges, toChangedEdges);
            }
            for (auto e : mNewK_1TrussEdges)
            {
                updateCandidatesP2_v1(mInsertEdge, e, mCandidateEdges, mP2ofCandidateEdges, mNewK_1TrussEdges, toChangedEdges);
            }
            for (auto e : toChangedEdges)
            {
                if (e.second > 0)
                {
                    mCandidateEdgeQueue.push({mP1ofCandidateEdges[e.first] - mP2ofCandidateEdges[e.first], e.first});
                    // 第一种情况引入的新边
                    if (mCandidateEdges.find(e.first) == mCandidateEdges.end())
                    {
                        mCandidateEdges.insert(e.first);
                        int mP1 = calP1(e.first), mP2 = calP2(e.first);
                        mP1ofCandidateEdges[e.first] = mP1;
                        mP2ofCandidateEdges[e.first] = mP2;
                    }
                    
                }
            }
        };

        int alreadyIsInKTruss = 2; // 新插入的，一个端点是w，一个端点在k_truss_community中，且truss值>=k 的边个数

        USetEdge mCandidateEdges;
        priority_queue<pair<int, Edge>> mCandidateEdgeQueue;
        UMapEdgeInt mP1ofCandidateEdges;
        UMapEdgeInt mP2ofCandidateEdges;

        for (auto mEdge : mNotInKTrussEdges)
        {
            getCandidates(mEdge, mCandidateEdges, mCandidateEdgeQueue, mP1ofCandidateEdges, mP2ofCandidateEdges);
        }
        while (true)
        {
            // 种子边已被加入k_truss_community，终止
            if (tau[w][mSeedu] >= k && tau[w][mSeedv] >= k)
            {
                break;
            }
            // 所有侯选边中，p1-p2最大的那条边作为本轮的插入边
            Edge mInsertEdge;
            while (!mCandidateEdgeQueue.empty())
            {
                Edge mEdge = mCandidateEdgeQueue.top().second;
                mCandidateEdgeQueue.pop();

                if (mCandidateEdges.find(mEdge) != mCandidateEdges.end())
                {
                    mInsertEdge = mEdge;
                    break;
                }
            }
            
            USetEdge mNewKtrussEdges, mNewK_1TrussEdges;
            allTrussUpdateWithInsertion(mInsertEdge, mRecoverTau, mNewKtrussEdges, mNewK_1TrussEdges);
            mAllInsertedEdges.insert(mInsertEdge);

            // 插入边后，将该边从侯选边中删去
            mCandidateEdges.erase(mInsertEdge);

            updateCandidates(mInsertEdge, mNotInKTrussEdges, mCandidateEdges, mCandidateEdgeQueue, mP1ofCandidateEdges, mP2ofCandidateEdges, mNewKtrussEdges, mNewK_1TrussEdges);

            // 插入边集合已被处理，不再处理
            if (mAllInsertedEdges.size() > 0 && visitEdgeSet(mAllInsertedEdges))
            {
                mIsVisited = true;
               
                break;
            }

            // 扩展中>最小插入边，终止
            if (newEdges.size() > 0 && newEdges.size() <= mAllInsertedEdges.size())
            {
                mTerminateFlag = true;
                
                break;
            }
        }

        //  恢复
        for (auto it : mRecoverTau)
        {
            int u = it.first.u, v = it.first.v, tau_uv = it.second;
            tau[u][v] = tau_uv;
            tau[v][u] = tau_uv;
        }
        for (auto e : mAllInsertedEdges)
        {
            int u = e.u, v = e.v;
            graph[u].erase(v);
            graph[v].erase(u);
            tau[u].erase(v);
            tau[v].erase(u);
        }

        // 已被访问
        if (mIsVisited || mTerminateFlag)
            continue;
       
        for (auto e : mAllInsertedEdges)
        {
            if (mNotInKTrussEdges.find(e) != mNotInKTrussEdges.end())
            {
                mAllInsertedEdges.erase(e);
            }
        }
       
        // 更新newEdges
        if (newEdges.size() == 0 || newEdges.size() > mAllInsertedEdges.size())
        {
            newEdges = vector<Edge>(mAllInsertedEdges.begin(), mAllInsertedEdges.end());
        }
        if (newEdges.size() == 1)
        {
            break;
        }
    }
    // out.close();
}

void myAlgorithm::mDebug()
{

}

USetEdge myAlgorithm::calculateDependentEdges(Edge e, USetEdge &newComers)
{
    clock_t start, end;
    start = clock();
    USetEdge infectedEdges;
    int u0 = e.u, v0 = e.v;
    // 改变e0
    graph[u0].erase(v0);
    graph[v0].erase(u0);

    tau[u0].erase(v0);
    tau[v0].erase(u0);

    infectedEdges.insert(makeEdge(u0, v0));
  
    // 将与e0形成三角形的边，truss值<=tau[eo][v0]的边加入queue
    USetEdge inL;  // 所有边只入一次队列
    queue<Edge> q; // 记录可能truss值减少1的边
    vector<int> neighbors = getNeighbors(u0, v0);

    for (auto w : neighbors)
    {
        int minn = min(tau[w][u0], tau[w][v0]);
        if (minn == k)
        {
            if (tau[w][u0] == minn && newComers.find(makeEdge(w, u0)) != newComers.end())
            {
                q.push(makeEdge(w, u0));
                inL.insert(makeEdge(w, u0));
            }
            if (tau[w][v0] == minn && newComers.find(makeEdge(w, v0)) != newComers.end())
            {
                q.push(makeEdge(w, v0));
                inL.insert(makeEdge(w, v0));
            }
        }
    }

    UMapEdgeInt s; // support
    int tk = k;

    while (!q.empty())
    {
        int u = q.front().u;
        int v = q.front().v;
        q.pop();
        s[makeEdge(u, v)] = 0;
        vector<int> neighbors = getNeighbors(u, v);
        for (auto w : neighbors)
        {
            if (tau[w][v] >= tk && tau[w][u] >= tk)
            {
                s[makeEdge(u, v)]++;
            }
        }
        queue<Edge> supQueue;
        if (s[makeEdge(u, v)] < tk - 2)
        {
            supQueue.push(makeEdge(u, v));
        }

        while (!supQueue.empty())
        {
            int x = supQueue.front().u;
            int y = supQueue.front().v;
            supQueue.pop();
            tau[x][y] = tk - 1;
            tau[y][x] = tk - 1;
            // cout<<"truss reduce: "<<x<<" "<<y<<" "<<tau[x][y]<<endl;
            infectedEdges.insert(makeEdge(x, y));

            s[makeEdge(x, y)] = -1;
            vector<int> neighbors = getNeighbors(x, y);
            for (auto z : neighbors)
            {
                if (tau[x][z] < tk || tau[y][z] < tk)
                    continue;
                if (tau[x][z] == tk && newComers.find(makeEdge(x, z)) != newComers.end())
                {
                    if (s.find(makeEdge(x, z)) != s.end() && s[makeEdge(x, z)] != -1)
                    {
                        s[makeEdge(x, z)]--;
                        if (s[makeEdge(x, z)] == tk - 3)
                        {
                            supQueue.push(makeEdge(x, z));
                        }
                    }
                    if (inL.find(makeEdge(x, z)) == inL.end())
                    {
                        inL.insert(makeEdge(x, z));
                        q.push(makeEdge(x, z));
                    }
                }
                if (tau[y][z] == tk && newComers.find(makeEdge(y, z)) != newComers.end())
                {
                    if (s.find(makeEdge(y, z)) != s.end() && s[makeEdge(y, z)] != -1)
                    {
                        s[makeEdge(y, z)]--;
                        if (s[makeEdge(y, z)] == tk - 3)
                        {
                            supQueue.push(makeEdge(y, z));
                        }
                    }
                    if (inL.find(makeEdge(y, z)) == inL.end())
                    {
                        inL.insert(makeEdge(y, z));
                        q.push(makeEdge(y, z));
                    }
                }
            }
        }
    }

    end = clock();
    return infectedEdges;
}
/*
description:
    variation of the  k truss maximization algorithm
    only works when intimacy of w is 1
*/
void myAlgorithm::gtm()
{

    newEdges.clear();
    vector<Edge> temp = selectCandidates();
    USetEdge candidates(temp.begin(), temp.end());

    // 当侯选边全部插入而算法还没结束的时候提前退出
    // cout<<"candidates size: "<<candidates.size()<<endl;

    while (candidates.size())
    {
        if (checkWhetherTimeOut())
            return;
        Edge selectedEdge;
        int newComersNum = 0;
        for (auto candidate : candidates)
        {
            int u = candidate.u;
            int v = candidate.v;

            vector<Edge> newTau = kTrussUpdateWithInsertion(candidate);
            if (newTau.size() >= newComersNum)
            {
                selectedEdge = candidate;
                newComersNum = newTau.size();
            }

            // 还原
            for (int j = 0; j < newTau.size(); j++)
            {
                int u = newTau[j].u;
                int v = newTau[j].v;
                tau[u][v] -= 1;
                tau[v][u] -= 1;
            }
            tau[u].erase(v);
            tau[v].erase(u);
            graph[u].erase(v);
            graph[v].erase(u);
        }
        // cout<<candidates.size()<<endl;
        // 将本轮贡献最大的侯选边插入图中
        newEdges.push_back(selectedEdge);
        // cout<<"selected edge: "<<selectedEdge.u<<" "<<selectedEdge.v<<endl;

        allTrussUpdateWithInsertion(selectedEdge);
        candidates.erase(selectedEdge);
        bool isWInKTruss = false;
        vector<Edge> k_truss_community_new = queryProcess(k);
        // check whether w is in k_truss_community_new
        for (auto edge : k_truss_community_new)
        {
            if (edge.u == w || edge.v == w)
            {
                isWInKTruss = true;
                break;
            }
        }
        if (isWInKTruss)
            break;
    }
    // cout<<"budget: "<<budget<<endl;
    // check whether w is in k_truss_community_new
    bool isWInKTruss = false;
    vector<Edge> k_truss_community_new = queryProcess(k);
    for (auto edge : k_truss_community_new)
    {
        if (edge.u == w || edge.v == w)
        {
            isWInKTruss = true;
            break;
        }
    }

    // recover
    for (auto e : newEdges)
    {
        allTrussUpdateWithDeletion(e);
    }
    if (!isWInKTruss)
    {
        newEdges.clear();
    }
    // cout<<"newEdges size: "<<newEdges.size()<<endl;
}
bool myAlgorithm::checkWhetherTriangleConnected(vector<int> &outofFollowers)
{
    UMapEdgeInt visited;
    for (auto v : outofFollowers)
    {
        if (visited[makeEdge(w, v)])
            continue;
        visited[makeEdge(w, v)] = 1;

        queue<Edge> Q;
        Q.push(makeEdge(w, v));
        while (!Q.empty())
        {
            int x = Q.front().u;
            int y = Q.front().v;
            Q.pop();
            if (k_truss_community.find(makeEdge(x, y)) != k_truss_community.end())
            {
                return true;
            }
            vector<int> neighbors = getNeighbors(x, y);
            for (int i = 0; i < neighbors.size(); i++)
            {
                int z = neighbors[i];
                if (tau[x][z] >= k && tau[y][z] >= k)
                {
                    if (visited.find(makeEdge(x, z)) == visited.end())
                    {
                        Q.push(makeEdge(x, z));
                        visited[makeEdge(x, z)] = 1;
                    }
                    if (visited.find(makeEdge(y, z)) == visited.end())
                    {
                        Q.push(makeEdge(y, z));
                        visited[makeEdge(y, z)] = 1;
                    }
                }
            }
        }
    }
    return false;
}
USetEdge myAlgorithm::calculateDependentEdges(Edge e, USetEdge &newComers, vector<TRI> &mRecoverTau)
{
    // Timer timer("calculateDependentEdges");
    // cout<<"newComers: "<<endl;
  
    clock_t start, end;
    start = clock();

    USetEdge infectedEdges;
    int u0 = e.u, v0 = e.v;
    // 改变e0
    graph[u0].erase(v0);
    graph[v0].erase(u0);
    mRecoverTau.push_back(TRI{u0, v0, tau[u0][v0]});
    tau[u0].erase(v0);
    tau[v0].erase(u0);
    infectedEdges.insert(makeEdge(u0, v0));

    // 将与e0形成三角形的边，truss值<=tau[eo][v0]的边加入queue
    USetEdge inL;  // 所有边只入一次队列
    queue<Edge> q; // 记录可能truss值减少1的边
    vector<int> neighbors = getNeighbors(u0, v0);

    for (auto w : neighbors)
    {
        int minn = min(tau[w][u0], tau[w][v0]);
        if (minn == k)
        {
            if (tau[w][u0] == minn && newComers.find(makeEdge(w, u0)) != newComers.end())
            {
                q.push(makeEdge(w, u0));
                inL.insert(makeEdge(w, u0));
            }
            if (tau[w][v0] == minn && newComers.find(makeEdge(w, v0)) != newComers.end())
            {
                q.push(makeEdge(w, v0));
                inL.insert(makeEdge(w, v0));
            }
        }
    }

    UMapEdgeInt s; // support
    int tk = k;

    while (!q.empty())
    {
        int u = q.front().u;
        int v = q.front().v;
        q.pop();
        s[makeEdge(u, v)] = 0;
        vector<int> neighbors = getNeighbors(u, v);
        for (auto w : neighbors)
        {
            if (tau[w][v] >= tk && tau[w][u] >= tk)
            {
                s[makeEdge(u, v)]++;
            }
        }
        queue<Edge> supQueue;
        if (s[makeEdge(u, v)] < tk - 2)
        {
            supQueue.push(makeEdge(u, v));
        }

        while (!supQueue.empty())
        {
            int x = supQueue.front().u;
            int y = supQueue.front().v;
            supQueue.pop();
            tau[x][y] = tk - 1;
            tau[y][x] = tk - 1;
            infectedEdges.insert(makeEdge(x, y));
            mRecoverTau.push_back(TRI{x, y, k});
            s[makeEdge(x, y)] = -1;
            vector<int> neighbors = getNeighbors(x, y);
            for (auto z : neighbors)
            {
                if (tau[x][z] < tk || tau[y][z] < tk)
                    continue;
                if (tau[x][z] == tk && newComers.find(makeEdge(x, z)) != newComers.end())
                {
                    if (s.find(makeEdge(x, z)) != s.end() && s[makeEdge(x, z)] != -1)
                    {
                        s[makeEdge(x, z)]--;
                        if (s[makeEdge(x, z)] == tk - 3)
                        {
                            supQueue.push(makeEdge(x, z));
                        }
                    }
                    if (inL.find(makeEdge(x, z)) == inL.end())
                    {
                        inL.insert(makeEdge(x, z));
                        q.push(makeEdge(x, z));
                    }
                }
                if (tau[y][z] == tk && newComers.find(makeEdge(y, z)) != newComers.end())
                {
                    if (s.find(makeEdge(y, z)) != s.end() && s[makeEdge(y, z)] != -1)
                    {
                        s[makeEdge(y, z)]--;
                        if (s[makeEdge(y, z)] == tk - 3)
                        {
                            supQueue.push(makeEdge(y, z));
                        }
                    }
                    if (inL.find(makeEdge(y, z)) == inL.end())
                    {
                        inL.insert(makeEdge(y, z));
                        q.push(makeEdge(y, z));
                    }
                }
            }
        }
    }

    end = clock();
    // cout<<"c1 c2: "<<c1<<" "<<c2<<endl;
    // if (double(end - start) / CLOCKS_PER_SEC * 1000 > 200)
    // cout << "边删除时间新: " << double(end - start) / CLOCKS_PER_SEC << endl;

    return infectedEdges;
}
vector<Edge> myAlgorithm::queryProcess(int k)
{
    UMapEdgeInt visited;

    int v = (*k_truss_community.begin()).u;
    int q = (*k_truss_community.begin()).v;
    Edge et = makeEdge(q, v);
    vector<Edge> community;
    if (tau[q][v] >= k && visited.find(et) == visited.end())
    {
        queue<Edge> Q;
        Q.push(et);
        visited[et] = 1;
        while (!Q.empty())
        {
            int x = Q.front().u;
            int y = Q.front().v;
            Q.pop();
            community.push_back({x, y});
            vector<int> neighbors = getNeighbors(x, y);
            for (int i = 0; i < neighbors.size(); i++)
            {
                int z = neighbors[i];
                if (tau[x][z] >= k && tau[y][z] >= k)
                {
                    if (visited.find(makeEdge(x, z)) == visited.end())
                    {
                        Q.push(makeEdge(x, z));
                        visited[makeEdge(x, z)] = 1;
                    }
                    if (visited.find(makeEdge(y, z)) == visited.end())
                    {
                        Q.push(makeEdge(y, z));
                        visited[makeEdge(y, z)] = 1;
                    }
                }
            }
        }
    }
    return community;
}

int myAlgorithm::getNewEdgesSize()
{
    return newEdges.size();
}

vector<Edge> myAlgorithm::selectCandidates()
{
    vector<Edge> candidate;
    vector<Edge> Ck_1 = queryProcess(k - 1);

    UMapEdgeBool isInCandidate;
    vector<Edge> L_k_1;
    for (auto e : Ck_1)
    {
        int u = e.u, v = e.v;
        if (tau[u][v] > k - 1)
            continue;
        // check  whether edge(v,u) is unstable
        int sup = 0;
        vector<int> neighbors = getNeighbors(u, v);
        for (auto ww : neighbors)
        {
            if (tau[ww][u] >= k - 1 && tau[ww][v] >= k - 1)
                sup++;
        }
        if (sup != k - 3)
            continue;
        for (auto iter : graph[u])
        {
            int ww = iter.first;
            if (ww == v)
                continue;
            if (graph[ww].find(v) == graph[ww].end() && tau[u][ww] >= k - 1 && isInCandidate.find(makeEdge(v, ww)) == isInCandidate.end())
            {
                candidate.push_back(makeEdge(ww, v));
                isInCandidate[makeEdge(ww, v)] = true;
            }
        }
        for (auto iter : graph[v])
        {
            int ww = iter.first;
            if (ww == u)
                continue;
            if (graph[ww].find(u) == graph[ww].end() && tau[v][ww] >= k - 1 && isInCandidate.find(makeEdge(u, ww)) == isInCandidate.end())
            {
                candidate.push_back(makeEdge(ww, u));
                isInCandidate[makeEdge(ww, u)] = true;
            }
        }
    }

    // printEdges("candidate", candidate);
    // cout << "hi " << candidate.size() << endl;
    // 根据权重筛选侯选边
    vector<Edge> candidate2;
    for (int i = 0; i < candidate.size(); i++)
    {
        int u = candidate[i].u;
        int v = candidate[i].v;
        int lambda = 0;
        vector<int> neighbors = getNeighbors(u, v);
        for (int j = 0; j < neighbors.size(); j++)
        {
            int ww = neighbors[j];
            if (tau[ww][u] >= k - 1 && tau[ww][v] >= k - 1)
            {
                lambda++;
            }
        }
        if (lambda >= k - 2)
        {
            candidate2.push_back(candidate[i]);
        }
    }
    // printEdges("candidate", candidate2);
    return candidate2;
}

bool myAlgorithm::checkWhetherTimeOut()
{
    if (terminateFlag.load())
    {
        return true;
    }
    return false;
}

void myAlgorithm::printNewCommunity(string &algorithm_name)
{
    fstream outFile("newCommunity_" + algorithm_name + ".txt", ios::out);
    vector<Edge> Ck_1 = queryProcess(k);
    cout << Ck_1.size() << endl;
    for (auto e : Ck_1)
    {
        if (k_truss_community.find(e) == k_truss_community.end())
        {
            outFile << e.u << " " << e.v << endl;
        }
    }
    fstream outFile2("newEdges_" + algorithm_name + ".txt", ios::out);
    for (auto e : newEdges)
    {
        outFile2 << e.u << " " << e.v << endl;
    }
   
    outFile.close();
    outFile2.close();
    // outFile3.close();
}

void myAlgorithm::verify(string algorithm_name)
{
    cout << "----------------------------------Verification----------------------------------" << endl;

    vector<Edge> edges = newEdges;
    cout << "number of inserted edges: " << edges.size() << endl;
    printEdges("they are ", edges);
    // printEdges("inserted edges:", edges);
    // 保存副本
    vector<unordered_map<int, int>> graph1 = this->graph, tau1 = this->tau;
    // 将edges插入图中
    for (int i = 0; i < edges.size(); i++)
    {
        Edge e = edges[i];
        int u = e.u;
        int v = e.v;
        allTrussUpdateWithInsertion(e);
    }

    // printVector("k_truss_vertices: ",k_truss_vertices);
    int f = 0;
    vector<Edge> ans = queryProcess(k);
    for (auto e : ans)
    {
        if (e.u == w || e.v == w)
        {
            f = 1;
            break;
        }
    }
    DEBUG_ASSERT(f == 1, "error!");

    cout << "correct! w is in k_truss_community" << endl;

    //  printNewCommunity(algorithm_name);

    // 还原图结构
    graph = graph1;
    tau = tau1;
}
void myAlgorithm::allTrussUpdateWithDeletion(Edge e0)
{
    clock_t start, end;
    start = clock();

    int u0 = e0.u, v0 = e0.v;
    // cout << "delete u,v: " << u0 << " " << v0 << endl;/
    int kMax = tau[u0][v0];
    // 改变e0
    graph[u0].erase(v0);
    graph[v0].erase(u0);
    tau[u0].erase(v0);
    tau[v0].erase(u0);

    vector<vector<Edge>> L(kMax + 1);

    // 将与e0形成三角形的边，truss值<=tau[eo][v0]的边加入queue
    USetEdge inL; // 所有边只入一次队列
    vector<int> neighbors = getNeighbors(u0, v0);

    for (auto w : neighbors)
    {
        int minn = min(tau[w][u0], tau[w][v0]);
        if (minn <= kMax)
        {
            if (tau[w][u0] == minn)
            {
                L[minn].push_back(makeEdge(w, u0));
                inL.insert(makeEdge(w, u0));
            }
            if (tau[w][v0] == minn)
            {
                L[minn].push_back(makeEdge(w, v0));
                inL.insert(makeEdge(w, v0));
            }
        }
    }

    UMapEdgeInt s; // support
    for (int tk = 3; tk <= kMax; tk++)
    {
        queue<Edge> q; // 记录可能truss值减少1的边
        for (auto e : L[tk])
            q.push(e);
        while (!q.empty())
        {
            int u = q.front().u;
            int v = q.front().v;
            q.pop();
            s[makeEdge(u, v)] = 0;
            vector<int> neighbors = getNeighbors(u, v);
            for (auto w : neighbors)
            {
                if (tau[w][v] >= tk && tau[w][u] >= tk)
                {
                    s[makeEdge(u, v)]++;
                }
            }
            queue<Edge> supQueue;
            if (s[makeEdge(u, v)] < tk - 2)
            {
                supQueue.push(makeEdge(u, v));
            }

            while (!supQueue.empty())
            {
                int x = supQueue.front().u;
                int y = supQueue.front().v;
                supQueue.pop();
                // cout << "truss reduce: " << x << " " << y << " " << tau[x][y] << endl;

                tau[x][y]--;
                tau[y][x]--;
                s[makeEdge(x, y)] = -1;
                vector<int> neighbors = getNeighbors(x, y);
                for (auto z : neighbors)
                {
                    if (tau[x][z] < tk || tau[y][z] < tk)
                        continue;
                    if (tau[x][z] == tk)
                    {
                        if (s.find(makeEdge(x, z)) != s.end() && s[makeEdge(x, z)] != -1)
                        {
                            s[makeEdge(x, z)]--;
                            if (s[makeEdge(x, z)] == tk - 3)
                            {
                                supQueue.push(makeEdge(x, z));
                            }
                        }
                        if (inL.find(makeEdge(x, z)) == inL.end())
                        {
                            inL.insert(makeEdge(x, z));
                            q.push(makeEdge(x, z));
                        }
                    }
                    if (tau[y][z] == tk)
                    {
                        if (s.find(makeEdge(y, z)) != s.end() && s[makeEdge(y, z)] != -1)
                        {
                            s[makeEdge(y, z)]--;
                            if (s[makeEdge(y, z)] == tk - 3)
                            {
                                supQueue.push(makeEdge(y, z));
                            }
                        }
                        if (inL.find(makeEdge(y, z)) == inL.end())
                        {
                            inL.insert(makeEdge(y, z));
                            q.push(makeEdge(y, z));
                        }
                    }
                }
            }
        }
    }
    end = clock();
    // cout<<"c1 c2: "<<c1<<" "<<c2<<endl;
    // if (double(end - start) / CLOCKS_PER_SEC * 1000 > 200)
    //     cout << "边删除时间: " << double(end - start) / CLOCKS_PER_SEC << endl;
    return;
}
/*
功能：模拟e插入，得到所有truss值改变的边，并改变了e0的graph和影响边的tau
输出：插入边
输出：所有truss值变化的边
*/
void myAlgorithm::allTrussUpdateWithInsertion(Edge e0)
{
    clock_t start, end;
    start = clock();

    // printGraph();
    // printTau();
    int k1 = 2, k2 = 3; // e0 插入后的truss上下界
    calbound(e0, k1, k2);
    int u = e0.u;
    int v = e0.v;

    // 插入e0
    tau[u][v] = k1; // 记得改回来
    tau[v][u] = k1;
    graph[u][v] = 0;
    graph[v][u] = 0;
    // cout<<k1<<" "<<k2<<endl;

    int kmax = k2 - 1;

    // 相邻的放入L
    vector<vector<Edge>> L(kmax + 1);
    USetEdge inQueue;
    vector<int> wList = getNeighbors(u, v);

    // e0 的truss值也可能更新
    if (k1 == kmax)
    {
        L[k1].push_back(makeEdge(u, v));
    }
    // cout<<wList.size()<<endl;
    for (int i = 0; i < wList.size(); i++)
    {
        int w = wList[i];
        int tk = min(tau[w][u], tau[w][v]);
        if (tk <= kmax)
        {
            if (tau[w][u] == tk && inQueue.find(makeEdge(w, u)) == inQueue.end()) // this,mp?
            {
                L[tk].push_back(makeEdge(w, u));
                inQueue.insert(makeEdge(w, u));
            }
            if (tau[w][v] == tk && inQueue.find(makeEdge(w, v)) == inQueue.end())
            {
                L[tk].push_back(makeEdge(w, v));
                inQueue.insert(makeEdge(w, v));
            }
        }
    }
    // cout<<"hi2"<<endl;

    // 边连通的放入L
    for (int tk = kmax; tk >= 2; tk--)
    {
        queue<Edge> Q;
        UMapEdgeInt s;

        UMapEdgeInt mMovedToBin; // 记录是否移动到bin
        // cout<<"tk: "<<tk<<endl;
        for (int i = 0; i < L[tk].size(); i++)
        {
            Q.push(L[tk][i]);
        }
        while (!Q.empty())
        {

            int x = Q.front().u;
            int y = Q.front().v;
            Q.pop();
            // cout<<"x y "<<x<<" "<<y<<endl;

            vector<int> wList = getNeighbors(x, y);
            vector<Edge> tempEdges;

            for (int i = 0; i < wList.size(); i++)
            {
                int z = wList[i];
                if ((tau[z][x] < tk) || (tau[z][y] < tk))
                    continue;
                if (tau[z][x] == tk && mMovedToBin[makeEdge(z, x)] == 2)
                    continue;
                if (tau[z][y] == tk && mMovedToBin[makeEdge(z, y)] == 2)
                    continue;
                s[{x, y}]++;
                if ((tau[z][x] == tk) && inQueue.find(makeEdge(z, x)) == inQueue.end())
                {
                    tempEdges.push_back(makeEdge(z, x));
                }
                if ((tau[z][y] == tk) && inQueue.find(makeEdge(z, y)) == inQueue.end())
                {
                    tempEdges.push_back(makeEdge(z, y));
                }
            }
            mMovedToBin[makeEdge(x, y)] = 1;
            // cout<<"x,y,s "<<x<<" "<<y<<" "<<s[{x, y}]<<endl;
            if (s[{x, y}] > tk - 2)
            {

                for (auto edge : tempEdges)
                {
                    if (inQueue.find(edge) == inQueue.end())
                    {
                        Q.push(edge);
                        inQueue.insert(edge);
                    }
                }
            }
            else
            {
                queue<Edge> Q1; // mToUpdatedEges里面s值<=tk-2的边
                Q1.push(makeEdge(x, y));

                while (!Q1.empty())
                {
                    int x1 = Q1.front().u, y1 = Q1.front().v;
                    mMovedToBin[makeEdge(x1, y1)] = 2;
                    // cout<<"x1 y1 "<<x1<<" "<<y1<<endl;
                    Q1.pop();
                    vector<int> wList = getNeighbors(x1, y1);
                    for (int i = 0; i < wList.size(); i++)
                    {
                        int z1 = wList[i];
                        if (tau[x1][z1] < tk || tau[y1][z1] < tk)
                            continue;
                        if (tau[x1][z1] == tk && mMovedToBin[makeEdge(x1, z1)] == 2)
                            continue;
                        if (tau[y1][z1] == tk && mMovedToBin[makeEdge(y1, z1)] == 2)
                            continue;
                        if (tau[x1][z1] == tk)
                        {
                            if (mMovedToBin[makeEdge(x1, z1)] == 1)
                            {
                                s[makeEdge(x1, z1)] = s[makeEdge(x1, z1)] - 1;
                                if (s[makeEdge(x1, z1)] == tk - 2)
                                {
                                    Q1.push(makeEdge(x1, z1));
                                }
                            }
                        }
                        if (tau[y1][z1] == tk)
                        {
                            if (mMovedToBin[makeEdge(y1, z1)] == 1)
                            {
                                s[makeEdge(y1, z1)] = s[makeEdge(y1, z1)] - 1;
                                if (s[makeEdge(y1, z1)] == tk - 2)
                                {
                                    Q1.push(makeEdge(y1, z1));
                                }
                            }
                        }
                    }
                }
            }
        }
        // cout << "tk: " << tk << " "
        //      << inQueue.size: " <<inQueue.size() << endl;
        for (auto edge : inQueue)
        {
            if (mMovedToBin[edge] == 1)
            {
                int u = edge.u;
                int v = edge.v;
                tau[u][v] += 1;
                tau[v][u] += 1;

              
            }

           
        }

       
    }
    end = clock();
    // if (double(end - start) / CLOCKS_PER_SEC * 1000 > 200)
    //     cout << "边插入时间: " << double(end - start) / CLOCKS_PER_SEC << endl;
}

void myAlgorithm::allTrussUpdateWithInsertion(Edge e0, UMapEdgeInt &mRecoverTau)
{
    clock_t start, end;
    start = clock();

    // printGraph();
    // printTau();
    int k1 = 2, k2 = 3; // e0 插入后的truss上下界
    calbound(e0, k1, k2);
    int u = e0.u;
    int v = e0.v;

    // 插入e0
    tau[u][v] = k1; // 记得改回来
    tau[v][u] = k1;
    graph[u][v] = 0;
    graph[v][u] = 0;
    // cout<<k1<<" "<<k2<<endl;

    int kmax = k2 - 1;

    // 相邻的放入L
    vector<vector<Edge>> L(kmax + 1);
    USetEdge inQueue;
    vector<int> wList = getNeighbors(u, v);

    // e0 的truss值也可能更新
    if (k1 == kmax)
    {
        L[k1].push_back(makeEdge(u, v));
    }
    // cout<<wList.size()<<endl;
    for (int i = 0; i < wList.size(); i++)
    {
        int w = wList[i];
        int tk = min(tau[w][u], tau[w][v]);
        if (tk <= kmax)
        {
            if (tau[w][u] == tk && inQueue.find(makeEdge(w, u)) == inQueue.end()) // this,mp?
            {
                L[tk].push_back(makeEdge(w, u));
                inQueue.insert(makeEdge(w, u));
            }
            if (tau[w][v] == tk && inQueue.find(makeEdge(w, v)) == inQueue.end())
            {
                L[tk].push_back(makeEdge(w, v));
                inQueue.insert(makeEdge(w, v));
            }
        }
    }
    // cout<<"hi2"<<endl;

    // 边连通的放入L
    for (int tk = kmax; tk >= 2; tk--)
    {
        queue<Edge> Q;
        UMapEdgeInt s;

        UMapEdgeInt mMovedToBin; // 记录是否移动到bin
        // cout<<"tk: "<<tk<<endl;
        for (int i = 0; i < L[tk].size(); i++)
        {
            Q.push(L[tk][i]);
        }
        while (!Q.empty())
        {

            int x = Q.front().u;
            int y = Q.front().v;
            Q.pop();
            // cout<<"x y "<<x<<" "<<y<<endl;

            vector<int> wList = getNeighbors(x, y);
            vector<Edge> tempEdges;

            for (int i = 0; i < wList.size(); i++)
            {
                int z = wList[i];
                if ((tau[z][x] < tk) || (tau[z][y] < tk))
                    continue;
                if (tau[z][x] == tk && mMovedToBin[makeEdge(z, x)] == 2)
                    continue;
                if (tau[z][y] == tk && mMovedToBin[makeEdge(z, y)] == 2)
                    continue;
                s[{x, y}]++;
                if ((tau[z][x] == tk) && inQueue.find(makeEdge(z, x)) == inQueue.end())
                {
                    tempEdges.push_back(makeEdge(z, x));
                }
                if ((tau[z][y] == tk) && inQueue.find(makeEdge(z, y)) == inQueue.end())
                {
                    tempEdges.push_back(makeEdge(z, y));
                }
            }
            mMovedToBin[makeEdge(x, y)] = 1;
            // cout<<"x,y,s "<<x<<" "<<y<<" "<<s[{x, y}]<<endl;
            if (s[{x, y}] > tk - 2)
            {

                for (auto edge : tempEdges)
                {
                    if (inQueue.find(edge) == inQueue.end())
                    {
                        Q.push(edge);
                        inQueue.insert(edge);
                    }
                }
            }
            else
            {
                queue<Edge> Q1; // mToUpdatedEges里面s值<=tk-2的边
                Q1.push(makeEdge(x, y));

                while (!Q1.empty())
                {
                    int x1 = Q1.front().u, y1 = Q1.front().v;
                    mMovedToBin[makeEdge(x1, y1)] = 2;
                    // cout<<"x1 y1 "<<x1<<" "<<y1<<endl;
                    Q1.pop();
                    vector<int> wList = getNeighbors(x1, y1);
                    for (int i = 0; i < wList.size(); i++)
                    {
                        int z1 = wList[i];
                        if (tau[x1][z1] < tk || tau[y1][z1] < tk)
                            continue;
                        if (tau[x1][z1] == tk && mMovedToBin[makeEdge(x1, z1)] == 2)
                            continue;
                        if (tau[y1][z1] == tk && mMovedToBin[makeEdge(y1, z1)] == 2)
                            continue;
                        if (tau[x1][z1] == tk)
                        {
                            if (mMovedToBin[makeEdge(x1, z1)] == 1)
                            {
                                s[makeEdge(x1, z1)] = s[makeEdge(x1, z1)] - 1;
                                if (s[makeEdge(x1, z1)] == tk - 2)
                                {
                                    Q1.push(makeEdge(x1, z1));
                                }
                            }
                        }
                        if (tau[y1][z1] == tk)
                        {
                            if (mMovedToBin[makeEdge(y1, z1)] == 1)
                            {
                                s[makeEdge(y1, z1)] = s[makeEdge(y1, z1)] - 1;
                                if (s[makeEdge(y1, z1)] == tk - 2)
                                {
                                    Q1.push(makeEdge(y1, z1));
                                }
                            }
                        }
                    }
                }
            }
        }
        // cout << "tk: " << tk << " "
        //      << inQueue.size: " <<inQueue.size() << endl;
        for (auto edge : inQueue)
        {
            if (mMovedToBin[edge] == 1)
            {
                int u = edge.u;
                int v = edge.v;
                if (mRecoverTau.find(makeEdge(u, v)) == mRecoverTau.end())
                {
                    mRecoverTau[makeEdge(u, v)] = tau[u][v];
                }
                tau[u][v] += 1;
                tau[v][u] += 1;

                // cout << "truss increase: " << u << " " << v << " " << tau[u][v] << endl;
            }

            // cout<<"u "<<u<<" v "<<v<<endl;
        }

        // printEdges(str, newTau);
    }
    end = clock();
    // if (double(end - start) / CLOCKS_PER_SEC * 1000 > 200)
    //     cout << "边插入时间: " << double(end - start) / CLOCKS_PER_SEC << endl;
}
void myAlgorithm::allTrussUpdateWithInsertion(Edge e0, UMapEdgeInt &mRecoverTau, USetEdge &newKtrussEdges, USetEdge &newK_1trussEdges)
{
    clock_t start, end;
    start = clock();

    // printGraph();
    // printTau();
    int k1 = 2, k2 = 3; // e0 插入后的truss上下界
    calbound(e0, k1, k2);
    int u = e0.u;
    int v = e0.v;
    // 插入e0
    tau[u][v] = k1; // 记得改回来
    tau[v][u] = k1;
    graph[u][v] = 0;
    graph[v][u] = 0;
    // cout<<k1<<" "<<k2<<endl;

    int kmax = k2 - 1;

    // 相邻的放入L
    vector<vector<Edge>> L(kmax + 1);
    USetEdge inQueue;
    vector<int> wList = getNeighbors(u, v);

    // e0 的truss值也可能更新
    if (k1 == kmax)
    {
        L[k1].push_back(makeEdge(u, v));
    }
    // cout<<wList.size()<<endl;
    for (int i = 0; i < wList.size(); i++)
    {
        int w = wList[i];
        int tk = min(tau[w][u], tau[w][v]);
        if (tk <= kmax)
        {
            if (tau[w][u] == tk && inQueue.find(makeEdge(w, u)) == inQueue.end()) // this,mp?
            {
                L[tk].push_back(makeEdge(w, u));
                inQueue.insert(makeEdge(w, u));
            }
            if (tau[w][v] == tk && inQueue.find(makeEdge(w, v)) == inQueue.end())
            {
                L[tk].push_back(makeEdge(w, v));
                inQueue.insert(makeEdge(w, v));
            }
        }
    }
    // cout<<"hi2"<<endl;

    // 边连通的放入L
    for (int tk = kmax; tk >= 2; tk--)
    {
        queue<Edge> Q;
        UMapEdgeInt s;

        UMapEdgeInt mMovedToBin; // 记录是否移动到bin
        // cout<<"tk: "<<tk<<endl;
        for (int i = 0; i < L[tk].size(); i++)
        {
            Q.push(L[tk][i]);
        }
        while (!Q.empty())
        {

            int x = Q.front().u;
            int y = Q.front().v;
            Q.pop();
            // cout<<"x y "<<x<<" "<<y<<endl;

            vector<int> wList = getNeighbors(x, y);
            vector<Edge> tempEdges;

            for (int i = 0; i < wList.size(); i++)
            {
                int z = wList[i];
                if ((tau[z][x] < tk) || (tau[z][y] < tk))
                    continue;
                if (tau[z][x] == tk && mMovedToBin[makeEdge(z, x)] == 2)
                    continue;
                if (tau[z][y] == tk && mMovedToBin[makeEdge(z, y)] == 2)
                    continue;
                s[{x, y}]++;
                if ((tau[z][x] == tk) && inQueue.find(makeEdge(z, x)) == inQueue.end())
                {
                    tempEdges.push_back(makeEdge(z, x));
                }
                if ((tau[z][y] == tk) && inQueue.find(makeEdge(z, y)) == inQueue.end())
                {
                    tempEdges.push_back(makeEdge(z, y));
                }
            }
            mMovedToBin[makeEdge(x, y)] = 1;
            // cout<<"x,y,s "<<x<<" "<<y<<" "<<s[{x, y}]<<endl;
            if (s[{x, y}] > tk - 2)
            {

                for (auto edge : tempEdges)
                {
                    if (inQueue.find(edge) == inQueue.end())
                    {
                        Q.push(edge);
                        inQueue.insert(edge);
                    }
                }
            }
            else
            {
                queue<Edge> Q1; // mToUpdatedEges里面s值<=tk-2的边
                Q1.push(makeEdge(x, y));

                while (!Q1.empty())
                {
                    int x1 = Q1.front().u, y1 = Q1.front().v;
                    mMovedToBin[makeEdge(x1, y1)] = 2;
                    // cout<<"x1 y1 "<<x1<<" "<<y1<<endl;
                    Q1.pop();
                    vector<int> wList = getNeighbors(x1, y1);
                    for (int i = 0; i < wList.size(); i++)
                    {
                        int z1 = wList[i];
                        if (tau[x1][z1] < tk || tau[y1][z1] < tk)
                            continue;
                        if (tau[x1][z1] == tk && mMovedToBin[makeEdge(x1, z1)] == 2)
                            continue;
                        if (tau[y1][z1] == tk && mMovedToBin[makeEdge(y1, z1)] == 2)
                            continue;
                        if (tau[x1][z1] == tk)
                        {
                            if (mMovedToBin[makeEdge(x1, z1)] == 1)
                            {
                                s[makeEdge(x1, z1)] = s[makeEdge(x1, z1)] - 1;
                                if (s[makeEdge(x1, z1)] == tk - 2)
                                {
                                    Q1.push(makeEdge(x1, z1));
                                }
                            }
                        }
                        if (tau[y1][z1] == tk)
                        {
                            if (mMovedToBin[makeEdge(y1, z1)] == 1)
                            {
                                s[makeEdge(y1, z1)] = s[makeEdge(y1, z1)] - 1;
                                if (s[makeEdge(y1, z1)] == tk - 2)
                                {
                                    Q1.push(makeEdge(y1, z1));
                                }
                            }
                        }
                    }
                }
            }
        }
        // cout << "tk: " << tk << " "
        //      << inQueue.size: " <<inQueue.size() << endl;
        for (auto edge : inQueue)
        {
            if (mMovedToBin[edge] == 1)
            {
                int u = edge.u;
                int v = edge.v;
                if (mRecoverTau.find(makeEdge(u, v)) == mRecoverTau.end())
                {
                    mRecoverTau[makeEdge(u, v)] = tau[u][v];
                }
                tau[u][v] += 1;
                tau[v][u] += 1;
                if (tau[u][v] == k)
                {
                    newKtrussEdges.insert(makeEdge(u, v));
                }
                if (tau[u][v] == k - 1)
                {
                    newK_1trussEdges.insert(makeEdge(u, v));
                }

                // cout << "truss increase: " << u << " " << v << " " << tau[u][v] << endl;
            }

            // cout<<"u "<<u<<" v "<<v<<endl;
        }

        // printEdges(str, newTau);
    }
    end = clock();
    // if (double(end - start) / CLOCKS_PER_SEC * 1000 > 200)
    //     cout << "边插入时间: " << double(end - start) / CLOCKS_PER_SEC << endl;
}
vector<Edge> myAlgorithm::kTrussUpdateWithInsertion(Edge e0)
{
    vector<Edge> trussUpdatedEdges;
    clock_t start, end;
    start = clock();

    int k1 = 2, k2 = 3; // e0 插入后的truss上下界
    calbound(e0, k1, k2);
    int u = e0.u;
    int v = e0.v;

    // 插入e0
    tau[u][v] = k1; // 记得改回来
    tau[v][u] = k1;
    graph[u][v] = 0;
    graph[v][u] = 0;
    // cout<<k1<<" "<<k2<<endl;

    int kmax = k2 - 1;

    // 相邻的放入L
    vector<Edge> L;
    USetEdge inQueue;
    vector<int> wList = getNeighbors(u, v);

    // e0 的truss值也可能更新
    if (k1 == kmax)
    {
        L.push_back(makeEdge(u, v));
    }
    // cout<<wList.size()<<endl;
    for (int i = 0; i < wList.size(); i++)
    {
        int w = wList[i];
        int tk = min(tau[w][u], tau[w][v]);
        if (tk <= kmax && tk == k - 1)
        {
            if (tau[w][u] == tk && inQueue.find(makeEdge(w, u)) == inQueue.end()) // this,mp?
            {
                L.push_back(makeEdge(w, u));
                inQueue.insert(makeEdge(w, u));
            }
            if (tau[w][v] == tk && inQueue.find(makeEdge(w, v)) == inQueue.end())
            {
                L.push_back(makeEdge(w, v));
                inQueue.insert(makeEdge(w, v));
            }
        }
    }
    // cout<<"hi2"<<endl;

    // 边连通的放入L
    int tk = k - 1;
    queue<Edge> Q;
    UMapEdgeInt s;

    UMapEdgeInt mMovedToBin; // 记录是否移动到bin
    // cout<<"tk: "<<tk<<endl;
    for (int i = 0; i < L.size(); i++)
    {
        Q.push(L[i]);
    }
    while (!Q.empty())
    {

        int x = Q.front().u;
        int y = Q.front().v;
        Q.pop();
        // cout<<"x y "<<x<<" "<<y<<endl;

        vector<int> wList = getNeighbors(x, y);
        vector<Edge> tempEdges;

        for (int i = 0; i < wList.size(); i++)
        {
            int z = wList[i];
            if ((tau[z][x] < tk) || (tau[z][y] < tk))
                continue;
            if (tau[z][x] == tk && mMovedToBin[makeEdge(z, x)] == 2)
                continue;
            if (tau[z][y] == tk && mMovedToBin[makeEdge(z, y)] == 2)
                continue;
            s[{x, y}]++;
            if ((tau[z][x] == tk) && inQueue.find(makeEdge(z, x)) == inQueue.end())
            {
                tempEdges.push_back(makeEdge(z, x));
            }
            if ((tau[z][y] == tk) && inQueue.find(makeEdge(z, y)) == inQueue.end())
            {
                tempEdges.push_back(makeEdge(z, y));
            }
        }
        mMovedToBin[makeEdge(x, y)] = 1;
        // cout<<"x,y,s "<<x<<" "<<y<<" "<<s[{x, y}]<<endl;
        if (s[{x, y}] > tk - 2)
        {

            for (auto edge : tempEdges)
            {
                if (inQueue.find(edge) == inQueue.end())
                {
                    Q.push(edge);
                    inQueue.insert(edge);
                }
            }
        }
        else
        {
            queue<Edge> Q1; // mToUpdatedEges里面s值<=tk-2的边
            Q1.push(makeEdge(x, y));

            while (!Q1.empty())
            {
                int x1 = Q1.front().u, y1 = Q1.front().v;
                mMovedToBin[makeEdge(x1, y1)] = 2;
                // cout<<"x1 y1 "<<x1<<" "<<y1<<endl;
                Q1.pop();
                vector<int> wList = getNeighbors(x1, y1);
                for (int i = 0; i < wList.size(); i++)
                {
                    int z1 = wList[i];
                    if (tau[x1][z1] < tk || tau[y1][z1] < tk)
                        continue;
                    if (tau[x1][z1] == tk && mMovedToBin[makeEdge(x1, z1)] == 2)
                        continue;
                    if (tau[y1][z1] == tk && mMovedToBin[makeEdge(y1, z1)] == 2)
                        continue;
                    if (tau[x1][z1] == tk)
                    {
                        if (mMovedToBin[makeEdge(x1, z1)] == 1)
                        {
                            s[makeEdge(x1, z1)] = s[makeEdge(x1, z1)] - 1;
                            if (s[makeEdge(x1, z1)] == tk - 2)
                            {
                                Q1.push(makeEdge(x1, z1));
                            }
                        }
                    }
                    if (tau[y1][z1] == tk)
                    {
                        if (mMovedToBin[makeEdge(y1, z1)] == 1)
                        {
                            s[makeEdge(y1, z1)] = s[makeEdge(y1, z1)] - 1;
                            if (s[makeEdge(y1, z1)] == tk - 2)
                            {
                                Q1.push(makeEdge(y1, z1));
                            }
                        }
                    }
                }
            }
        }
        for (auto edge : inQueue)
        {
            if (mMovedToBin[edge] == 1)
            {
                int u = edge.u;
                int v = edge.v;
                tau[u][v] += 1;
                tau[v][u] += 1;
                trussUpdatedEdges.push_back(edge);
            }
        }
    }
    end = clock();
    return trussUpdatedEdges;
}
void myAlgorithm::printEdges(string message, vector<Edge> &edges)
{
    cout << message << endl;
    for (int i = 0; i < edges.size(); i++)
    {
        cout << " ( " << edges[i].u << " , " << edges[i].v << " ) " << endl;
    }
}
void myAlgorithm::printTau()
{
    cout << "tau:" << endl;
    for (int i = 0; i < n; i++)
    {
        for (unordered_map<int, int>::iterator it = tau[i].begin(); it != tau[i].end(); it++)
        {
            cout << i << " " << it->first << " " << it->second << endl;
        }
    }
}
void myAlgorithm::calbound(Edge e0, int &k1, int &k2)
{
    int u = e0.u;
    int v = e0.v;
    vector<int> wList = getNeighbors(u, v);
    map<int, int, greater<int>> level; // level[i]表示e0的另外两边较小的truss值=i的个数
    for (int i = 0; i < wList.size(); i++)
    {
        int w = wList[i];
        int t = min(tau[w][u], tau[w][v]);
        if (level.find(t) == level.end())
        {
            level[t] = 1;
        }
        else
        {
            level[t]++;
        }
    }
    int cnt = 0, f1 = 0, f2 = 0;
    int k = level.begin()->first;
    while (k >= 2)
    {
        cnt += level[k];
        if (cnt >= k - 2)
        {
            k1 = k;
            break;
        }
        k--;
    }
    k = level.begin()->first;
    cnt = 0;
    while (k >= 2)
    {
        cnt += level[k];
        if (cnt >= k - 1)
        {
            k2 = k + 1;
            break;
        }
        k--;
    }
}
vector<int> myAlgorithm::getNeighbors(int u, int v)
{
    if (graph[u].size() > graph[v].size())
        swap(u, v);
    vector<int> t;
    for (unordered_map<int, int>::iterator iter = graph[u].begin(); iter != graph[u].end(); iter++)
    {
        int w = iter->first;
        if (graph[v].find(w) != graph[v].end())
        {
            t.push_back(w);
        }
    }
    return t;
}


size_t myAlgorithm::getSetSizeInBytes(const set<set<Edge>> &s)
{
    size_t totalSize = sizeof(s); // Base size of the set structure
    for (const auto &edgeSet : s)
    {
        totalSize += sizeof(edgeSet); // Size of each edgeSet
        for (const auto &edge : edgeSet)
        {
            totalSize += sizeof(edge); // Size of each Edge (assuming it is a POD type)
        }
    }
    return totalSize;
}
