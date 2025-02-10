#pragma once
#include <condition_variable>
#include <atomic>
#include <chrono>
#include <memory>
#include <unordered_map>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <cstring>
#include <iostream>

using namespace std;
#define DEBUG_MODEL 1
#if DEBUG_MODEL == 1

#define LOG(fmt, args...)    \
    do                       \
    {                        \
        printf(fmt, ##args); \
    } while (0)
#else

#define LOG(...) // No operation if DEBUG_MODEL is not defined

#endif

typedef tuple<int, int, int> TRI;
typedef pair<int, int> PII;
typedef unordered_map<int, int> UMapIntInt;
typedef unsigned long long ULL;
struct Timer
{
    chrono::time_point<chrono::high_resolution_clock> start,end;
    string messages;
    Timer(string messages = "")
    {
        this->messages = messages;
        start = chrono::high_resolution_clock::now();
    }
    ~Timer()
    {
        end = chrono::high_resolution_clock::now();
        chrono::duration<double> duration = end - start;
        cout << "Time used: for "<< messages << " is " << duration.count() << " seconds" << endl;
    }
};
struct Edge
{
    int u, v;
    friend bool operator<(Edge e1, Edge e2)
    {
        return e1.u < e2.u || (e1.u == e2.u && e1.v < e2.v);
    }
    // Equality operator
    bool operator==(const Edge &other) const
    {
        return u == other.u && v == other.v;
    }
};
struct EdgeHasher
{
    size_t operator()(const Edge &e) const
    {
        return std::hash<int>{}(e.u) ^ (std::hash<int>{}(e.v) << 1);
    }
};
typedef unordered_map<Edge, int, EdgeHasher> UMapEdgeInt;
typedef unordered_map<Edge, bool, EdgeHasher> UMapEdgeBool;
typedef unordered_set<Edge, EdgeHasher> USetEdge;

inline void orderTwoVertex(int &u, int &v)
{
    if (u > v)
    {
        swap(u, v);
    }
}
inline Edge makeEdge(int u, int v)
{
    if (u > v)
        swap(u, v);
    return Edge{u, v};
}
// usedin global skip-vertex algorithm to find whether exists a parent string
// #define M_COMMUNITY_SIZE 1000 // 仅包含 '0' 和 '1'
struct Trie
{
    vector<vector<int>> nex;
    int cnt; // 节点计数
    Trie(int mDataSetEdgeNum) : cnt(0), nex(mDataSetEdgeNum, vector<int>(2, -1))
    {
    }
    int getSize()
    {
        int cnt1 = cnt == -1 ? nex.size() : cnt;
        // return cnt1*sizeof(int)*2;
        return cnt;
    }

    // 插入二进制字符串
    void insert(const string &s, int l)
    {
        // cout<<"insert: "<<s<<endl;
        // cout<<cnt<<" "<<l<<" "<<nex.size()<<endl;
        if (cnt + l >= nex.size())
        {
            cnt = -1; // 超出容量限制
            return;
        }

        int p = 0; // 当前节点
        for (int i = 0; i < l; i++)
        {
            int c = s[i] - '0'; // 将字符转换为 0 或 1
            if (nex[p][c] == -1)
                nex[p][c] = ++cnt; // 如果没有对应的子节点，则创建
            // cout<<"insert: "<<c<<" "<<p<<" "<<nex[p][0]<<" "<<nex[p][1]<<endl;
            p = nex[p][c]; // 移动到子节点
        }
        // cout<<"insert end: "<<cnt<<endl;
        // exist[p] = true; // 标记该节点为完整字符串
    }
    bool isFull()
    {
        return cnt == -1;
    }
    // 检查是否存在父串
    bool isParentString(const string &query, int l)
    {
        return isParentStringHelper(0, query, l, 0); // 从根节点开始递归检查
    }
    // check if exitst the same string
    bool isSameString(const string &query, int l)
    {
        // cout<<"query "<<query<<" "<<l<<endl;
        return isSameStringHelper(0, query, l, 0); // 从根节点开始递归检查
    }

private:
    bool isSameStringHelper(int p, const string &query, int l, int index)
    {
        // if(p!=-1)
        // cout<<"isSameStringHelper "<<p<<" "<<query[index]<<" "<<nex[p][0]<<" "<<nex[p][1]<<endl;
        if (p == -1)
            return false; // 如果已经到达根节点，说明不存在父串
        if (index == l)
            return true;            // 如果已经到达查询字符串的末尾，检查当前节点
        int c = query[index] - '0'; // 将字符转换为 0 或 1
        return (isSameStringHelper(nex[p][c], query, l, index + 1));
    }
    bool isParentStringHelper(int p, const string &query, int l, int index)
    {
        if (index == l)
            return true; // 如果已经到达查询字符串的末尾，检查当前节点
        if (p == -1)
            return false;           // 如果已经到达根节点，说明不存在父串
        int c = query[index] - '0'; // 将字符转换为 0 或 1

        // 如果当前字符是 '0'，我们可以选择忽略当前节点
        if (c == 0)
        {
            // 选择继续查找 '0' 或 '1'
            return isParentStringHelper(nex[p][0], query, l, index + 1) ||
                   isParentStringHelper(nex[p][1], query, l, index + 1);
        }
        else
        { // 当前字符是 '1'
            // 只能查找 '1' 的分支
            return (nex[p][1] != 0 && isParentStringHelper(nex[p][1], query, l, index + 1));
        }
    }
};

class DynamicBitset
{
private:
    std::vector<ULL> bits; // 存储位的数组
    ULL size;              // bitset的大小（位数）

    // 助手函数，计算所需的整数数量
    static ULL calculateCapacity(ULL bitsetSize)
    {
        return (bitsetSize + sizeof(ULL) * 8 - 1) / (sizeof(ULL) * 8);
    }

public:
    DynamicBitset(ULL numBits) : size(numBits)
    {
        bits.resize(calculateCapacity(numBits), 0);
    }
    DynamicBitset() {}

    // 直接访问位
    bool operator[](ULL index) const
    {
        if (index >= size)
            throw std::out_of_range("Index out of range");
        return (bits[index / (sizeof(ULL) * 8)] >> (index % (sizeof(ULL) * 8))) & 1;
    }

    // 设置位
    void set(ULL index, bool value)
    {
        if (index >= size)
            throw std::out_of_range("Index out of range");

        if (value)
        {
            // cout << "set " << index << " " << value << endl;
            // print();
            // cout<<index / (sizeof(ULL) * 8)<<" "<<(1 << (index % (sizeof(ULL) * 8)))<<endl;
            bits[index / (sizeof(ULL) * 8)] |= ((ULL)1 << (index % (sizeof(ULL) * 8)));
            // print();
        }
        else
        {
            // cout<<"set "<<index<<" "<<value<<endl;
            bits[index / (sizeof(ULL) * 8)] &= ~((ULL)1 << (index % (sizeof(ULL) * 8)));
            // print();
        }
    }

    // 获取大小
    ULL bitSize() const
    {
        return size;
    }
    int ullSize() const
    {
        return bits.size();
    }

    // 重置位
    void reset(ULL index)
    {
        if (index >= size)
            throw std::out_of_range("Index out of range");
        set(index, false);
    }

    // 清空所有位
    void clear()
    {
        std::fill(bits.begin(), bits.end(), 0);
    }

    // 打印位集（debug用途）
    void print() const
    {
        for (ULL i = 0; i < size; ++i)
        {
            std::cout << ((*this)[i] ? '1' : '0');
        }
        std::cout << std::endl;
    }
    // Function to find the position of the longest common prefix between two bitsets
    int longestCommonPrefix(const DynamicBitset &other) const
    {
        // cout << "longestCommonPrefix" << endl;
        // cout << bitSize() << " " << other.bitSize() << endl;
        // other.print();
        // Iterate through the bits of both bitsets
        for (ULL i = 0; i < this->ullSize() && i < other.ullSize(); ++i)
        {
            // Check if the bits at the current position are equal
            if ((this->bits[i] ^ other.bits[i]) == 0)
            {
                continue; // Continue if they are equal
            }

            // If they are not equal, check for the position of the first differing bit
            for (int j = 0; j < (sizeof(ULL) * 8); ++j)
            {
                if (j >= this->bitSize() || j >= other.bitSize())
                    return min(this->bitSize(), other.bitSize());
                if ((1 & ((this->bits[i]) >> j)) != (1 & ((other.bits[i]) >> j)))
                {
                    // cout << i * sizeof(ULL) * 8 + j << endl;
                    return i * sizeof(ULL) * 8 + j; // Return the position
                }
            }
        }
        return min(this->bitSize(), other.bitSize()); // Return result if no differing bits are found
    }
    // =
    DynamicBitset &operator=(const DynamicBitset &other)
    {
        if (this != &other)
        { // 自赋值检查
            size_t otherSize = other.bitSize();
            this->resize(otherSize);
            try
            {
                // 拷贝位数据
                for (size_t i = 0; i < otherSize; ++i)
                {
                    this->set(i, other[i]);
                }
            }
            catch (...)
            {
                // 发生异常时,恢复对象到原始状态
                this->resize(this->size);
                throw;
            }
        }
        return *this;
    }

    // subbits
    DynamicBitset subBits(ULL start, ULL end) const
    {
        DynamicBitset result(end - start);
        end = min(end, size);
        if (start < end)
        {
            for (ULL i = start; i < end; i++)
            {
                // cout<<"subBits "<<i-start<<" "<<(*this)[i]<<endl;
                result.set(i - start, (*this)[i]);
            }
        }
        // cout << "subBits " << start << " " << end << endl;
        // result.print();
        // print();
        return result;
    }

    // resize 方法
    void resize(ULL newSize)
    {
        bits.resize(calculateCapacity(newSize));
        size = newSize;
    }
};

class RadixTreeNode
{
public:
    RadixTreeNode(const DynamicBitset &bits)
        : bits(bits) {}
    RadixTreeNode() {}

    unordered_set<shared_ptr<RadixTreeNode>> children;
    DynamicBitset bits;
};

class RadixTree
{
public:
    RadixTree() : root(make_shared<RadixTreeNode>()), size(sizeof(RadixTreeNode)) {};

    virtual ~RadixTree() = default;

    RadixTree(const RadixTree &) = delete;
    RadixTree &operator=(const RadixTree &) = delete;

    void insert(DynamicBitset &bits)
    {
        // Timer timer("tree insert");
        insert_helper(bits, root);
        // LOG("insert_helper memory size: %d\n", size);
        // print();
    }

    bool search(const DynamicBitset &bits)
    {
        // Timer timer("tree search");

        bool flag = search_helper(bits, root);
        // LOG("search result: %d\n", flag);
        // print();
        // if(flag)
        // {
        //     cout<<"search result: true"<<endl;
        //      bits.print();
        // }

        return flag;
    }
    int getSize()
    {
        return size;
    }
    void print()
    {
        print_helper(root, 0);
    }

private:
    shared_ptr<RadixTreeNode> root;
    ULL size;
    void insert_helper(DynamicBitset &bits, std::shared_ptr<RadixTreeNode> root)
    {
        std::shared_ptr<RadixTreeNode> current = root;

        while (true)
        {
            if (!current || bits.bitSize() == 0)
            {
                return; // 边界条件检查
            }

            if (current->children.empty())
            {
                try
                {
                    auto newNode = std::make_shared<RadixTreeNode>(bits);
                    current->children.insert(newNode);
                    size += sizeof(RadixTreeNode);
                }
                catch (const std::bad_alloc &)
                {
                    // 内存分配失败的错误处理
                    LOG("Memory allocation failed for new node.");
                    return;
                }
                return;
            }

            std::shared_ptr<RadixTreeNode> &matchingNode = current;
            int maxPrefixLen = 0;

            for (const auto &child : current->children)
            {
                int prefixLen = child->bits.longestCommonPrefix(bits);
                if (prefixLen > 0)
                {
                    maxPrefixLen = prefixLen;
                    matchingNode = child;
                    break;
                }
            }
            // LOG("maxPrefixLen: %d\n", maxPrefixLen);
            if (maxPrefixLen > 0)
            {
                if (maxPrefixLen < matchingNode->bits.bitSize())
                {
                    try
                    {
                        auto newNode = std::make_shared<RadixTreeNode>(
                            matchingNode->bits.subBits(maxPrefixLen, matchingNode->bits.bitSize()));
                        size += sizeof(RadixTreeNode);
                        matchingNode->bits = matchingNode->bits.subBits(0, maxPrefixLen);
                        // matchingNode->bits.subBits(0, maxPrefixLen).print();
                        // matchingNode->bits.print();

                        matchingNode->children.swap(newNode->children);
                        matchingNode->children.insert(newNode);

                        if (maxPrefixLen < bits.bitSize())
                        {
                            auto newNode2 = std::make_shared<RadixTreeNode>(
                                bits.subBits(maxPrefixLen, bits.bitSize()));
                            size += sizeof(RadixTreeNode);
                            matchingNode->children.insert(newNode2);
                        }
                        //     cout << "1" << endl;
                        //     print();
                    }
                    catch (const std::bad_alloc &)
                    {
                        // 内存分配失败的错误处理
                        LOG("Memory allocation failed during node splitting.");
                        return;
                    }
                }
                else if (maxPrefixLen < bits.bitSize())
                {
                    bits = bits.subBits(maxPrefixLen, bits.bitSize());
                    current = matchingNode;
                    continue;
                }
                return;
            }

            try
            {
                auto newNode = std::make_shared<RadixTreeNode>(bits);
                current->children.insert(newNode);
                size += sizeof(RadixTreeNode);
            }
            catch (const std::bad_alloc &)
            {
                // 内存分配失败的错误处理
                LOG("Memory allocation failed for new node insertion.");
                return;
            }
            return;
        }
    }

    void print_helper(shared_ptr<RadixTreeNode> node, int depth)
    {
        cout << "depth " << depth << endl;
        node->bits.print();
        // Recursively print the children with updated prefix
        for (const auto &child : node->children)
        {
            print_helper(child, depth + 1);
        }
    }
    bool search_helper(const DynamicBitset &bits, shared_ptr<RadixTreeNode> node)
    {
        for (auto current : node->children)
        {
            int i = current->bits.longestCommonPrefix(bits);
            // cout << "longestCommonPrefix" << endl;
            // current->bits.print();
            // bits.print();
            // cout<<"i "<<i<<endl;
            if (i != 0)
            {
                // 情况一：当前节点的内容与字符串完全匹配，根据是否为完整单词判断结果
                if (i == bits.bitSize() && i == current->bits.bitSize())
                {
                    return true;
                }
                // 情况二：当前节点的内容是字符串的前缀
                else if (i == current->bits.bitSize())
                {
                    return search_helper(bits.subBits(i, bits.bitSize()), current);
                }
                // 情况三：字符串的内容是当前节点的前缀，直接返回错误
                else
                {
                    return false;
                }
            }
        }
        // 没有找到
        return false;
    }
};
