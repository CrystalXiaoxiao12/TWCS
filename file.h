#pragma once

#include <iostream>
#include "myGraph.h"
#include "common.h"
#include<fstream>
#include<unordered_set>

using namespace std;
void file_ReadGraph(myGraph &G, string dataset_file_name);
void file_ReadKTrussCommunity(USetEdge &k_truss_community, string k_truss_community_file_name);