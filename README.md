# WhyCom (TWCS)
## Introduction
This repository contains the code used in our paper: **Truss-based Why-not Community Search**. In this paper, we investigate a new problem of truss-based why-not community search. Given a $k$-truss community $\mathcal{C}$ in a graph $G$ and a why-not vertex $w \notin \mathcal{C}$, the goal is to insert the minimum number of new edges into $G$ to ensure that $w$ becomes part of the $k$-truss community.
## Datasets
We use ten real-world graphs, which are from the Stanford Large Network Dataset Collection\footnote{http://snap.stanford.edu/data/}.
## Algorithms
We implement two algorithms: SIM* (GlobalPlus), EXA* (LocalPlus). All algorithms are implemented in C++. All experiments are conducted on Linux Server with 3.4GHz eight-core CPU and
131GB main memory.
- decomp.h, decomp.cpp: truss decomposition algorithm
- algorithm.h, algorithm.cpp: EXA* and SIM* algorithm
- main.cpp: main function to run the algorithm
- graph.h, graph.cpp: graph data structure
## Usage
### compile and run
`make`
`./main.out dataSetFile k communityFile whyNotVertex algorithmName`
### parameters
1. input
- dataSet file path
- truss value k
- community file path
- why-not vertex
- algorithm name
2. output
- number of newly inserted edges
- running time
- verification result
### example
Run the following command to run the algorithm on Graph.example with k=5, community file Targetcommunity.example, why-not vertex 8, and algorithm GlobalPlus:
`make`
`./main.out Graph.example 5 Targetcommunity.example 8 GlobalPlus`
Output is as follows:
2 0.00337109
----------------------------------Verification----------------------------------

number of inserted edges: 2

they are 

 ( 4 , 8 ) 
 ( 5 , 8 ) 


correct! w is in k_truss_community

