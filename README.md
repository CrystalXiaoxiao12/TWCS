## input
- dataSet file path
- truss value k
- community file path
- why-not vertex
- algorithm name
## output
- number of newly inserted edges
- running time
- verification result
## example
### input

run  `make`

`./main.out Graph.example 5 Targetcommunity.example 8 GlobalPlus`
### output
2 0.00337109


----------------------------------Verification----------------------------------

number of inserted edges: 2

they are 
 ( 4 , 8 ) 
 ( 5 , 8 ) 
 
correct! w is in k_truss_community
