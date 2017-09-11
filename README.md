# TEAR
Time efficient energy aware routing in software defined networking

Description

Part1 : Tradeoff between energy saving and load balancing

TEAR.cpp is for the trade-off experiment which shows that adjusting the value of gamma and beta can control energy saving and load balancing.

The example data is tradeoff_big.zip.

Part2 : Performance experiment

TEAR_for_compare.cpp is for comparing ORPEAR and original shortest path to see four metrics performance : 
1. processing time  2. power saving   3. The number of link state change  4. The number of rule installation cost

The example data is big_topo2.zip


Execution

unzip the data ( tradeoff_big.zip , big_topo2.zip )

For tradeoff experiment :

>> g++ -std=c++11 -g -o TEAR TEAR.cpp

For performance experiment :

>> g++ -std=c++11 -g -o TEAR_for_compare TEAR_for_compare.cpp

Note:

Both ORPEAR and SP algorithm can be chosen in TEAR.cpp,TEAR_for_compare.cpp for executing.
The example data is just the part of evaluation in the paper, if you want the complete data
you can contact jimmyw86878@gmail.com.

Thank you!



