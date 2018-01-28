# sspp-algorithms
This repository contains code of several heuristic search algorithms for goal oriented Markov decision processes. The following algorithms are implemented: value iteration (vi), ilao*, hdp, lrtdp and fvi.Two classic MDP domains are given. The racetrack and the wetfloor.
All algorithms compute epsilon-consistent policies, with the exception of fvi where efficient bounds are available to compute epsilon-optimal policies.
The code successfully compiles with g++-6.The compile both domains simple do $make.
An example of how to call the algorithms on a given instance is given in usage.sh under racetrack and wetfloor folders.

