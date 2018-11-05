# Analog maxSAT solver

Many real-life optimization problems can be formulated in Boolean logic as MaxSAT, a class of problems where the task is finding Boolean assignments to variables satisfying the maximum number of logical constraints. Since MaxSAT is NP-hard, no algorithm is known to
efficiently solve these problems. Here we present a continuous-time analog solver for MaxSAT and show that the scaling of the escape rate, an invariant of the solver’s dynamics, can predict the maximum number of satisfiable constraints, often well before finding the
optimal assignment. Simulating the solver, we illustrate its performance on MaxSAT competition problems, then apply it to two-color Ramsey number R(m, m) problems. Although it finds colorings without monochromatic 5-cliques of complete graphs on N ≤ 42 vertices, the
best coloring for N = 43 has two monochromatic 5-cliques, supporting the conjecture that R(5, 5) = 43. This approach shows the potential of continuous-time analog dynamical systems as algorithms for discrete optimization.

For more details on how the system works please study:
B Molnár, F Molnár, M Varga, Z Toroczkai & M Ercsey-Ravasz
**A continuous-time Max-SAT solver with high analog performance**
*Nature Communications, DOI: 10.1038/s41467-018-07327-2, 2018*

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### What's in the box?

* makefile
* cfg folder - containing the configuration files
* cnf folder - in this folder you have to store the SAT problems in DIMACS cnf format
* includes - containing the header files of the project
* sh - contains the bash script for generating the config files based on the data in the cnf folder
* src - contains the source code

### Installing

Navigate into the project folder and type the command
```
make
```
This will build the project and create the executable file in the folder of the project named
```
analogsat
```

## Running analogsat

1. Copy your custom SAT problem files into the cnf folder. Please make sure that they are in DIMACS cnf format (see pre uploaded files for example).
2. Run generateCFG.sh from the sh folder using the relative path from the cnf root folder the folder containing the cnf file.
Eg.
```
./generateCFG.sh 2016/ms_random/highgirth/4sat
```
This will generate the cfg file into the cfg folder.
3. Run analogsat
```
./analogsat <CFG filename> <tmax> <Nprobb>
```
   where
   * CFG filename - config file filename without extension from cfg folder CFG file must include the names of the problems in DIMACS cnf format
   * tmax - maximum analog time to simulate (0 < tmax <= 150). We recommend 35 or 50 depending on the size and complextiy of the problem.
   * Nprobb - initial number of trajectories (1000 <= Nprobb <= 2000000)
   Eg.
```
./analogsat 2016.ms_random.highgirth.3sat 50 100000
```
4. The results will be generated in the results, solutions and logs folders.
   * results - final decision about the the minimum number of unsatisfied clauses, the predictions for the minimum and the trajectories needed to find that minimum, running time, etc.
   * solutions - solutions to the SAT problem satisfying each reached level of satisfied clauses
   * logs - real-time performance of the system, exporting all the necessary data regarding the fitting, aproximations, running times, etc.

## Built With

* C++
* [Numerical Recipies](http://numerical.recipes/)

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

* **Botond MOLNÁR** - *Initial work* - [molnarb14](https://github.com/molnarb14)
* **Ferenc MOLNÁR**
* **Melinda VARGA**
* **Zoltán TOROCZKAI**
* **Mária-Magdolna ERCSEY-RAVASZ**
See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

We thank to X.S. Hu, S. Datta, Z. Néda, E. Regan, X. Yin, and S. Kharel for useful comments and discussions. We also thank Brendan McKay for his corrections and comments on known Ramsey colorings for R(5, 5). This work was supported in part by grants of the Romanian National Authority for Scientific Research and Innovation CNCS/CCCDI-UEFISCDI, project nrs. PN-III-P2-2.1-BG-2016-0252 (MER), PN-III-P1-
1.1-TE-2016-1457 (MER, BM), COFUND-FLAGERA II-CORTICITY (MER), the GSCE-30260-2015 Grant for Supporting Excellent Research of the Babeş-Bolyai University (BM, MER). It was also supported in part by the National Science Foundation under Grants CCF-1644368 and 1640081, and by the Nanoelectronics Research Corporation, a wholly-owned subsidiary of the Semiconductor Research Corporation, through Extremely
Energy Efficient Collective Electronics, an SRC-NRI Nanoelectronics Research Initiative under Research Task ID 2698.004 (Z.T., M.V., and F.M.).

