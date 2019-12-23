package bnb.types;

import bnb.Configuration;

//TYPE:
//Node in Branch and Bound tree (solution, branched y for this node, all branched y until this node)
public class BnBNode extends Configuration {
    public int[] current_y = new int[4]; // tgli; index for y that was branched for this node
    public int[][][][] branched_y; // tgli; 0 = not branched, 1 = branched to zero, 2 = branched to one
    public SolutionScenario[][] sols; // solution for scenarios: [G[T-1]][L[T-1][g]]


    public BnBNode() {
        sols = new SolutionScenario[groupsPerPeriod[timePeriods - 1]][];
    }

    public BnBNode(int T, int[] G) {
        sols = new SolutionScenario[G[T - 1]][];
    }
}