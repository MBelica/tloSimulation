package bnb.types;

//TYPE:
//Solution in one scenario for all time stages
public class SolutionScenario {
    public double result; // objective value: double
    public double[][][] xval; // solution x: double[T][I][J]
    public double[][] yval; // solution y: double[T][I]
    public boolean solvable;
} //END TYPE: SolutionScenario