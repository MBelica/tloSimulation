package bnb.types;

//TYPE:
//Solution for original problem (for all nodes in Scenario Tree)
public class SolutionOriginal {
    public double result; // objective value: double
    public double[][][][][] xval; // solution x: double[T][G][L][I][J]
    public double[][][][] yval; // solution y: double[T][G][L][I]
} //END TYPE: SolutionOriginal