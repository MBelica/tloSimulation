package bnb.types;

//TYPE:
//Node in Scenario tree (index, probability, problem parameters)
public class NodeScenario {
    public int t; //time stage index
    public int g; //group index
    public int l; //leaf index
    public double P; //probability
    public double[][] alpha; //coefficients in objective function for x: double[I][J]
    public double[] beta; //coefficients in objective function for y: double[I];
    public double[] beq; //right-hand side values for constraints: double[J];
} //END TYPE: NodeScenario