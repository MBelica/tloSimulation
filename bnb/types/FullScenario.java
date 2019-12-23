package bnb.types;

import bnb.Configuration;

//TYPE:
//Full path in Scenario tree (index, probability, problem parameters per time stage)
public class FullScenario extends Configuration {

    public int[][] index; // 3 indices per time stage (t,g,l); i.e. one node per time stage
    public double P_s; // probability
    public double[][][] alpha_s; // coefficients in objective function for x: double[T][I][J]
    public double[][] beta_s; // coefficients in objective function for y: double[T][I]
    public double[][] beq_s; // right-hand side values for constraints: double[T][J]


    public FullScenario() {
        index = new int[timePeriods][3];
        alpha_s = new double[timePeriods][][];
        beta_s = new double[timePeriods][];
        beq_s = new double[timePeriods][];
    }

    public FullScenario(int t) {
        index = new int[t][3];
        alpha_s = new double[t][][];
        beta_s = new double[t][];
        beq_s = new double[t][];
    }
} //END TYPE: FullScenario