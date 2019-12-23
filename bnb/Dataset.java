package bnb;

/*DATA:
 * Class structures and input BnBUtil are defined
 * BnBUtil: Scenario Tree (in Node notation AND Scenario notation) with respective problem-specific parameters
 * classes: Node in Scenario, Path in Scenario, Solution for Scenario, Solution for Original Problem, BnB node, Fathoming result
 * connection between both notations (between scenario and node notation, i.e. original and splitting variable representation)
 */

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.*;

import bnb.types.FullScenario;
import bnb.types.NodeScenario;

import static bnb.Dataset.Core.generateData;

/**
 * DataSet Object for the Scenario Trees
 */
public class Dataset extends Configuration {

    //Scenario Tree
    public final int T; // no. time periods
    public final int[] G; // no. groups with same predecessor per time period t in T
    public final int[][] L; // no. leaves per time period t in T within group g in G

    //Problem specific
    public final int I; // no. facilities (dimension for decision variables)
    public final int J; // no. customers (dimension for decision variables)

    public NodeScenario[][][] DataTree;
    public FullScenario[][] ScenarioTree;

    public Dataset() {
        T = timePeriods;
        G = groupsPerPeriod;
        L = leavesPerPeriod;

        I = facilities;
        J = customers;

        initializeTree(T, G, L);
        loadSetup();
    }

    private void loadSetup() {

        DataTree = generateData(DataTree, T, G, L, I, J);

        createScenarioTree();
    }

    /**
     * initializeTree
     * Scenario tree consisting of multiple nodes
     *
     * @param T time periods to be initialized
     * @param G groups with same predecessor per time period
     * @param L leaves per time period within group
     */
    private void initializeTree(int T, int[] G, int[][] L) {
        //Define return variable (Scenario Tree)
        DataTree = new NodeScenario[T][][];
        for (int i = 0; i < T; i++) {
            DataTree[i] = new NodeScenario[G[i]][];
            for (int g = 0; g < G[i]; g++) {
                DataTree[i][g] = new NodeScenario[L[i][g]];
                for (int l = 0; l < L[i][g]; l++) {
                    DataTree[i][g][l] = new NodeScenario();
                } //for l
            } //for g
        } //for t

        //Initialize Node indices: t,g,l
        for (int i = 0; i < T; i++) {
            for (int j = 0; j < G[i]; j++) {
                for (int k = 0; k < L[i][j]; k++) {
                    DataTree[i][j][k].t = i;
                    DataTree[i][j][k].g = j;
                    DataTree[i][j][k].l = k;
                } //for l
            } //for g
        } //for t
    } //END DATA: DataTree

    /**
     * createScenarioTree
     * Scenario tree consisting of multiple scenarios (computed from nodes notation)
     */
    private void createScenarioTree() {

        //Define return variable (Scenario Tree)
        FullScenario[][] S = new FullScenario[G[T - 1]][];
        for (int g = 0; g < G[T - 1]; g++) {
            S[g] = new FullScenario[L[T - 1][g]];
            for (int l = 0; l < L[T - 1][g]; l++) {
                S[g][l] = new FullScenario(T);
            } //for l
        } //for g

        //Compute problem parameters from nodes notation
        for (int g = 0; g < G[T - 1]; g++) {
            for (int l = 0; l < L[T - 1][g]; l++) {
                int tp = T - 1;
                int gp = g;
                int lp = l;
                //Compute parameters for all predecessors
                while (tp >= 0) {
                    //Index
                    S[g][l].index[tp][0] = tp;
                    S[g][l].index[tp][1] = gp;
                    S[g][l].index[tp][2] = lp;
                    //Parameters
                    S[g][l].alpha_s[tp] = DataTree[tp][gp][lp].alpha;
                    S[g][l].beta_s[tp] = DataTree[tp][gp][lp].beta;
                    S[g][l].beq_s[tp] = DataTree[tp][gp][lp].beq;
                    //Compute Predecessor
                    int[] predecessor = Predecessor(tp, gp, L);
                    tp = predecessor[0];
                    gp = predecessor[1];
                    lp = predecessor[2];
                } //while(tp>=0)
                //Probability (same for all time periods)
                S[g][l].P_s = DataTree[T - 1][g][l].P;
            } //for l
        } //for g
        ScenarioTree = S;
    } //END DATA: Scenario tree

    /**
     * Predecessor
     * Compute indices of predecessor of specific node in scenario tree
     *
     * @param t time period
     * @param g group with same predecessor in time period t
     * @param L DataTree to be examined
     */
    public static int[] Predecessor(int t, int g, int[][] L) {
        int[] tgl = new int[3];
        int h = 0;
        while (g + 1 - L[t][h] > 0) {
            g = g - L[t][h];
            h = h + 1;
        }
        tgl[0] = t - 1;
        tgl[1] = h;
        tgl[2] = g;
        return tgl;
    } //END METHOD: Predecessor

    /**
     * numberScenarioWithNode
     * matrix that counts no. scenarios containing a specific nodes (for all nodes)
     *
     * @return int[][][] count of scenarios containing a specific nodes
     */
    public int[][][] numberScenarioWithNode() {

        //Define return value (binary)
        int[][][] n = new int[T][][];
        for (int t = 0; t < T; t++) {
            n[t] = new int[G[t]][];
            for (int g = 0; g < G[t]; g++) {
                n[t][g] = new int[L[t][g]];
            } //for g
        } //for t

        for (int g = 0; g < G[T - 1]; g++) {
            for (int l = 0; l < L[T - 1][g]; l++) {
                int tp = T - 1;
                int gp = g;
                int lp = l;
                while (tp >= 0) {
                    //count number of scenarios containing tp,gp,lp
                    n[tp][gp][lp] += 1;
                    //compute predecessor
                    int[] predecessor = Predecessor(tp, gp, L);
                    tp = predecessor[0];
                    gp = predecessor[1];
                    lp = predecessor[2];
                } //while(tp>=0)
            } //for l
        } //for g

        return n;
    } //END METHOD: numberScenarioWithNode

    /**
     * scenarioContainingNode
     * list of all scenarios that contain one specific node
     *
     * @param n_t time period to be examined
     * @param n_g groups with same predecessor to be examined
     * @param n_l leaves per time period within group to be examined
     * @param n tree to be examined
     * @return int[][]
     */
    public int[][] scenarioContainingNode(int n_t, int n_g, int n_l, int[][][] n) {
        //define return value
        int nn = n[n_t][n_g][n_l];
        int[][] s_vt = new int[nn][2]; //s_vt[0]: index g[T-1], s_vt[1]: index l[T-1][g]
        //fill list of scenarios
        int i = 0;
        for (int g = 0; g < G[T - 1]; g++) {
            for (int l = 0; l < L[T - 1][g]; l++) {
                int gp = g;
                int tp = T - 1;
                int lp = l;
                while (n_t != tp) {
                    int[] pre = Predecessor(tp, gp, L);
                    tp = pre[0];
                    gp = pre[1];
                    lp = pre[2];
                } //while(n_t != tp)
                //if indices match
                if (n_g == gp && n_l == lp) {
                    s_vt[i][0] = g;
                    s_vt[i][1] = l;
                    i = i + 1;
                } //end if
            } //for l
        } //for g

        return s_vt;
    } //END METHOD: scenarioContainingNode


    /**
     * DataSet Object for the Scenario Trees
     */
    public static class Core {

        /**
         * generateData
         * Calls function to generate data and fills the NodeScenario within DataTree,
         * which has to have been initialised first
         *
         * @param DataTree scenario to be filled with random data
         * @param T time periods in game
         * @param G groups with same predecessor per time period
         * @param L leaves per time period per group
         * @param I no. facilities
         * @param J no. customers
         * @return NodeScenario[][][]
         */
        static NodeScenario[][][] generateData(NodeScenario[][][] DataTree, int T, int[] G, int[][] L, int I, int J) {

            for (int i = 0; i < T; i++) {

                int[] noSubNode = new int[T];
                double[] nodeProbDistribution;

                for (int j = 0; j < G[i]; j++) noSubNode[i] += L[i][j];
                nodeProbDistribution = randomUProbDistribution(probDecimals, noSubNode[i]);


                for (int k = 0; k < G[i]; k++) {
                    for (int l = 0; l < L[i][k]; l++) {
                        DataTree[i][k][l].P   = nodeProbDistribution[l];

                        DataTree[i][k][l].alpha = new double[I][J];
                        if (k == 0) {
                            DataTree[i][k][l].beta = prepareArrayOfRandoms(I, betaParam1, betaParam2, betaDecimals, betaMethod);
                            DataTree[i][k][l].beq  = prepareArrayOfRandoms(J, beqParam1,  beqParam2, beqDecimals, beqMethod);
                            for (int m = 0; m < I; m++) DataTree[i][k][l].alpha[m] = prepareArrayOfRandoms(J, alphaParam1, alphaParam2, alphaDecimals, alphaMethod);
                        } else {
                            DataTree[i][k][l].beta  = DataTree[i][0][l].beta;
                            DataTree[i][k][l].beq   = DataTree[i][0][l].beq;
                            DataTree[i][k][l].alpha = DataTree[i][0][l].alpha;
                        } //if k==0
                    } //for l
                } //for k
            } //for i
            return DataTree;
        }

        /**
         * createArrayOfRandoms
         * prepares generating random array by deciding whether
         * last digit before decimals is supposed to be 0 or random
         *
         * @param size length of array to be created
         * @param param1 for uniform: minimal value each entry can have, for gaussian: mean
         * @param param2 for uniform: maximal value each entry can have, for gaussian: standard deviation
         * @param method which probability distribution is used for generating numbers
         * @param decimals number of decimals
         * @return double[]
         */
        private static double[] prepareArrayOfRandoms (int size, double param1, double param2, int decimals, String method) {
            double[] result;

            if ( (method.equals("uniform")) && (param2 <= param1 || decimals < 0) ) {
                throw new IllegalArgumentException("Parameters for uniform distribution illegal");
            }
            else if ( (method.equals("gaussian")) && (param2 <= 0 || decimals < 0) ) {
                throw new IllegalArgumentException("Parameters for gaussian distribution illegal");
            }
            else result = createArrayOfRandoms(size, param1, param2, decimals, method);

            return result;
        }

        /**
         * randomArrayOfDecimals
         * generates an array of random doubles with given size and given decimals
         * used "Random.nextInt() instead of Math.random() as it's both more efficient and less biased
         *
         * @param size length of array to be created
         * @param param1 for uniform: minimal value each entry can have, for gaussian: mean
         * @param param2 for uniform: maximal value each entry can have, for gaussian: standard deviation
         * @param decimals no. of decimals for each entry
         * @return double[]
         */
        private static double[] createArrayOfRandoms(int size, double param1, double param2, int decimals, String method) {
            Random r = new Random();
            double[] result = new double[size];

            for (int i = 0; i < size; i++) {
               switch (method) {
                   case "uniform": // param1 = min, param2 = max
                       result[i] = (r.nextInt( (int)((param2 - param1) * ((int) Math.pow(10, decimals))) )  / Math.pow(10, decimals) ) + param1;
                       break;
                   case "gaussian": // param1 = mean, param2 = standard deviation
                       result[i] = round(((r.nextGaussian() * param2) + param1), decimals);
                       break;
                   default: // assuming uniform: param1 = min, param2 = max
                       result[i] = (r.nextInt( (int)((param2 - param1) * ((int) Math.pow(10, decimals))) )  / Math.pow(10, decimals) ) + param1;
               }
            } //for i
            return result;
        }

        /**
         * round
         * the standard rounding method breaks down badly in corner cases with either
         * a very high number of decimal places or large integer part. Hence, this updated version is used
         *
         * @param number to be rounded
         * @param decimals places used during rounding
         * @return double
         */
        private static double round(double number, int decimals) {
            if (decimals < 0) throw new IllegalArgumentException();

            BigDecimal bd = new BigDecimal(number);
            bd = bd.setScale(decimals, RoundingMode.HALF_UP);
            return bd.doubleValue();
        }

        /**
         * randomUProbDistribution
         * generates an array of random doubles summing up to one
         * with a given number of decimals and using a uniform distribution
         *
         * @param decimals no. of decimals for each entry
         * @param numberOfDraws no. of probabilities to be generated
         * @return double[]
         */
        private static double[] randomUProbDistribution(int decimals, int numberOfDraws) {
            int targetSum = (int) Math.pow(10, decimals);
            Random r = new Random();
            List < Integer > load = new ArrayList < > ();

            //random numbers
            int sum = 0;
            for (int i = 0; i < numberOfDraws; i++) {
                int next = r.nextInt(targetSum) + 1;
                load.add(next);
                sum += next;
            } //for i

            //scale to the desired target sum
            double scale = 1d * targetSum / sum;
            sum = 0;
            for (int i = 0; i < numberOfDraws; i++) {
                load.set(i, (int)(load.get(i) * scale));
                sum += load.get(i);
            }

            //take rounding issues into account
            while (sum++ < targetSum) {
                int i = r.nextInt(numberOfDraws);
                load.set(i, load.get(i) + 1);
            } //while (sum++ < targetSum)

            double[] result = new double[numberOfDraws];
            for (int l = 0; l < numberOfDraws; l++) result[l] = load.get(l) / Math.pow(10, decimals);

            return result;
        }

        /**
         * saveDataAsFile
         * generates a file saving amongst other things the
         * NodeScenario[][][] DataTree of the given instance in a file
         *
         * @param data DataSet containing Scenario to save
         */
        public static void saveDataAsFile (Dataset data, String filePath) {

            NodeScenario[][][] DataTree = data.DataTree;

            StringBuilder probDist = new StringBuilder();
            for (int i = 0; i < data.T; i++) {
                probDist.append("(");
                for (int j = 0; j < data.G[i]; j++) {
                    for (int k = 0; k < data.L[i][j]; k++) {
                        probDist.append(data.DataTree[i][j][k].P);
                        if (k < (data.L[i][j]-1)) probDist.append(",");
                    }
                }
                probDist.append(")");
                if (i < ( data.T-1)) probDist.append(",");
            }

            StringBuilder leaves = new StringBuilder();
            for (int i = 0; i < data.L.length; i++) {
                leaves.append(Arrays.toString(data.L[i]));
                if (i < (data.L.length - 1)) leaves.append(",");
            }

            Date dt = new Date();
            StringBuilder builder = new StringBuilder();
            builder.append("DataTree;;").append(dt).append(" \n \n")
                    .append(";Facilities:;").append(data.I).append(";;Probability Distribution:;;(").append(probDist).append(") \n")
                    .append(";Customers:;").append(data.J).append(";;Groups with same Predecessor.:;;").append(Arrays.toString(data.G)).append(" \n")
                    .append(";Time Periods:;").append(data.T).append(";;Leaves per Time Period:;;[").append(leaves).append("] \n \n \n");

            for (int i = 0; i < DataTree.length; i++) {
                int j = 0; // the entries are for all j the same
                builder.append("Period ").append(i + 1).append("\n ;fixed costs f_i^1  \n");
                for (int k = 0; k < DataTree[i][j].length; k++) {
                    builder.append(";;Scenario ").append(k + 1).append(";");
                    for (int l = 0; l < DataTree[i][j][k].beta.length; l++) {

                        builder.append(" ").append(DataTree[i][j][k].beta[l]);
                        if (l < DataTree[i][j][k].beta.length - 1) builder.append(";");
                    }
                    builder.append("\n");
                }
                builder.append("\n \n ;c_{ij}^1  \n");
                for (int k = 0; k < DataTree[i][j].length; k++) {
                    builder.append(";;Scenario ").append(k + 1).append(";");
                    for (int l = 0; l < DataTree[i][j][k].alpha.length; l++) {
                        for (int m = 0; m < DataTree[i][j][k].alpha[l].length; m++) {

                            builder.append(" ").append(DataTree[i][j][k].alpha[l][m]);
                            if (l < DataTree[i][j][k].alpha[l].length - 1) builder.append(";");
                        }
                        builder.append("\n ;;;");
                    }
                    builder.append("\n");
                }

                builder.append("\n ;d_j^1  \n");
                for (int k = 0; k < DataTree[i][j].length; k++) {
                    builder.append(";;Scenario ").append(k + 1).append(";");
                    for (int l = 0; l < DataTree[i][j][k].beq.length; l++) {

                        builder.append(" ").append(DataTree[i][j][k].beq[l]);
                        if (l < DataTree[i][j][k].beq.length - 1) builder.append(";");
                    }
                    builder.append("\n");
                }
                builder.append("\n \n \n");

            }
            try {
                File file = new File(filePath);

                if (!file.exists()) //noinspection ResultOfMethodCallIgnored
                    file.createNewFile();

                try (BufferedWriter writer = new BufferedWriter(new FileWriter(filePath))) {
                    writer.write(builder.toString()); //save the string representation of the board
                    writer.close();
                }
            } catch (IOException e) {
                throw new IllegalArgumentException(e);
            }
        }
    }
}