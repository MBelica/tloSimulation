package bnb;

public class Configuration {

    //General
    public static String dataFilePath       = "dataTree.csv";           // path and filename used for saveDataAsFile

    //Scenario Tree
    public static int timePeriods           = 3;					    // no. time periods
    public static int[] groupsPerPeriod     = {1,2,4};				    // no. groups with same predecessor per time period t in T
    public static int[][] leavesPerPeriod   = { {2}, {2,2},{2,2,2,2} };	// no. leaves per time period t in T within group g in G

    //Problem specific
    public static int facilities            = 4;		                // no. facilities (dimension for decision variables)
    public static int customers             = 5;			            // no. customers (dimension for decision variables)

    //Program specific
    static int probDecimals                 = 5;                        // no. of decimals for probability distribution

    /**
     * For beta, beq and alpha the following parameters are used
     *   - method: which probability distribution is used for generating numbers, currently implemented:
     *      + "uniform"
     *      + "gaussian"
     *   - param1: depends on method used
     *      + for uniform distribution: minimal value each entry can have,
     *      + for gaussian distribution: mean
     *   - param2:
     *      + for uniform distribution: maximal value each entry can have,
     *      + for gaussian distribution: standard deviation
     *   - decimals: number of decimals
     */
    static String betaMethod                = "gaussian";
    static double betaParam1                = 75.0;
    static double betaParam2                = 10.0;
    static int betaDecimals                 = 2;

    static String beqMethod                 = "gaussian";
    static double beqParam1                 = 150.0;
    static double beqParam2                 = 25.0;
    static int beqDecimals                  = 2;

    static String alphaMethod               = "uniform";
    static double alphaParam1               = 0.0;
    static double alphaParam2               = 3.0;
    static int alphaDecimals                = 1;
}
