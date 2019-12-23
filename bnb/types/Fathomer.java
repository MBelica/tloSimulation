package bnb.types;

//TYPE:
//result of fathoming method (boolean fathomed, upper bound, best solution)
public class Fathomer {
    public boolean fathom; // fathom this node?
    public double F_U; // current best objective value
    public SolutionOriginal sol = new SolutionOriginal(); // current best solution
} //END TYPE: Fathomer