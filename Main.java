import java.util.*;

import bnb.Configuration;
import bnb.Dataset;
import bnb.types.*;

import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

import static bnb.Dataset.Core.saveDataAsFile;


public class Main {

	static Dataset data = new Dataset();

	//import Data
	static int I     = data.I;
	static int J     = data.J;

	static int T     = data.T;
	static int[] G   = data.G;
	static int[][] L = data.L;

	// MAIN: BnB framework

	public static void main(String[] args) {
		final long timeStart = System.currentTimeMillis();
		//compute matrix scenarios with specific node
		int[][][] n = data.numberScenarioWithNode( );
		//define best solution
		SolutionOriginal bestSol = new SolutionOriginal();
		
		//Initialize origin node
		BnBNode OriginNode = new BnBNode(T, G);
		OriginNode.branched_y = new int[T][][][];
		for(int t=0; t<T; t++){
			OriginNode.branched_y[t] = new int[G[t]][][];
			for(int g=0; g<G[t]; g++){
				OriginNode.branched_y[t][g] = new int[L[t][g]][I];
				for(int l=0; l<L[t][g]; l++){
					for(int i=0; i<I; i++){
						OriginNode.branched_y[t][g][l][i] = 0;
					}//for i
				}//for l
			}//for g
		}//for t
		
		//Solve Q_s^{SV-NA} for all scenarios s
		SolutionScenario[][] resultScenario = new SolutionScenario[G[T-1]][];
		for(int g=0; g<G[T-1]; g++){
			resultScenario[g] = new SolutionScenario[L[T-1][g]];
			for(int l=0; l<L[T-1][g]; l++){
				resultScenario[g][l] = splittingVariable(data.ScenarioTree[g][l], OriginNode.branched_y);
			}//for l
		}//for g
						
		//Compute averages and round them to obtain solution for original problem
		SolutionOriginal result_orig    = avgSolution(resultScenario); 			//no objective value yet
		SolutionOriginal input_rounded  = roundY(result_orig); 					//no objective value yet
		SolutionOriginal result_rounded = originalModel(input_rounded.yval);
		bestSol = result_rounded;
		
		//Initialize upper bound F^U
		double F_U = result_rounded.result;
		OriginNode.sols = resultScenario;
				
		//Branch node k=0
		BnBNode[] ChildNodes = new BnBNode[2];
		ChildNodes = branch(OriginNode);

		
		
		//Start BnB Node List with two children of Origin Node
		Collection <BnBNode> BnBNodes = new HashSet<BnBNode>();
		BnBNodes.add(ChildNodes[0]);
		BnBNodes.add(ChildNodes[1]);

		//Branch and Bound
		while(BnBNodes.isEmpty() == false){
			
			//Choose next node
			BnBNode node = BnBNodes.iterator().next();
			//Compute new solutions for all scenarios containing branched y
			int[][] s_vt = data.scenarioContainingNode(node.current_y[0], node.current_y[1], node.current_y[2], n);
			int si=0;
			boolean callfathom = true;
		
			while(si<s_vt.length){
				int g = s_vt[si][0];
				int l = s_vt[si][1];
				node.sols[g][l] = splittingVariable(data.ScenarioTree[g][l],node.branched_y);
				//break if one problem is not solvable
				if(node.sols[g][l].solvable == false){
					callfathom = false;
					break;
				}
				si = si + 1;
			}
			double s_result = 0;
			
			for(int g=0; g<G[T-1]; g++){
				for(int l=0; l<L[T-1][g]; l++){
					s_result = s_result + node.sols[g][l].result;
				}
			}
			
			//Check if node must be fathomed
			if(callfathom == true){
				Fathomer fathom = fathom(node, F_U, s_result);
				F_U = fathom.F_U;
				if(fathom.fathom == false){
					ChildNodes = branch(node);
					BnBNodes.add(ChildNodes[0]);
					BnBNodes.add(ChildNodes[1]);
				}
				bestSol = fathom.sol;
			}
						
			//Remove node from list
			BnBNodes.remove(node);
		}
		final long timeEnd = System.currentTimeMillis();
		System.out.println("Optimal objective value:   " + F_U);
		System.out.println("Laufzeit:  " + (timeEnd-timeStart)+"  millisek.");

		saveDataAsFile(data, Configuration.dataFilePath);
	}
	
	//solve Q with fixed y
	public static SolutionOriginal originalModel(double[][][][] y){
		SolutionOriginal result = new SolutionOriginal();
		result.yval = y;
		try{

			IloCplex cplex = new IloCplex();
			cplex.setOut(null);
			
			//DVAR
			IloNumVar[][][][][] x = new IloNumVar[T][][][][];
			for(int t=0; t<T; t++){
				x[t] = new IloNumVar[G[t]][][][];
				for(int g=0; g<G[t]; g++){
					x[t][g] = new IloNumVar[L[t][g]][I][J];
					for(int l=0; l<L[t][g]; l++){
						for (int i=0; i<I;i++){
								x[t][g][l][i] = cplex.numVarArray(J, 0, Double.MAX_VALUE);
						}//for i
					}//for l
				}
			}
			
			//OBJ
			IloLinearNumExpr[][][] exprobj = new IloLinearNumExpr[T][][];
			IloLinearNumExpr summedobj = cplex.linearNumExpr();
			for(int t=0; t<T; t++){
				exprobj[t] = new IloLinearNumExpr[G[t]][];
				for(int g=0; g<G[t]; g++){
					exprobj[t][g] = new IloLinearNumExpr[L[t][g]];
					for(int l=0; l<L[t][g]; l++){
						//for all nodes in Tree:
						exprobj[t][g][l] = cplex.linearNumExpr();
						for (int i=0; i<I;i++){
							for(int j=0; j<J; j++){
								exprobj[t][g][l].addTerm(data.DataTree[t][g][l].P*data.DataTree[t][g][l].alpha[i][j], x[t][g][l][i][j]);
							}
						}
						summedobj.add(exprobj[t][g][l]);
						
					}//end l = L[t][g]
				}
			}
			
			cplex.addMinimize(summedobj); //Min

			
			//S.T.
			for(int t=0; t<T; t++){
				for(int g=0; g<G[t]; g++){
					for(int l=0; l<L[t][g]; l++){
						
						//x <= y;
						for(int i=0; i<I; i++){
							for(int j=0; j<J; j++){
								cplex.addLe(x[t][g][l][i][j], data.DataTree[t][g][l].beq[j]*y[t][g][l][i]);
							}
						}
						
						//sum(i) x = d; forall j
						for(int j=0; j<J; j++){
							cplex.addEq(cplex.sum(x[t][g][l][0][j],x[t][g][l][1][j],x[t][g][l][2][j],x[t][g][l][3][j]), data.DataTree[t][g][l].beq[j]);
						}
					}
				}
			}
			
			//Solve
			cplex.solve();
			double objval = cplex.getObjValue();
			
			for(int t=0; t<T; t++){
				for(int g=0; g<G[t]; g++){
					for(int l=0; l<L[t][g]; l++){
						//for all nodes in Tree:
						for (int i=0; i<I;i++){
							objval = objval + data.DataTree[t][g][l].P*data.DataTree[t][g][l].beta[i]*y[t][g][l][i];
						}						
					}//end l = L[t][g]
				}//end g = G[t]
			}//end t = T
			
			
			result.result = objval;
			//result.xval
			result.xval = new double[T][][][][];
			for(int t=0; t<T; t++){
				result.xval[t] = new double[G[t]][][][];
				for(int g=0; g<G[t]; g++){
					result.xval[t][g] = new double[L[t][g]][I][J];
					for(int l=0; l<L[t][g]; l++){	
						for(int i=0; i<I; i++){
							for(int j=0; j<J; j++){
								result.xval[t][g][l][i][j] = cplex.getValue(x[t][g][l][i][j]);
							}
						}
					}
				}
			}
			
			return result;
		}
		catch(IloException exc){
			System.err.println("Concert exception caught:" + exc);
			result.result = 0;
			result.xval = null;
			return result;
		}
	}
	
	//solve Q_s^SV-NA for a scenario s
	public static SolutionScenario splittingVariable(FullScenario s, int[][][][] branched_y){
		SolutionScenario result = new SolutionScenario();
		
		try{
			
			IloCplex cplex = new IloCplex();
			cplex.setOut(null);
			
			//DVAR
			IloNumVar[][][] x = new IloNumVar[T][I][J];
			IloIntVar[][] y = new IloIntVar[T][I];
			
			for(int t=0; t<T; t++){
				for (int i=0; i<I;i++){
					x[t][i] = cplex.numVarArray(J, 0, Double.MAX_VALUE);
				}
				y[t] = cplex.boolVarArray(I);
			}

			//OBJ
			IloLinearNumExpr exprobj = cplex.linearNumExpr();
			
			for(int t=0; t<T; t++){
				for(int i=0; i<I; i++){
					for(int j=0; j<J; j++){
						exprobj.addTerm(s.P_s*s.alpha_s[t][i][j], x[t][i][j]);
					}
					exprobj.addTerm(s.P_s*s.beta_s[t][i], y[t][i]);
				}
			}
			
			cplex.addMinimize(exprobj); //Min

			//S.T.
			
			//fix branched y
			int tt = T-1;
			int gt = s.index[T-1][1];
			int lt = s.index[T-1][2];
			
			while(tt>=0){
				for(int i=0; i<I; i++){
					if(branched_y[tt][gt][lt][i]==1){
						cplex.addEq(y[tt][i], 0);
					}else if(branched_y[tt][gt][lt][i]==2){
						cplex.addEq(y[tt][i], 1);
					}
				}
				
				int[] pre = Dataset.Predecessor(tt, gt, L);
				tt = pre[0];
				gt = pre[1];
				lt = pre[2];
				
			}
			
			
			for(int t=0; t<T; t++){
				
				//x[t][i][j] <= y[t][i]
				for(int i=0; i<I; i++){
					for(int j=0; j<J; j++){
						cplex.addLe(cplex.diff(x[t][i][j], cplex.prod(s.beq_s[t][j], y[t][i])), 0);
					}
				}
				
				//sum(i) x = d; forall j
				for(int j=0; j<J; j++){
					cplex.addEq(cplex.sum(x[t][0][j],x[t][1][j],x[t][2][j],x[t][3][j]), s.beq_s[t][j]);
				}
				
			}
			
			//time period overlapping constraints
			for(int t=1; t<T; t++){
				for(int i=0; i<I; i++){
					cplex.addGe(cplex.diff(y[t][i], y[t-1][i]), 0);
				}
			}
			
			//Solve
			result.solvable = cplex.solve();
			double objval = cplex.getObjValue(); //Zielsetzung
			
			//DVAR
			result.xval = new double[T][I][J];
			result.yval = new double[T][I];
			
			for(int t=0; t<T; t++){
				for(int i=0; i<I; i++){
					System.arraycopy(cplex.getValues(x[t][i]), 0, result.xval[t][i], 0,J-1);
					result.yval[t][i] = Math.round(cplex.getValue(y[t][i]));
				}
				
			}
			
			result.result = objval;

			return result;
		}
		catch(IloException exc){
			System.err.println("Concert exception caught:" + exc);
			result.result = 0;
			result.xval = null;
			result.yval = null;
			return result;
		}
	}

	//weighted average solution
	public static SolutionOriginal avgSolution(SolutionScenario[][] sol_s){
		
		SolutionOriginal avgSolution = new SolutionOriginal();
		avgSolution.xval = new double[T][][][][];
		avgSolution.yval = new double[T][][][];
		
		for(int t=0; t<T; t++){
			avgSolution.xval[t] = new double[G[t]][][][];
			avgSolution.yval[t] = new double[G[t]][][];
			for(int g=0; g<G[t]; g++){
				avgSolution.xval[t][g] = new double[L[t][g]][I][J];
				avgSolution.yval[t][g] = new double[L[t][g]][I];
			}
		}
		
		int[][][] n_matrix = data.numberScenarioWithNode( );
		
		for(int t=0; t<T; t++){
			for(int g=0; g<G[t]; g++){
				for(int l=0; l<L[t][g]; l++){

					int[][] s_vt = data.scenarioContainingNode(t, g, l, n_matrix);
					
					int s_index = s_vt.length;
					
					//P_summed: could be replaced by dataTree[t][g][l].P
					double P_summed = 0;
					for(int s=0; s<s_index; s++){
						P_summed += data.ScenarioTree[s_vt[s][0]][s_vt[s][1]].P_s;
					}
					
			
					for(int s=0; s<s_index; s++){
						
						for(int i=0; i<I; i++){
							for(int j=0; j<J; j++){
								avgSolution.xval[t][g][l][i][j] += (data.ScenarioTree[s_vt[s][0]][s_vt[s][1]].P_s*sol_s[s_vt[s][0]][s_vt[s][1]].xval[t][i][j])/P_summed;
							}
							avgSolution.yval[t][g][l][i] += (data.ScenarioTree[s_vt[s][0]][s_vt[s][1]].P_s*sol_s[s_vt[s][0]][s_vt[s][1]].yval[t][i])/P_summed;
						}
						
					}
					
				}
			}
		}
		
		return avgSolution;
	}
	
	//round y
	public static SolutionOriginal roundY (SolutionOriginal s){
		
		//simple rounding for period 1
		for(int g=0; g<G[0]; g++){
			for(int l=0; l<L[0][g]; l++){
				for(int i=0; i<I; i++){
					s.yval[0][g][l][i] = Math.ceil(s.yval[0][g][l][i]);
				}
			}
		}
		
		
		for(int t=1; t<T; t++){
			for(int g=0; g<G[t]; g++){
				for(int l=0; l<L[t][g]; l++){
					for(int i=0; i<I; i++){
						s.yval[t][g][l][i] = Math.ceil(s.yval[t][g][l][i]);
						int[] pre = Dataset.Predecessor(t, g, L);
						if(s.yval[t][g][l][i] < s.yval[pre[0]][pre[1]][pre[2]][i]){
							s.yval[t][g][l][i] = 1;
						}
					}
				}
			}
		}
		
		return s;
	}

	//branch node k
	public static BnBNode[] branch(BnBNode parentNode){
		BnBNode[] childNodes = new BnBNode[2];
		childNodes[0] = new BnBNode(T, G);
		childNodes[1] = new BnBNode(T, G);
		childNodes[0].branched_y = parentNode.branched_y;
		childNodes[1].branched_y = parentNode.branched_y;
		
		double dist_smallest = 0.5;
		SolutionOriginal S_orig = avgSolution(parentNode.sols);
		
		for(int t=0; t<T; t++){
			for(int g=0; g<G[t]; g++){
				for(int l=0; l<L[t][g]; l++){
					for(int i=0; i<I; i++){
						if(parentNode.branched_y[t][g][l][i] == 0){
							if(Math.abs(S_orig.yval[t][g][l][i]-0.5) < dist_smallest){
								
								dist_smallest = Math.abs(S_orig.yval[t][g][l][i]-0.5);
	
								childNodes[0].current_y[0] = t;
								childNodes[0].current_y[1] = g;
								childNodes[0].current_y[2] = l;
								childNodes[0].current_y[3] = i;
								
								childNodes[1].current_y[0] = t;
								childNodes[1].current_y[1] = g;
								childNodes[1].current_y[2] = l;
								childNodes[1].current_y[3] = i;								
							}
						}	
					}
				}
			}
		}
		
		//branch/fix all y
		childNodes[0].sols = parentNode.sols;
		childNodes[1].sols = parentNode.sols;
		
		childNodes[0].branched_y[childNodes[0].current_y[0]][childNodes[0].current_y[1]][childNodes[0].current_y[2]][childNodes[0].current_y[3]] = 1;
		childNodes[1].branched_y[childNodes[1].current_y[0]][childNodes[1].current_y[1]][childNodes[1].current_y[2]][childNodes[1].current_y[3]] = 2;
		
		return childNodes;
	}
	
	
	public static Fathomer fathom(BnBNode node, double F_U, double s_result){
		Fathomer F = new Fathomer();
		F.fathom = false;
		F.F_U = F_U;
		
		
		if(s_result >= F_U){
			F.fathom = true;
		}else if(nonant(node)[1]==true){
			if(nonant(node)[0]==true){
				F.F_U = s_result;
				F.sol = avgSolution(node.sols); //stimmt das?
				F.fathom = true;
			}else{
				SolutionOriginal sol_vt = avgSolution(node.sols);
				SolutionOriginal sol_fixed = originalModel(sol_vt.yval);
				F.sol = sol_fixed;
				if(sol_fixed.result < F_U){
					F.F_U = sol_fixed.result;
				}
				F.fathom = false;
			}
		} else{
			SolutionOriginal sol_vt = avgSolution(node.sols);
			SolutionOriginal sol_rounded = roundY(sol_vt);
			SolutionOriginal complete_sol_rounded = originalModel(sol_rounded.yval);
			F.sol = complete_sol_rounded;
			if(complete_sol_rounded.result < F_U){
				F.F_U = complete_sol_rounded.result;
			}
			F.fathom = false;
		}
		
		return F;
	}
	
	public static boolean[] nonant(BnBNode node){
		boolean[] nonant = new boolean[2]; //for x and for y
		nonant[1] = true;
		
		int[][][] n_matrix = data.numberScenarioWithNode( );

		SolutionOriginal sol_avg = avgSolution(node.sols);
		
		//check non-anticipativity for y
		int t=0;
		while(t<T){
			for(int g=0; g<G[t]; g++){
				for(int l=0; l<L[t][g]; l++){
					for(int i=0; i<I; i++){
						if(sol_avg.yval[t][g][l][i]<1&&sol_avg.yval[t][g][l][i]>0){
							nonant[1]=false;
							break;
						}
					}
				}
			}
			t = t+1;
		}
		
		//check non-anticipativity for x
		t=0;
		while(t<T){
			for(int g=0; g<G[t]; g++){
				for(int l=0; l<L[t][g]; l++){
					int[][] s_vt = data.scenarioContainingNode(t, g, l, n_matrix);
					int s_index = s_vt.length;

					for(int s=0; s<s_index; s++){
						for (int i = 0; i < I; i++) {
							for (int j = 0; j < J; j++) {
								if(sol_avg.xval[t][g][l][i][j] != node.sols[s_vt[s][0]][s_vt[s][1]].xval[t][i][j]){
									nonant[0]=false;
									break;
								}
							}
						}
					}
		
				}
			}
			t = t+1;
		}
		
		
		return nonant;
	}
	
	
}
