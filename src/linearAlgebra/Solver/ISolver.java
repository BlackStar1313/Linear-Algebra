package linearAlgebra.Solver;

import linearAlgebra.Matrix.Matrix;
import linearAlgebra.vectors.Vector;

public interface ISolver {
	
	boolean isNonSingular();
	
	/**
	 * Solves the linear system of equations Ax = b.
	 * 
	 * @param b	a vector.
	 * @return
	 */
	Vector Solve(Vector b);
}
