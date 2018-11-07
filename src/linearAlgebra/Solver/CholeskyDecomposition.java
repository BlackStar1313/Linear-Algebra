package linearAlgebra.Solver;

import java.io.Serializable;

import linearAlgebra.Exception.MatrixException;
import linearAlgebra.Exception.VectorException;
import linearAlgebra.Matrix.Matrix;
import linearAlgebra.vectors.Vector;
import static java.lang.Math.*;
import static linearAlgebra.utils.MathsUtils.*;

public final class CholeskyDecomposition implements ISolver, Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -8883111599884051447L;

	/**
	 * An array to store the decomposition
	 */
	public final Matrix L;
	
	private final Matrix A;
	
	
	/**
	 * Cholesky algorithm for symmetric and positive definite matrix
	 * 
	 * @param A
	 */
	public CholeskyDecomposition(final Matrix A) {
		
		if(!A.isSquare()) {
			throw new MatrixException("The matrix to be decomposed is not square!!");
		}else if(!A.isSymmetric()) {
			throw new MatrixException("The matrix to be decomposed is not symmetric!!");
		}
		
		this.A = A;
		
		this.L = new Matrix(A.rows, A.cols);
		
		decompose();
		
	}
	
	private void decompose() {
		
		for(int i = 0; i < A.rows; i++) {
			
			//diagonal element of the matrix A
			double diag = A.data[i][i];
			
			for(int k = 0; k <= i; k++) {
				double sum = 0D;
				for(int j = 0; j < k; j++) {
					sum += (L.data[i][j] * L.data[k][j]);
				}
				
				if(i == k) {
					L.data[i][i] = sqrt(diag - sum);
				}else {
					L.data[i][k] = (A.data[i][k] - sum) / L.data[k][k];
				}
			}
			
			if(L.data[i][i] <= 0D) {
				throw new MatrixException("The matrix to be decomposed is not positive definite!!");
			}
		}
	}
	
	/**
	 * Check whether A == LL<sup>T</sup>
	 * 
	 * @return true if A == LL<sup>T</sup>, otherwise false.
	 */
	public boolean check(final Matrix A) {
		
		Matrix temp = L.Multiply(L.Transpose());
		
		for(int i = 0; i < A.rows; i++) {
			for(int j = 0; j < A.cols; j++) {
				if(A.data[i][j] != temp.data[i][j]) {
					return false;
				}
			}
		}
		return true;
	}

	@Override
	public boolean isNonSingular() {
		return L.Determinant() != 0;
	}

	@Override
	public Vector Solve(Vector b) {
		
		if(b.length != L.rows) {
			throw new VectorException("Your vector length does not correspond to the lenth of the matrix!!");
		}
		
		/**
		 * First solve L y = b;
		 */
		Vector y = forwardSub(L, b);
		
		/**
		 * Solve L<sup>T</sup>x = y;
		 */
		Vector x = backSub(L.Transpose(), y);
		
		return x;
	}
	
	@Override
	public String toString() {
		return L.toString();
	}

}
