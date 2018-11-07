package linearAlgebra.Solver;

import java.io.Serializable;

import linearAlgebra.Exception.MatrixException;
import linearAlgebra.Matrix.Matrix;
import linearAlgebra.vectors.Vector;


/**
 * LU Decomposition.
 * <P>
 * For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n unit
 * lower triangular matrix L, an n-by-n upper triangular matrix U, and a
 * permutation vector piv of length m so that A(piv,:) = L*U. If m < n, then L
 * is m-by-m and U is m-by-n.
 * <P>
 * The LU decompostion with pivoting always exists, even if the matrix is
 * singular, so the constructor will never fail. The primary use of the LU
 * decomposition is in the solution of square systems of simultaneous linear
 * equations. This will fail if isNonsingular() returns false.
 */
public final class LUDecomposition implements ISolver, Serializable {

	private static final long serialVersionUID = 1L;

	/** Entries of LU decomposition. */
	private Matrix LU;

	/** Pivot permutation associated with LU decomposition. */
	private final int[] pivot;

	/** Cached value of L. */
	private Matrix L;

	/** Cached value of U. */
	private Matrix U;

	/** Cached value of P. */
	private Matrix P;

	/**
	 * Row and column dimensions, and pivot sign.
	 * 
	 * @serial row dimension.
	 * @serial column dimension.
	 * @serial pivot sign.
	 */
	private final int mRows, nCols;
	private int pivotSign;

	public LUDecomposition(final Matrix matrix) {

		if(!matrix.isSquare()) {
			throw new MatrixException("The matrix is not squared!!");
		}

		this.LU = matrix.deepCopy();
		this.mRows = matrix.rows;
		this.nCols = matrix.cols;

		// Initialize permutation array and parity
		this.pivot = new int[mRows];
		for (int i = 0; i < mRows; i++) {
			pivot[i] = i;
		}
		pivotSign = 1;

		// Loop over columns
		for(int col = 0; col < nCols; col++) {





			// Divide the lower elements by the "winning" diagonal elt.
			final double luDiag = LU.data[col][col];
			for (int row = col + 1; row < mRows; row++) {
				LU.data[row][col] /= luDiag;
			}
		}
	}
	
	public Matrix getL() {
		return null;
	}
	
	public Matrix getU() {
		return null;
	}

	@Override
	public boolean isNonSingular() {
		for(int i = 0; i < nCols; i++) {
			if(LU.data[i][i] == 0) {
				return false;
			}
		}
		return true;
	}

	@Override
	public Vector Solve(Vector b) {
		// TODO Auto-generated method stub
		return null;
	}

}
