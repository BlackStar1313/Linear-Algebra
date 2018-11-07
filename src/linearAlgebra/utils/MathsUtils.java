package linearAlgebra.utils;

import static java.lang.Math.*;
import java.util.*;

import linearAlgebra.Exception.MatrixException;
import linearAlgebra.Matrix.Matrix;
import linearAlgebra.vectors.Vector;

public final class MathsUtils {

	/** sqrt(a^2 + b^2) without under/overflow. **/
	public static double hypot(double a, double b) {
		double r = 0D;
		if (abs(a) > abs(b)) {
			r = b/a;
			r = abs(a) * sqrt(1+r*r);
		} else if (b != 0) {
			r = a/b;
			r = abs(b)* sqrt(1+r*r);
		} else {
			r = 0D;
		}
		return r;
	}
	
	/**
	 * Generate a random number between minimum number(inclusive) and maximum number(inclusive).
	 * 
	 * @param min		The minimum number.
	 * @param max		The maximum number.
	 * @return			return the random number in range of min and max. Otherwise it returns -1 if the min is greater than the max.
	 */
	public static int Random(int min, int max) throws IllegalStateException{
		int randomInt = 0;
		
		if(min > max) {
			throw new IllegalStateException("the max number must be grater than the min");
		}else if(min == max) {
			return min;
		}else {
			Random rand = new Random();
			randomInt = rand.nextInt((max - min) + 1) + min;
		}
		return randomInt;
	}
	
	/**
	 * Solves for the vector y such that <code>Ay = b</code>
	 * @param A		an upper triangular matrix.
	 * @param b		vector which length is equal to the rows in A.
	 * @return		y such that <code>Ay = b</code>
	 */
	public static Vector backSub(Matrix A, Vector b) {
		
		if(!A.isTriangularSuperior()) {
			throw new MatrixException("Your matrix is not in an upper triangular form!!");
		}
		
		int n = A.rows;
		
		Vector y = b.deepCopy();
		
		for(int i = n-1; i >= 0; i--) {
			
			double b_i = y.components[i];
			
			for(int j = i+1; j < n; j++) {
				b_i -= (A.data[i][j] * y.components[j]);
			}
			
			//Occurs when A[i][i] = 0
			if(Double.isInfinite(b_i)) {
				b_i = 0D;
			}
			
			b_i /= A.data[i][i];
			y.components[i] = b_i;
		}
		
		return y;
	}
	
	/**
	 * Solves for the vector x such that <code>Ax = y</code>
	 * @param A		A lower triangular matrix.
	 * @param vect	A vector which length is equal tot the rows in A.
	 * @return		x such that <code>Ax = y</code>
	 */
	public static Vector forwardSub(Matrix A, Vector y) {
		
		if(!A.isTriangularInferior()) {
			throw new MatrixException("Your matrix is not in an lower triangular form!!");
		}
		int n = A.rows;
		
		Vector x = y.deepCopy();
		
		for(int i = 0; i < n; i++) {
			
			double x_i = x.components[i];
			
			for(int j = 0; j < i; j++) {
				x_i -= (A.data[i][j] * x.components[j]);
			}
			
			x_i /= A.data[i][i];
			x.components[i] = x_i;
		}
		
		return x;
	}
}
