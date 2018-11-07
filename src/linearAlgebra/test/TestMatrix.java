package linearAlgebra.test;

import linearAlgebra.Exception.MatrixException;
import linearAlgebra.Matrix.Matrix;
import linearAlgebra.Solver.CholeskyDecomposition;
import linearAlgebra.vectors.Vector;

public class TestMatrix {

	public static void main(String[] args) {
		
		double[][] M1 = {
				{1, 1, 0, 0, 0, 0, -5200},
				{0, 0, 1, 1, 0, 0, -7400},
				{0, 0, 0, 0, 1, 1, -3900},
				{-1, 0, 1, 0, 0, 1, 0},
				{0, 1, 0, -1, 0, 1, 0},
				{0, 1, 1, 0, -1, 0, 0}
		};
		
		double[][] A = { 
				{2, -2, -3},
				{-2, 5, 4},
				{-3, 4, 5}
				/*{1, 1, 1, 1},
				{1, 5, 5, 5},
				{1, 5, 14, 14},
				{1, 5, 14, 15}*/
		};
		
		double[][] M2 = {
				{2, 0, 0},
				{17, -2, 0},
				{4, 3, 20}
		};
		
		double[][] M3 = {
				{1, -1, 4},
				{2, 1, 1},
				{0, 1, 3}
		};
		Matrix a = new Matrix(A);
		CholeskyDecomposition chol = new CholeskyDecomposition(a);
		System.out.println(chol);
		Vector b = new Vector(7, -12, -12);
		System.out.println(chol.Solve(b));
		
		a = null;
		b = null;
		//System.out.println("matrice triangulaire inférieur? : \n" + a.isTriangularSuperior());
		//System.out.println("est : \n" + a.Inverse());
		/*A.GaussianElimination();
		System.out.println(A + "\nrank = " + A.rank());*/
		//System.out.println(A.norm1());
		/*Matrix B = new Matrix(M2);
		Matrix C = null;
		try {
			C = A.Multiply(B);
		} catch (MatrixException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		System.out.println("original matrix \n" + C);
		
		Matrix d = C.deepCopy();
		
		System.out.println("\nduplicated matrix \n" + d);
		
		System.out.println(d.MultiplyByScalar(2));
		System.out.println("original matrix \n" + C);*/
		
		/*try {
			System.out.println(A.Multiply(B));
		} catch (MatrixException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}*/
	}

}
