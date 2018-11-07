package linearAlgebra.Matrix;

import java.util.Scanner;

import static java.lang.Math.*;
import static linearAlgebra.utils.MathsUtils.*;

import linearAlgebra.Exception.MatrixException;
import linearAlgebra.vectors.Vector;

public final class Matrix {

	private static final int DEFAULT_ROWS = 2;
	private static final int DEFAULT_COLUMNS = 2;
	
	public final double[][] data;
	public final int rows;
	public final int cols;
	
	//Different Constructor
	/**
	 * Constructor with a matrix of 2 x 2;
	 * @throws MatrixException
	 */
	public Matrix(){
		this(DEFAULT_ROWS, DEFAULT_COLUMNS);//we call the constructor just right under this one..
	}
	
	/**
	 * Constructor to create a new matrix of row and col.
	 * 
	 * @param row	The number of rows.
	 * @param col	The number of columns.
	 * @throws MatrixException 		throws an exception if we try to create a matrix smaller than 1X1.
	 */
	public Matrix(int row, int col){
		if(row < 1 || col < 1) {
			throw new IndexOutOfBoundsException("You cannot create a matrix smaller than 1 X 1");
		}else {
			this.rows = row;
			this.cols = col;
			this.data = new double[this.rows][this.cols];
		}
	}
	
	/**
	 * Constructor to create a new matrix of row and columns filled with a scalar value.
	 * @param row		The number of rows.
	 * @param col		The number of columns.
	 * @param scalarValue		The scalar value.
	 * @throws MatrixException
	 */
	public Matrix(int row, int col, double scalarValue){
		this(row, col);
		setMatrix(scalarValue);
	}
	
	
	/**
	 * Constructor for creating a new matrix from the one in argument.
	 * @param matrix		the matrix.
	 */
	public Matrix(double[][] matrix){
		this.data = matrix;
		this.rows = matrix.length;
		this.cols = matrix[0].length;
	}
	
	/**
	 * Copy Constructor
	 * @param anotherMatrix
	 */
	public Matrix(Matrix anotherMatrix) {
		this.data = anotherMatrix.data;
		this.rows = anotherMatrix.rows;
		this.cols = anotherMatrix.cols;
	}
	
	/*public double get(int row, int col){
		return data[row][col];
	}
	
	//Setters...
	public void set(int row, int col, double value){
		data[row][col] = value;
	}*/
	
	
	/**
	 * We redefine the function equals.
	 */
	public boolean equals(Object o) {
		//self-check
		if(this == o) {
			return true;
		//null check and type check and cast. Here we don't use instanceof because we only want to check if the object is a matrix and not an inherited object.
		}else if(o == null || getClass() != o.getClass()) {
			return false;
		}else {
			Matrix matrix = (Matrix)o;
			return (rows == matrix.rows) && (cols == matrix.cols) && (this.data == matrix.data);
		}
	}
	
	/**
	 * Check if the current matrix is diagonal.
	 * That is to say if (i != j) and if matrix(i,j) = 0. 
	 * @return	true if the current matrix is diagonal, otherwise false.
	 */
	public boolean isDiagonal(){
		Matrix matrix = this;
		
			if(!matrix.isSquare()) {
				throw new MatrixException("Your matrix is not squared!");
			}else {
				for(int i = 0; i < matrix.rows; i++) {
					for(int j = 0; j < matrix.cols; j++) {
						if(i != j && matrix.data[i][j] != 0) {
							return false;
						}
					}
				}
			}
		return true;
	}
	
	/**
	 * Check if the current matrix is triangular inferior. 
	 * That is to say if (i > j) and if matrix(i,j) = 0.
	 * @return
	 */
	public boolean isTriangularSuperior(){
		Matrix matrix = this;
		
			if(!matrix.isSquare()) {
				throw new MatrixException("Your matrix is not squared!");
			}else {
			
				for(int i = 0; i < matrix.rows; i++) {
					for(int j = 0; j < i; j++) {
						if(matrix.data[i][j] != 0) {
							return false;
						}
					}
				}
			}
		return true;
	}
	
	/**
	 * Check if the matrix is triangular inferior. 
	 * That is to say if (i < j) and if matrix(i,j) = 0.
	 * @return		true if the current matrix is triangular inferior, otherwise false.
	 */
	public boolean isTriangularInferior(){
		Matrix matrix = this;
		
			if(!matrix.isSquare()) {
				throw new MatrixException("Your matrix is not squared!");
			}else {
				for(int i = 0; i < matrix.rows; i++) {
					for(int j = i+1; j < matrix.cols; j++) {
						if(matrix.data[i][j] != 0) {
							return false;
						}
					}
				}
			}
		return true;
	}
	
	/**
	 * Check if the current matrix is symmetric.
	 * That is to say if matrix(i,j) == matrix(j,i).
	 * @return		true if the current matrix is symmetric, otherwise false.
	 */
	public boolean isSymmetric(){
		Matrix matrix = this;
		
			if(!matrix.isSquare()) {
				throw new MatrixException("Your matrix is not squared!");
			}else {
				for(int i = 0; i < matrix.rows; i++) {
					for(int j = 0; j < matrix.cols; j++) {
						if(matrix.data[i][j] !=  matrix.data[j][i]) {
							return false;
						}
					}
				}
			}
		return true;
	}
	
	/**
	 * Check if the current is anti-symmetric.
	 * That is to say if matrix(i,j) == -matrix(j,i).
	 * @return		true if the current is anti-symmetric, otherwise false.
	 */
	public boolean isAntiSymmetric(){
		Matrix matrix = this;
		
			if(!matrix.isSquare()) {
				throw new MatrixException("Your matrix is not squared!");
			}else {
			
				for(int i = 0; i < matrix.rows; i++) {
					for(int j = 0; j < matrix.cols; j++) {
						if(matrix.data[i][j] !=  -matrix.data[j][i]) {
							return false;
						}
					}
				}
			}
		return true;
	}
	
	/**
	 * Check if this matrix size equals to a given a matrix B. 
	 * @param B			the Matrix B
	 * @return			true if both of them have the same size, otherwise false.
	 */
	public boolean isSameDimensionsAs(final Matrix B) {
		if(this == B) {
			return true;
		}else if(B == null) {
			return false;
		}
		
		return (this != null && (this.cols == B.cols) && (this.rows == B.rows));
	}
	
	
	/**
	 * Check whether the multiplication of two matrix is possible that is to say if the number of columns for a given matrix A
	 * is the same number of rows of a given matrix B.
	 * @param B		the other matrix to check to
	 * @return		true if the multiplication is doable, otherwise false.
	 */
	public boolean isMultiplicationPossible(final Matrix B) {
		if(this == B) {
			return true;
		}else if(B == null) {
			return false;
		}
		
		return (this != null && this.cols == B.rows);
	}
	
	
	/**
	 * Check if the matrix is squared.
	 * @param data	the matrix to check.
	 * @return			true if the matrix is squared, otherwise false.
	 */
	public boolean isSquare(){
		Matrix matrix = this;
		return (matrix.rows == matrix.cols);
	}
	
	/**
	 * Check if the current matrix is invertible.
	 * @param matrix	the current matrix.
	 * @return			true if its determinant is non-zero, otherwise false.
	 */
	public boolean isInvertible(final Matrix matrix){
		return (this.Determinant() != 0);
	}
	
	/**
	 * The determinant is a useful value that can be computed from the elements of a square matrix. The determinant of a matrix A is denoted <i>det(A), det A, or |A|</i>.
	 * <n>It can be viewed as the scaling factor of the transformation described by the matrix.
	 * If the current matrix is 1X1 dimension then its determinant is just itself.
	 * If the current matrix is 2X2 dimensions then its determinant is (ad-bc).
	 * if the current matrix is triangular superior or inferior then its determinant is the product of numbers on its diagonal.
	 * otherwise we compute its determinant from its column.
	 * @param data  the current matrix.
	 * @return			the determinant of the matrix.
	 * @throws MatrixException		throws an exception if the matrix is not squared.
	 */
	public double Determinant(){
		Matrix matrix = this;
		double det = 0.0;
		
		if(!matrix.isSquare()){
			throw new MatrixException("Your matrix is not squared!");
		}else {
			if(matrix.rows == 1){
				return matrix.data[0][0];
			}else if(matrix.rows == 2){
				return (matrix.data[0][0] * matrix.data[1][1]) - ((matrix.data[0][1] * matrix.data[1][0]));
			}else if(matrix.isTriangularInferior() || matrix.isTriangularSuperior()){
				det = 1.0;
				for(int i = 0; i < matrix.rows; i++) {
					det *= matrix.data[i][i];
				}
			}else{
				for(int col = 0; col < matrix.rows; col++){
					Matrix minor = Minor(col, 0);
					det += Math.pow(-1, col) * matrix.data[col][0] * minor.Determinant();
				}
			}
		}
		return det;
	}
	
	/**
	 * Calculates the trace of the current matrix if that is squared.
	 * That is to say the sum of all element within its diagonal.
	 * @return		the trace of the current matrix.
	 * @throws MatrixException
	 */
	public double trace(){
		double trace = 0.0;
		Matrix matrix = this;
		
		if(!matrix.isSquare()) {
			throw new MatrixException("Your matrix is not squared!");
		}else {
			for(int i = 0; i < matrix.rows; i++) {
				trace += matrix.data[i][i];
			}
		}
		return trace;
	}
	
	public double norm1() {
		double maxSum = 0D;
		Matrix matrix = this;
		
		for(int j = 0; j < this.cols; j++) {
			double currentSum = 0D;
			for(int i = 0; i < matrix.rows; i++) {
				currentSum += (abs(matrix.data[i][j]));
			}
			if(maxSum < currentSum) {
				maxSum = currentSum;
			}
		}
		return maxSum;
	}
	
	/**
	 * Compute the infinity norm which is the maximum row sum.
	 * 
	 * @return the maximum row sum of this matrix.
	 */
	public double normInf() {
		
		double maxSum = 0D;
		Matrix matrix = this;
		
		for(int i = 0; i < this.rows; i++) {
			double currentSum = 0D;
			for(int j = 0; j < this.cols; j++) {
				currentSum += (abs(matrix.data[i][j]));
			}
			
			if(maxSum < currentSum) {
				maxSum = currentSum;
			}
		}
		return maxSum;
	}
	
	
	/**
	 * Executes the swap of two rows in bidimensional array.
	 * @param first		The first element.
	 * @param second	The second element.
	 */
	private void swapRows(final int first, final int second){
		Matrix matrix = this;
		double[] temp = matrix.data[first];
		matrix.data[first] = matrix.data[second];
		matrix.data[second] = temp;
	}
	
	/**
	 * Swap two elements in bi-dimensional array.
	 * 
	 * @param row1		The first row.
	 * @param col1		The first column.
	 * @param row2		The second row.
	 * @param col2		The second column.
	 */
	private void swapElementAt(final int row1, final int col1,  final int row2, final int col2) {
		Matrix matrix = this;
		double temp = data[row1][col1];
		data[row1][col1] = data[row2][col2];
		data[row2][col2] = temp;
	}
	
	
	/**
	 * Here we set the asked matrix.
	 * @return  a matrix in which we have just filled.
	 */
	public Matrix Fill(){
		Matrix fill = this;
		Scanner scan = new Scanner(System.in);
		
		for(int i = 0; i < fill.rows; i++){
			for(int j = 0; j < fill.cols; j++){
				System.out.print("matrice[" + (i + 1) + "]" + "[" + (j + 1) + "]" + " = ");
				double v = scan.nextDouble();
				fill.data[i][j] = v;
			}
		}
		
		return fill;
	}
	
	/**
	 * Create a new matrix in which rows and columns have changed.
	 * That is to say A(i,j) = A(j,i) which means A = A<sup>T</sup>.
	 * @return The transposed Matrix.
	 */
	public Matrix Transpose(){
		Matrix matrix = this;
		Matrix transposed = new Matrix(this.cols, this.rows);
		for(int i = 0; i < transposed.rows; i++){
			for(int j = 0; j < transposed.cols; j++){
				transposed.data[i][j] =  matrix.data[j][i];
			}
		}
		return transposed;
	}
	
	/**
	 * Perform the calculation of the current matrix to the power n.
	 * Which means A<sup>n</sup>.
	 * By default A<sup>0</sup> equals to I which is the <tt><i>matrix identity</i></tt>.
	 * @param power		the power n.
	 * @return			A<sup>n</sup>.
	 */
	public Matrix toPower(final int power){
		Matrix matrix = this;
		if(power < 0) {
			throw new MatrixException("The power n is smaller than 0");
		}else if(power == 0) {
			matrix = identity(matrix.rows);
		}else {
			for(int k = 1; k < power; k++) {
				matrix = Multiply(matrix);
			}
		}
		return matrix;
	}
	
	
	public Matrix Minor(int removedRow, int removedCol){
		Matrix subMatrix = new Matrix(this.rows - 1, this.cols - 1);
		Matrix matrix = this;
		
			if((removedRow < 0 || removedCol < 0) || (removedRow >= this.rows || removedCol >=this.cols)){
				throw new IndexOutOfBoundsException("The index of the removed row or the removed column is out of bound !");
			}else {
			
				for(int i = 0; i < this.rows; i++){
					int row = 0;
					int col = 0;
					for(int j = 0; (i != removedRow) && j < this.cols; j++){
						if(j != removedCol){
							row = (i < removedRow) ? i : i - 1;
							col = (j < removedCol) ? j : j - 1;
							subMatrix.data[row][col] = matrix.data[i][j];
						}
					}
				}
			}
		return subMatrix;
	}
	
	/**
	 * Compute the inverse of this matrix with the sub-matrix methods.
	 * @return
	 * @throws MatrixException
	 */
	public Matrix Inverse(){
		Matrix matrix = this;
		Matrix inverse = null;
		
		if(!isInvertible(this)){//"this" means the current matrix!
			inverse = matrix;
			throw new MatrixException("The determinant of the matrix equals to zero!");
		}else {
			//creating a new matrix.
			inverse = new Matrix(matrix.rows, matrix.cols);
			
			//minors and cofactors.
			for(int i = 0; i < matrix.rows; i++) {
				for(int j = 0; j < matrix.cols; j++) {
					Matrix minor = Minor(i , j);
					inverse.data[i][j] = Math.pow(-1, i + j) * minor.Determinant();
				}
			}
			
			//adjugate and determinant..
			double det = 1 / matrix.Determinant();
			inverse = inverse.Transpose();
			inverse.MultiplyByScalar(det);
		}
		
		return inverse;
	}
	
	/**
	 * Computer the multiplication of two given matrices. That is to say C = A * B.
	 * 
	 * @param B		the ohter matrix
	 * @return		a new matrix from the multiplication of A and B-
	 * @throws MatrixException
	 */
	public Matrix Multiply(final Matrix B){
		Matrix A = this;
		Matrix C = new Matrix(A.rows, B.rows);
	
		if(!isMultiplicationPossible(B)){
			throw new MatrixException("The matrices are not the same dimension !");
		}else {
			for(int i = 0; i < A.rows; i++){
				for(int j = 0; j < B.cols; j++){
					for(int k = 0; k < A.cols; k++){
						C.data[i][j] += A.data[i][k] * B.data[k][j];
					}
				}
			}
		}
		return C;
	}
	
	/**
	 * 
	 * @param B
	 * @return
	 * @throws MatrixException
	 */
	public Matrix Subtract(Matrix B) {
		Matrix A = this;
		Matrix C = new Matrix(A.rows, A.cols);
		
		if(!isSameDimensionsAs(B)) {
			throw new MatrixException("The column's length of A and B are not the same dimension!");
		}else {
			for(int i = 0; i < A.rows; i++) {
				for(int j = 0; j < A.cols; j++) {
					C.data[i][j] = A.data[i][j] - B.data[i][j];
				}
			}
		}
		return C;
	}
	
	/**
	 * C = A + B.
	 * @param B
	 * @return
	 * @throws MatrixException
	 */
	public Matrix plus(final Matrix B){
		Matrix A = this;
		Matrix C = new Matrix(A.rows, A.cols);
		
		if(!isSameDimensionsAs(B)) {
			throw new MatrixException("The column's length of A and the row's length of B are not the same dimension!");
		}else {
			for(int i = 0; i < A.rows; i++) {
				for(int j = 0; j < A.cols; j++) {
					C.data[i][j] = A.data[i][j] + B.data[i][j];
				}
			}
		}
		return C;
	}
	
	/**
	 * Multiply a matrix by a vector.
	 * @param vector	the vector 
	 * @return			a vector.
	 */
	public Vector times(Vector vector){
		Matrix matrix = this;
		Vector newVector = new Vector(vector.length);
		
		if(matrix.cols != vector.length) {
			throw new MatrixException("The matrix does not have the same size as the vector!");
		}else {
			for(int i = 0; i < matrix.rows; i++){
				for(int j = 0; j < matrix.cols; j++){
					newVector.components[i] += (matrix.data[i][j] * vector.components[j]); 
				}
			}
		}
		return newVector;
	}
	
	/**
	 * Create a distinct copy of this matrix.
	 * 
	 * @return a distinct copy of this matrix.
	 */
	public Matrix deepCopy() {
		Matrix duplicated = new Matrix(this.rows, this.cols);
		for(int i = 0; i < duplicated.rows; i++) {
			for(int j = 0; j < duplicated.cols; j++) {
				duplicated.data[i][j] = this.data[i][j];
			}
		}
		return duplicated;
	}
	
	/**
	 * Create a copy of this matrix.
	 * 
	 * @return a copy of this matrix.
	 */
	public Matrix copy() {
		Matrix duplicated = new Matrix(this);
		return duplicated;
	}
	
	/**
	 * Creates a diagonal Matrix
	 * 
	 * @param values
	 * @throws MatrixException
	 */
	public static Matrix diag(final double... values){
		Matrix d = null;
		
		if(values == null || values.length == 0) {
			throw new IllegalArgumentException("We cannot create a diagonal matrix with no values!!");
		}else {
			int row = values.length, col = values.length;
			d = new Matrix(row, col);
			for(int i = 0; i < values.length; i++) {
				d.data[i][i] = values[i];
			}
		}
		
		return d;
	}
	
	
	/**
	 * Create the matrix identity of a given size n.
	 * 
	 * @param n		the number of both rows and cols.
	 * @return		An n x n matrix with ones on the diagonal and zeros elsewhere
	 */
	public static Matrix identity(int n){
		Matrix id = new Matrix(n, n);
		for(int k = 0; k < id.rows; k++) {
			id.data[k][k] = 1.0;
		}
		return id;
	}
	
	/**
	 * Create a random matrix in which digits are set randomly.
	 * 
	 * @param rows				the number of rows.
	 * @param cols				the number of columns.
	 * @param minNUmber			the minimum number
	 * @param maxNumber			the maximum number
	 * @param rows				the number of rows of the matrix.
	 * @param cols				the number of columns of the matrix.
	 * @return					a random matrix that digits are set randomly in range of minNumber and MaxNumber.
	 */
	public static Matrix generateRandomMatrix(final int minNumber, final int maxNumber, final int rows, final int cols){
		Matrix randomMatrix = new Matrix(rows, cols);
		
		for(int j = 0; j < randomMatrix.rows; j++){
			for(int k = 0; k < randomMatrix.cols; k++){
				randomMatrix.data[j][k] =  Random(minNumber, maxNumber);
			}
		}
		return randomMatrix;
	}
	
	
	/**
	 * Fill the current matrix with a scalar value.
	 * @param value		the scalar value.
	 */
	public void setMatrix(double value) {
		for(int i = 0; i < this.rows; i++) {
			for(int j = 0; j < this.cols; j++) {
				data[i][j] = value;
			}
		}
	}
	
	/**
	 * A = scalar * A;
	 * 
	 * @param scalar
	 * @return
	 */
	public void MultiplyByScalar(final double scalar){
		Matrix matrix = this;
		
		for(int i = 0; i < matrix.rows; i++) {
			for(int j = 0; j < matrix.cols; j++) {
				matrix.data[i][j] *= scalar;
			}
		}
	}
	
	
	/**
	 * Gaussian elimination (also known as row reduction) is an algorithm for solving systems of linear equations.
	 * It is usually understood as a sequence of operations performed on the corresponding matrix of coefficients. 
	 * This method can also be used to find the rank of a matrix, to calculate the determinant of a matrix, and to calculate the inverse of an invertible square matrix. 
	 * @return			the reduced row echelon form of the matrix.
	 */
	public void GaussianElimination(){
		Matrix matrix = this;
		int r = 0;
		
		for(int col = 0; col < min(matrix.rows, matrix.cols); col++) {
			int k = r;
			for(int row = r + 1; row < matrix.rows; row++) {
				if(abs(matrix.data[row][col]) > abs(matrix.data[k][col])) {
					k = row;
				}
			}
		
			//We want to avoid NAN -> Not A Number which occurs when we try to divide by zero!!!
			if(abs(matrix.data[k][col]) != 0) {
				double divisor = 1.0 / matrix.data[k][col];
				for(int c = 0; c < matrix.cols; c++) {
					matrix.data[k][c] *= divisor;
				}
				
				//if k equals to r, there is no need to switch rows....
				if(k != r)matrix.swapRows(k, r);
				
				for(int i = 0; i < matrix.rows; i++) {
					if(i != r) {
						double simplify = matrix.data[i][col];
						for(int j = 0; j < matrix.cols; j++) {
							matrix.data[i][j] -= (simplify * matrix.data[r][j]);
						}//fin de la boucle de j
					}//fin de if(i!=r)...
				}//fin de la boucle de i....
				r++;
			}//fin de if(M(k, col) != 0)...
		}//fin de la première boucle de col...
	}
	
	/**
	 *  the rank of a matrix A is the dimension of the vector space generated (or spanned) by its columns.
	 *  This corresponds to the maximal number of linearly independent columns of A.
	 *  In other words, the rank is calculated from the reduced row echelon form.
	 *  That is to say the number of pivots (or basic columns) and also the number of non-zero rows.
	 * @return		the rank of the matrix.
	 */
	public int rank() {
		Matrix matrix = this;
		int rank = 0;
		double pivot = 1.0;
		
		for(int i = 0; i < matrix.rows; i++) {
			boolean firstPivot = false;
			int indexOfPivot = -1;
			for(int j = 0; j < matrix.cols && !firstPivot; j++) {
				if(matrix.data[i][j] == pivot) {
					indexOfPivot = j;
					firstPivot = true;
				}
			}
			
			//If we have found the pivot then start counting the zero..otherwise increment the very first for loop.
			if(firstPivot) {
				int rowIndex = 0;
				int count = 0;
				
				while(rowIndex < matrix.rows) {
					if(matrix.data[rowIndex][indexOfPivot] != pivot) {
						if(matrix.data[rowIndex][indexOfPivot] == 0) {
							count++;
						}
					}
					rowIndex++;
				}
				
				if(count == matrix.rows - 1) {
					rank++;
				}
			}
		}
		return rank;
	}
	
	public String toString(){
		String text = "";
		
		for(int i = 0; i < this.rows; i++){
			for(int j = 0; j < this.cols; j++){
				text += data[i][j] + "\t";
			}
			text += "\n";
		}
		return text;
	}
}
