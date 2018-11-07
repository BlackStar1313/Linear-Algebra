package linearAlgebra.vectors;
import linearAlgebra.Exception.VectorException;

public final class Vector{

	public final double[] components;
	public final int length;


	public Vector(int length) {
		this.length = length;
		this.components = new double[this.length];
	}

	public Vector(double...value) {
		if(value == null || value.length == 0)
			throw new IllegalArgumentException("You must provide some value to the vector!");

		length = value.length;
		components = new double[this.length];

		for(int i = 0; i < length; i++)
			components[i] = value[i];
	}

	public Vector(Vector v) {
		this.components = v.components;
		this.length = v.length;
	}


	public double dot(Vector B) throws VectorException{
		Vector A = this;
		double value = 0.0;

		if(!A.isSameDimension(B)) {
			throw new VectorException("The dimensions are not equal!");
		}else {
			for(int i = 0; i < A.length; i++)
				value += (A.components[i] * B.components[i]);
		}
		return value;
	}

	/**
	 * 
	 * @param B
	 * @return
	 */
	public Vector plus(Vector B) {
		Vector A = this;
		Vector C = null;

		if(!A.isSameDimension(B)) {
			throw new VectorException("The dimensions are not equal!");
		}else {
			C = new Vector(A.length);
			for(int i = 0; i < A.length; i++)
				C.components[i] = A.components[i] + B.components[i];
		}
		return C;
	}

	/**
	 * 
	 * @param B
	 * @return
	 */
	public Vector minus(Vector B) {
		Vector A = this;
		Vector C = null;

		if(!A.isSameDimension(B)) {
			throw new VectorException("The dimensions are not equal!");
		}else {
			C = new Vector(A.length);
			for(int i = 0; i < A.length; i++)
				C.components[i] = A.components[i] - B.components[i];
		}
		return C;
	}

	public Vector times(double factor) {
		Vector c = new Vector(this.length);
		for (int i = 0; i < c.length; i++)
			c.components[i] *= factor;
		return c;
	}


	/**
	 * Create a distinct/defensive copy so that client can't alter our copy of data[].
	 * 
	 * @return A defensive/distinct copy.
	 */
	public Vector deepCopy() {
		Vector vectorToCopyFrom = this;
		Vector duplicated = new Vector(vectorToCopyFrom.length);

		for(int i = 0; i < duplicated.length; i++) {
			duplicated.components[i] = vectorToCopyFrom.components[i];
		}

		return duplicated;
	}

	/**
	 * Create a copy of this vector.
	 * 
	 * @return A copy of this vector.
	 */
	public Vector copy() {
		Vector copy = new Vector(this);
		return copy;
	}

	/**
	 * 
	 * @param vector
	 * @return
	 */
	public boolean isSameDimension(Vector vector) {
		return (this.length == vector.length);
	}

	public String toString() {
		String text = "";
		for(double number : components) {
			text += number + "\n";
		}
		return text;
	}
}
