package edu.scripps.yates.pcq.compare.model;

public class TTestMatrix {
	private final MyTTest[][] data;

	public TTestMatrix(int numrows, int numcolumns) {
		data = new MyTTest[numrows][numcolumns];
	}

	public MyTTest get(int row, int col) {
		return data[row][col];
	}

	public void set(int row, int col, MyTTest ttest) {
		data[row][col] = ttest;
	}

	public int nrows() {
		return data.length;
	}

	public int ncols() {
		return data[0].length;
	}
}
