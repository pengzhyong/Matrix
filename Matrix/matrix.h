#pragma once
#include <stdio.h>
#include <iostream>
class Matrix
{
	friend std::ostream& operator<<(std::ostream& os, const Matrix& matr);
public:
	explicit Matrix(int row = 0, int col = 0, double initValue = 0);
	Matrix(const Matrix&);
	~Matrix();
	int rows() const { return mRows; } ;
	int cols() const { return mCols; } ;
	void print();
	Matrix multi(const Matrix& B);
	void LU(Matrix& L, Matrix& U);
	Matrix& operator=(const Matrix& B);
	Matrix T();
	Matrix inv();
	double det();

private:
	double** mat;
	int mRows;
	int mCols;
	Matrix& dnTriInv();
	Matrix& upTriInv();
};

std::ostream& operator<<(std::ostream& os, const Matrix& matr);

