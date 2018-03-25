#pragma once
#include <stdio.h>
#include <iostream>
class Matrix
{
	friend std::ostream& operator<<(std::ostream& os, const Matrix& matr);
	friend Matrix operator*(double a, const Matrix& A);
	friend Matrix operator*(const Matrix& A, double a);

public:
	explicit Matrix(int row = 0, int col = 0, double initValue = 0);
	Matrix(const Matrix&);
	~Matrix();
	int rows() const { return mRows; } ;
	int cols() const { return mCols; } ;
	void print();
	Matrix multi(const Matrix& B);
	void LU(Matrix& L, Matrix& U);
	Matrix T();
	Matrix inv();//LU·Ö½âÇóÄæ
	double det();
	Matrix dot(const Matrix& B);

	Matrix& operator=(const Matrix& B);
	Matrix operator*(const Matrix& B);
	double* operator[](int i);
	Matrix operator+(const Matrix& B);
	Matrix operator-(const Matrix& B);
	bool operator==(const Matrix& B);
	bool operator!=(const Matrix& B);

private:
	double** mat;
	int mRows;
	int mCols;
	Matrix& dnTriInv();
	Matrix& upTriInv();
};

std::ostream& operator<<(std::ostream& os, const Matrix& matr);
Matrix operator*(double a, const Matrix& A);
Matrix operator*(const Matrix& A, double a);

