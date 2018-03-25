#include "matrix.h"
#define EPS 0.000001
#define zero(x) x < EPS && x > -EPS 

Matrix::Matrix(int row, int col, double initValue) :mRows(row), mCols(col)
{
	mat = new double*[mRows];
	for (int i = 0; i < mRows; i++)
		mat[i] = new double[mCols];
	for (int i = 0; i < mRows; i++)
	{
		for (int j = 0; j < mCols; j++)
		{
			mat[i][j] = initValue;
		}
	}
}

Matrix::Matrix(const Matrix& B)
{
	mRows = B.rows();
	mCols = B.cols();
	mat = new double*[mRows];
	for (int i = 0; i < mRows; i++)
		mat[i] = new double[mCols];
	for (int i = 0; i < mRows; i++)
	{
		for (int j = 0; j < mCols; j++)
		{
			mat[i][j] = B.mat[i][j];
		}
	}
}

Matrix::~Matrix()
{
	for (int i = 0; i < mRows; i++)
		delete[] mat[i];
	delete[] mat;
}

void Matrix::print()
{
	for (int i = 0; i < mRows; i++)
	{
		for (int j = 0; j < mCols; j++)
		{
			printf("%f ", mat[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

Matrix Matrix::multi(const Matrix& B)
{
	Matrix C(mRows, B.cols());
	int rowsA = this->mRows;
	int colsA = this->mCols;
	int rowsB = B.rows();
	int colsB = B.cols();
	for (int ra = 0; ra < rowsA; ra++)
	{
		for (int ca = 0; ca < B.cols(); ca++)
		{
			C.mat[ra][ca] = 0;
			for (int index = 0; index < colsA; index++)
			{
				C.mat[ra][ca] += this->mat[ra][index] * B.mat[index][ca];
			}
		}
	}
	return C;
}

Matrix Matrix::T()
{
	Matrix tmpMat(mCols, mRows, 0);
	for (int r = 0; r < mCols; r++)
	{
		for (int c = 0; c < mRows; c++)
		{
			tmpMat.mat[r][c] = this->mat[c][r];
		}
	}
	return tmpMat;
}

double Matrix::det()
{
	Matrix triMat = *this;
	int signFlag = 1;
	for (int r = 0; r < mRows - 1; r++)
	{
		if (triMat.mat[r][r] == 0)//如果对角线上元素为0,则把非0元素列交换到此列
		{
			int c = r + 1;
			bool flag = false;//flag为false表示未找到不为0的列
			for (; c < mCols; c++)
			{
				if (!zero(triMat.mat[r][c]))
				{
					flag = true;
					break;
				}
			}
			if (!flag)
				return 0;
			//交换两列
			for (int row = 0; row < mRows; row++)
			{
				double tmp = triMat.mat[row][r];
				triMat.mat[row][r] = triMat.mat[row][c];
				triMat.mat[row][c] = tmp;
			}
			signFlag = -signFlag;
		}
		for (int row = r + 1; row < mRows; row++)
		{
			double div = triMat.mat[row][r] / triMat.mat[r][r];
			for (int col = r; col < mCols; col++)
			{
				triMat.mat[row][col] -= div * triMat.mat[r][col];
			}
		}

	}
	double deta = 1;
	for (int i = 0; i < mRows; i++)
	{
		deta *= triMat.mat[i][i];
	}
	return deta;
}

void Matrix::LU(Matrix& L, Matrix& U)
{
	if (mRows != mCols || mRows == 0 || mCols == 0)
	{
		printf("Unable to decomposition !\n");
		return;
	}
	if (L.mRows != mRows || L.mCols != mCols || U.mRows != mRows || U.mCols != mCols)
	{
		printf("Dimisioin not match !");
		return;
	}
	for (int i = 0; i < mRows; i++)
	{
		L.mat[i][0] = mat[i][0] / mat[0][0];
		U.mat[0][i] = mat[0][i];
		L.mat[i][i] = 1;
	}
	for (int pos = 1; pos < mRows; pos++)
	{
		//先计算U
		for (int c = pos; c < mCols; c++)
		{
			double sum = 0;
			for (int index = 0; index < pos; index++)
			{
				sum += L.mat[pos][index] * U.mat[index][c];
			}
			U.mat[pos][c] = mat[pos][c] - sum;
		}
		//再计算L
		if (pos < mRows)
		{
			for (int r = pos + 1; r < mRows; r++)
			{
				double sum = 0;
				for (int index = 0; index < pos; index++)
				{
					sum += L.mat[r][index] * U.mat[index][pos];
				}
				L.mat[r][pos] = (mat[r][pos] - sum) / U.mat[pos][pos];
			}
		}
		
	}

}

Matrix& Matrix::dnTriInv()
{
	double eps = 0.0000001;
    for(int i = 0; i < mRows; i++)
	{
		for (int j = i + 1; j < mCols; j++)
		{
			if ((mat[i][j]) < -eps || mat[i][j] > eps)
			{
				printf("This is not a requied Triangle Matrix !");
				return *this;
			}
		}
	}
	Matrix L = *this;
	Matrix prevMat = *this;
	for (int i = 0; i < L.rows(); i++)
		L.mat[i][i] = 1 / L.mat[i][i];
	for (int r = 1; r < L.rows(); r++)
	{
		for (int c = r - 1; c >=0; c--)
		{
			double sum = 0;
			for (int index = c; index < r; index++)
			{
				sum += prevMat.mat[r][index] * L.mat[index][c];
			}
			L.mat[r][c] = -L.mat[r][r] * sum;
		}
	}
	*this = L;
	return *this;
}

Matrix& Matrix::upTriInv()
{
	*this =  this->T();
	*this = dnTriInv();
	*this = this->T();
	return *this;
}

Matrix Matrix::inv()
{
	double deta = this->det();
	if (zero(deta))//如果行列式为0，说明矩阵为奇异矩阵，无法求逆
	{
		printf("Singular Matrix, unable to get inverse !");
		return *this;
	}
	Matrix L(mRows, mCols);
	Matrix U(mRows, mCols);
	Matrix Multi(mRows, mCols);
	this->LU(L, U);
	L = L.dnTriInv();
	U = U.upTriInv();
	Multi = U.multi(L);
	
	return Multi;
}

Matrix& Matrix::operator=(const Matrix& B)
{
	if (B.rows() != mRows || B.cols() != mCols)
	{
		printf("Dimishion not match in operator '=' !");
		return *this;
	}
	for (int r = 0; r < mRows; r++)
	{
		for (int c = 0; c < mCols; c++)
		{
			this->mat[r][c] = B.mat[r][c];
		}
	}
	return *this;
}

std::ostream& operator<<(std::ostream& os, const Matrix& matr)
{
	for (int r = 0; r < matr.mRows; r++)
	{
		for (int c = 0; c < matr.mCols; c++)
		{
			os << matr.mat[r][c] << " ";
		}
		os << std::endl;
	}
	return os;
}
