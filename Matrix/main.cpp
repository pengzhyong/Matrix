#include <iostream>
#include "matrix.h"
using namespace std;

int main()
{
	Matrix A(3, 4, 2);
	cout << A << endl;
	Matrix B(4, 3, 1);
	Matrix C = A*B;
	cout << C;
	Matrix D = C * 10;
	cout << D;
	cout << D[2][2] << endl;
	D[1][0] = -10;
	cout << D << endl;
	C = A.dot(B);
	cout << C << endl;
	system("pause");
	return 0;
}