#ifndef MATRIX
#define MATRIX

#include <iostream>
#include <vector>

using namespace std;

class Matrix{

  private:
    vector< vector<int> > m;
    int nfil, ncol;

  public:

    Matrix();

    Matrix(int f, int c);

    void rellenar();

    void resize(int f, int c);

    int getF();

    int getC();

    vector<int> & operator[](int a);

    vector<int> operator[](int a) const;

    Matrix operator*(const Matrix &b);

    void print();

    Matrix t(); // Traspuesta

    int product_sum(Matrix b); // Operacion <M, b>

};

#endif
