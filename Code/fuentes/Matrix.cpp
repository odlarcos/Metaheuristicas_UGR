
#include "Matrix.h"

Matrix::Matrix(){
  nfil=0;
  ncol=0;
}

Matrix::Matrix(int f, int c){
  resize(f,c);
}

void Matrix::rellenar(){
  for(int i=0; i<nfil; i++)
    for(int j=0; j<ncol; j++)
      m[i][j]=i+j;
}

void Matrix::resize(int f, int c){
  nfil = f;
  ncol = c;
  m.resize(f);
  for(int i=0; i<f; i++)
    m[i].resize(c,0);
}

int Matrix::getF(){
  return nfil;
}

int Matrix::getC(){
  return ncol;
}

vector<int> & Matrix::operator[](int a){
  return m[a];
}

vector<int> Matrix::operator[](int a) const{
  return m[a];
}

Matrix Matrix::operator*(const Matrix &b){
  if( this->ncol != b.nfil ){
    exit(-1);
    cerr<<" No se pueden multiplicar"<<endl;
  }else{
    Matrix r(this->nfil, b.ncol);

    for(int i=0; i < r.nfil; i++)
      for(int j=0; j < r.ncol; j++)
        for(int k=0; k < b.nfil; k++)
          r.m[i][j] += m[i][k] * b.m[k][j];

      return r;
  }
}

void Matrix::print(){
  for(int i=0; i<nfil; i++){
    for(int j=0; j<ncol; j++)
      cout << m[i][j]<< " ";
    cout<<endl;
  }
}

Matrix Matrix::t(){
  Matrix t(ncol,nfil);
  for(int i=0; i<nfil; i++)
    for(int j=0; j<ncol; j++){
      t.m[j][i]=m[i][j];
    }
    return t;
}

int Matrix::product_sum(Matrix b){

  if( ncol != b.ncol || nfil != b.nfil){
    exit(-1);
    cerr<<" No tienen las mismas dimensiones"<<endl;

  }else{
    int result = 0;

    for(int i=0; i < nfil; i++)
      for(int j=0; j < ncol; j++)
        result += m[i][j] * b.m[i][j];

    return result;
  }
}
