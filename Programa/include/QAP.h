#ifndef _QAP
#define _QAP

#include "Matrix.h"
#include <iostream>
#include <string>

using namespace std;

class QAP{
  private:
    Matrix D;
    Matrix F;
    int size;
    //vector<int> S;
  public:
    QAP();
    void load(string fichero);
    int Coste(vector<int> S);
    vector<int> solucion_Greedy();
    vector<int> solucion_BL(int seed);
      vector<int> LocalSearch(vector<int> v, int maxIt);
      int BL_diferenciaCoste(const vector<int> &v ,int i, int j);
      vector<int> jump (vector<bool> & DBL, vector<int> result);
    pair< vector<int>, int> solucion_GeneticoAGG(int seed, int cross, int meme, int tam);
      void generarPoblacion(vector< pair< vector<int>, int> > &poblacion, int tam_poblacion);
      vector<int> seleccionarPadres(const vector< pair< vector<int>, int> > &poblacion, int num_padres);
      void calcularFitness(vector< pair< vector<int>, int> > &poblacion);
      vector< pair< vector<int>, int>  > cruzarPoblacion(const vector< pair< vector<int>, int>  > &poblacion, vector<int> indices_padres, double probabilidad, int cross);
      pair< vector<int>, int > crucePosicion(vector<int> padre1, vector<int> padre2);
      pair<  pair<vector<int>, int>, pair<vector<int>, int> > cruceOrden(vector<int> padre1, vector<int> padre2);
      int searchBestSolution(const vector< pair< vector<int>, int>  > &poblacion);
      int searchWorstSolution(const vector< pair< vector<int>, int>  > &poblacion);
      pair<int,int> search2WorstPositions(const vector< pair< vector<int>, int>  > &poblacion);
      void Mutate( vector< pair< vector<int>, int>  > &poblacion, double Pm, bool AGE);
      void reemplazoAGE(vector< pair< vector<int>, int>  > &poblacion, vector< pair< vector<int>, int>  > hijos);
    pair< vector<int>, int> solucion_GeneticoAGE(int seed, int cross, int tam);

    void AM1(vector< pair< vector<int>, int>  > &poblacion);
    void AM2(vector< pair< vector<int>, int>  > &poblacion);
    void AM3(vector< pair< vector<int>, int>  > &poblacion);


    vector<int> randomSolution();
    pair< vector<int>, int> ES(int seed);
    pair<int,int> returnRandomPositions(int tam);
    vector<int> swap2Positions(vector<int> s, int i, int j);

    pair< vector<int>, int> BMB(int seed);

    pair< vector<int>, int> GRASP(int seed);
    vector<int> greedyGRASP();
    int costeGRASP(pair<int,int> par, vector< pair<int,int> > asignaciones_realizadas);
    void EliminarPosicion( vector<int> &array, int pos);

    pair< vector<int>, int> ILS(int seed);
    vector<int> mutacionILS(vector<int> sol);

    pair< vector<int>, int> ILS_ES(int seed);
    pair< vector<int>, int> ES(vector<int> sol);

    Matrix getF();
    Matrix getD();
};

#endif

/*
La busqueda local tiene que poder partir de una solucion inicial
Debes poder pasarle por parametro el numero de iteraciones 50000 predefinidas

mascara de bits vector<bool> inicializado a false
*/
