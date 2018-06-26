#include <string>
#include "QAP.h"
#include "Matrix.h"
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <random>

QAP::QAP(){
}

// Cargar problema
void QAP::load(string fichero){

  ifstream fe(fichero) ;

  fe >> size ;

  F.resize(size, size);
  D.resize(size, size);

  for (int i = 0 ; i < size ; i++) {
    for (int j = 0 ; j < size ; j++) {
      fe >> F[i][j] ;
    }
  }

  for (int i = 0 ; i < size ; i++) {
    for (int j = 0 ; j < size ; j++) {
      fe >> D[i][j] ;
    }
  }

  fe.close();
}

// Calcula fitness de una solución
int QAP::Coste(vector<int> S){
  int coste=0;
  Matrix Mp(size, size);

  for(int i=0; i<S.size(); i++)
    Mp[i][S[i]] = 1;

  Matrix MpT = Mp.t();
  Matrix MatrixAux = Mp*D*MpT;
  coste = F.product_sum(MatrixAux);
  return coste;
}

Matrix QAP::getF(){
  return F;
}

Matrix QAP::getD(){
  return D;
}

vector<int> QAP::solucion_Greedy(){

  int sum;
  // Par valor-posicion
  pair<int,int> par;
  vector< pair<int, int> > f;
  vector< pair<int, int> > d;

  for(int i=0; i<size; i++){
    sum=0;
    for(int j=0; j<size; j++)
      sum += F[i][j];
    par.first=sum; par.second=i;
    f.push_back(par);
  }

  for(int i=0; i<size; i++){
    sum=0;
    for(int j=0; j<size; j++)
      sum += D[i][j];
    par.first=sum; par.second=i;
    d.push_back(par);
  }

  sort(f.begin(),f.end());
  sort(d.begin(),d.end());

  vector<int> resultado;
  resultado.resize(size);

  int tam = size-1;
  for(int i=0; i<size; i++){
    resultado[f[tam].second] = d[i].second;
    tam--;
  }

  return resultado;
}

// Diferencia entre 2 soluciones factorizado
int QAP::BL_diferenciaCoste(const vector<int> &v ,int r, int s){

  int dif=0;

  for(int k=0; k < v.size(); k++){

    if( k != r && k != s){
      dif += ( F[r][k] * (D[v[s]][v[k]] - D[v[r]][v[k]]) )+( F[s][k] * (D[v[r]][v[k]] - D[v[s]][v[k]]) )+
            ( F[k][r] * (D[v[k]][v[s]] - D[v[k]][v[r]]) )+( F[k][s] * (D[v[k]][v[r]] - D[v[k]][v[s]]) );

    }
  }

  return dif;
}

int myrandom (int i) { return std::rand()%i;}

vector<int> QAP::solucion_BL(int seed){

  vector<int> v;

  for(int i=0; i<size; i++)
    v.push_back(i);

  if(seed == 0)
    random_shuffle(v.begin(), v.end());
  else{
    srand(seed);
    random_shuffle(v.begin(), v.end(),myrandom);
  }

  return LocalSearch(v,50000);
}

vector<int> QAP::LocalSearch(vector<int> v, int maxIt){

  vector<bool> DBL (size,false);
  int i = 0;
  vector<int> new_result;
  bool improve = true ;

  int max_iter = maxIt;
  while (i < max_iter && improve) {
    new_result = jump (DBL,  v);
    improve = (new_result != v);
    v = new_result;
    i++;
  }
  return v;
}

vector<int> QAP::jump (vector<bool> & DBL, vector<int> result) {
  for (int i = 0 ; i < size ; i++){
    // Si ha sido probado antes nos lo saltamos
    if (DBL[i]) continue;
    for (int j = 0 ; j < size ; j++){
      // Comprobamos intercambios consigo mismo
      if (j == i) continue ;
      // Intercambiamos los valores de la solucion
      // Si produce una mejora
      //cout << factorize_function (result, i, j) << endl ;
      if (BL_diferenciaCoste(result, i, j) < 0){
        // Reactivamos el candidato y devolvemos la nueva solucion
        swap (result[i], result[j]);
        DBL[j] = false;
        return result;
      }
    }
    // Unidad que no produce mejoras
    DBL[i] = true ;
  }
  return result;
}

// GENÉTICOS
void printPoblacion(int size, const vector< pair< vector<int>, int>  > &poblacion){

  for(int i=0; i < poblacion.size(); i++){
      cout<<"POS: "<<i<<" ->";
      for(int j=0; j<size; j++)
        cout<<(poblacion[i].first)[j]<<", ";

      cout<<"\nCOSTE: "<<poblacion[i].second<<endl;

      cout<<"-----------------"<<endl;
    }
}

pair< vector<int>, int> QAP::solucion_GeneticoAGG(int seed, int cross, int meme, int tam){

  int tam_poblacion = tam;
  int evaluaciones = 0;
  int max_fitness = 50000;
  int generacion = 1;

  srand(seed);
  pair< vector<int>, int> solucion;
  int bestParentPosition, worstChildPosition;
  vector< pair< vector<int>, int>  > poblacion, hijos; // solucion-coste
  vector<int> indices_padres;

  generarPoblacion(poblacion, tam_poblacion);

  evaluaciones += tam_poblacion;

  while(evaluaciones < max_fitness){


    if( generacion == 10 && meme != 0){

      switch(meme){

        case 1:
          AM1(poblacion);
          evaluaciones += tam_poblacion;
          break;
        case 2:
          AM2(poblacion);
          evaluaciones += tam_poblacion * 0.1;
          break;
        case 3:
          AM3(poblacion);
          evaluaciones += tam_poblacion * 0.1;
          break;

      }
      generacion = 1;
    }

    // Selecciono los padres
    indices_padres = seleccionarPadres(poblacion, poblacion.size()); //De cada dos padres saldrá un crío

    // Genero hijos
    hijos = cruzarPoblacion(poblacion, indices_padres, 0.7, cross);

    evaluaciones += 0.7*indices_padres.size();

    // Genero mutaciones
    Mutate(hijos, 0.001, false);

    // Mantengo elitismo
    bestParentPosition = searchBestSolution(poblacion);
    worstChildPosition = searchWorstSolution(hijos);

    hijos[worstChildPosition] = poblacion[bestParentPosition];

    // Reemplazo
    poblacion = hijos;

    generacion++;
  }

  bestParentPosition = searchBestSolution(poblacion);
  solucion = poblacion[bestParentPosition];
  return solucion;
}


pair< vector<int>, int> QAP::solucion_GeneticoAGE(int seed, int cross, int tam){

  int tam_poblacion = tam;
  int evaluaciones = 0;
  int max_fitness = 50000;

  srand(seed);
  pair< vector<int>, int> solucion;
  int bestParentPosition;
  vector< pair< vector<int>, int>  > poblacion, hijos; // solucion-coste
  vector<int> indices_padres;

  generarPoblacion(poblacion, tam_poblacion);
  evaluaciones += tam_poblacion;

  int iteracion = 0;
  int generacion = 1;
  while(evaluaciones < max_fitness){


    // Selecciono los padres
    indices_padres = seleccionarPadres(poblacion, 2); //De cada dos padres saldrá un crío

    // Genero hijos
    hijos = cruzarPoblacion(poblacion, indices_padres, 1, cross);
    evaluaciones += indices_padres.size();

    // Genero mutaciones
    Mutate(hijos, 0.001, true);

    // Reemplazo
    reemplazoAGE(poblacion, hijos);

  }

  bestParentPosition = searchBestSolution(poblacion);
  solucion = poblacion[bestParentPosition];
  return solucion;
}

void QAP::generarPoblacion(vector< pair< vector<int>, int> > &poblacion, int tam_poblacion){

  for(int i = 0; i < tam_poblacion; i++){
    vector<int> v;

    for(int j=0; j<size; j++)
      v.push_back(j);

    random_shuffle(v.begin(), v.end(), myrandom);

    pair< vector<int>, int > par = make_pair(v, Coste(v));
    poblacion.push_back(par);
  }
}

void QAP::calcularFitness(vector< pair< vector<int>, int> > &poblacion){

  for(int i=0; i<poblacion.size(); i++)
    poblacion[i].second = Coste(poblacion[i].first);

}


vector<int> QAP::seleccionarPadres(const vector< pair< vector<int>, int> > &poblacion, int num_padres){
  vector<int> padres;
  int p1,p2,tam_poblacion=poblacion.size();

  for(int i=0; i < num_padres; i++){

    p1 = rand() % tam_poblacion;
    p2 = rand() % tam_poblacion;

    if(poblacion[p1].second < poblacion[p2].second)
      padres.push_back(p1);
    else
      padres.push_back(p2);
  }

  return padres;
}

vector< pair< vector<int>, int>  > QAP::cruzarPoblacion(const vector< pair< vector<int>, int>  > &poblacion, vector<int> indices_padres, double probabilidad, int cross){
  vector< pair< vector<int>, int>  > hijos;
  pair<  pair<vector<int>, int>, pair<vector<int>, int> > hijosOrden;
  pair< vector<int>, int> hijo1, hijo2;

  int total_parejas = indices_padres.size()/2;
  int maxIndex = 2 *( probabilidad * total_parejas );

  for(int i = 0; i < maxIndex; i+=2){

    if(cross == 0){
      hijo1 = crucePosicion(poblacion[indices_padres[i]].first, poblacion[indices_padres[i+1]].first);
      hijo2 = crucePosicion(poblacion[indices_padres[i]].first, poblacion[indices_padres[i+1]].first);
    }else{
      hijosOrden = cruceOrden(poblacion[indices_padres[i]].first, poblacion[indices_padres[i+1]].first);
      hijo1 = hijosOrden.first;
      hijo2 = hijosOrden.second;
    }
    hijos.push_back(hijo1);
    hijos.push_back(hijo2);
  }

  for(int i=maxIndex+1; i < indices_padres.size(); i++){
    hijos.push_back(poblacion[indices_padres[i]]);
  }

  return hijos;
}

pair<  vector<int>, int > QAP::crucePosicion(vector<int> padre1, vector<int> padre2){
  vector<int> hijo;
  vector<int> valores_distintos;
  vector<bool> posiciones_libres(size,true);
  hijo.resize(size);

  for(int i=0; i < size; i++){

    if(padre1[i] == padre2[i]){
      hijo[i] = padre1[i];
      posiciones_libres[i] = false;
    }else{
      valores_distintos.push_back(padre1[i]);
    }
  }

  random_shuffle(valores_distintos.begin(), valores_distintos.end(), myrandom);

  for(int i = 0, contador_valores_distintos=0; i < size; i++){

    if(posiciones_libres[i]){
      hijo[i] = valores_distintos[contador_valores_distintos];
      contador_valores_distintos++;
    }
  }

  return make_pair(hijo, Coste(hijo));
}

// Cambiar para que devuelva el coste y eso
pair<  pair<vector<int>, int>, pair<vector<int>, int> > QAP::cruceOrden(vector<int> padre1, vector<int> padre2){

  vector<int> child1, child2;
  vector<int> centerValuesC1(size,-1), centerValuesC2(size,-1);
  child1.resize(size);
  child2.resize(size);

  // Generar puntos de corte (Primer/Ultimo elemento del centro)
  int random = rand()%(size*size);
  int firstMiddle = random / size;
  int lastMiddle = random % size;

  if(firstMiddle>lastMiddle) swap(firstMiddle,lastMiddle);

  // Relleno vectores
  // Escribo los del centro
  for(int i=firstMiddle; i<=lastMiddle; i++){
    child1[i] = padre2[i];
    centerValuesC1[ padre2[i] ] = padre1[i]; // Ese valor esta pillado por el centro
    child2[i] = padre1[i];
    centerValuesC2[ padre1[i] ] = padre2[i];
  }

  int aux;
  // Escribo extremos salvando incongruencias
  for(int i=0; i<size; i++){

    if(i == firstMiddle){
      i = lastMiddle;
      continue;
    }

    aux = padre1[i];
    while(centerValuesC1[aux] != -1){
      aux = centerValuesC1[aux];
    }
    child1[i] = aux;

    aux = padre2[i];
    while(centerValuesC2[aux] != -1){
      aux = centerValuesC2[aux];
    }
    child2[i] = aux;
  }

  return make_pair( make_pair(child1, Coste(child1)), make_pair(child2, Coste(child2)));
}

int QAP::searchBestSolution(const vector< pair< vector<int>, int>  > &poblacion){
  int indexOfBest, bestCost = std::numeric_limits<int>::max();

  for(int i=0; i < poblacion.size(); i++){

    if(poblacion[i].second < bestCost){
      bestCost = poblacion[i].second;
      indexOfBest = i;
    }
  }
  return indexOfBest;
}

int QAP::searchWorstSolution(const vector< pair< vector<int>, int>  > &poblacion){
  int indexOfWorst, worstCost = -1;

  for(int i=0; i < poblacion.size(); i++){

    if(poblacion[i].second > worstCost){
      worstCost = poblacion[i].second;
      indexOfWorst = i;
    }
  }
  return indexOfWorst;
}

pair<int,int> QAP::search2WorstPositions(const vector< pair< vector<int>, int>  > &poblacion){
  int indexOfWorst, worstCost = -1;
  int indexOfSecondWorst, secondWorstCost = -1;

  for(int i=0; i < poblacion.size(); i++){

    if(poblacion[i].second > secondWorstCost){

      if( poblacion[i].second > worstCost){

        secondWorstCost = worstCost;
        indexOfSecondWorst = indexOfWorst;

        worstCost = poblacion[i].second;
        indexOfWorst = i;

      }else{
        secondWorstCost = poblacion[i].second;
        indexOfSecondWorst = i;
      }
    }
  }
  return make_pair(indexOfWorst, indexOfSecondWorst);
}

void QAP::Mutate( vector< pair< vector<int>, int>  > &poblacion, double Pm, bool AGE){

  int numberOfChromosomes = poblacion.size();
  int numberOfGenes = size;

  int numberOfMutations;
  int differenceFitness, aux;
  int randomChromosome, randomGene, i, j;

  int random, P;
  if( AGE ){
    P = (1/Pm);
    numberOfMutations = 0;
    for(int i=0; i<poblacion.size(); i++){
      random = rand() % P;
      if(random < size)
        numberOfMutations++;
    }
  }else{

    numberOfMutations = numberOfChromosomes*numberOfGenes*Pm;
  }

  for(int i=0; i<numberOfMutations; i++){

    randomChromosome = rand()%numberOfChromosomes;

    do{
      randomGene = rand()%(numberOfGenes*numberOfGenes);
      i = randomGene / numberOfGenes;
      j = randomGene % numberOfGenes;
    }while(j == i); // Probabilidad de que sean iguales (1/size)

    differenceFitness = BL_diferenciaCoste(poblacion[randomChromosome].first,i,j);
    aux = (poblacion[randomChromosome].first)[i];
    (poblacion[randomChromosome].first)[i] = (poblacion[randomChromosome].first)[j];
    (poblacion[randomChromosome].first)[j] = aux;

    poblacion[randomChromosome].second += differenceFitness;
  }

}

void QAP::reemplazoAGE(vector< pair< vector<int>, int>  > &poblacion, vector< pair< vector<int>, int>  > hijos){

  pair<int,int> worst2Index;

  // Hijos[1] el mejor
  if(hijos[0].second < hijos[1].second)
    swap(hijos[0], hijos[1]);
  bool firstMatch = false, secondMatch = false;
  worst2Index = search2WorstPositions(poblacion); //WORST/SECONDWORST
  // worst <-> secondBest
  if( poblacion[worst2Index.first].second > hijos[0].second)
    firstMatch = true;
  // secondWorst <-> Best
  if( poblacion[worst2Index.second].second > hijos[1].second)
    secondMatch = true;

  if( firstMatch && secondMatch){ // Si se cumplen ambos
    poblacion[worst2Index.first] = hijos[0];
    poblacion[worst2Index.second] = hijos[1];
  }else if( !firstMatch && !secondMatch){ // si no se ha cumplido ninguna, compruebo si el mejor tiene menos coste que el peor
    if(poblacion[worst2Index.first].second > hijos[1].second)
      poblacion[worst2Index.first] = hijos[1];
  }else
    poblacion[worst2Index.first] = hijos[1];

}

void QAP::AM1(vector< pair< vector<int>, int>  > &poblacion){

  vector<int> new_solution;
  for(int i=0; i<poblacion.size(); i++){
    new_solution = LocalSearch(poblacion[i].first, 400);
    poblacion[i].first = new_solution;
    poblacion[i].second = Coste(new_solution);
  }
}

void QAP::AM2(vector< pair< vector<int>, int>  > &poblacion){
  double Pls = 0.1;
  vector<int> indices;
  for(int i=0; i < poblacion.size(); i++)
    indices.push_back(i);

  random_shuffle(indices.begin(), indices.end(), myrandom);

  int maxIndex = poblacion.size()*Pls;

  vector<int> new_solution;
  for(int i=0; i<maxIndex; i++){
    new_solution = LocalSearch(poblacion[indices[i]].first, 400);
    poblacion[indices[i]].first = new_solution;
    poblacion[indices[i]].second = Coste(new_solution);
  }

}

bool cmp(pair<int, int> i, pair<int, int> j){
  return( i.second < j.second);
}
void QAP::AM3(vector< pair< vector<int>, int>  > &poblacion){

  double Pls = 0.1;
  vector< pair<int, int> > indices_coste;
  pair<int, int> par;
  for(int i=0; i < poblacion.size(); i++){
    par.first = i;
    par.second = poblacion[i].second;

    indices_coste.push_back(par);
  }

  sort(indices_coste.begin(), indices_coste.end(), cmp);

  int maxIndex = poblacion.size()*Pls;

  vector<int> new_solution;
  for(int i=0; i<maxIndex; i++){
    new_solution = LocalSearch(poblacion[indices_coste[i].first].first, 400);
    poblacion[indices_coste[i].first].first = new_solution;
    poblacion[indices_coste[i].first].second = Coste(new_solution);
  }

}

//========================
//***** PRACTICA 3 *******
//========================

vector<int> QAP::randomSolution(){
  vector<int> v;

  for(int i=0; i<size; i++)
    v.push_back(i);

  //srand(seed);
  random_shuffle(v.begin(), v.end(),myrandom);

  return v;
}

pair<int,int> QAP::returnRandomPositions(int tam){
  int i,j,random;
  do{
    random = rand()%(tam*tam);
    i = random / tam;
    j = random % tam;
  }while(j == i);
  return make_pair(i,j);
}

vector<int> QAP::swap2Positions(vector<int> s, int i, int j){
  vector<int> r = s;
  swap(r[i], r[j]);
  return r;
}

double calculoBeta(int T_0, int T_f, int M){
  return((T_0-T_f)/(M*T_0*T_f));
}

pair< vector<int>, int> QAP::ES(int seed){

  default_random_engine re;
  uniform_real_distribution<double> unif(0,1);
  srand(seed);

  pair< vector<int>, int> s;
  pair< vector<int>, int> bestSolution;
  pair< vector<int>, int> r;
  pair<int, int> randomPositions;

  int maxNeighbors = 10*size;
  int MaxSuccesses = size;
  int maxIter = 50000/(maxNeighbors);

  s.first = randomSolution();
  s.second = Coste(s.first);
  bestSolution = s;

  int costDifference;
  double T_0 = (s.second*0.3)/(-log(0.3));
  double T_k = T_0;
  double T_f = 0.001;
  int it = 0, neighbors, successes;
  do{
    successes=0;
    neighbors=0;
    do{
      randomPositions = returnRandomPositions(size);
      r.first = swap2Positions(s.first, randomPositions.first, randomPositions.second);
      costDifference =  BL_diferenciaCoste(s.first, randomPositions.first, randomPositions.second);
      r.second = s.second + costDifference;
      if(costDifference < 0 || ( unif(re)<=(exp(-(costDifference)/T_k)) )){
        successes++;
        s = r;
        if(s.second < bestSolution.second){
          bestSolution = s;
        }
      }
      neighbors++;
    }while(neighbors < maxNeighbors && successes < MaxSuccesses);
    it++;
    T_k = T_k/(1+(calculoBeta(T_0, T_f, maxIter)*T_k));
    //T_k = T_k*0.95;
  }while(it < maxIter && successes != 0);

  return bestSolution;
}

pair< vector<int>, int> QAP::BMB(int seed){

  srand(seed);

  int maxIter = 25;
  vector<int> v, bestSolution;
  int cost, bestCost = std::numeric_limits<int>::max();

  for(int i=0; i < maxIter; i++){
    v = randomSolution();
    v = LocalSearch(v, 50000);
    cost = Coste(v);
    if(cost < bestCost){
      bestSolution = v;
      bestCost = cost;
    }
  }
  return make_pair(bestSolution, bestCost);
}

bool lowToHigh (pair<int,int> i, pair<int,int> j) {
  return (i.first < j.first);
}
bool HighToLow (pair<int,int> i, pair<int,int> j) {
  return (i.first > j.first);
}
bool segunCosto(pair< pair<int,int>, int> i, pair< pair<int,int>, int> j){
  return (i.second < j.second);
}

// Genero Soluciones Greedy Aleatorizadas
vector<int> QAP::greedyGRASP(){

  int sum;
  // Par valor-posicion
  pair<int,int> par;
  // Creo solucion
  vector<int> solucion;
  solucion.resize(size);
  // Creo las listas LCu y LCl
  vector< pair<int, int> > f; // LCu
  vector< pair<int, int> > d; // LCf


//  ===========================
//  ======== ETAPA 1 ==========
//  ===========================

//  ========= Inicio LCs ======

  // F = LCunidades
  for(int i=0; i<size; i++){
    sum=0;
    for(int j=0; j<size; j++)
      sum += F[i][j];
    par.first=sum; par.second=i;
    f.push_back(par);
  }
  sort(f.begin(),f.end(), HighToLow);
  // D = LClocalizaciones
  for(int i=0; i<size; i++){
    sum=0;
    for(int j=0; j<size; j++)
      sum += D[i][j];
    par.first=sum; par.second=i;
    d.push_back(par);
  }
  sort(d.begin(),d.end(), lowToHigh);

  //  ========= Inicio LRCs ========

  vector< pair<int, int> > lrcu; // F
  vector< pair<int, int> > lrcl; // D
  double alfa = 0.3;

  // Calculo umbrales
  int umbral_u = f[0].first - alfa*(f[0].first - f[size-1].first);
  int umbral_l = d[0].first + alfa*(d[size-1].first - d[0].first);

/*
  cout<<"\nF: ";
  for(int i=0; i<f.size(); i++)
    cout<<f[i].first<<" ";

  cout<<"\nD: ";
  for(int i=0; i<d.size(); i++)
    cout<<d[i].first<<" ";
*/

  // Inicio listas

  for(int i=0; i < size; i++){
    if(f[i].first > umbral_u)
      lrcu.push_back(f[i]);

    if(d[i].first < umbral_l)
      lrcl.push_back(d[i]);
  }

  // Por si solo se inserta uno
  if(lrcu.size() < 2 && size > 1){
    if(lrcu.empty()){
      lrcu.push_back(f[0]);
    }
    lrcu.push_back(f[1]);
  }if(lrcl.size() < 2 && size > 1){
    if(lrcl.empty()){
      lrcl.push_back(d[0]);
    }
    lrcl.push_back(d[1]);
}

/*
  cout<<"\nLRCu: ";
  for(int i=0; i<lrcu.size(); i++)
    cout<<lrcu[i].first<<" ";

  for(int i=0; i<lrcl.size(); i++)
    cout<<lrcl[i].first<<" ";
*/
//  cout<<endl;

  // Escoger dos posiciones aleatorias distintas
  pair<int, int> lrcu_positions = returnRandomPositions(lrcu.size());
  pair<int, int> lrcl_positions = returnRandomPositions(lrcl.size());

  // Inserto en la solución las parejas obtenidas (1-1,2-2)
  solucion[lrcu[lrcu_positions.first].second] = lrcl[lrcl_positions.first].second;
  solucion[lrcu[lrcu_positions.second].second] = lrcl[lrcl_positions.second].second;

  // Elimino las unidades - localizaciones usadas
  if(lrcu_positions.first < lrcu_positions.second)
    swap(lrcu_positions.first , lrcu_positions.second);
  if(lrcl_positions.first < lrcl_positions.second)
    swap(lrcl_positions.first , lrcl_positions.second);

  f.erase( f.begin() + lrcu_positions.first);
  f.erase( f.begin() + lrcu_positions.second);
  d.erase( d.begin() + lrcl_positions.first);
  d.erase( d.begin() + lrcl_positions.second);

  // Incluyo ambas asignaciones en asignaciones_realizadas
  vector< pair<int,int> > asignaciones_realizadas;
  par = make_pair(lrcu[lrcu_positions.first].second, lrcl[lrcl_positions.first].second);
  asignaciones_realizadas.push_back(par);
  par = make_pair(lrcu[lrcu_positions.second].second, lrcl[lrcl_positions.second].second);
  asignaciones_realizadas.push_back(par);

//  ===========================
//  ======== ETAPA 2 ==========
//  ===========================

  // Almaceno unidades/localizaciones sin usar
  // size-2 porque ya he quitado 2 elementos
  vector<int> unidadesSinUsar;
  vector<int> localizacionesSinUsar;
  for(int i=0; i<size-2; i++){
    unidadesSinUsar.push_back(f[i].second);
    localizacionesSinUsar.push_back(d[i].second);
  }

  int random;
  vector< pair< pair<int,int>, int> > LC;
  vector< pair< pair<int,int>, int> > LRC;

  while(!unidadesSinUsar.empty()){

    LC.clear();
    LRC.clear();
    // Creo pairs unidad-localizacion
    for(int i=0; i<unidadesSinUsar.size(); i++)
      for(int j=0; j<localizacionesSinUsar.size(); j++){
        par = make_pair(unidadesSinUsar[i], localizacionesSinUsar[j]);
        LC.push_back(make_pair(par, costeGRASP(par, asignaciones_realizadas)));
      }

/*
    cout<<"\nUNIDADES SIN USAR: ";
    for(int i=0; i<unidadesSinUsar.size(); i++)
      cout<<unidadesSinUsar[i]<<" ";
    cout<<"\nLOCALIZ SIN USAR: ";
    for(int i=0; i<localizacionesSinUsar.size(); i++)
      cout<<localizacionesSinUsar[i]<<" ";
*/

    sort(LC.begin(), LC.end(), segunCosto);

    int umbral = LC[0].second + alfa*(LC[LC.size()-1].second - LC[0].second);

    if(localizacionesSinUsar.size() > 1){
      for(int i=0; i < LC.size(); i++){
        if(LC[i].second <= umbral)
          LRC.push_back(LC[i]);
      }
      random = rand()%LRC.size();
      solucion[LRC[random].first.first] = LRC[random].first.second;
      asignaciones_realizadas.push_back(LRC[random].first);
    }else{
      random = 0;
      LRC.push_back(LC[random]);
      solucion[LRC[random].first.first] = LRC[random].first.second;
      asignaciones_realizadas.push_back(LRC[random].first);
    }

    EliminarPosicion(unidadesSinUsar, LRC[random].first.first);
    EliminarPosicion(localizacionesSinUsar, LRC[random].first.second);

  }
  return solucion;
}

int QAP::costeGRASP(pair<int,int> par, vector< pair<int,int> > asignaciones_realizadas){

  int coste=0;
  for(int i=0; i<asignaciones_realizadas.size(); i++){
    coste += F[par.first][asignaciones_realizadas[i].first] *
          D[par.second][asignaciones_realizadas[i].second];
  }
  return coste;
}

void QAP::EliminarPosicion( vector<int> &array, int pos){

  for(int i=0; i<array.size(); i++){
    if(array[i] == pos){
      array.erase(array.begin()+i);
      break;
    }
  }
}

pair< vector<int>, int> QAP::GRASP(int seed){
  srand(seed);
  pair< vector<int>, int> bestSolution, solution;
  bestSolution.second = std::numeric_limits<int>::max();

  for(int i=0; i < 25; i++){
    solution.first = greedyGRASP();
    solution.first = LocalSearch(solution.first, 50000);
    solution.second = Coste(solution.first);
    if(solution.second < bestSolution.second)
      bestSolution = solution;
  }
  return bestSolution;
}

pair< vector<int>, int> QAP::ILS(int seed){

  srand(seed);
  int maxIter = 24;
  pair<vector<int>, int> bestSolution;
  bestSolution.first = randomSolution();
  bestSolution.second = Coste(bestSolution.first);

  pair<vector<int>, int> solution;
  solution.first = LocalSearch(bestSolution.first, 50000);
  solution.second = Coste(solution.first);

  int it = 0;
  while(it < maxIter){

    if(solution.second < bestSolution.second){
      bestSolution = solution;
    }
    solution = bestSolution;

    solution.first = mutacionILS(solution.first);
    solution.first = LocalSearch(solution.first, 50000);
    solution.second = Coste(solution.first);

    it++;
  }
  return solution;
}

vector<int> QAP::mutacionILS(vector<int> sol){
  int t = size/4;
  int random = rand() % (size-t);

  random_shuffle(sol.begin()+random, sol.begin()+random+t);
  return sol;
}

pair< vector<int>, int> QAP::ILS_ES(int seed){

  srand(seed);
  int maxIter = 24;
  pair<vector<int>, int> bestSolution;
  bestSolution.first = randomSolution();
  bestSolution.second = std::numeric_limits<int>::max();

  pair<vector<int>, int> solution;
  solution = ES(bestSolution.first);

  int it = 0;
  while(it < maxIter){

    if(solution.second < bestSolution.second){
      bestSolution = solution;
    }
    solution = bestSolution;
    solution.first = mutacionILS(solution.first);
    solution = ES(solution.first);

    it++;
  }
  return solution;
}

pair< vector<int>, int> QAP::ES(vector<int> sol){

  default_random_engine re;
  uniform_real_distribution<double> unif(0,1);

  pair< vector<int>, int> s;
  pair< vector<int>, int> bestSolution;
  pair< vector<int>, int> r;
  pair<int, int> randomPositions;

  int maxNeighbors = 10*size;
  int MaxSuccesses = size;
  int maxIter = 50000/(maxNeighbors);

  s.first = sol;
  s.second = Coste(s.first);
  bestSolution = s;

  int costDifference;
  double T_0 = (s.second*0.3)/(-log(0.3));
  double T_k = T_0;
  double T_f = 0.001;
  int it = 0, neighbors, successes;
  do{
    successes=0;
    neighbors=0;
    do{
      randomPositions = returnRandomPositions(size);
      r.first = swap2Positions(s.first, randomPositions.first, randomPositions.second);
      costDifference =  BL_diferenciaCoste(s.first, randomPositions.first, randomPositions.second);
      r.second = s.second + costDifference;
      if(costDifference < 0 || ( unif(re)<=(exp(-(costDifference)/T_k)) )){
        successes++;
        s = r;
        if(s.second < bestSolution.second){
          bestSolution = s;
        }
      }
      neighbors++;
    }while(neighbors < maxNeighbors && successes < MaxSuccesses);
    it++;
    T_k = T_k/(1+(calculoBeta(T_0, T_f, maxIter)*T_k));
  }while(it < maxIter && successes != 0);

  return bestSolution;
}
