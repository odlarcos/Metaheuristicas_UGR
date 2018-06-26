
#include <iostream>
#include <fstream>
#include <string>
#include "Matrix.h"
#include "QAP.h"
#include <vector>
#include <ctime>

using namespace std;

int main(int argc, char *argv[]){

  unsigned t0, t1;
  double time;
  QAP problema;
  string ruta = "./Instancias_QAP/";
  int seed = 0;
  int cruce, meme;

  // Fijar semilla
  if(argc == 2)
    seed = atoi(argv[1]);

  vector<string> archivos;
  archivos.push_back("chr22a.dat");
  archivos.push_back("chr22b.dat");
  archivos.push_back("chr25a.dat");
  archivos.push_back("esc128.dat");
  archivos.push_back("had20.dat");
  archivos.push_back("lipa60b.dat");
  archivos.push_back("lipa80b.dat");
  archivos.push_back("nug28.dat");
  archivos.push_back("sko81.dat");
  archivos.push_back("sko90.dat");
  archivos.push_back("sko100a.dat");
  archivos.push_back("sko100f.dat");
  archivos.push_back("tai100a.dat");
  archivos.push_back("tai100b.dat");
  archivos.push_back("tai150b.dat");
  archivos.push_back("tai256c.dat");
  archivos.push_back("tho40.dat");
  archivos.push_back("tho150.dat");
  archivos.push_back("wil50.dat");
  archivos.push_back("wil100.dat");

  pair< vector<int>, int> solucion;
  vector<int> sol;
  string archivo;

  // Genético y memético
  cout<<"AGG-Posicion: "<<endl;
  for(int i=0; i<archivos.size(); i++){
    archivo = archivos[i];

    problema.load(ruta+archivo);

    cruce = 0;
    meme = 0;

    t0 = clock();
    solucion = problema.solucion_GeneticoAGG(seed, cruce, meme, 50);
    t1 = clock();
    time = (double(t1-t0)/CLOCKS_PER_SEC);

    cout<<i<<"\t";
    cout<<(solucion.second)<<"\t";
    cout<<time<<"\t";

    cout<<"Sol=";
    for(int i=0; i<(solucion.first).size(); i++){
      cout<<(solucion.first)[i]<<",";
    }
    cout<<endl;
  }

  cout<<"AGG-PMX: "<<endl;
  for(int i=0; i<archivos.size(); i++){
    archivo = archivos[i];

    problema.load(ruta+archivo);

    cruce = 1;
    meme = 0;

    t0 = clock();
    solucion = problema.solucion_GeneticoAGG(seed, cruce, meme, 50);
    t1 = clock();
    time = (double(t1-t0)/CLOCKS_PER_SEC);

    cout<<i<<"\t";
    cout<<(solucion.second)<<"\t";
    cout<<time<<"\t";

    cout<<"Sol=";
    for(int i=0; i<(solucion.first).size(); i++){
      cout<<(solucion.first)[i]<<",";
    }
    cout<<endl;
  }


cout<<"AGE-Posicion: "<<endl;
  for(int i=0; i<archivos.size(); i++){
    archivo = archivos[i];

    problema.load(ruta+archivo);

    cruce = 0;
    meme = 0;

    t0 = clock();
    solucion = problema.solucion_GeneticoAGE(seed, cruce, 50);
    t1 = clock();
    time = (double(t1-t0)/CLOCKS_PER_SEC);

    cout<<i<<"\t";
    cout<<(solucion.second)<<"\t";
    cout<<time<<"\t";

    cout<<"Sol=";
    for(int i=0; i<(solucion.first).size(); i++){
      cout<<(solucion.first)[i]<<",";
    }
    cout<<endl;
  }

  cout<<"\nAGE-PMX: "<<endl;
  for(int i=0; i<archivos.size(); i++){
    archivo = archivos[i];

    problema.load(ruta+archivo);

    cruce = 1;
    meme = 0;

    t0 = clock();
    solucion = problema.solucion_GeneticoAGE(seed, cruce, 50);
    t1 = clock();
    time = (double(t1-t0)/CLOCKS_PER_SEC);

    cout<<i<<"\t";
    cout<<(solucion.second)<<"\t";
    cout<<time<<"\t";

    cout<<"Sol=";
    for(int i=0; i<(solucion.first).size(); i++){
      cout<<(solucion.first)[i]<<",";
    }
    cout<<endl;
  }

  cout<<"\nAM-PMX-1: "<<endl;
  for(int i=0; i<archivos.size(); i++){
    archivo = archivos[i];
    //cout<<endl<<"Archivo: "<<archivo<<endl;
    problema.load(ruta+archivo);

    cruce = 1;
    meme = 1;

    t0 = clock();
    solucion = problema.solucion_GeneticoAGG(seed, cruce, meme, 10);
    t1 = clock();
    time = (double(t1-t0)/CLOCKS_PER_SEC);

    cout<<i<<"\t";
    cout<<(solucion.second)<<"\t";
    cout<<time<<"\t";

    cout<<"Sol=";
    for(int i=0; i<(solucion.first).size(); i++){
      cout<<(solucion.first)[i]<<",";
    }
    cout<<endl;
  }

  cout<<"\nAM-PMX-2: "<<endl;
  for(int i=0; i<archivos.size(); i++){
    archivo = archivos[i];
    //cout<<endl<<"Archivo: "<<archivo<<endl;
    problema.load(ruta+archivo);

    cruce = 1;
    meme = 2;

    t0 = clock();
    solucion = problema.solucion_GeneticoAGG(seed, cruce, meme, 10);
    t1 = clock();
    time = (double(t1-t0)/CLOCKS_PER_SEC);

    cout<<i<<"\t";
    cout<<(solucion.second)<<"\t";
    cout<<time<<"\t";

    cout<<"Sol=";
    for(int i=0; i<(solucion.first).size(); i++){
      cout<<(solucion.first)[i]<<",";
    }
    cout<<endl;
  }

  cout<<"\nAM-PMX-3: "<<endl;
  for(int i=0; i<archivos.size(); i++){
    archivo = archivos[i];
    problema.load(ruta+archivo);

    cruce = 1;
    meme = 3;

    t0 = clock();
    solucion = problema.solucion_GeneticoAGG(seed, cruce, meme, 10);
    t1 = clock();
    time = (double(t1-t0)/CLOCKS_PER_SEC);

    cout<<i<<"\t";
    cout<<(solucion.second)<<"\t";
    cout<<time<<"\t";

    cout<<"Sol=";
    for(int i=0; i<(solucion.first).size(); i++){
      cout<<(solucion.first)[i]<<",";
    }
    cout<<endl;
  }

  // Basados en Trayectorias

  cout<<"\nES: "<<endl;
  for(int i=0; i<archivos.size(); i++){
    archivo = archivos[i];

    problema.load(ruta+archivo);


    t0 = clock();
    solucion = problema.ES(seed);
    t1 = clock();
    time = (double(t1-t0)/CLOCKS_PER_SEC);

    cout<<i<<"\t";
    cout<<(solucion.second)<<"\t";
    cout<<time<<"\t";

    cout<<"Sol=";
    for(int i=0; i<(solucion.first).size(); i++){
      cout<<(solucion.first)[i]<<",";
    }
    cout<<endl;
  }

  cout<<"\nBMB: "<<endl;
  for(int i=0; i<archivos.size(); i++){
    archivo = archivos[i];

    problema.load(ruta+archivo);

    t0 = clock();
    solucion = problema.BMB(seed);
    t1 = clock();
    time = (double(t1-t0)/CLOCKS_PER_SEC);

    cout<<i<<"\t";
    cout<<(solucion.second)<<"\t";
    cout<<time<<"\t";

    cout<<"Sol=";
    for(int i=0; i<(solucion.first).size(); i++){
      cout<<(solucion.first)[i]<<",";
    }
    cout<<endl;
  }

  cout<<"\nGRASP: "<<endl;
  for(int i=0; i<archivos.size(); i++){
    archivo = archivos[i];

    problema.load(ruta+archivo);

    t0 = clock();
    solucion = problema.GRASP(seed);
    t1 = clock();
    time = (double(t1-t0)/CLOCKS_PER_SEC);

    cout<<i<<"\t";
    cout<<(solucion.second)<<"\t";
    cout<<time<<"\t";

    cout<<"Sol=";
    for(int i=0; i<(solucion.first).size(); i++){
      cout<<(solucion.first)[i]<<",";
    }
    cout<<endl;
  }

  cout<<"\nILS: "<<endl;
  for(int i=0; i<archivos.size(); i++){
    archivo = archivos[i];

    problema.load(ruta+archivo);

    t0 = clock();
    solucion = problema.ILS(seed);
    t1 = clock();
    time = (double(t1-t0)/CLOCKS_PER_SEC);

    cout<<i<<"\t";
    cout<<(solucion.second)<<"\t";
    cout<<time<<"\t";

    cout<<"Sol=";
    for(int i=0; i<(solucion.first).size(); i++){
      cout<<(solucion.first)[i]<<",";
    }
    cout<<endl;
  }

  cout<<"\nILS-ES: "<<endl;
  for(int i=0; i<archivos.size(); i++){
    archivo = archivos[i];

    problema.load(ruta+archivo);

    t0 = clock();
    solucion = problema.ILS_ES(seed);
    t1 = clock();
    time = (double(t1-t0)/CLOCKS_PER_SEC);

    cout<<i<<"\t";
    cout<<(solucion.second)<<"\t";
    cout<<time<<"\t";

    cout<<"Sol=";
    for(int i=0; i<(solucion.first).size(); i++){
      cout<<(solucion.first)[i]<<",";
    }
    cout<<endl;
  }

}
