# CONSTANTES SIMBÓLICAS
BIN = bin
SRC = fuentes
OBJ = obj
LIB = lib
INCLUDE = include
CC = g++


all: $(BIN)/main

$(BIN)/main : $(OBJ)/main.o $(OBJ)/Matrix.o $(OBJ)/QAP.o
	$(CC) -o $(BIN)/main $(OBJ)/main.o $(OBJ)/Matrix.o $(OBJ)/QAP.o  -std=c++14 -O2

$(OBJ)/main.o : $(SRC)/main.cpp $(INCLUDE)/Matrix.h $(INCLUDE)/QAP.h
	$(CC) -c -o $(OBJ)/main.o $(SRC)/main.cpp -I$(INCLUDE) -std=c++14 -O2

$(OBJ)/Matrix.o : $(SRC)/Matrix.cpp $(INCLUDE)/Matrix.h
	$(CC) -c -o $(OBJ)/Matrix.o $(SRC)/Matrix.cpp -I$(INCLUDE) -std=c++14 -O2

$(OBJ)/QAP.o : $(SRC)/QAP.cpp $(INCLUDE)/QAP.h $(INCLUDE)/Matrix.h
	$(CC) -c -o $(OBJ)/QAP.o $(SRC)/QAP.cpp -I$(INCLUDE) -std=c++14 -O2


# Limpieza
clean:
	-rm $(BIN)/*
	-rm $(OBJ)/*.*

tar:
	tar -cvf practica1.tar ./
