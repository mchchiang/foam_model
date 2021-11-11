// model.h

#ifndef MODEL_H
#define MODEL_H

#include "cell.h"
#include "dump.h"

typedef struct Model {
  int lx;
  int ly;
  int numOfCells;
  double dt;
  double thickness;
  double cahnHilliardCoeff;
  double volumePenaltyCoeff;
  double repulsionCoeff;
  double cellRadius;
  int cellLx;
  int cellLy;
  double* cellXCM;
  double* cellYCM;
  int* cellXBoundCount;
  int* cellYBoundCount;
  Cell** cells;

  // Reduction fields
  double** totalFieldSq;
  
  // Dumps
  Dump** dumps;
  int ndumps;
} Model;

Model* createModel(int lx, int ly,int ncells);
void deleteModel(Model* model);
void initCellsFromFile(Model* model, char* cmFile, char* shapeFile,
		       unsigned long seed);
void run(Model* model, int nsteps);
void output(Model* model, int step);

void updateTotalField(Model* model);
void updateCellProperties(Cell* cell, int step);
void updateCellField(Model* model, Cell* cell, int cellIndex, int step);
void updateCellCM(Model* model, Cell* cell, int cellIndex);

#endif
