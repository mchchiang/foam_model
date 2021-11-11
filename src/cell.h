// cell.h

#ifndef CELL_H
#define CELL_H

#include "random.h"
#include "cell_radius.h"

typedef struct {
  double** field[2]; // phase field
  int setIndex;
  int getIndex;
  int lx; // x size of lattice
  int ly; // y size of lattice
  int x; // x pos of the cell relative to main lattice
  int y; // y pos of the cell relative to main lattice
  int type;
  Random* random;
  double xcm; // x centre of mass in cell's frame
  double ycm; // y centre of mass in cell's frame
  double drx; // Change in x centre of mass
  double dry; // Change in y centre of mass
  double deltaXCM;
  double deltaYCM;
  double volume; // Total volume of the cell
  double radius; // Ideal radius of the cell
  CellRadiusKernel* radiusKernel; // Control the radius of the cell
  double incell;
} Cell;

Cell* createCell(int x, int y, int lx, int ly, double incell, 
		 unsigned long seed);
void deleteCell(Cell* cell);
void setRadiusKernel(Cell* cell, char* kernelType, char* kernelArgs);
void setField(Cell* cell, double** field);
void updateRadius(Cell* cell, int step);
void updateVolume(Cell* cell);
void updateCM(Cell* cell);
void shiftCoordinates(Cell* cell, int xShift, int yShift);
void calculateCM(Cell* cell, double* xcm, double* ycm);
void startUpdateCellField(Cell* cell);
void endUpdateCellField(Cell* cell);

#endif
