// phase_field_model.c

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "array.h"
#include "cell.h"
#include "model.h"
#include "dump.h"
#include "random.h"
#include "constant.h"
#include "util.h"

Model* createModel(int lx, int ly, int ncells) {
  Model* model =  malloc(sizeof *model);
  model->lx = lx;
  model->ly = ly;
  model->numOfCells = ncells;
  model->dt = 0.01;
  model->cahnHilliardCoeff = 1.0;
  model->volumePenaltyCoeff = 1.0;
  model->repulsionCoeff = 1.0;
  model->thickness = 1.0;
  model->cellRadius = 1.0;
  model->cellLx = 1;
  model->cellLy = 1;
  model->cells = malloc(sizeof *model->cells * ncells);
  model->totalFieldSq = create2DDoubleArray(model->lx, model->ly);
  model->cellXCM = create1DDoubleArray(ncells);
  model->cellYCM = create1DDoubleArray(ncells);
  model->cellXBoundCount = create1DIntArray(ncells);
  model->cellYBoundCount = create1DIntArray(ncells);
  model->dumps = NULL;
  model->ndumps = 0;
  return model; 
}

void deleteModel(Model* model) {
  for (int i = 0; i < model->numOfCells; i++) {
    if (model->cells[i] != NULL) {
      deleteCell(model->cells[i]);
    }
  }
  free(model->cells);
  free(model->totalFieldSq);
  free(model->cellXCM);
  free(model->cellYCM);
  free(model->cellXBoundCount);
  free(model->cellYBoundCount);
  free(model);
}

void initCellsFromFile(Model* model, char* cellFile,
		       char* shapeFile, unsigned long seed) {
  FILE* fcell = fopen(cellFile, "r");
  if (fcell == NULL) {
    printf("Problem with opening the centre of mass file!\n");
    return;
  }

  FILE* fshape = fopen(shapeFile, "r");
  if (fshape == NULL) {
    printf("Problem with opening the shape file!\n");
    return;
  }

  char line [1000];
  int nargs, x, y;
  double val;
  int clx = model->cellLx;
  int cly = model->cellLy;

  // Allocate memory for the template field
  double** field = create2DDoubleArray(clx, cly);
  // Read the template field from the shape file
  while (fgets(line, sizeof(line), fshape) != NULL) {
    nargs = sscanf(line, "%d %d %lf", &x, &y, &val);
    if (nargs == 3 && x >= 0 && x < clx && y >= 0 && y < cly) {
      field[x][y] = val;
    }
  }
  
  int index;
  double xcm, ycm;
  int count = 0;
  Cell* cell;
  char kernelArgs[80];
  while (fgets(line, sizeof(line), fcell) != NULL) {
    nargs = sscanf(line, "%d %lf %lf", &index, &xcm, &ycm);
    if (nargs == 3) {
      x = (int) round(xcm-clx/2.0);
      y = (int) round(ycm-cly/2.0);
      cell = createCell(x, y, clx, cly, 0.5, index+seed);
      model->cells[index] = cell;
      sprintf(kernelArgs, "%.5f", model->cellRadius);
      setRadiusKernel(cell, "const", kernelArgs);
      setField(cell, field);
      count++;
    } else {
      printf("ERROR: not enough arguments supplied for a cell\n");
      exit(1);
    }
  }

  if (count != model->numOfCells) {
    printf("ERROR: not all cells initialised!\n");
    exit(1);
  }
  free(field);
  fclose(fcell);
  fclose(fshape);
}

void run(Model* model, int nsteps) {
  output(model, 0); // Output initial state of the model
  for (int step = 1; step <= nsteps; step++) {
    //printf("Step: %d\n", step);
    //printf("Total field ...\n");
    updateTotalField(model);

#pragma omp parallel for default(none) shared(model, step) schedule(static)
    for (int m = 0; m < model->numOfCells; m++) {
      Cell* cell = model->cells[m];
      //printf("%d: Cell properties ...\n", m);
      updateCellProperties(cell, step);
      //printf("%d: Cell field ...\n", m);
      updateCellField(model, cell, m, step);
      //printf("%d: Cell cm ...\n", m);
      updateCellCM(model, cell, m);
    }
    output(model, step);
  }
}

void output(Model* model, int step) {
  if (step % 1000 == 0) {
    printf("Step %d\n", step);
  }
  for (int i = 0; i < model->ndumps; i++) {
    dumpOutput(model->dumps[i], model, step);
  }
}

void updateTotalField(Model* model) {
  // Reset all global fields
#pragma omp parallel for default(none) shared(model) schedule(static)
  for (int i = 0; i < model->lx; i++) {
    for (int j = 0; j < model->ly; j++) {
      model->totalFieldSq[i][j] = 0.0;
    }
  }

  Cell* cell;
  int clx, cly, x, y, cx, cy;
  double phi;
  double** cellField;
  
  for (int m = 0; m < model->numOfCells; m++) {
    cell = model->cells[m];
    clx = cell->lx;
    cly = cell->ly;
    cx = cell->x;
    cy = cell->y;
    cellField = cell->field[cell->getIndex];
#pragma omp parallel for default(none)			\
  shared(model, cell, clx, cly, cx, cy, cellField)	\
  private(x, y, phi) schedule(static)
    for (int i = 0; i < clx; i++) {
      x = iwrap(model->lx, cx+i);
      for (int j = 0; j < cly; j++) {
	y = iwrap(model->ly, cy+j);
	phi = cellField[i][j];
	model->totalFieldSq[x][y] += phi*phi;
      }
    }
  }
}

void updateCellProperties(Cell* cell, int step) {
  updateRadius(cell, step);
  updateVolume(cell);
}

void updateCellField(Model* model, Cell* cell, int m, int step) {
  startUpdateCellField(cell);
  int clx = cell->lx;
  int cly = cell->ly;
  int cx = cell->x;
  int cy = cell->y;
  int set = cell->setIndex;
  int get = cell->getIndex;
  int x, y; // Lab frame coordinates of a lattice element
  int iu, id, ju, jd; // Nearest neighbours
  double** cellField = cell->field[get];
  double thicknessSq = model->thickness * model->thickness;
  double phi, cahnHilliard, volumeConstraint, repulsion;
  double scaledVolume = cell->volume / (PF_PI * cell->radius * cell->radius);
  
  // Apply fixed (Dirichlet) boundary conditions (phi = 0 at boundaries)
  // i and j are coordinates in the cell's own reference frame
  for (int i = 2; i < clx-2; i++) {
    iu = iup(clx, i);
    id = idown(clx, i);
    x = iwrap(model->lx, cx+i);
    for (int j = 2; j < cly-2; j++) {
      ju = iup(cly, j);
      jd = idown(cly, j);
      y = iwrap(model->ly, cy+j);
      phi = cellField[i][j];
      
      // Cahn-Hilliard term
      cahnHilliard =  2.0 * model->cahnHilliardCoeff * 
	(2.0 * phi * (phi - 1.0) * (phi - 0.5) - 
	 thicknessSq * laplacian(i, j, iu, id, ju, jd, cellField));

      // Volume term
      volumeConstraint = 4.0 * model->volumePenaltyCoeff * phi * 
	(scaledVolume - 1.0);
      
      // Repulsion term
      repulsion = 4.0 * model->repulsionCoeff * phi * 
	(model->totalFieldSq[x][y] - phi*phi);

      // Update cell field
      cell->field[set][i][j] = cell->field[get][i][j] - model->dt *
	(cahnHilliard + volumeConstraint + repulsion);
    }
  }
  endUpdateCellField(cell);
}

void updateCellCM(Model* model, Cell* cell, int cellIndex) {
  updateCM(cell);
  int ix, iy;
  double x, y, cx, cy;
  cx = cell->x;
  cy = cell->y;
  x = cx+cell->xcm;
  y = cy+cell->ycm;
  ix = (int) floor(x / model->lx);
  iy = (int) floor(y / model->ly);
  model->cellXCM[cellIndex] = x - ix * model->lx;
  model->cellYCM[cellIndex] = y - iy * model->ly;
  model->cellXBoundCount[cellIndex] = ix;
  model->cellYBoundCount[cellIndex] = iy;
}
