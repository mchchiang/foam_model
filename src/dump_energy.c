// dump_energy.c
// Dump the cm of the cells

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "array.h"
#include "dump.h"
#include "model.h"
#include "cell.h"
#include "util.h"
#include "constant.h"

typedef struct EnergyDump {
  Dump super; // Base struct must be the first element
  bool overwrite;
  int lx;
  int ly;
  double** phi2Field; // Store the sum of phi^2 field
  double** phi4Field; // Store the sum of phi^4 field
} EnergyDump;

void energyOutput(EnergyDump* dump, Model* model, int step) {
  FILE* f;
  f = fopen(dump->super.filename, "a");
  
  // Reset fields
  for (int i = 0; i < dump->lx; i++) {
    for (int j = 0; j < dump->ly; j++) {
      dump->phi2Field[i][j] = 0.0;
      dump->phi4Field[i][j] = 0.0;
    }
  }

  // Compute the auxillary field (sum of phi^2)
  Cell* cell;
  double phi, phi2;
  int clx, cly, x, y, cx, cy;
  double** cellField;
  for (int m = 0; m < model->numOfCells; m++) {
    cell = model->cells[m];
    clx = cell->lx;
    cly = cell->ly;
    cx = cell->x;
    cy = cell->y;
    cellField = cell->field[cell->getIndex];
    for (int i = 0; i < clx; i++) {
      x = iwrap(model->lx, cx+i);
      for (int j = 0; j < cly; j++) {
	y = iwrap(model->ly, cy+j);
	phi = cellField[i][j];
	phi2 = phi*phi;
	dump->phi2Field[x][y] += phi2;
	dump->phi4Field[x][y] += phi2*phi2;
      }
    }
  }

  // Compute the free energy
  int iu, iuu, id, idd, ju, juu, jd, jdd;
  double dphi, dphi2, phi2Sum;
  double gphix, gphiy;
  double cellVolume, piR2, dV;
  double cahnHilliard = 0.0;
  double volumeConstraint = 0.0;
  double repulsion = 0.0;
  double thicknessSq = model->thickness * model->thickness;
  for (int m = 0; m < model->numOfCells; m++) {
    cell = model->cells[m];
    clx = cell->lx;
    cly = cell->ly;
    cx = cell->x;
    cy = cell->y;
    cellField = cell->field[cell->getIndex];
    cellVolume = 0.0;
    piR2 = PF_PI * cell->radius * cell->radius;
    for (int i = 2; i < clx-2; i++) {
      iu = iup(clx, i);
      iuu = iup(clx, iu);
      id = idown(clx, i);
      idd = idown(clx, id);
      x = iwrap(model->lx, cx+i);
      for (int j = 2; j < cly-2; j++) {
	ju = iup(cly, j);
	juu = iup(cly, ju);
	jd = idown(cly, j);
	jdd = idown(cly, jd);
	y = iwrap(model->ly, cy+j);
	
	phi = cellField[i][j];
	phi2 = phi*phi;
	dphi = phi-1.0;
	dphi2 = dphi*dphi;
	
	// Cahn-Hilliard term
	gphix = grad4(i, j, iuu, iu, id, idd, 0, cellField);
	gphiy = grad4(i, j, juu, ju, jd, jdd, 1, cellField);
	cahnHilliard += model->cahnHilliardCoeff * 
	  (phi2 * dphi2 + thicknessSq * (gphix * gphix + gphiy * gphiy));
	
	// Calculate cell volume
	cellVolume += phi2;
	
	// Repulsion term
	phi2Sum = dump->phi2Field[x][y];
	repulsion += model->repulsionCoeff * 
	  (phi2Sum*phi2Sum - dump->phi4Field[x][y]);
      }
    }
    dV = 1.0 - cellVolume / piR2;
    volumeConstraint += model->volumePenaltyCoeff * piR2 * dV * dV;
  }
  //  energy = cahnHilliard + volumeConst + repulsion;
  double energy = cahnHilliard + volumeConstraint + repulsion;
  fprintf(f, "%d %.5f %.5f %.5f %.5f\n", step, cahnHilliard, volumeConstraint, 
	  repulsion, energy);
  fclose(f);
}

void deleteEnergyDump(EnergyDump* dump) {
  free(dump->phi2Field);
  free(dump->phi4Field);
  free(dump);
}

DumpFuncs energyDumpFuncs =
  {
   .output = (void (*)(Dump*, Model*, int)) &energyOutput,
   .destroy = (void (*)(Dump*)) &deleteEnergyDump
  };

Dump* createEnergyDump(char* filename, int lx, int ly, int printInc, 
		       bool overwrite) {
  EnergyDump* dump = malloc(sizeof *dump);
  setDump(&dump->super, dump, filename, printInc, &energyDumpFuncs);
  dump->overwrite = overwrite;
  dump->lx = lx;
  dump->ly = ly;
  dump->phi2Field = create2DDoubleArray(lx, ly);
  dump->phi4Field = create2DDoubleArray(lx, ly);
  return (Dump*) dump;
}
