// cell_radius.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cell_radius.h"
#include "constant.h"

double computeCellRadius(CellRadiusKernel* kernel, int step) {
  return kernel->funcs->computeCellRadius(kernel, step);
}

void deleteCellRadiusKernel(CellRadiusKernel* kernel) {
  kernel->funcs->destroy(kernel);
}

void setCellRadiusKernel(CellRadiusKernel* kernel, void* derived,
			 CellRadiusFuncs* funcs) {
  kernel->derived = derived;
  kernel->funcs = funcs;
}

// For constant cell radius
typedef struct ConstCellRadiusKernel {
  CellRadiusKernel super;
  double radius;
} ConstCellRadiusKernel;

double constCellRadius(ConstCellRadiusKernel* kernel, int step) {
  return kernel->radius;
}

void deleteConstCellRadiusKernel(ConstCellRadiusKernel* kernel) {
  free(kernel);
}

CellRadiusFuncs constCellRadiusFuncs = 
  {
    .computeCellRadius = (double(*)(CellRadiusKernel*, int)) &constCellRadius,
    .destroy = (void(*)(CellRadiusKernel*)) &deleteConstCellRadiusKernel
  };

CellRadiusKernel* createConstCellRadiusKernel(double r) {
  ConstCellRadiusKernel* kernel = malloc(sizeof *kernel);
  setCellRadiusKernel(&kernel->super, kernel, &constCellRadiusFuncs);
  kernel->radius = r;
  return (CellRadiusKernel*) kernel;
}

// For cosine cell radius
// i.e., R = 0.5*(R_max-R_min)*[cos(2*pi*t/period+shift)+1.0]+R_min
typedef struct CosineCellRadiusKernel {
  CellRadiusKernel super;
  double radiusMin;
  double radiusMax;
  double freq;
  double shift;
} CosineCellRadiusKernel;

double cosineCellRadius(CosineCellRadiusKernel* kernel, int step) {
  return 0.5 * (kernel->radiusMax - kernel->radiusMin) * 
    (1.0 - cos(kernel->freq * step + kernel->shift)) + kernel->radiusMin;
}

void deleteCosineCellRadiusKernel(CosineCellRadiusKernel* kernel) {
  free(kernel);
}

CellRadiusFuncs cosineCellRadiusFuncs = 
  {
    .computeCellRadius = (double(*)(CellRadiusKernel*, int)) &cosineCellRadius,
    .destroy = (void(*)(CellRadiusKernel*)) &deleteCosineCellRadiusKernel
  };

CellRadiusKernel* createCosineCellRadiusKernel(double rmin, double rmax,
					       double period, double shift) {
  CosineCellRadiusKernel* kernel = malloc(sizeof *kernel);
  setCellRadiusKernel(&kernel->super, kernel, &cosineCellRadiusFuncs);
  kernel->radiusMin = rmin;
  kernel->radiusMax = rmax;
  kernel->freq = 2.0 * PF_PI / period;
  kernel->shift = 2.0 * PF_PI * shift;
  return (CellRadiusKernel*) kernel;  
};
