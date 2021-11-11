// cell_radius.h

#ifndef CELL_RADIUS_H
#define CELL_RADIUS_H

struct CellRadiusKernel;

typedef struct CellRadiusFuncs {
  double (*computeCellRadius)(struct CellRadiusKernel* kernel, int step);
  void (*destroy)(struct CellRadiusKernel* kernel);
} CellRadiusFuncs;

typedef struct CellRadiusKernel {
  CellRadiusFuncs* funcs;
  void* derived;
} CellRadiusKernel;

double computeCellRadius(CellRadiusKernel* kernel, int step);
void deleteCellRadiusKernel(CellRadiusKernel* kernel);
void setCellRadiusKernel(CellRadiusKernel* kernel, void* derived,
			 CellRadiusFuncs* funcs);

// For constant cell radius
CellRadiusKernel* createConstCellRadiusKernel(double r);

// For cosine cell radius
// i.e., R = 0.5*(R_max-R_min)*[cos(2*pi*t/period+shift)+1.0]+R_min
CellRadiusKernel* createCosineCellRadiusKernel(double rmin, double rmax,
					       double period, double shift);
#endif
