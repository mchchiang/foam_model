// run_model.c

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "dump.h"
#include "model.h"
#include "cell_radius.h"
#include "constant.h"

int main (int argc, char* argv[]) {
  
  if (argc != 2) {
    printf("Usage: run_model paramsFile\n");
    return 1;
  }

  int argi = 0;
  char* paramsFile = argv[++argi];
  
  printf("Reading parameters ...\n");

  FILE* fparams = fopen(paramsFile, "r");
  if (fparams == NULL) {
    printf("Error in opening parameter file!\n");
    return 1;
  }

  char line [PF_DIR_SIZE], dumpMode [PF_DIR_SIZE];
  char cellFile [PF_DIR_SIZE], shapeFile [PF_DIR_SIZE], dumpFile [PF_DIR_SIZE];
  double dt, radius, thickness, cahnHilliardCoeff, 
    volumePenaltyCoeff, repulsionCoeff;
  int cellLx = -1;
  int cellLy = -1;
  int lx = -1;
  int ly = -1;
  int nequil, nsteps, ncells;
  unsigned long seed;
  int nparams = 0;

  int ndumps = 0;
  int nedumps = 0;
  int nradii = 0;

  while (fgets(line, sizeof(line), fparams) != NULL) {
    nparams += sscanf(line, "cahnHilliardCoeff = %lf", &cahnHilliardCoeff);
    nparams += sscanf(line, "volumePenaltyCoeff = %lf", &volumePenaltyCoeff);
    nparams += sscanf(line, "repulsionCoeff = %lf", &repulsionCoeff);
    nparams += sscanf(line, "thickness = %lf", &thickness);
    nparams += sscanf(line, "cellRadius = %lf", &radius);
    nparams += sscanf(line, "cellLx = %d", &cellLx);
    nparams += sscanf(line, "cellLy = %d", &cellLy);
    nparams += sscanf(line, "lx = %d", &lx);
    nparams += sscanf(line, "ly = %d", &ly);
    nparams += sscanf(line, "nsteps = %d", &nsteps);
    nparams += sscanf(line, "nequil = %d", &nequil);
    nparams += sscanf(line, "ncells = %d", &ncells);
    nparams += sscanf(line, "dt = %lf", &dt);
    nparams += sscanf(line, "cell_file = %s", cellFile);
    nparams += sscanf(line, "shape_file = %s", shapeFile);
    nparams += sscanf(line, "seed = %ld", &seed);
    
    // Count number of dumps 
    if (strstr(line, "dump_") != NULL) {
      if (strstr(line, "equil") != NULL) {
	nedumps++;
      } else if (strstr(line, "main") != NULL) {
	ndumps++;
      }
      // Count number of radius fixes
    } else if (strstr(line, "fix_radius") != NULL) {
      nradii++;
    }
  }
  
  fclose(fparams);
  
  if (nparams != 16) {
    printf("Not enough parameters specified!\n");
    return 1;
  }  
  printf("Read parameters:\n");
  printf("lx = %d\n", lx);
  printf("ly = %d\n", ly);
  printf("cellLx = %d\n", cellLx);
  printf("cellLy = %d\n", cellLy);
  printf("ncells = %d\n", ncells);
  printf("cahnHilliardCoeff = %.5f\n", cahnHilliardCoeff);
  printf("volumePenaltyCoeff = %.5f\n", volumePenaltyCoeff);
  printf("repulsionCoeff = %.5f\n", repulsionCoeff);
  printf("thickness = %.5f\n", thickness);
  printf("radius = %.5f\n", radius);
  printf("dt = %.5f\n", dt);
  printf("seed = %ld\n", seed);
  printf("cell_file = %s\n", cellFile);
  printf("shape_file = %s\n", shapeFile);  
  printf("nsteps = %d\n", nsteps);
  printf("nequil = %d\n", nequil);
  printf("nedumps = %d\n", nedumps);
  printf("ndumps = %d\n", ndumps);
  printf("nradii = %d\n", nradii);
					    
  // Read dumps and fixes
  int idump = 0;
  int iedump = 0;
  int iradius = 0;
  int printInc, overwrite, cellIndex;
  
#if PF_HAS_ARMA
  int fieldScale, kernelLength, sgolayDegree, sgolayLength;
  double kernelSigma;
#endif
  Dump** equilDumps = malloc(sizeof *equilDumps * nedumps);
  Dump** dumps = malloc(sizeof *dumps * ndumps);

  int nchar;
  int* radiusKernelCellIndex = malloc(sizeof *radiusKernelCellIndex * nradii);
  char** radiusKernelType = malloc(sizeof *radiusKernelType * nradii);
  char** radiusKernelArgs = malloc(sizeof *radiusKernelArgs * nradii);
  for (int i = 0; i < nradii; i++) {
    radiusKernelType[i] = malloc(sizeof *radiusKernelType[i] * 50);
    radiusKernelArgs[i] = malloc(sizeof *radiusKernelArgs[i] * 1000);
  }
  
  fparams = fopen(paramsFile, "r");
  while (fgets(line, sizeof(line), fparams) != NULL) {
    // Read dumps
    // CM dump
    if (sscanf(line, "dump_cm %d %d %s %s", 
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0) {
	equilDumps[iedump] = createCMDump(dumpFile, printInc, overwrite);
	iedump++;
      } else if (strcmp(dumpMode, "main") == 0) {
	dumps[idump] = createCMDump(dumpFile, printInc, overwrite);
	idump++; 
      }
    }
    // Bulk CM dump
    if (sscanf(line, "dump_bulk_cm %d %s %s",
	       &printInc, dumpMode, dumpFile) == 3) {
      if (strcmp(dumpMode, "equil") == 0) {
	equilDumps[iedump] = createBulkCMDump(dumpFile, printInc);
	iedump++;
      } else if (strcmp(dumpMode, "main") == 0) {
	dumps[idump] = createBulkCMDump(dumpFile, printInc);
	idump++; 
      }
    }
    // Gyration dump
    if (sscanf(line, "dump_gyr %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0) {
	equilDumps[iedump] = 
	  createGyrationDump(dumpFile, printInc, overwrite);
	iedump++;
      } else if (strcmp(dumpMode, "main") == 0) {
	dumps[idump] = createGyrationDump(dumpFile, printInc, overwrite);
	idump++;
      }
    }
    // Gyration field dump
    if (sscanf(line, "dump_gyr_field %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0) {
	equilDumps[iedump] = 
	  createGyrationFieldDump(dumpFile, printInc, overwrite);
	iedump++;
      } else if (strcmp(dumpMode, "main") == 0) {
	dumps[idump] = createGyrationFieldDump(dumpFile, printInc, overwrite);
	idump++;
      }
    }
    // Deform dump
    if (sscanf(line, "dump_deform %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0) {
	equilDumps[iedump] = 
	  createDeformDump(dumpFile, printInc, overwrite);
	iedump++;
      } else if (strcmp(dumpMode, "main") == 0) {
	dumps[idump] = createDeformDump(dumpFile, printInc, overwrite);
	idump++;
      }
    }
    // Deform field dump
    if (sscanf(line, "dump_deform_field %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0) {
	equilDumps[iedump] = 
	  createDeformFieldDump(dumpFile, printInc, overwrite);
	iedump++;
      } else if (strcmp(dumpMode, "main") == 0) {
	dumps[idump] = createDeformFieldDump(dumpFile, printInc, overwrite);
	idump++;
      }
    }
    // Total field dump
    if (sscanf(line, "dump_field %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0) {
	equilDumps[iedump] = createFieldDump(dumpFile, printInc, overwrite);
	iedump++;
      } else if (strcmp(dumpMode, "main") == 0) {
	dumps[idump] = createFieldDump(dumpFile, printInc, overwrite);
	idump++;
      }
    }
    // Total index field dump
    if (sscanf(line, "dump_index_field %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      if (strcmp(dumpMode, "equil") == 0) {
	equilDumps[iedump] =
	  createIndexFieldDump(dumpFile, printInc, overwrite);
	iedump++;
      } else if (strcmp(dumpMode, "main") == 0) {
	dumps[idump] = createIndexFieldDump(dumpFile, printInc, overwrite);
	idump++;
      }
    }    
    // Individual cell field dump
    if (sscanf(line, "dump_cell_field %d %d %d %s %s",
	       &cellIndex, &printInc, &overwrite, dumpMode, dumpFile) == 5) {
      if (cellIndex < 0 || cellIndex >= ncells) {
	printf("ERROR: cell index out of bounds\n");
	exit(1);
      }
      if (strcmp(dumpMode, "equil") == 0) {
	equilDumps[iedump] =
	  createCellFieldDump(dumpFile, cellIndex, printInc, overwrite);
	iedump++;
      } else if (strcmp(dumpMode, "main") == 0) {
	dumps[idump] =
	  createCellFieldDump(dumpFile, cellIndex, printInc, overwrite);
	idump++;
      }
    }
    // Neighbour dump
    if (sscanf(line, "dump_neighbour %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      // Only created the dump when the field size is known,
      // as it is needed for creating the neighbour analysers
      if (lx > 0 && ly > 0) {
	if (strcmp(dumpMode, "equil") == 0) {
	  equilDumps[iedump] =
	    createNeighbourDump(dumpFile, lx, ly, printInc, overwrite);
	  iedump++;
	} else if (strcmp(dumpMode, "main") == 0) {
	  dumps[idump] =
	    createNeighbourDump(dumpFile, lx, ly, printInc, overwrite);
	  idump++;
	}
      }
    }
    // Energy dump
    if (sscanf(line, "dump_energy %d %d %s %s",
	       &printInc, &overwrite, dumpMode, dumpFile) == 4) {
      // Only created the dump when the field size is known,
      // as it is needed for creating the neighbour analysers
      if (lx > 0 && ly > 0) {
	if (strcmp(dumpMode, "equil") == 0) {
	  equilDumps[iedump] =
	    createEnergyDump(dumpFile, lx, ly, printInc, overwrite);
	  iedump++;
	} else if (strcmp(dumpMode, "main") == 0) {
	  dumps[idump] =
	    createEnergyDump(dumpFile, lx, ly, printInc, overwrite);
	  idump++;
	}
      }
    }
#if PF_HAS_ARMA
    // Shape dump
    if (sscanf(line, "dump_shape %d %d %lf %d %d %d %d %s %s",
	       &fieldScale, &kernelLength, &kernelSigma,
	       &sgolayDegree, &sgolayLength,
	       &printInc, &overwrite, dumpMode, dumpFile) == 9) {
      // Only created the dump when cell field size and phi0 are known,
      // as they are needed for creating the shape analysers
      if (cellLx > 0 && cellLy > 0) {
	if (strcmp(dumpMode, "equil") == 0) {
	  equilDumps[iedump] =
	    createShapeDump(dumpFile, fieldScale, cellLx, cellLy,
			    kernelLength, kernelSigma, sgolayDegree,
			    sgolayLength, 0.5, printInc, overwrite);
	  iedump++;
	} else if (strcmp(dumpMode, "main") == 0) {
	  dumps[idump] =
	    createShapeDump(dumpFile, fieldScale, cellLx, cellLy,
			    kernelLength, kernelSigma, sgolayDegree,
			    sgolayLength, 0.5, printInc, overwrite);
	  idump++;
	}
      }
    }
#endif
    
    // Read cell radius kernels
    if (sscanf(line, "fix_radius %d %s %n", &cellIndex,
	       radiusKernelType[iradius], &nchar) == 2) {
      if (cellIndex < 0 || cellIndex >= ncells) {
	printf("ERROR: cell index out of bounds\n");
	exit(1);
      }
      radiusKernelCellIndex[iradius] = cellIndex;
      strcpy(radiusKernelArgs[iradius], line + nchar);
      iradius++;
    } 
  } // Close loop over reading parameters
  
  printf("Initialising model ...\n");
  
  Model* model = createModel(lx, ly, ncells);
  model->cahnHilliardCoeff = cahnHilliardCoeff;
  model->volumePenaltyCoeff = volumePenaltyCoeff;
  model->repulsionCoeff = repulsionCoeff;
  model->thickness = thickness;
  model->cellRadius = radius;
  model->cellLx = cellLx;
  model->cellLy = cellLy;
  model->ndumps = nedumps;
  model->dumps = equilDumps;
  
  initCellsFromFile(model, cellFile, shapeFile, seed);
  
  printf("Done initialisation.\n");
  
  model->dt = dt;

  printf("Doing equilibration run ...\n");

#ifdef _OPENMP
  double start, end, duration;
  start = omp_get_wtime();
#endif
  run(model, nequil);

#ifdef _OPENMP
  end = omp_get_wtime();
  duration = end-start;
  printf("Time taken (sec): %.5f\n", duration);
  printf("\n");
#endif

  model->ndumps = ndumps;
  model->dumps = dumps;
  
  // Fix cell radius kernel
  for (int i = 0; i < nradii; i++) {
    cellIndex = radiusKernelCellIndex[i];
    setRadiusKernel(model->cells[cellIndex], radiusKernelType[i],
		    radiusKernelArgs[i]);
  }

  printf("Doing main simulation run ...\n");

#ifdef _OPENMP
  start = omp_get_wtime();
#endif

  run(model, nsteps);

#ifdef _OPENMP
  end = omp_get_wtime();
  duration = end-start;
  printf("Time taken (sec): %.5f\n", duration);
#endif
  
  // Clean up
  deleteModel(model);

  for (int i = 0; i < nedumps; i++) {
    deleteDump(equilDumps[i]);
  }
  free(equilDumps);

  for (int i = 0; i < ndumps; i++) {
    deleteDump(dumps[i]);
  }
  free(dumps);

  for (int i = 0; i < nradii; i++) {
    free(radiusKernelType[i]);
    free(radiusKernelArgs[i]);
  }
  free(radiusKernelCellIndex);
  free(radiusKernelType);
  free(radiusKernelArgs);
}
