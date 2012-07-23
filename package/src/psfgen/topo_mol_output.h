#ifndef TOPO_MOL_OUTPUT_H
#define TOPO_MOL_OUTPUT_H

#include <stdio.h>
#include "topo_mol.h"



int topo_mol_write_translate(topo_mol *mol, FILE *file, void *, void (*print_msg)(void *, const char *));

int topo_mol_write_top(topo_mol *mol, FILE *file);






#endif

