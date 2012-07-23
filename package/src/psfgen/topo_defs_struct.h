#ifndef TOPO_DEFS_STRUCT_H
#define TOPO_DEFS_STRUCT_H

#include "memarena.h"
#include "hasharray.h"
#include <vector>
#include <string>

#define NAMEMAXLEN 8
#define NAMETOOLONG(X) ( strlen(X) >= NAMEMAXLEN )

struct topo_defs_type_t {
	char name[NAMEMAXLEN];
	char element[NAMEMAXLEN];
	int id;
	double mass;
};

struct topo_defs_atom_t {
	struct topo_defs_atom_t *next;
	char name[NAMEMAXLEN];
	char type[NAMEMAXLEN];
	double charge;
	int res, rel;

#ifdef SWIG
	%rename(_del) del;
#endif

	int del;
};

struct topo_defs_bond_t {
	struct topo_defs_bond_t *next;
	char atom1[NAMEMAXLEN];
	char atom2[NAMEMAXLEN];
	int res1, rel1;
	int res2, rel2;

#ifdef SWIG
	%rename(_del) del;
#endif
	int del;
};

struct topo_defs_angle_t {
	struct topo_defs_angle_t *next;
	char atom1[NAMEMAXLEN];
	char atom2[NAMEMAXLEN];
	char atom3[NAMEMAXLEN];
	int res1, rel1;
	int res2, rel2;
	int res3, rel3;

#ifdef SWIG
	%rename(_del) del;
#endif
	int del;
};

struct topo_defs_dihedral_t {
	struct topo_defs_dihedral_t *next;
	char atom1[NAMEMAXLEN];
	char atom2[NAMEMAXLEN];
	char atom3[NAMEMAXLEN];
	char atom4[NAMEMAXLEN];
	int res1, rel1;
	int res2, rel2;
	int res3, rel3;
	int res4, rel4;

#ifdef SWIG
	%rename(_del) del;
#endif
	int del;
};

struct topo_defs_improper_t {
	struct topo_defs_improper_t *next;
	char atom1[NAMEMAXLEN];
	char atom2[NAMEMAXLEN];
	char atom3[NAMEMAXLEN];
	char atom4[NAMEMAXLEN];
	int res1, rel1;
	int res2, rel2;
	int res3, rel3;
	int res4, rel4;

#ifdef SWIG
	%rename(_del) del;
#endif
	int del;
};

struct topo_defs_cmap_t {
	struct topo_defs_cmap_t *next;
	char atoml[8][NAMEMAXLEN];
	int resl[8], rell[8];

#ifdef SWIG
	%rename(_del) del;
#endif
	int del;
};

struct topo_defs_conformation_t {
	struct topo_defs_conformation_t *next;
	char atom1[NAMEMAXLEN];
	char atom2[NAMEMAXLEN];
	char atom3[NAMEMAXLEN];
	char atom4[NAMEMAXLEN];
	int res1, rel1;
	int res2, rel2;
	int res3, rel3;
	int res4, rel4;
	int improper;
	double dist12, angle123, dihedral, angle234, dist34;

#ifdef SWIG
	%rename(_del) del;
#endif
	int del;

};

#ifdef SWIG
%feature("python:slot", "tp_str",  functype="reprfunc") topo_defs_residue::__str__;
%feature("python:slot", "tp_repr",  functype="reprfunc") topo_defs_residue::__str__;
#endif

struct topo_defs_residue {
	char name[NAMEMAXLEN];
	int patch;
	topo_defs_atom_t *atoms;
	topo_defs_bond_t *bonds;
	topo_defs_angle_t *angles;
	topo_defs_dihedral_t *dihedrals;
	topo_defs_improper_t *impropers;
	topo_defs_cmap_t *cmaps;
	topo_defs_conformation_t *conformations;
	char pfirst[NAMEMAXLEN];
	char plast[NAMEMAXLEN];
    std::string __str__();
};

struct topo_defs_topofile_t {
	/*   struct topo_defs_topofile_t *next; */
	char filename[256];
};


#ifdef SWIG
%feature("python:slot", "tp_str",  functype="reprfunc") topo_defs::__str__;
%feature("python:slot", "tp_repr",  functype="reprfunc") topo_defs::__str__;
#endif

struct topo_defs {
	int auto_angles;
	int auto_dihedrals;
	int cmaps_present;
	char pfirst[NAMEMAXLEN];
	char plast[NAMEMAXLEN];

#ifndef SWIG
	topo_defs_topofile_t *topo_array;
	hasharray *topo_hash;

	topo_defs_type_t *type_array;
	hasharray *type_hash;

	topo_defs_residue *residue_array;
	hasharray *residue_hash;
	topo_defs_residue *buildres;
	int buildres_no_errors;
	memarena *arena;
#endif

    /**
     * creates a topo_defs object
     */
    topo_defs();

    ~topo_defs();

	/**
	 * checks if this defs has a residue with the given residue name.
	 * if one exists, returns 0. If one does not exist, a new one is created
	 * and inserted into the residue array and the buildres member points to
	 * this new residue.
	 *
	 * The new residue is intialized to all zeros except patch is set to the given
	 * patch value.
	 *
	 * buildres is set to 0 if either a residue already exists or fails to be created.
	 *
	 * @return 0 on success,
	 * -2 if name too long,
	 * -4 if the name fails to insert into hash array
	 */
	int build_residue(const char *rname, int patch);

	int read_topology(const char* fname, bool all_caps=false);

    std::vector<topo_defs_residue*> get_residues();
    std::vector<topo_defs_residue*> get_patches();
    std::vector<std::string> get_topology_filenames();

    std::string __str__();
};

#endif

