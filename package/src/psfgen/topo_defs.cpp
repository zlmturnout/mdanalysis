#include "topo_defs.h"
#include "topo_mol.h"
#include "charmm_parse_topo_defs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdexcept>
#include <sstream>

using namespace std;

topo_defs::topo_defs() {

	topo_defs *defs = this;

    defs->auto_angles = 0;
    defs->auto_dihedrals = 0;
    defs->cmaps_present = 0;
    strcpy(defs->pfirst,"");
    strcpy(defs->plast,"");
    defs->buildres = 0;
    defs->buildres_no_errors = 0;
    defs->topo_hash = hasharray_create(
        (void**) &(defs->topo_array), sizeof(topo_defs_topofile_t));
    defs->type_hash = hasharray_create(
        (void**) &(defs->type_array), sizeof(topo_defs_type_t));
    defs->residue_hash = hasharray_create(
        (void**) &(defs->residue_array), sizeof(topo_defs_residue));
    defs->arena = memarena_create();
    
    if ( ! defs->type_hash || ! defs->residue_hash ||
         ! defs->arena || ! defs->topo_hash ||
         defs->build_residue("NONE",1) ||
         defs->build_residue("None",1) ||
         defs->build_residue("none",1) ) {
        //TODO bad  
    }
    topo_defs_end(defs);
}

topo_defs::~topo_defs() {
    topo_defs *defs = this;
	int i,n;
	struct topo_defs_atom_t *a, *a2;
	struct topo_defs_bond_t *b, *b2;
	struct topo_defs_angle_t *an, *an2;
	struct topo_defs_dihedral_t *di, *di2;
	struct topo_defs_improper_t *im, *im2;
	struct topo_defs_cmap_t *cm, *cm2;
	struct topo_defs_conformation_t *c, *c2;

	hasharray_destroy(defs->topo_hash);
	hasharray_destroy(defs->type_hash);
	n = hasharray_count(defs->residue_hash);
	for ( i=0; i<n; ++i ) {
		a = defs->residue_array[i].atoms;
		while ( a ) {
			a2 = a->next;
			free((void*)a);
			a = a2;
		}
		b = defs->residue_array[i].bonds;
		while ( b ) {
			b2 = b->next;
			free((void*)b);
			b = b2;
		}
		an = defs->residue_array[i].angles;
		while ( an ) {
			an2 = an->next;
			free((void*)an);
			an = an2;
		}
		di = defs->residue_array[i].dihedrals;
		while ( di ) {
			di2 = di->next;
			free((void*)di);
			di = di2;
		}
		im = defs->residue_array[i].impropers;
		while ( im ) {
			im2 = im->next;
			free((void*)im);
			im = im2;
		}
		cm = defs->residue_array[i].cmaps;
		while ( cm ) {
			cm2 = cm->next;
			free((void*)cm);
			cm = cm2;
		}
		c = defs->residue_array[i].conformations;
		while ( c ) {
			c2 = c->next;
			free((void*)c);
			c = c2;
		}
	}
	hasharray_destroy(defs->residue_hash);
	memarena_destroy(defs->arena);
}


topo_defs * topo_defs_create(void) {
    return new topo_defs();
}

void topo_defs_destroy(topo_defs *defs) {
    delete defs;
}


/* internal method */
void topo_defs_log_error(topo_defs *defs, const char *msg) {
	PSFGEN_ERROR(msg);
}

void topo_defs_auto_angles(topo_defs *defs, int autogen) {
	if ( defs ) defs->auto_angles = ! ! autogen;
}

void topo_defs_auto_dihedrals(topo_defs *defs, int autogen) {
	if ( defs ) defs->auto_dihedrals = ! ! autogen;
}

int topo_defs_type(topo_defs *defs, const char *atype, const char *element, double mass, int id) {
	int i;
	topo_defs_type_t *newitem;
	char errmsg[64 + NAMEMAXLEN];
	if ( ! defs ) return -1;
	if ( NAMETOOLONG(atype) ) return -2;
	if ( NAMETOOLONG(element) ) return -3;
	if ( ( i = hasharray_index(defs->type_hash,atype) ) != HASHARRAY_FAIL ) {
		sprintf(errmsg,"duplicate type key %s",atype);
		topo_defs_log_error(defs,errmsg);
		newitem = &defs->type_array[i];
	} else {
		i = hasharray_insert(defs->type_hash,atype);
		if ( i == HASHARRAY_FAIL ) return -4;
		newitem = &defs->type_array[i];
		strcpy(newitem->name,atype);
	}
	newitem->id = id;
	strcpy(newitem->element,element);
	newitem->mass = mass;
	return 0;
}

int topo_defs::build_residue(const char *rname, int patch) {
	topo_defs *defs = this;
	int i;
	topo_defs_residue *newitem;
	char errmsg[64 + NAMEMAXLEN];
	if ( ! defs ) return -1;
	defs->buildres = 0;
	defs->buildres_no_errors = 0;
	if ( NAMETOOLONG(rname) ) return -2;
	if ( ( i = hasharray_index(defs->residue_hash,rname) ) != HASHARRAY_FAIL ) {
		sprintf(errmsg,"duplicate residue key %s will be ignored",rname);
		topo_defs_log_error(defs,errmsg);
		/* newitem = &defs->residue_array[i]; */
		defs->buildres_no_errors = 1;
		return 0;
	} else {
		i = hasharray_insert(defs->residue_hash,rname);
		if ( i == HASHARRAY_FAIL ) return -4;
		newitem = &defs->residue_array[i];
		strcpy(newitem->name,rname);
	}
	newitem->patch = patch;
	newitem->atoms = 0;
	newitem->bonds = 0;
	newitem->angles = 0;
	newitem->dihedrals = 0;
	newitem->impropers = 0;
	newitem->cmaps = 0;
	newitem->conformations = 0;
	strcpy(newitem->pfirst,defs->pfirst);
	strcpy(newitem->plast,defs->plast);
	defs->buildres = newitem;
	return 0;
}

int topo_defs_end(topo_defs *defs) {
	if ( ! defs ) return -1;
	defs->buildres = 0;
	defs->buildres_no_errors = 0;
	return 0;
}

int topo_defs_atom(topo_defs *defs, const char *rname, int del,
		const char *aname, int ares, int arel,
		const char *atype, double charge) {
	topo_defs_atom_t *newitem;
	if ( ! defs ) return -1;
	if ( ! defs->buildres ) {
		if ( defs->buildres_no_errors ) return 0;
		topo_defs_log_error(defs,"no residue in progress for atom");
		return -1;
	}
	if ( NAMETOOLONG(aname) ) return -2;
	if ( NAMETOOLONG(atype) ) return -3;
	if ( ares && ! defs->buildres->patch ) return -4;
	if ( arel && ! defs->buildres->patch ) return -4;
	if ( del && ! defs->buildres->patch ) return -5;
	newitem = (topo_defs_atom_t*) malloc(sizeof(topo_defs_atom_t));
	if ( ! newitem )  return -6;
	newitem->res = ares;
	newitem->rel = arel;
	newitem->del = del;
	newitem->charge = charge;
	strcpy(newitem->name,aname);
	strcpy(newitem->type,atype);
	newitem->next = defs->buildres->atoms;
	defs->buildres->atoms = newitem;
	return 0;
}

int topo_defs_bond(topo_defs *defs, const char *rname, int del,
		const char *a1name, int a1res, int a1rel,
		const char *a2name, int a2res, int a2rel) {
	topo_defs_bond_t *newitem;
	if ( ! defs ) return -1;
	if ( ! defs->buildres ) {
		if ( defs->buildres_no_errors ) return 0;
		topo_defs_log_error(defs,"no residue in progress for bond");
		return -1;
	}
	if ( NAMETOOLONG(a1name) ) return -2;
	if ( NAMETOOLONG(a2name) ) return -3;
	if ( del && ! defs->buildres->patch ) return -4;
	if ( ( a1res || a2res ) && ! defs->buildres->patch ) return -4;
	newitem = (topo_defs_bond_t*) malloc(sizeof(topo_defs_bond_t));
	if ( ! newitem )  return -5;
	newitem->res1 = a1res;
	newitem->rel1 = a1rel;
	newitem->res2 = a2res;
	newitem->rel2 = a2rel;
	newitem->del = del;
	strcpy(newitem->atom1,a1name);
	strcpy(newitem->atom2,a2name);
	newitem->next = defs->buildres->bonds;
	defs->buildres->bonds = newitem;
	return 0;
}

int topo_defs_angle(topo_defs *defs, const char *rname, int del,
		const char *a1name, int a1res, int a1rel,
		const char *a2name, int a2res, int a2rel,
		const char *a3name, int a3res, int a3rel) {
	topo_defs_angle_t *newitem;
	if ( ! defs ) return -1;
	if ( ! defs->buildres ) {
		if ( defs->buildres_no_errors ) return 0;
		topo_defs_log_error(defs,"no residue in progress for angle");
		return -1;
	}
	if ( NAMETOOLONG(a1name) ) return -2;
	if ( NAMETOOLONG(a2name) ) return -3;
	if ( NAMETOOLONG(a3name) ) return -4;
	if ( del && ! defs->buildres->patch ) return -5;
	if ( ( a1res || a2res || a3res ) && ! defs->buildres->patch ) return -6;
	newitem = (topo_defs_angle_t*) malloc(sizeof(topo_defs_angle_t));
	if ( ! newitem )  return -7;
	newitem->res1 = a1res;
	newitem->rel1 = a1rel;
	newitem->res2 = a2res;
	newitem->rel2 = a2rel;
	newitem->res3 = a3res;
	newitem->rel3 = a3rel;
	newitem->del = del;
	strcpy(newitem->atom1,a1name);
	strcpy(newitem->atom2,a2name);
	strcpy(newitem->atom3,a3name);
	newitem->next = defs->buildres->angles;
	defs->buildres->angles = newitem;
	return 0;
}

int topo_defs_dihedral(topo_defs *defs, const char *rname, int del,
		const char *a1name, int a1res, int a1rel,
		const char *a2name, int a2res, int a2rel,
		const char *a3name, int a3res, int a3rel,
		const char *a4name, int a4res, int a4rel) {
	topo_defs_dihedral_t *newitem;
	if ( ! defs ) return -1;
	if ( ! defs->buildres ) {
		if ( defs->buildres_no_errors ) return 0;
		topo_defs_log_error(defs,"no residue in progress for dihedral");
		return -1;
	}
	if ( NAMETOOLONG(a1name) ) return -2;
	if ( NAMETOOLONG(a2name) ) return -3;
	if ( NAMETOOLONG(a3name) ) return -4;
	if ( NAMETOOLONG(a4name) ) return -5;
	if ( del && ! defs->buildres->patch ) return -6;
	if ( ( a1res || a2res || a3res || a4res ) &&
			! defs->buildres->patch ) return -7;
	newitem = (topo_defs_dihedral_t*) malloc(sizeof(topo_defs_dihedral_t));
	if ( ! newitem )  return -8;
	newitem->res1 = a1res;
	newitem->rel1 = a1rel;
	newitem->res2 = a2res;
	newitem->rel2 = a2rel;
	newitem->res3 = a3res;
	newitem->rel3 = a3rel;
	newitem->res4 = a4res;
	newitem->rel4 = a4rel;
	newitem->del = del;
	strcpy(newitem->atom1,a1name);
	strcpy(newitem->atom2,a2name);
	strcpy(newitem->atom3,a3name);
	strcpy(newitem->atom4,a4name);
	newitem->next = defs->buildres->dihedrals;
	defs->buildres->dihedrals = newitem;
	return 0;
}

int topo_defs_improper(topo_defs *defs, const char *rname, int del,
		const char *a1name, int a1res, int a1rel,
		const char *a2name, int a2res, int a2rel,
		const char *a3name, int a3res, int a3rel,
		const char *a4name, int a4res, int a4rel) {
	topo_defs_improper_t *newitem;
	if ( ! defs ) return -1;
	if ( ! defs->buildres ) {
		if ( defs->buildres_no_errors ) return 0;
		topo_defs_log_error(defs,"no residue in progress for improper");
		return -1;
	}
	if ( NAMETOOLONG(a1name) ) return -2;
	if ( NAMETOOLONG(a2name) ) return -3;
	if ( NAMETOOLONG(a3name) ) return -4;
	if ( NAMETOOLONG(a4name) ) return -5;
	if ( del && ! defs->buildres->patch ) return -6;
	if ( ( a1res || a2res || a3res || a4res ) &&
			! defs->buildres->patch ) return -7;
	newitem = (topo_defs_improper_t*) malloc(sizeof(topo_defs_improper_t));
	if ( ! newitem )  return -8;
	newitem->res1 = a1res;
	newitem->rel1 = a1rel;
	newitem->res2 = a2res;
	newitem->rel2 = a2rel;
	newitem->res3 = a3res;
	newitem->rel3 = a3rel;
	newitem->res4 = a4res;
	newitem->rel4 = a4rel;
	newitem->del = del;
	strcpy(newitem->atom1,a1name);
	strcpy(newitem->atom2,a2name);
	strcpy(newitem->atom3,a3name);
	strcpy(newitem->atom4,a4name);
	newitem->next = defs->buildres->impropers;
	defs->buildres->impropers = newitem;
	return 0;
}

int topo_defs_cmap(topo_defs *defs, const char *rname, int del,
		const char* const anamel[8], const int aresl[8], const int arell[8]) {
	int i;
	topo_defs_cmap_t *newitem;
	if ( ! defs ) return -1;
	if ( ! defs->buildres ) {
		if ( defs->buildres_no_errors ) return 0;
		topo_defs_log_error(defs,"no residue in progress for cmap");
		return -1;
	}
	for ( i=0; i<8; ++i ) {
		if ( NAMETOOLONG(anamel[i]) ) return -2-i;
	}
	if ( del && ! defs->buildres->patch ) return -10;
	if ( ( aresl[0] || aresl[1] || aresl[2] || aresl[3] ||
			aresl[4] || aresl[5] || aresl[6] || aresl[7] ) &&
			! defs->buildres->patch ) return -11;
	newitem = (topo_defs_cmap_t*) malloc(sizeof(topo_defs_cmap_t));
	if ( ! newitem )  return -12;
	for ( i=0; i<8; ++i ) {
		newitem->resl[i] = aresl[i];
		newitem->rell[i] = arell[i];
		strcpy(newitem->atoml[i],anamel[i]);
	}
	newitem->del = del;
	newitem->next = defs->buildres->cmaps;
	defs->buildres->cmaps = newitem;
	if ( ! defs->cmaps_present ) {
		topo_defs_log_error(defs,"cross-term entries present in topology definitions");
	}
	defs->cmaps_present = 1;
	return 0;
}

int topo_defs_conformation(topo_defs *defs, const char *rname, int del,
		const char *a1name, int a1res, int a1rel,
		const char *a2name, int a2res, int a2rel,
		const char *a3name, int a3res, int a3rel,
		const char *a4name, int a4res, int a4rel,
		double dist12, double angle123, double dihedral, int improper,
		double angle234, double dist34) {
	topo_defs_conformation_t *newitem;
	if ( ! defs ) return -1;
	if ( ! defs->buildres ) {
		if ( defs->buildres_no_errors ) return 0;
		topo_defs_log_error(defs,"no residue in progress for conformation");
		return -1;
	}
	if ( NAMETOOLONG(a1name) ) return -2;
	if ( NAMETOOLONG(a2name) ) return -3;
	if ( NAMETOOLONG(a3name) ) return -4;
	if ( NAMETOOLONG(a4name) ) return -5;
	if ( del && ! defs->buildres->patch ) return -6;
	if ( ( a1res || a2res || a3res || a4res ) &&
			! defs->buildres->patch ) return -7;
	newitem = (topo_defs_conformation_t*)malloc(sizeof(topo_defs_conformation_t));
	if ( ! newitem )  return -8;
	newitem->res1 = a1res;
	newitem->rel1 = a1rel;
	newitem->res2 = a2res;
	newitem->rel2 = a2rel;
	newitem->res3 = a3res;
	newitem->rel3 = a3rel;
	newitem->res4 = a4res;
	newitem->rel4 = a4rel;
	newitem->del = del;
	newitem->improper = improper;
	newitem->dist12 = dist12;
	newitem->angle123 = angle123;
	newitem->dihedral = dihedral;
	newitem->angle234 = angle234;
	newitem->dist34 = dist34;
	strcpy(newitem->atom1,a1name);
	strcpy(newitem->atom2,a2name);
	strcpy(newitem->atom3,a3name);
	strcpy(newitem->atom4,a4name);
	newitem->next = defs->buildres->conformations;
	defs->buildres->conformations = newitem;
	return 0;
}


int topo_defs_default_patching_first(topo_defs *defs, const char *pname) {
	if ( ! defs ) return -1;
	if ( NAMETOOLONG(pname) ) return -2;
	strcpy(defs->pfirst,pname);
	return 0;
}

int topo_defs_default_patching_last(topo_defs *defs, const char *pname) {
	if ( ! defs ) return -1;
	if ( NAMETOOLONG(pname) ) return -2;
	strcpy(defs->plast,pname);
	return 0;
}

int topo_defs_patching_first(topo_defs *defs, const char *rname,
		const char *pname) {
	if ( ! defs ) return -1;
	if ( ! defs->buildres ) {
		if ( defs->buildres_no_errors ) return 0;
		topo_defs_log_error(defs,"no residue in progress for patching");
		return -1;
	}
	if ( NAMETOOLONG(pname) ) return -2;
	strcpy(defs->buildres->pfirst,pname);
	return 0;
}

int topo_defs_patching_last(topo_defs *defs, const char *rname,
		const char *pname) {
	if ( ! defs ) return -1;
	if ( ! defs->buildres ) {
		if ( defs->buildres_no_errors ) return 0;
		topo_defs_log_error(defs,"no residue in progress for patching");
		return -1;
	}
	if ( NAMETOOLONG(pname) ) return -2;
	strcpy(defs->buildres->plast,pname);
	return 0;
}

/* int topo_defs_add_topofile(topo_defs *defs, const char *filename) { */
/*   topo_defs_topofile_t **topofiles; */
/*   topo_defs_topofile_t *topofiletmp; */
/*   if ( ! defs ) return -1; */
/*   if ( strlen(filename)>=256 ) return -2; */
/*   topofiles = &(defs->topofiles); */
/*   topofiletmp = 0; */
/*   topofiletmp = memarena_alloc(defs->arena,sizeof(topo_defs_topofile_t)); */
/*   if ( ! topofiletmp ) return -3; */

/*   strcpy(topofiletmp->filename,filename); */

/*   printf("add_topo %i %s;\n", defs->ntopo, topofiletmp->filename);  */
/*   defs->ntopo++; */
/*   topofiletmp->next = *topofiles; */
/*   *topofiles = topofiletmp; */
/*   return 0; */
/* } */

int topo_defs_add_topofile(topo_defs *defs, const char *filename) {
	/*   topo_defs_topofile_t **topofiles; */
	/*   topo_defs_topofile_t *topofiletmp; */
	/*   if ( ! defs ) return -1; */
	/*   if ( strlen(filename)>=256 ) return -2; */
	/*   topofiles = &(defs->topofiles); */
	/*   topofiletmp = 0; */
	/*   topofiletmp = memarena_alloc(defs->arena,sizeof(topo_defs_topofile_t)); */
	/*   if ( ! topofiletmp ) return -3; */

	/*   strcpy(topofiletmp->filename,filename); */

	/*   printf("add_topo %i %s;\n", defs->ntopo, topofiletmp->filename);  */
	/*   defs->ntopo++; */
	/*   topofiletmp->next = *topofiles; */
	/*   *topofiles = topofiletmp; */
	/*   return 0; */

	int i;
	topo_defs_topofile_t *newitem;
	char errmsg[64 + NAMEMAXLEN];
	if ( ! defs ) return -1;
	if ( strlen(filename)>=256 ) return -2;
	if ( ( i = hasharray_index(defs->type_hash,filename) ) != HASHARRAY_FAIL ) {
		sprintf(errmsg,"duplicate type key %s",filename);
		topo_defs_log_error(defs,errmsg);
		newitem = &defs->topo_array[i];
	} else {
		i = hasharray_insert(defs->topo_hash,filename);
		if ( i == HASHARRAY_FAIL ) return -4;
		newitem = &defs->topo_array[i];
		strcpy(newitem->filename,filename);
	}
	return 0;
}



//int tcl_topology(ClientData data, Tcl_Interp *interp,
//					int argc, CONST84 char *argv[]) {
//  FILE *defs_file;
//  const char *filename;
//  char msg[2048];
//  int itopo,ntopo;
//  psfgen_data *psf = *(psfgen_data **)data;
//  PSFGEN_TEST_MOL(interp,psf);
//
//  if ( argc == 1 ) {
//    Tcl_SetResult(interp,"no topology file specified",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  if ( argc > 2 ) {
//    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }


//  if (argc == 2 && !strcasecmp(argv[1], "residues") ) {
//    psfgen_data *psf = *(psfgen_data **)data;
//    topo_defs *defs = psf->defs;
//    /* Return a list of the known residue definitions */
//    int n = hasharray_count(defs->residue_hash);
//    int i;
//    for (i=0; i<n; i++) {
//      if (!defs->residue_array[i].patch)
//        Tcl_AppendElement(interp, defs->residue_array[i].name);
//    }
//    return TCL_OK;
//  } 
std::vector<topo_defs_residue*> topo_defs::get_residues() {
    int n = hasharray_count(residue_hash);
    std::vector<topo_defs_residue*> result;

    for(int i = 0; i < n; i++) {
    	if (!residue_array[i].patch) {
    		result.push_back(&residue_array[i]);
    	}
    }

    return result;
}

// else if (argc == 2 && !strcasecmp(argv[1], "patches") ) {
//    psfgen_data *psf = *(psfgen_data **)data;
//    topo_defs *defs = psf->defs;
//    /* Return a list of the known residue definitions */
//    int n = hasharray_count(defs->residue_hash);
//    int i;
//    for (i=0; i<n; i++) {
//      if (defs->residue_array[i].patch)
//        Tcl_AppendElement(interp, defs->residue_array[i].name);
//    }
//    return TCL_OK;
//  }
std::vector<topo_defs_residue*> topo_defs::get_patches() {
    int n = hasharray_count(residue_hash);
    std::vector<topo_defs_residue*> result;

    for (int i=0; i<n; i++) {
    	if (residue_array[i].patch)
    		result.push_back(&residue_array[i]);
    }

    return result;
}
 

// else if (argc == 2 && !strcasecmp(argv[1], "list") ) {
//    psfgen_data *psf = *(psfgen_data **)data;
//    topo_defs *defs = psf->mol->defs;
//    topo_defs_topofile_t *topo;
//    ntopo = hasharray_count(defs->topo_hash);
//    for ( itopo=0; itopo<ntopo; ++itopo ) {
//      topo = &(defs->topo_array[itopo]);
//      Tcl_AppendElement(interp, topo->filename);
//    }
//    return TCL_OK;
//  }
std::vector<std::string> topo_defs::get_topology_filenames() {
    topo_defs_topofile_t *topo;
    int ntopo = hasharray_count(topo_hash);
    std::vector<std::string> result;

    for (int itopo=0; itopo<ntopo; ++itopo ) {
        topo = &(topo_array[itopo]);
        result.push_back(topo->filename);
    }

    return result;
}


//  filename = argv[1];
//  if ( ! ( defs_file = fopen(filename,"r") ) ) {
//    sprintf(msg,"ERROR: Unable to open topology file %s\n",filename);
//    Tcl_SetResult(interp,msg,TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  } else {
//    sprintf(msg,"reading topology file %s\n",filename);
//    newhandle_msg(interp,msg);
//    charmm_parse_topo_defs(psf->defs,defs_file,psf->all_caps,interp,newhandle_msg);
//    topo_defs_add_topofile(psf->defs, filename);
//    fclose(defs_file);
//  }
//  return TCL_OK;
int topo_defs::read_topology(const char* filename, bool all_caps) {
	char msg [256];
	int result = -1;
	FILE *defs_file = fopen(filename,"r");
    
	if ( ! defs_file) {
		snprintf(msg,sizeof(msg)/sizeof(char), "ERROR: Unable to open topology file %s\n",filename);
		return result;
    } else {
        snprintf(msg,sizeof(msg)/sizeof(char), "reading topology file %s\n",filename);
        psfgen_msg(0, msg);
        charmm_parse_topo_defs(this, defs_file,all_caps);
        topo_defs_add_topofile(this, filename);
        fclose(defs_file);
        return 0;
    }
}

std::string topo_defs::__str__() {
    stringstream ss;
    int n = hasharray_count(residue_hash);

    ss << "<topo_defs, res count:";
    ss << n;
    ss << ">";
    return ss.str();
}


//}

std::string topo_defs_residue::__str__() {
    stringstream ss;
    ss << "<topo_defs_residue";
    ss << ", name:\'" << name << "\'";
    ss << "\', patch:" << patch;
    ss << ", pfirst:\'" << pfirst << "\'";
    ss << ", plast:\'" << plast << "\'";
    ss << ">";
    return ss.str();
}







