#include "topo_defs.h"
#include "topo_mol.h"
#include "topo_mol_struct.h"
#include "topo_mol.h"
#include "pdb_file.h"
#include "topo_mol_output.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdexcept>
#include <sstream>
#include <algorithm>

int* translate =0;

static int topo_mol_write_pdb(topo_mol *mol, FILE *file);
static int topo_mol_write_psf(topo_mol *mol, FILE *file, int charmmfmt, int nocmap);


// from gromacs:
/* this MUST correspond to the
   t_interaction_function[F_NRE] in gmxlib/ifunc.c */
enum {
  F_BONDS,
  F_G96BONDS,
  F_MORSE,
  F_CUBICBONDS,
  F_CONNBONDS,
  F_HARMONIC,
  F_FENEBONDS,
  F_TABBONDS,
  F_TABBONDSNC,
  F_ANGLES,
  F_G96ANGLES,
  F_CROSS_BOND_BONDS,
  F_CROSS_BOND_ANGLES,
  F_UREY_BRADLEY,
  F_QUARTIC_ANGLES,
  F_TABANGLES,
  F_PDIHS,
  F_RBDIHS,
  F_FOURDIHS,
  F_IDIHS,
  F_PIDIHS,
  F_TABDIHS,
  F_LJ14,
  F_COUL14,
  F_LJC14_Q,
  F_LJC_PAIRS_NB,
  F_LJ,
  F_BHAM,
  F_LJ_LR,
  F_BHAM_LR,
  F_DISPCORR,
  F_COUL_SR,
  F_COUL_LR,
  F_RF_EXCL,
  F_COUL_RECIP,
  F_DPD,
  F_POLARIZATION,
  F_WATER_POL,
  F_THOLE_POL,
  F_POSRES,
  F_DISRES,
  F_DISRESVIOL,
  F_ORIRES,
  F_ORIRESDEV,
  F_ANGRES,
  F_ANGRESZ,
  F_DIHRES,
  F_DIHRESVIOL,
  F_CONSTR,
  F_CONSTRNC,
  F_SETTLE,
  F_VSITE2,
  F_VSITE3,
  F_VSITE3FD,
  F_VSITE3FAD,
  F_VSITE3OUT,
  F_VSITE4FD,
  F_VSITE4FDN,
  F_VSITEN,
  F_COM_PULL,
  F_EQM,
  F_EPOT,
  F_EKIN,
  F_ETOT,
  F_ECONSERVED,
  F_TEMP,
  F_PRES,
  F_DVDL,
  F_DKDL,
  F_DGDL_CON,
  F_NRE		/* This number is for the total number of energies	*/
};


typedef enum {
  d_defaults,
  d_atomtypes,
  d_bondtypes,
  d_constrainttypes,
  d_pairtypes,
  d_angletypes,
  d_dihedraltypes,
  d_nonbond_params,
  d_moleculetype,
  d_atoms,
  d_vsites2,
  d_vsites3,
  d_vsites4,
  d_vsitesn,
  d_bonds,
  d_exclusions,
  d_pairs,
  d_pairs_nb,
  d_angles,
  d_dihedrals,
  d_constraints,
  d_settles,
  d_polarization,
  d_water_polarization,
  d_thole_polarization,
  d_system,
  d_molecules,
  d_position_restraints,
  d_angle_restraints,
  d_angle_restraints_z,
  d_distance_restraints,
  d_orientation_restraints,
  d_dihedral_restraints,
  d_maxdir,
  d_invalid,
  d_none
} directive;


static char *ds[d_maxdir+1] = {
  "defaults",
  "atomtypes",
  "bondtypes",
  "constrainttypes",
  "pairtypes",
  "angletypes",
  "dihedraltypes",
  "nonbond_params",
  /* All the directives above can not appear after moleculetype */
  "moleculetype",
  "atoms",
  "virtual_sites2",
  "virtual_sites3",
  "virtual_sites4",
  "virtual_sitesn",
  "bonds",
  "exclusions",
  "pairs",
  "pairs_nb",
  "angles",
  "dihedrals",
  "constraints",
  "settles",
  "polarization",
  "water_polarization",
  "thole_polarization",
  "system",
  "molecules",
  "position_restraints",
  "angle_restraints",
  "angle_restraints_z",
  "distance_restraints",
  "orientation_restraints",
  "dihedral_restraints",
  "invalid"
  };

enum { ebtsBONDS, ebtsANGLES, ebtsPDIHS, ebtsIDIHS, ebtsEXCLS, ebtsNR };

///////////////////////////////////////////////////////////////
////////////////// FOR REFERENCE: ///////////////////////////
///////////////////////////////////////////////////////////////

  /* these bonded parameters will overwritten be when  *
   * there is a [ bondedtypes ] entry in the .rtp file */
//  bts[0] = 1; /* normal bonds     */
//  bts[1] = 1; /* normal angles    */
//  bts[2] = 1; /* normal dihedrals */
//  bts[3] = 2; /* normal impropers */


  /* Column 5 & 6 aren't really bonded types, but we include
   * them here to avoid introducing a new section:
   * Column 5: 1 means generate all dihedrals, 0 not.
   * Column 6: Number of bonded neighbors to exclude.
   * Coulmn 7: Generate 1,4 interactions between pairs of hydrogens
   * Column 8: Remove impropers over the same bond as a proper dihedral
  */
////////////////////////////////////////////////////////////////////////

typedef struct gpp_atomtype *t_atomtype;

#define PERTURBED(a) (((a).mB != (a).m) || ((a).qB != (a).q) || ((a).typeB != (a).type))

	    /*! \brief Double precision accuracy */
#define GMX_DOUBLE_EPS   1.11022302E-16

	    /*! \brief Maximum double precision value */
#define GMX_DOUBLE_MAX   1.79769313E+308

	    /*! \brief Minimum double precision value */
#define GMX_DOUBLE_MIN   2.22507386E-308

	    /*! \brief Single precision accuracy */
#define GMX_FLOAT_EPS    5.96046448E-08

	    /*! \brief Maximum single precision value */
#define GMX_FLOAT_MAX    3.40282347E+38

	    /*! \brief Minimum single precision value */
#define GMX_FLOAT_MIN    1.17549435E-38

typedef int     	atom_id;	/* To indicate an atoms id         */


/* check kernel/toppush.c when you change these numbers */
#define MAXATOMLIST	5
#define MAXFORCEPARAM	12
#define NR_RBDIHS	6
#define NR_FOURDIHS     4
#define MAXSLEN 32
#define STRLEN 4096
#define NOTSET -12345




typedef float real;


//typedef struct {
//  real m,q;		/* Mass and charge			*/
//  real 		mB,qB;		/* Mass and charge for Free Energy calc */
//  unsigned short type;		/* Atom type				*/
//  unsigned short typeB;		/* Atom type for Free Energy calc	*/
//  int           ptype;		/* Particle type			*/
//  int 		resnr;		/* Residue number			*/
//  int           atomnumber;     /* Atomic Number or NOTSET              */
//  char name[MAXSLEN];
//  char type_name[MAXSLEN];          /* type name                     */
//  unsigned char chain;          /* chain identifier                     */
//  char segname[MAXSLEN];          /* segment name                     */
//  char resname[MAXSLEN];          /* residue name                     */
//} t_atom;
//
//typedef struct {
//  atom_id    a[MAXATOMLIST];	/* The atom list (eg. bonds: particle	*/
//				/* i = a[0] (AI), j = a[1] (AJ))	*/
//  real 	     c[MAXFORCEPARAM];	/* Force parameters (eg. b0 = c[0])	*/
//  char       s[MAXSLEN];        /* A string (instead of parameters),    *
//				 * read from the .rtp file in pdb2gmx   */
//} t_param;
//
//typedef struct {
//  int		nr;		/* The number of bonds in this record 	*/
//  t_param 	*param;		/* Array of parameters (dim: nr)	*/
//} t_params;
//
//typedef struct {
//  int  nr;			/* Number of different groups		*/
//  int  *nm_ind;                 /* Index in the group names             */
//} t_grps;
//
//
//typedef struct {
//  int  type;                    /* PDB record name                      */
//  int  atomnr;                  /* PDB atom number                      */
//  char altloc;                  /* Alternate location indicator         */
//  char atomnm[6];               /* True atom name including spaces      */
//  char pdbresnr[6];             /* PDB res number                       */
//real occup;                   /* Occupancy                            */
//  real bfac;                  /* B-factor                             */
//int bb;
////  bool bAnisotropic;            /* (an)isotropic switch                 */
//  int  uij[6];                  /* Anisotropic B-factor                 */
//} t_pdbinfo;
//
//
//typedef struct {
//  int           nr;             /* Nr of atoms                          */
//  t_atom	*atom;		/* Array of atoms (dim: nr)		*/
//				/* The following entries will not 	*/
//				/* allways be used (nres==0)	 	*/
//  char		***atomname;	/* Array of pointers to atom name	*/
//				/* use: (*(atomname[i]))		*/
//  char		***atomtype;	/* Array of pointers to atom types	*/
//				/* use: (*(atomtype[i]))		*/
//  char		***atomtypeB;	/* Array of pointers to B atom types	*/
//				/* use: (*(atomtypeB[i]))		*/
//  int		nres;		/* Nr of residue names			*/
//  char		***resname; 	/* Array of pointers to residue names 	*/
//				/* use: (*(resname[i]))	       	*/
//  t_pdbinfo     *pdbinfo;       /* PDB Information, such as aniso. Bfac */
//} t_atoms;
//
//typedef struct {
//  int           nr;              /* number of atomtypes                     */
//  real         *radius;         /* GBSA radius for each atomtype        */
//  real         *vol;            /* GBSA efective volume for each atomtype   */
//  real         *surftens;       /* implicit solvent surftens for each atomtype */
//  int          *atomnumber;     /* Atomic number, used for QM/MM */
//} t_atomtypes;
//
//typedef struct {
//  int           nr;		/* The number of atomtypes		*/
//  t_atom	*atom;		/* Array of atoms			*/
//  char          ***atomname;	/* Names of the atomtypes		*/
//  t_param	*nb;		/* Nonbonded force default params	*/
//  int           *bondatomtype;  /* The bond_atomtype for each atomtype  */
//  real          *radius;        /* Radius for GBSA stuff                */
//  real          *vol;           /* Effective volume for GBSA            */
//  real          *surftens;      /* Surface tension with water, for GBSA */
//  int           *atomnumber;    /* Atomic number, used for QM/MM        */
//} gpp_atomtype;

//	typedef struct symbuf {
//	  int bufsize;
//	  char **buf;
//	  struct symbuf *next;
//	} t_symbuf;
//
//	typedef struct
//	{
//	  int      nr;
//	  t_symbuf *symbuf;
//	} t_symtab;
//
//	typedef struct
//	{
//	  char    *name;	/* the name of this function			*/
//	  char    *longname;    /* The name for printing etc.                   */
//	  int     nratoms;	/* nr of atoms needed for this function		*/
//	  int     nrfpA,nrfpB;  /* number of parameters for this function.      */
//	                        /* this corresponds to the number of params in  */
//	                        /* iparams struct! (see idef.h)                 */
//	  /* A and B are for normal and free energy components respectively.    */
//	  unsigned long   flags;        /* Flags (see above)                            */
//	  int     nrnb_ind;     /* index for nrnb (-1 if unknown)               */
//	  t_ifunc *ifunc;	/* the function it self				*/
//	} t_interaction_function;

#define NRFPA(ftype) (interaction_function[(ftype)].nrfpA)
#define NRFPB(ftype) (interaction_function[(ftype)].nrfpB)
#define NRFP(ftype)  (NRFPA(ftype)+NRFPB(ftype))
#define NRAL(ftype) (interaction_function[(ftype)].nratoms)



#define snew(ptr,nelem) (ptr)=save_calloc(#ptr,__FILE__,__LINE__,\
			(nelem),sizeof(*(ptr)))
#define srenew(ptr,nelem) (ptr)=save_realloc(#ptr,__FILE__,__LINE__,\
			(ptr),(nelem),sizeof(*(ptr)))
#define smalloc(ptr,size) (ptr)=save_malloc(#ptr,__FILE__,__LINE__,size)
#define scalloc(ptr,nelem,elsize)\
		(ptr)=save_calloc(#ptr,__FILE__,__LINE__,nelem,elsize)
#define srealloc(ptr,size) (ptr)=save_realloc(#ptr,__FILE__,__LINE__,\
			(ptr),size,1)
#define sfree(ptr) save_free(#ptr,__FILE__,__LINE__,(ptr))


int topo_mol_write_translate(topo_mol *mol, FILE *file, void *v, void (*print_msg)(void *, const char *)) {

	int iseg,nseg,ires,nres,atomid;
	int has_guessed_atoms = 0;
	double x,y,z,o,b;

	int ncgrps;
	topo_mol_segment_t *seg;
	topo_mol_residue_t *res;
	topo_mol_atom *atom;
	topo_mol_bond *bond;
	int nbonds;
	topo_mol_angle *angl;
	int nangls;
	topo_mol_dihedral *dihe;
	int ndihes;
	topo_mol_improper *impr;
	int nimprs;
	topo_mol_chargegroup *cgrp;
	topo_mol_cmap *cmap;
	int ncmaps;
	char buf[128];	int i;


	if ( ! mol ) return -1;

	atomid = 0;
	nbonds = 0;
	nangls = 0;
	ndihes = 0;
	nimprs = 0;
	ncgrps = 0;
	ncmaps = 0;
	int nrestotal = 0;
	nseg = hasharray_count(mol->segment_hash);
	for ( iseg=0; iseg<nseg; ++iseg ) {
		seg = mol->segment_array[iseg];
		if (! seg) continue;
		nres = hasharray_count(seg->residue_hash);
		nrestotal += nres;
		for ( ires=0; ires<nres; ++ires ) {
			res = &(seg->residue_array[ires]);
			for ( atom = res->atoms; atom; atom = atom->next ) {
				atom->atomid = ++atomid;
				for ( bond = atom->bonds; bond;
						bond = topo_mol_bond_next(bond,atom) ) {
					if ( bond->atom[0] == atom && ! bond->del ) {
						++nbonds;
					}
				}
				for ( angl = atom->angles; angl;
						angl = topo_mol_angle_next(angl,atom) ) {
					if ( angl->atom[0] == atom && ! angl->del ) {
						++nangls;
					}
				}
				for ( dihe = atom->dihedrals; dihe;
						dihe = topo_mol_dihedral_next(dihe,atom) ) {
					if ( dihe->atom[0] == atom && ! dihe->del ) {
						++ndihes;
					}
				}
				for ( impr = atom->impropers; impr;
						impr = topo_mol_improper_next(impr,atom) ) {
					if ( impr->atom[0] == atom && ! impr->del ) {
						++nimprs;
					}
				}
				for ( cgrp = atom->chargegroups; cgrp;
						cgrp = topo_mol_chargegroup_next(cgrp,atom) ) {
					if ( cgrp->atom == atom && ! cgrp->del ) {
						++ncgrps;
					}
				}


				for ( cmap = atom->cmaps; cmap;
						cmap = topo_mol_cmap_next(cmap,atom) ) {
					if ( cmap->atom[0] == atom && ! cmap->del ) {
						++ncmaps;
					}
				}
			}
		}
	}



	//////////




	// create array of charge group numbers.
	int* cgnr = (int*) malloc((atomid+1)*sizeof(int));
	int curcgnr = 0;

	for ( iseg=0; iseg<nseg; ++iseg ) {
		seg = mol->segment_array[iseg];
		if (! seg) continue;
		nres = hasharray_count(seg->residue_hash);

		for ( ires=0; ires<nres; ++ires ) {
			res = &(seg->residue_array[ires]);
			for ( atom = res->atoms; atom; atom = atom->next ) {
				for ( cgrp = atom->chargegroups; cgrp; cgrp = topo_mol_chargegroup_next(cgrp,atom) ) {
					if ( cgrp->atom == atom && ! cgrp->del ) {
						curcgnr++;
					}
				}
				cgnr[atom->atomid]= curcgnr;
			}
		}
	}

	// scan all bonds and detect a special pair -> overwrite charge group number (combine)
	for ( iseg=0; iseg<nseg; ++iseg ) {
		seg = mol->segment_array[iseg];
		if (! seg) continue;
		nres = hasharray_count(seg->residue_hash);

		for ( ires=0; ires<nres; ++ires ) {
			res = &(seg->residue_array[ires]);
			for ( atom = res->atoms; atom; atom = atom->next ) {
				for ( bond = atom->bonds; bond; bond = topo_mol_bond_next(bond,atom) ) {
					int bondfound=0;
					if ( bond->atom[0] == atom && ! bond->del ) {
						// known charge group problem:
						if ((strcmp(res->name,"AGLC")==0) || (strcmp(res->name,"BGLC")==0)) {
							// cellulose:
							// PRES GB14
							if ((strcmp(bond->atom[0]->name,"O1")==0) && (strcmp(bond->atom[1]->name,"C4")==0))
								bondfound = 1;
							if ((strcmp(bond->atom[0]->name,"C4")==0) && (strcmp(bond->atom[1]->name,"O1")==0))
								bondfound = 1;
							// PRES BETA
							if ((strcmp(bond->atom[0]->name,"C1")==0) && (strcmp(bond->atom[1]->name,"O1")==0))
								bondfound = 1;
							if ((strcmp(bond->atom[0]->name,"O1")==0) && (strcmp(bond->atom[1]->name,"C4")==0))
								bondfound = 1;
						} else {
							// lignin:
							// PRES BO4R , PRES BO4L
							if ((strcmp(bond->atom[0]->name,"C8")==0) && (strcmp(bond->atom[1]->name,"O4")==0))
								bondfound = 1;
							if ((strcmp(bond->atom[0]->name,"O4")==0) && (strcmp(bond->atom[1]->name,"C8")==0))
								bondfound = 1;
							// PRES  55
							if ((strcmp(bond->atom[0]->name,"C5")==0) && (strcmp(bond->atom[1]->name,"C5")==0))
								bondfound = 1;
							if ((strcmp(bond->atom[0]->name,"C5")==0) && (strcmp(bond->atom[1]->name,"C5")==0))
								bondfound = 1;
							// PRES AO4R , PRES AO4L
							if ((strcmp(bond->atom[0]->name,"C7")==0) && (strcmp(bond->atom[1]->name,"O4")==0))
								bondfound = 1;
							if ((strcmp(bond->atom[0]->name,"O4")==0) && (strcmp(bond->atom[1]->name,"C7")==0))
								bondfound = 1;
							// PRES B5R, B5L
							if ((strcmp(bond->atom[0]->name,"C7")==0) && (strcmp(bond->atom[1]->name,"O4")==0))
								bondfound = 1;
							if ((strcmp(bond->atom[0]->name,"O4")==0) && (strcmp(bond->atom[1]->name,"C7")==0))
								bondfound = 1;
							// PRES B5R, B5L
							if ((strcmp(bond->atom[0]->name,"C8")==0) && (strcmp(bond->atom[1]->name,"C8")==0))
								bondfound = 1;
							if ((strcmp(bond->atom[0]->name,"C5")==0) && (strcmp(bond->atom[1]->name,"C5")==0))
								bondfound = 1;
						}
						// now do the work:
						if (bondfound==0) continue;

						int cgnr_from=cgnr[bond->atom[0]->atomid];
						int cgnr_to=cgnr[bond->atom[1]->atomid];
						if (cgnr_from>cgnr_to) {
							for (i=0;i<atomid;i++) { if (cgnr[i+1]==cgnr_from) cgnr[i+1]=cgnr_to; }			
						} else {
							for (i=0;i<atomid;i++) { if (cgnr[i+1]==cgnr_to) cgnr[i+1]=cgnr_from; }			
						}
						//		        fprintf (file,"%5d %5d %5d \n",bond->atom[0]->atomid,bond->atom[1]->atomid,1);
					}
				}
			}
		}
	}



	// test for non-neutral charge groups
	int maxcgnr=0;
	int distinctchargegroups=0;
	int numcgfixes=0;
	// set maximum charge group number;
	for (i=0;i<atomid;i++) maxcgnr = (cgnr[i]>maxcgnr)? cgnr[i] : maxcgnr; 
	double* qchargegroups = (double*) calloc((maxcgnr+1),sizeof(double));
	int* cmpcg = (int*) calloc((maxcgnr+1),sizeof(int));

	for ( iseg=0; iseg<nseg; ++iseg ) {
		seg = mol->segment_array[iseg];
		if (! seg) continue;
		nres = hasharray_count(seg->residue_hash);
		for ( ires=0; ires<nres; ++ires ) {
			res = &(seg->residue_array[ires]);
			for ( atom = res->atoms; atom; atom = atom->next ) {
				qchargegroups[cgnr[atom->atomid]]+=atom->charge;
			}
		}
	}
	sprintf(buf,"trying to make non-neutral charge groups neutral");
	PSFGEN_INFO(buf);

	for (i=0;i<maxcgnr;i++) {
		if (fabs(qchargegroups[i]) > 4*GMX_DOUBLE_EPS) {
			sprintf(buf,"WARNING!!!!!!!!!!!!!!!!!!!");
			PSFGEN_WARN(buf);
			sprintf(buf,"charge group %d is not neutral: q=%f",i,qchargegroups[i]);
			PSFGEN_WARN(buf);
		}
	}

	// scans the array and count distinct chargegroups. uses temp. array cmpcg
	for ( iseg=0; iseg<nseg; ++iseg ) {
		seg = mol->segment_array[iseg];
		if (! seg) continue;
		nres = hasharray_count(seg->residue_hash);
		for ( ires=0; ires<nres; ++ires ) {
			res = &(seg->residue_array[ires]);
			for ( atom = res->atoms; atom; atom = atom->next ) {
				sprintf(buf,"current charge group: %d",cgnr[atom->atomid]);
				PSFGEN_INFO(buf);
				int found =0;
				for (i=0;i<maxcgnr;i++) {
					if (cmpcg[i]==cgnr[atom->atomid]) found=1;				
				}
				if (found==0) {
					sprintf(buf,"cmpcg[%d]=cgnr[%d]",distinctchargegroups,atom->atomid);
					PSFGEN_INFO(buf);
					cmpcg[distinctchargegroups]=cgnr[atom->atomid];
					distinctchargegroups++;
				}

				qchargegroups[cgnr[atom->atomid]]+=atom->charge;
			}
		}
	}

	int* translate = (int*) calloc(atomid+1,sizeof(int));

	int atomcgcounter=0;
	for (i=0;i<distinctchargegroups;i++) {	
		if (cmpcg[i]==0) continue; // ignore "charge group 0" since this is an empty field
		// current charge group number: cmpcg[i] 

		for ( iseg=0; iseg<nseg; ++iseg ) {
			seg = mol->segment_array[iseg];
			if (! seg) continue;
			nres = hasharray_count(seg->residue_hash);

			for ( ires=0; ires<nres; ++ires ) {
				res = &(seg->residue_array[ires]);
				for ( atom = res->atoms; atom; atom = atom->next ) {
					if (cmpcg[i]==cgnr[atom->atomid]) { // if atom has current charge group
						translate[atom->atomid]=++atomcgcounter;
					}
					// here we will a "translate table:"
					// translate[atom->atomid] gives the converted atomid for a sequence of atoms with a sort on charge number
				}
			}
		}

	}

	// print translate table:
	for ( iseg=0; iseg<nseg; ++iseg ) {
		seg = mol->segment_array[iseg];
		if (! seg) continue;
		nres = hasharray_count(seg->residue_hash);

		for ( ires=0; ires<nres; ++ires ) {
			res = &(seg->residue_array[ires]);
			for ( atom = res->atoms; atom; atom = atom->next ) {

				sprintf(buf,"ATOMID %d -> %d (cg:%d)",atom->atomid,translate[atom->atomid],cgnr[atom->atomid]);
				PSFGEN_INFO(buf);
			}
		}
	}


	////////




	for (i=1;i<atomid+1;i++) {
		fprintf(file,"%d\t%d\n",i,translate[i]);
	}

	free(translate);

	return 0;
}




int topo_mol_write_top(topo_mol *mol, FILE *file) {

	char buf[128];
	int i,j;
	int iseg,nseg,ires,nres,atomid;
	topo_mol_segment_t *seg;
	topo_mol_residue_t *res;
	topo_mol_atom *atom;
	topo_mol_bond *bond;
	int nbonds;
	topo_mol_angle *angl;
	int nangls;
	topo_mol_dihedral *dihe;
	int ndihes;
	topo_mol_improper *impr;
	int nimprs;
	topo_mol_chargegroup *cgrp;
	int ncgrps;
	topo_mol_cmap *cmap;
	int ncmaps;
	int numinline;
	int npres,ipres,ntopo,itopo;
	topo_defs_topofile_t *topo;
	topo_mol_patch_t *patch;
	topo_mol_patchres_t *patchres;
	char defpatch[10];
	fpos_t ntitle_pos, save_pos;
	const char *ntitle_fmt;
	int ntitle_count;
	strcpy(defpatch,"");

	if ( ! mol ) return -1;

	////////////////////////////////////////////////////////////////////////////////
	// Gather infos about the internal representation of the topology
	////////////////////////////////////////////////////////////////////////////////


	atomid = 0;
	nbonds = 0;
	nangls = 0;
	ndihes = 0;
	nimprs = 0;
	ncgrps = 0;
	ncmaps = 0;

	int nrestotal = 0;
	nseg = hasharray_count(mol->segment_hash);
	for ( iseg=0; iseg<nseg; ++iseg ) {
		seg = mol->segment_array[iseg];
		if (! seg) continue;
		nres = hasharray_count(seg->residue_hash);
		nrestotal += nres;
		for ( ires=0; ires<nres; ++ires ) {
			res = &(seg->residue_array[ires]);
			for ( atom = res->atoms; atom; atom = atom->next ) {
				atom->atomid = ++atomid;
				for ( bond = atom->bonds; bond;
						bond = topo_mol_bond_next(bond,atom) ) {
					if ( bond->atom[0] == atom && ! bond->del ) {
						++nbonds;
					}
				}
				for ( angl = atom->angles; angl;
						angl = topo_mol_angle_next(angl,atom) ) {
					if ( angl->atom[0] == atom && ! angl->del ) {
						++nangls;
					}
				}
				for ( dihe = atom->dihedrals; dihe;
						dihe = topo_mol_dihedral_next(dihe,atom) ) {
					if ( dihe->atom[0] == atom && ! dihe->del ) {
						++ndihes;
					}
				}
				for ( impr = atom->impropers; impr;
						impr = topo_mol_improper_next(impr,atom) ) {
					if ( impr->atom[0] == atom && ! impr->del ) {
						++nimprs;
					}
				}
				for ( cgrp = atom->chargegroups; cgrp;
						cgrp = topo_mol_chargegroup_next(cgrp,atom) ) {
					if ( cgrp->atom == atom && ! cgrp->del ) {
						++ncgrps;
					}
				}

				for ( cmap = atom->cmaps; cmap;
						cmap = topo_mol_cmap_next(cmap,atom) ) {
					if ( cmap->atom[0] == atom && ! cmap->del ) {
						++ncmaps;
					}
				}
			}
		}
	}

	sprintf(buf,"Initital statistics of topology (internal):",atomid);
	PSFGEN_INFO(buf);
	sprintf(buf,"total of %d segments",nseg);
	PSFGEN_INFO(buf);
	sprintf(buf,"total of %d residues",nrestotal);
	PSFGEN_INFO(buf);
	sprintf(buf,"total of %d atoms",atomid);
	PSFGEN_INFO(buf);
	sprintf(buf,"total of %d bonds",nbonds);
	PSFGEN_INFO(buf);
	sprintf(buf,"total of %d angles",nangls);
	PSFGEN_INFO(buf);
	sprintf(buf,"total of %d dihedrals",ndihes);
	PSFGEN_INFO(buf);
	sprintf(buf,"total of %d impropers",nimprs);
	PSFGEN_INFO(buf);
	sprintf(buf,"total of %d chargegroups",ncgrps);
	PSFGEN_INFO(buf);

	////////////////////////////////////////////////////////////////////////////////
	// Check for charge groups.
	// If none defined set the charge group ID to unique number
	////////////////////////////////////////////////////////////////////////////////

	if (ncgrps==0) {
		sprintf(buf,"WARNING: No charge groups defined. Will assign a chargegroup for each atom!");
		PSFGEN_INFO(buf);


		long atomic_chargroup_id=0;

		nseg = hasharray_count(mol->segment_hash);
		for ( iseg=0; iseg<nseg; ++iseg ) {
			seg = mol->segment_array[iseg];
			if (! seg) continue;
			nres = hasharray_count(seg->residue_hash);
			nrestotal += nres;
			for ( ires=0; ires<nres; ++ires ) {
				res = &(seg->residue_array[ires]);
				for ( atom = res->atoms; atom; atom = atom->next ) {
					topo_mol_chargegroup *group;

					group = (topo_mol_chargegroup*)memarena_alloc(mol->arena,sizeof(topo_mol_chargegroup));
					group->next = atom->chargegroups;
					group->atom = atom;
					group->del = 0;

					atom->chargegroups = group;
				}
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////////
	// scan all bonds and detect a special pair -> overwrite charge group number (combine)
	// this is a HARDCODED fix
	////////////////////////////////////////////////////////////////////////////////
	int* cgnr = (int*) malloc((atomid+1)*sizeof(int));
	int curcgnr = 0;

	for ( iseg=0; iseg<nseg; ++iseg ) {
		seg = mol->segment_array[iseg];
		if (! seg) continue;
		nres = hasharray_count(seg->residue_hash);

		for ( ires=0; ires<nres; ++ires ) {
			res = &(seg->residue_array[ires]);
			for ( atom = res->atoms; atom; atom = atom->next ) {
				for ( cgrp = atom->chargegroups; cgrp; cgrp = topo_mol_chargegroup_next(cgrp,atom) ) {
					if ( cgrp->atom == atom && ! cgrp->del ) {
						curcgnr++;
					}
				}
				cgnr[atom->atomid]= curcgnr;
				sprintf(buf,"curcgnr: %d",curcgnr);
				PSFGEN_INFO(buf);

			}
		}
	}

	//	for ( iseg=0; iseg<nseg; ++iseg ) {
	//		seg = mol->segment_array[iseg];
	//		if (! seg) continue;
	//		nres = hasharray_count(seg->residue_hash);
	//
	//		for ( ires=0; ires<nres; ++ires ) {
	//			res = &(seg->residue_array[ires]);
	//			for ( atom = res->atoms; atom; atom = atom->next ) {
	//				for ( bond = atom->bonds; bond; bond = topo_mol_bond_next(bond,atom) ) {
	//					int bondfound=0;
	//					if ( bond->atom[0] == atom && ! bond->del ) {
	//						// known charge group problem:
	//						if ((strcmp(res->name,"AGLC")==0) || (strcmp(res->name,"BGLC")==0)) {
	//						// cellulose:
	//						// PRES GB14
	//						if ((strcmp(bond->atom[0]->name,"O1")==0) && (strcmp(bond->atom[1]->name,"C4")==0))
	//							bondfound = 1;
	//						if ((strcmp(bond->atom[0]->name,"C4")==0) && (strcmp(bond->atom[1]->name,"O1")==0))
	//							bondfound = 1;
	//						// PRES BETA
	//						if ((strcmp(bond->atom[0]->name,"C1")==0) && (strcmp(bond->atom[1]->name,"O1")==0))
	//							bondfound = 1;
	//						if ((strcmp(bond->atom[0]->name,"O1")==0) && (strcmp(bond->atom[1]->name,"C4")==0))
	//							bondfound = 1;
	//						} else {
	//						// lignin:
	//						// PRES BO4R , PRES BO4L
	//						if ((strcmp(bond->atom[0]->name,"C8")==0) && (strcmp(bond->atom[1]->name,"O4")==0))
	//							bondfound = 1;
	//						if ((strcmp(bond->atom[0]->name,"O4")==0) && (strcmp(bond->atom[1]->name,"C8")==0))
	//							bondfound = 1;
	//						// PRES  55
	//						if ((strcmp(bond->atom[0]->name,"C5")==0) && (strcmp(bond->atom[1]->name,"C5")==0))
	//							bondfound = 1;
	//						if ((strcmp(bond->atom[0]->name,"C5")==0) && (strcmp(bond->atom[1]->name,"C5")==0))
	//							bondfound = 1;
	//						// PRES AO4R , PRES AO4L
	//						if ((strcmp(bond->atom[0]->name,"C7")==0) && (strcmp(bond->atom[1]->name,"O4")==0))
	//							bondfound = 1;
	//						if ((strcmp(bond->atom[0]->name,"O4")==0) && (strcmp(bond->atom[1]->name,"C7")==0))
	//							bondfound = 1;
	//						// PRES B5R, B5L
	//						if ((strcmp(bond->atom[0]->name,"C7")==0) && (strcmp(bond->atom[1]->name,"O4")==0))
	//							bondfound = 1;
	//						if ((strcmp(bond->atom[0]->name,"O4")==0) && (strcmp(bond->atom[1]->name,"C7")==0))
	//							bondfound = 1;
	//						// PRES B5R, B5L
	//						if ((strcmp(bond->atom[0]->name,"C8")==0) && (strcmp(bond->atom[1]->name,"C8")==0))
	//							bondfound = 1;
	//						if ((strcmp(bond->atom[0]->name,"C5")==0) && (strcmp(bond->atom[1]->name,"C5")==0))
	//							bondfound = 1;
	//						}
	//						// now do the work:
	//						if (bondfound==0) continue;
	//
	//						int cgnr_from=cgnr[bond->atom[0]->atomid];
	//						int cgnr_to=cgnr[bond->atom[1]->atomid];
	//						if (cgnr_from>cgnr_to) {
	//							for (i=0;i<atomid;i++) { if (cgnr[i+1]==cgnr_from) cgnr[i+1]=cgnr_to; }
	//						} else {
	//							for (i=0;i<atomid;i++) { if (cgnr[i+1]==cgnr_to) cgnr[i+1]=cgnr_from; }
	//						}
	//	//		        fprintf (file,"%5d %5d %5d \n",bond->atom[0]->atomid,bond->atom[1]->atomid,1);
	//					}
	//				}
	//			}
	//		}
	//	}

	////////////////////////////////////////////////////////////////////////////////
	// test for non-neutral charge groups
	////////////////////////////////////////////////////////////////////////////////


	// 
	int maxcgnr=0;
	int distinctchargegroups=0;
	int numcgfixes=0;
	// set maximum charge group number;
	for (i=0;i<atomid;i++) maxcgnr = (cgnr[i]>maxcgnr)? cgnr[i] : maxcgnr; 
	double* qchargegroups = (double*) calloc((maxcgnr+1),sizeof(double));
	int* cmpcg = (int*) calloc((maxcgnr+1),sizeof(int));

	for ( iseg=0; iseg<nseg; ++iseg ) {
		seg = mol->segment_array[iseg];
		if (! seg) continue;
		nres = hasharray_count(seg->residue_hash);
		for ( ires=0; ires<nres; ++ires ) {
			res = &(seg->residue_array[ires]);
			for ( atom = res->atoms; atom; atom = atom->next ) {
				qchargegroups[cgnr[atom->atomid]]+=atom->charge;
			}
		}
	}
	sprintf(buf,"trying to make non-neutral charge groups neutral");
	PSFGEN_INFO(buf);

	for (i=0;i<maxcgnr;i++) {
		if (fabs(qchargegroups[i]) > 4*GMX_DOUBLE_EPS) {
			sprintf(buf,"WARNING!!!!!!!!!!!!!!!!!!!");
			PSFGEN_INFO(buf);
			sprintf(buf,"charge group %d is not neutral: q=%f",i,qchargegroups[i]);
			PSFGEN_INFO(buf);
		}
	}

	////////////////////////////////////////////////////////////////////////////////
	// scans the array and count distinct chargegroups. uses temp. array cmpcg
	////////////////////////////////////////////////////////////////////////////////

	// 
	for ( iseg=0; iseg<nseg; ++iseg ) {
		seg = mol->segment_array[iseg];
		if (! seg) continue;
		nres = hasharray_count(seg->residue_hash);
		for ( ires=0; ires<nres; ++ires ) {
			res = &(seg->residue_array[ires]);
			for ( atom = res->atoms; atom; atom = atom->next ) {
				sprintf(buf,"current charge group: %d",cgnr[atom->atomid]); PSFGEN_INFO(buf);					
				int found =0;
				for (i=0;i<maxcgnr;i++) {
					if (cmpcg[i]==cgnr[atom->atomid]) found=1;				
				}
				if (found==0) {
					sprintf(buf,"cmpcg[%d]=cgnr[%d]",distinctchargegroups,atom->atomid); PSFGEN_INFO(buf);
					cmpcg[distinctchargegroups]=cgnr[atom->atomid];
					distinctchargegroups++;
				}

				qchargegroups[cgnr[atom->atomid]]+=atom->charge;
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////////
	// make charge groups consecutive:
	// prepare a translate table which lists atoms in the order of charge groups
	// instead of residue
	// use this translate table when writing the atoms
	////////////////////////////////////////////////////////////////////////////////

	// now prepare a TRANSLATE table. ie. we have to reorganize the order of the atoms written 
	// to the top file in a way that atoms in the same charge group are consequtive 

	// whooohooo ineffienency!
	// sort our cmpcg
	// 	int tmp;
	// 	for (i=0; i<(maxcgnr+1)-1; i++) {
	// 	  for (j=0; j<(maxcgnr+1)-1-i; j++)
	// 		// make 0 goto the end of the array....
	// 		if (cmpcg[j+1]=0) continue;
	// 	    if (cmpcg[j+1] < cmpcg[j]) {  /* compare the two neighbors */
	// 	      tmp = cmpcg[j];         /* swap a[j] and a[j+1]      */
	// 	      cmpcg[j] = cmpcg[j+1];
	// 	      cmpcg[j+1] = tmp;
	// 	  }
	// 	}

	int* translate = (int*) calloc(atomid+1,sizeof(int));

	int atomcgcounter=0;
	for (i=0;i<distinctchargegroups;i++) {	
		if (cmpcg[i]==0) continue; // ignore "charge group 0" since this is an empty field
		// current charge group number: cmpcg[i] 

		for ( iseg=0; iseg<nseg; ++iseg ) {
			seg = mol->segment_array[iseg];
			if (! seg) continue;
			nres = hasharray_count(seg->residue_hash);

			for ( ires=0; ires<nres; ++ires ) {
				res = &(seg->residue_array[ires]);
				for ( atom = res->atoms; atom; atom = atom->next ) {
					if (cmpcg[i]==cgnr[atom->atomid]) { // if atom has current charge group
						translate[atom->atomid]=++atomcgcounter;
					}
					// here we will a "translate table:"
					// translate[atom->atomid] gives the converted atomid for a sequence of atoms with a sort on charge number
				}
			}
		}

	}

	// print translate table:
	for ( iseg=0; iseg<nseg; ++iseg ) {
		seg = mol->segment_array[iseg];
		if (! seg) continue;
		nres = hasharray_count(seg->residue_hash);

		for ( ires=0; ires<nres; ++ires ) {
			res = &(seg->residue_array[ires]);
			for ( atom = res->atoms; atom; atom = atom->next ) {

				sprintf(buf,"ATOMID %d -> %d (cg:%d)",atom->atomid,translate[atom->atomid],cgnr[atom->atomid]); PSFGEN_INFO(buf);
			}
		}
	}

	// general header for charmm force field
	fprintf(file,"; Include forcefield parameters\n");
	fprintf(file,"#include \"ffcharmm27.itp\"\n");
	fprintf(file,"\n");

	int atomentries=0;
	int bondentries=0;
	int pairentries=0;
	int angleentries=0;	
	int dihedralentries=0;		
	int improperentries=0;		
	int chargegroupiter=0;

	for ( iseg=0; iseg<nseg; ++iseg ) {
		seg = mol->segment_array[iseg];
		if (! seg) continue;
		nres = hasharray_count(seg->residue_hash);

		/* top writing starts here! */
		// this section assumes the following data structures exist:
		// FILE *out // filepointer, initalized
		// char *pr  // optional, for print_top_posre
		// char *molname // optional, if not set Protein will be used
		// t_atoms *at // list of atoms
		// int bts[] // ???
		// t_params plist[] // ???
		// t_excls excls[] // ???
		// t_atomtype atype // list-type containing atom types (pointer)
		// int *cgnr // vector of charge group numbers, in sync with atomlist at
		// int nrexcl // ???

		char* molname = seg->segid;
		int nrexcl=3; // 1-3 interactions are excluded from nonbonded
		//  if (at && atype && cgnr) {


		fprintf(file,"[ moleculetype ]\n");
		fprintf(file,"; %-15s %5s\n","Name","nrexcl");
		fprintf(file,"%-15s %5d\n\n",molname?molname:"UnknownMolecule",nrexcl);

		fprintf(file,"[ atoms ]\n");
		fprintf(file,"; %4s %10s %6s %7s%6s %6s %10s %10s %6s %10s %10s\n",
				"nr","type","resnr","residue","atom","cgnr","charge","mass","typeB","chargeB","massB");

		//#define NOREORDER
		double qtot  = 0;

		for (i=1;i<atomid+1;i++) {
			for ( ires=0; ires<nres; ++ires ) {
				res = &(seg->residue_array[ires]);
				chargegroupiter++;
				for ( atom = res->atoms; atom; atom = atom->next ) {
					//#ifdef NOREORDER
					//				if (atom->atomid==i)
					//#else
					if (translate[atom->atomid]==i) {
						//#endif

						atomentries++;
						fprintf(file,"%6d %10s %6d %6s %6s %6d %10g %10g",
								translate[atom->atomid], // nr
								atom->type,
								1,
								//#ifdef NOREORDER
								//					ires+1, // resnr
								//#else
								//					1,
								//#endif
								res->name, // residue
								atom->name, // atom
								//			chargegroupiter, // charge group!
								cgnr[atom->atomid],	//charge group
								atom->charge, // charge
								atom->mass); // mass

						qtot += (double)atom->charge;
						if ( fabs(qtot) < 4*GMX_DOUBLE_EPS ) qtot=0;
						fprintf(file,"   ; qtot %.4g\n",qtot);
					}
				}
			}
		}
	}

	fprintf(file,"\n");

	fprintf(file,"[ bonds ]\n");
	fprintf(file,";  ai    aj funct            c0            c1            c2            c3\n");
	for ( ires=0; ires<nres; ++ires ) {
		res = &(seg->residue_array[ires]);
		for ( atom = res->atoms; atom; atom = atom->next ) {
			for ( bond = atom->bonds; bond; bond = topo_mol_bond_next(bond,atom) ) {
				if ( bond->atom[0] == atom && ! bond->del ) {
					bondentries++;
					fprintf (file,"%5d %5d %5d \n",translate[bond->atom[0]->atomid],translate[bond->atom[1]->atomid],1);
				}
			}	
		}
	}


	fprintf(file,"\n");

	//	top_mol_bond_t* mybonds = malloc(2*nbonds*sizeof(top_mol_bond_t))
	//
	//	/* duplicate all bonds to allow topological search */
	//	for ( iseg=0; iseg<nseg; ++iseg ) {
	//		seg = mol->segment_array[iseg];
	//			if (! seg) continue;
	//			nres = hasharray_count(seg->residue_hash);
	//			for ( ires=0; ires<nres; ++ires ) {
	//			res = &(seg->residue_array[ires]);
	//			for ( atom = res->atoms; atom; atom = atom->next ) {
	//				for ( bond = atom->bonds; bond;
	//			         bond = topo_mol_bond_next(bond,atom) ) {
	//					if ( bond->atom[0] == atom && ! bond->del ) {
	//
	//					}
	//				}
	//			}
	//		}
	//	}

	//	fprintf(file,"[ pairs13 ]\n");
	//	fprintf(file,";  ai    aj funct            c0            c1            c2            c3\n");



	int duplatomid[100]; // static array to test for topological duplicates (case of cyclic rings)
	int duplcounter=0;
	int npairs = 0;
	//	int i;




	// first scan all possible 1-3 pairs and store them in a temporary array
	// this we need as excludes for 1-4 pairs, since 5-ring topologies can produce false 1-4 entries!
	// pairs13 is array on int with length 2*(atomid+1)
	//	int pair13duplatomid[100];
	//	int pair13duplcounter=0;
	int pair13counter=0;
	int* pairs13 = (int*) calloc(12*(atomid+1),sizeof(int)); // just make it big enough...

	for ( ires=0; ires<nres; ++ires ) {
		res = &(seg->residue_array[ires]);
		for ( atom = res->atoms; atom; atom = atom->next ) {
			//			for (i=0;i<100;i++) pair13duplatomid[i]=-1;
			//			pair13duplcounter=0;

			topo_mol_atom *atom2,*atom3;
			topo_mol_bond *bond2;

			for ( bond = atom->bonds; bond; bond = topo_mol_bond_next(bond,atom) ) {

				if ( ! bond->del ) {
					if (bond->atom[0] == atom) atom2 = bond->atom[1]; else atom2 = bond->atom[0];
					if (atom==atom2) continue;

					for ( bond2 = atom2->bonds; bond2; bond2 = topo_mol_bond_next(bond2,atom2) ) {

						if ( ! bond2->del ) {
							if (bond2->atom[0] == atom2) atom3 = bond2->atom[1]; else atom3 = bond2->atom[0];
							if (atom==atom3 || atom2==atom3) continue;


							if (atom->atomid<atom3->atomid) {
								//							int pair13duplcounter=0;
								//							for (i=0;i<100;i++) {
								//	sprintf(buf,"writing bond %d %d",bond->atom[0]->atomid,bond->atom[1]->atomid);
								//						sprintf(buf,"writing pair13 %d with atomid: %d",pair13counter,atomid);
								//						PSFGEN_INFO(buf);

								pairs13[2*pair13counter]=atom->atomid;
								pairs13[2*pair13counter+1]=atom3->atomid;									
								pair13counter++;								
								//							}
								//   							fprintf (file,"%5d %5d %5d \n",atom->atomid,atom3->atomid,1);

								//if (atom3->atomid==pair13duplatomid[i]) pair13duplcounter=1; }
								//							if ((pair13duplcounter==0)) {


								//									pair13duplatomid[pair13duplcounter++] = atom3->atomid;
								//							}	else { sprintf(buf,"WARNING. FOUND DUPLICATE FOR PAIR13"); PSFGEN_INFO(buf);}
							}
						}

					}		
				}
			}
		}
	}
	//	fprintf(file,"\n");


	fprintf(file,"[ pairs ]\n");
	fprintf(file,";  ai    aj funct            c0            c1            c2            c3\n");

	for ( ires=0; ires<nres; ++ires ) {
		res = &(seg->residue_array[ires]);
		for ( atom = res->atoms; atom; atom = atom->next ) {
			for (i=0;i<100;i++) duplatomid[i]=-1;
			duplcounter=0;

			topo_mol_atom *atom2,*atom3,*atom4;
			topo_mol_bond *bond2,*bond3;

			for ( bond = atom->bonds; bond; bond = topo_mol_bond_next(bond,atom) ) {

				if ( ! bond->del ) {
					if (bond->atom[0] == atom) atom2 = bond->atom[1]; else atom2 = bond->atom[0];
					if (atom==atom2) continue;

					for ( bond2 = atom2->bonds; bond2; bond2 = topo_mol_bond_next(bond2,atom2) ) {

						if ( ! bond2->del ) {
							if (bond2->atom[0] == atom2) atom3 = bond2->atom[1]; else atom3 = bond2->atom[0];
							if (atom==atom3 || atom2==atom3) continue;

							for ( bond3 = atom3->bonds; bond3; bond3 = topo_mol_bond_next(bond3,atom3) ) {

								if ( ! bond3->del ) {
									if (bond3->atom[0] == atom3) atom4 = bond3->atom[1]; else atom4 = bond3->atom[0];
									if (atom==atom4 || atom2==atom4 || atom3==atom4) continue;

									//						fprintf(file,"third bond is  %d %d (%s,%s)\n",atom3->atomid,atom4->atomid,atom3->name,atom4->name);

									if (atom->atomid<atom4->atomid) {

										// test for appearance in pairs13.
										int p13found=0;
										for (i=0;i<pair13counter;i++) {
											if (pairs13[2*i]==atom->atomid && pairs13[2*i+1]==atom4->atomid) p13found=1;
										}
										if (p13found==1) continue;

										int dupl=0;
										for (i=0;i<100;i++) {
											if (atom4->atomid==duplatomid[i]) dupl=1;
										}
										if ((dupl==0)) {
											pairentries++;
											// translate atomid
											fprintf (file,"%5d %5d %5d \n",translate[atom->atomid],translate[atom4->atomid],1);
											npairs++;
											duplatomid[duplcounter++] = atom4->atomid;
										}
									}

								}
							}		
						}
					}

				}
			}
		}
	}

	fprintf(file,"\n");

	fprintf(file,"[ angles ]");
	fprintf(file,";  ai    aj    ak funct            c0            c1            c2            c3\n");
	for ( ires=0; ires<nres; ++ires ) {
		res = &(seg->residue_array[ires]);
		for ( atom = res->atoms; atom; atom = atom->next ) {
			for ( angl = atom->angles; angl; angl = topo_mol_angle_next(angl,atom) ) {
				if ( angl->atom[0] == atom && ! angl->del ) {
					angleentries++;
					// translate atomid					
					fprintf (file,"%5d %5d %5d %5d \n",translate[angl->atom[0]->atomid],translate[angl->atom[1]->atomid],translate[angl->atom[2]->atomid],5);
				}
			}
		}
	}

	fprintf(file,"\n");

	fprintf(file,"[ dihedrals ]\n");
	fprintf(file,";  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5\n");
	for ( ires=0; ires<nres; ++ires ) {
		res = &(seg->residue_array[ires]);
		for ( atom = res->atoms; atom; atom = atom->next ) {
			for ( dihe = atom->dihedrals; dihe; dihe = topo_mol_dihedral_next(dihe,atom) ) {
				if ( dihe->atom[0] == atom && ! dihe->del ) {
					dihedralentries++;
					fprintf (file,"%5d %5d  %5d  %5d %5d \n",translate[dihe->atom[0]->atomid],translate[dihe->atom[1]->atomid],translate[dihe->atom[2]->atomid],translate[dihe->atom[3]->atomid],9);
				}
			}
		}
	}
	fprintf(file,"; impropers (there's actually no improper section: function type 2 instead of 9(charmm, dihedral) )\n");
	for ( ires=0; ires<nres; ++ires ) {
		res = &(seg->residue_array[ires]);
		for ( atom = res->atoms; atom; atom = atom->next ) {
			for ( impr = atom->impropers; impr; impr = topo_mol_improper_next(impr,atom) ) {
				if ( impr->atom[0] == atom && ! impr->del ) {
					improperentries++;
					fprintf (file,"%5d %5d  %5d  %5d %5d \n",translate[impr->atom[0]->atomid],translate[impr->atom[1]->atomid],translate[impr->atom[2]->atomid],translate[impr->atom[3]->atomid],2);
				}
			}
		}
	}




	sprintf(buf,"Actually written:");
	PSFGEN_INFO(buf);
	sprintf(buf,"total of %d atoms",atomentries);
	PSFGEN_INFO(buf);
	sprintf(buf,"total of %d bonds",bondentries);
	PSFGEN_INFO(buf);
	sprintf(buf,"total of %d angles",angleentries);
	PSFGEN_INFO(buf);
	sprintf(buf,"total of %d dihedrals",dihedralentries);
	PSFGEN_INFO(buf);
	sprintf(buf,"total of %d impropers",improperentries);
	PSFGEN_INFO(buf);
	sprintf(buf,"total of %d chargegroups",distinctchargegroups);
	PSFGEN_INFO(buf);
	sprintf(buf,"total of %d pairs",pairentries);
	PSFGEN_INFO(buf);


	fprintf(file,"; edit the following sections as you need it:");

	fprintf(file,"; Include Position restraint file                  \n");
	fprintf(file,"#ifdef POSRES                                      \n");
	fprintf(file,"#include \"posre.itp\"                             \n");
	fprintf(file,"#endif                                             \n");
	fprintf(file,"                                                   \n");
	fprintf(file,"; Include water topology                           \n");
	fprintf(file,"#include \"tip3p.itp\"                             \n");
	fprintf(file,"                                                   \n");
	fprintf(file,"#ifdef POSRES_WATER                                \n");
	fprintf(file,"; Position restraint for each water oxygen         \n");
	fprintf(file,"[ position_restraints ]                            \n");
	fprintf(file,";  i funct       fcx        fcy        fcz         \n");
	fprintf(file,"   1    1       1000       1000       1000         \n");
	fprintf(file,"#endif                                             \n");
	fprintf(file,"                                                   \n");
	fprintf(file,"; Include generic topology for ions                \n");
	fprintf(file,"#include \"ions.itp\"                              \n");
	fprintf(file,"                                                   \n");

	fprintf(file,"[ system ]                                         \n");
	fprintf(file,"; Name                                             \n");
	fprintf(file,";Protein                                           \n");
	fprintf(file,"                                                   \n");

	fprintf(file,"[ molecules ]                                      \n");
	fprintf(file,"; Compound        #mols                            \n");
	fprintf(file,";define number of solvent atoms here!              \n");
	fprintf(file,";SOL		    44722                                \n");
	for ( iseg=0; iseg<nseg; ++iseg ) {
		seg = mol->segment_array[iseg];
		if (! seg) continue;
		nres = hasharray_count(seg->residue_hash);
		if (nres>0) {
			fprintf(file,"%-15s\t 1\n",seg->segid);
		}
	}


	free(cmpcg);
	free(qchargegroups);

	fprintf(file,"\n");	

	return 0;
	/* end of function here */

}


int topo_mol_write_pdb(topo_mol *mol, FILE *file) {

	char buf[128];
	int iseg,nseg,ires,nres,atomid;
	int has_guessed_atoms = 0;
	double x,y,z,o,b;
	topo_mol_segment_t *seg;
	topo_mol_residue_t *res;
	topo_mol_atom *atom;

	if ( ! mol ) return -1;

	write_pdb_remark(file,"original generated coordinate pdb file");

	atomid = 0;
	nseg = hasharray_count(mol->segment_hash);
	for ( iseg=0; iseg<nseg; ++iseg ) {
		seg = mol->segment_array[iseg];
		if (! seg) continue;

		if ( strlen(seg->segid) > 4 ) {
			sprintf(buf,
					"truncating segid %s to 4 characters allowed by pdb format",
					seg->segid);
			PSFGEN_WARN(buf);
		}

		nres = hasharray_count(seg->residue_hash);
		for ( ires=0; ires<nres; ++ires ) {
			res = &(seg->residue_array[ires]);
			for ( atom = res->atoms; atom; atom = atom->next ) {
				/* Paranoid: make sure x,y,z,o are set. */
				x = y = z = 0.0; o = -1.0;
				++atomid;
				switch ( atom->xyz_state ) {
				case TOPO_MOL_XYZ_SET:
					x = atom->x;  y = atom->y;  z = atom->z;  o = 1.0;
					break;
				case TOPO_MOL_XYZ_GUESS:
				case TOPO_MOL_XYZ_BADGUESS:
					x = atom->x;  y = atom->y;  z = atom->z;  o = 0.0;
					has_guessed_atoms = 1;
					break;
				default:
					PSFGEN_ERROR("Internal error, atom has invalid state.");
					PSFGEN_ERROR("Treating as void.");
					/* Yes, fall through */
				case TOPO_MOL_XYZ_VOID:
					x = y = z = 0.0;  o = -1.0;
					break;
				}
				b = atom->partition;
				write_pdb_atom(file,atomid,atom->name,res->name,atoi(res->resid),
						"",(float)x,(float)y,(float)z,(float)o,(float)b,res->chain,
						seg->segid,atom->element);
			}
		}
	}

	write_pdb_end(file);
	if (has_guessed_atoms) {
		PSFGEN_INFO("Atoms with guessed coordinates will have occupancy of 0.0.");
	}
	return 0;
}

int topo_mol_write_psf(topo_mol *mol, FILE *file, int charmmfmt, int nocmap) {

	char buf[128];
	int iseg,nseg,ires,nres,atomid;
	topo_mol_segment_t *seg;
	topo_mol_residue_t *res;
	topo_mol_atom *atom;
	topo_mol_bond *bond;
	int nbonds;
	topo_mol_angle *angl;
	int nangls;
	topo_mol_dihedral *dihe;
	int ndihes;
	topo_mol_improper *impr;
	int nimprs;
	topo_mol_cmap *cmap;
	int ncmaps;
	int numinline;
	int npres,ipres,ntopo,itopo;
	topo_defs_topofile_t *topo;
	topo_mol_patch_t *patch;
	topo_mol_patchres_t *patchres;
	char defpatch[10];
	fpos_t ntitle_pos, save_pos;
	const char *ntitle_fmt;
	int ntitle_count;
	strcpy(defpatch,"");

	if ( ! mol ) return -1;

	atomid = 0;
	nbonds = 0;
	nangls = 0;
	ndihes = 0;
	nimprs = 0;
	ncmaps = 0;
	nseg = hasharray_count(mol->segment_hash);
	for ( iseg=0; iseg<nseg; ++iseg ) {
		seg = mol->segment_array[iseg];
		if (! seg) continue;

		if ( strlen(seg->segid) > 4 ) {
			sprintf(buf,
					"segid %s is longer than 4 characters allowed by psf format",
					seg->segid);
			PSFGEN_ERROR(buf);
			return -2;
		}

		nres = hasharray_count(seg->residue_hash);
		for ( ires=0; ires<nres; ++ires ) {
			res = &(seg->residue_array[ires]);
			for ( atom = res->atoms; atom; atom = atom->next ) {
				atom->atomid = ++atomid;
				for ( bond = atom->bonds; bond;
						bond = topo_mol_bond_next(bond,atom) ) {
					if ( bond->atom[0] == atom && ! bond->del ) {
						++nbonds;
					}
				}
				for ( angl = atom->angles; angl;
						angl = topo_mol_angle_next(angl,atom) ) {
					if ( angl->atom[0] == atom && ! angl->del ) {
						++nangls;
					}
				}
				for ( dihe = atom->dihedrals; dihe;
						dihe = topo_mol_dihedral_next(dihe,atom) ) {
					if ( dihe->atom[0] == atom && ! dihe->del ) {
						++ndihes;
					}
				}
				for ( impr = atom->impropers; impr;
						impr = topo_mol_improper_next(impr,atom) ) {
					if ( impr->atom[0] == atom && ! impr->del ) {
						++nimprs;
					}
				}
				for ( cmap = atom->cmaps; cmap;
						cmap = topo_mol_cmap_next(cmap,atom) ) {
					if ( cmap->atom[0] == atom && ! cmap->del ) {
						++ncmaps;
					}
				}
			}
		}
	}
	sprintf(buf,"total of %d atoms",atomid);
	PSFGEN_INFO(buf);
	sprintf(buf,"total of %d bonds",nbonds);
	PSFGEN_INFO(buf);
	sprintf(buf,"total of %d angles",nangls);
	PSFGEN_INFO(buf);
	sprintf(buf,"total of %d dihedrals",ndihes);
	PSFGEN_INFO(buf);
	sprintf(buf,"total of %d impropers",nimprs);
	PSFGEN_INFO(buf);

	ntitle_fmt = "PSF\n\n%8d !NTITLE\n";
	if ( nocmap ) {
		sprintf(buf,"total of %d cross-terms (not written to file)",ncmaps);
	} else {
		sprintf(buf,"total of %d cross-terms",ncmaps);
		if ( ncmaps ) {
			ntitle_fmt = "PSF CMAP\n\n%8d !NTITLE\n";
		} else {
			nocmap = 1;
		}
	}
	PSFGEN_INFO(buf);

	fgetpos(file,&ntitle_pos);
	fprintf(file,ntitle_fmt,1);
	ntitle_count = 1;
	if ( charmmfmt )
		fprintf(file," REMARKS %s\n","original generated structure charmm psf file");
	else
		fprintf(file," REMARKS %s\n","original generated structure x-plor psf file");

	if (mol->npatch) {
		ntitle_count++;
		fprintf(file," REMARKS %i patches were applied to the molecule.\n", mol->npatch);
	}

	ntopo = hasharray_count(mol->defs->topo_hash);
	for ( itopo=0; itopo<ntopo; ++itopo ) {
		topo = &(mol->defs->topo_array[itopo]);
		ntitle_count++;
		fprintf(file," REMARKS topology %s \n", topo->filename);
	}

	nseg = hasharray_count(mol->segment_hash);
	for ( iseg=0; iseg<nseg; ++iseg ) {
		char angles[20], diheds[20];
		seg = mol->segment_array[iseg];
		if (! seg) continue;
		strcpy(angles,"none");
		strcpy(diheds,"");
		if (seg->auto_angles)    strcpy(angles,"angles");
		if (seg->auto_dihedrals) strcpy(diheds,"dihedrals");
		ntitle_count++;
		fprintf(file," REMARKS segment %s { first %s; last %s; auto %s %s }\n", seg->segid, seg->pfirst, seg->plast, angles, diheds);
	}

	for ( patch = mol->patches; patch; patch = patch->next ) {
		strcpy(defpatch,"");
		if (patch->deflt) strcpy(defpatch,"default");
		npres = patch->npres;
		ipres = 0;
		for ( patchres = patch->patchresids; patchres; patchres = patchres->next ) {
			/* Test the existence of segid:resid for the patch */
			if (!topo_mol_validate_patchres(mol,patch->pname,patchres->segid, patchres->resid)) {
				break;
			}
		}
		if ( patchres ) continue;

		for ( patchres = patch->patchresids; patchres; patchres = patchres->next ) {
			if (ipres==0) {
				ntitle_count++;
				fprintf(file," REMARKS %spatch %s ", defpatch, patch->pname);
			}
			if (ipres>0 && !ipres%6) {
				ntitle_count++;
				fprintf(file,"\n REMARKS patch ---- ");
			}
			fprintf(file,"%s:%s  ", patchres->segid, patchres->resid);
			if (ipres==npres-1) fprintf(file,"\n");
			ipres++;
		}
	}
	fprintf(file,"\n");
	fgetpos(file,&save_pos);
	fsetpos(file,&ntitle_pos);
	fprintf(file,ntitle_fmt,ntitle_count);
	fsetpos(file,&save_pos);

	fprintf(file,"%8d !NATOM\n",atomid);
	for ( iseg=0; iseg<nseg; ++iseg ) {
		seg = mol->segment_array[iseg];
		if (! seg) continue;
		nres = hasharray_count(seg->residue_hash);
		for ( ires=0; ires<nres; ++ires ) {
			res = &(seg->residue_array[ires]);
			if ( charmmfmt ) for ( atom = res->atoms; atom; atom = atom->next ) {
				int idef,type_id;
				idef = hasharray_index(mol->defs->type_hash,atom->type);
				if ( idef == HASHARRAY_FAIL ) {
					sprintf(buf,"unknown atom type %s",atom->type);
					PSFGEN_INFO(buf);
					return -3;
				}
				type_id = mol->defs->type_array[idef].id;
				fprintf(file,"%8d %-4s %-4s %-4s %-4s %4d %10.6f     %9.4f  %10d\n",
						atom->atomid, seg->segid,res->resid,res->name,
						atom->name,type_id,atom->charge,atom->mass,0);
			} else for ( atom = res->atoms; atom; atom = atom->next ) {
				fprintf(file,"%8d %-4s %-4s %-4s %-4s %-4s %10.6f     %9.4f  %10d\n",
						atom->atomid, seg->segid,res->resid,res->name,
						atom->name,atom->type,atom->charge,atom->mass,0);
			}
		}
	}
	fprintf(file,"\n");

	fprintf(file,"%8d !NBOND: bonds\n",nbonds);
	numinline = 0;
	for ( iseg=0; iseg<nseg; ++iseg ) {
		seg = mol->segment_array[iseg];
		if (! seg) continue;
		nres = hasharray_count(seg->residue_hash);
		for ( ires=0; ires<nres; ++ires ) {
			res = &(seg->residue_array[ires]);
			for ( atom = res->atoms; atom; atom = atom->next ) {
				for ( bond = atom->bonds; bond;
						bond = topo_mol_bond_next(bond,atom) ) {
					if ( bond->atom[0] == atom && ! bond->del ) {
						if ( numinline == 4 ) { fprintf(file,"\n");  numinline = 0; }
						fprintf(file," %7d %7d",atom->atomid,bond->atom[1]->atomid);
						++numinline;
					}
				}
			}
		}
	}
	fprintf(file,"\n\n");

	fprintf(file,"%8d !NTHETA: angles\n",nangls);
	numinline = 0;
	for ( iseg=0; iseg<nseg; ++iseg ) {
		seg = mol->segment_array[iseg];
		if (! seg) continue;
		nres = hasharray_count(seg->residue_hash);
		for ( ires=0; ires<nres; ++ires ) {
			res = &(seg->residue_array[ires]);
			for ( atom = res->atoms; atom; atom = atom->next ) {
				for ( angl = atom->angles; angl;
						angl = topo_mol_angle_next(angl,atom) ) {
					if ( angl->atom[0] == atom && ! angl->del ) {
						if ( numinline == 3 ) { fprintf(file,"\n");  numinline = 0; }
						fprintf(file," %7d %7d %7d",atom->atomid,
								angl->atom[1]->atomid,angl->atom[2]->atomid);
						++numinline;
					}
				}
			}
		}
	}
	fprintf(file,"\n\n");

	fprintf(file,"%8d !NPHI: dihedrals\n",ndihes);
	numinline = 0;
	for ( iseg=0; iseg<nseg; ++iseg ) {
		seg = mol->segment_array[iseg];
		if (! seg) continue;
		nres = hasharray_count(seg->residue_hash);
		for ( ires=0; ires<nres; ++ires ) {
			res = &(seg->residue_array[ires]);
			for ( atom = res->atoms; atom; atom = atom->next ) {
				for ( dihe = atom->dihedrals; dihe;
						dihe = topo_mol_dihedral_next(dihe,atom) ) {
					if ( dihe->atom[0] == atom && ! dihe->del ) {
						if ( numinline == 2 ) { fprintf(file,"\n");  numinline = 0; }
						fprintf(file," %7d %7d %7d %7d",atom->atomid,
								dihe->atom[1]->atomid,dihe->atom[2]->atomid,
								dihe->atom[3]->atomid);
						++numinline;
					}
				}
			}
		}
	}
	fprintf(file,"\n\n");

	fprintf(file,"%8d !NIMPHI: impropers\n",nimprs);
	numinline = 0;
	for ( iseg=0; iseg<nseg; ++iseg ) {
		seg = mol->segment_array[iseg];
		if (! seg) continue;
		nres = hasharray_count(seg->residue_hash);
		for ( ires=0; ires<nres; ++ires ) {
			res = &(seg->residue_array[ires]);
			for ( atom = res->atoms; atom; atom = atom->next ) {
				for ( impr = atom->impropers; impr;
						impr = topo_mol_improper_next(impr,atom) ) {
					if ( impr->atom[0] == atom && ! impr->del ) {
						if ( numinline == 2 ) { fprintf(file,"\n");  numinline = 0; }
						fprintf(file," %7d %7d %7d %7d",atom->atomid,
								impr->atom[1]->atomid,impr->atom[2]->atomid,
								impr->atom[3]->atomid);
						++numinline;
					}
				}
			}
		}
	}
	fprintf(file,"\n\n");

	fprintf(file,"%8d !NDON: donors\n\n\n",0);
	fprintf(file,"%8d !NACC: acceptors\n\n\n",0);
	fprintf(file,"%8d !NNB\n\n",0);
	/* Pad with zeros, one for every atom */
	{
		int i, fullrows;
		fullrows = atomid/8;
		for (i=0; i<fullrows; ++i)
			fprintf(file, "%8d%8d%8d%8d%8d%8d%8d%8d\n",0,0,0,0,0,0,0,0);
		for (i=atomid - fullrows*8; i; --i)
			fprintf(file, "%8d",0);
	}
	fprintf(file,"\n\n");

	fprintf(file,"%8d %7d !NGRP\n%8d%8d%8d\n\n",1,0,0,0,0);

	if ( ! nocmap ) {
		fprintf(file,"%8d !NCRTERM: cross-terms\n",ncmaps);
		for ( iseg=0; iseg<nseg; ++iseg ) {
			seg = mol->segment_array[iseg];
			if (! seg) continue;
			nres = hasharray_count(seg->residue_hash);
			for ( ires=0; ires<nres; ++ires ) {
				res = &(seg->residue_array[ires]);
				for ( atom = res->atoms; atom; atom = atom->next ) {
					for ( cmap = atom->cmaps; cmap;
							cmap = topo_mol_cmap_next(cmap,atom) ) {
						if ( cmap->atom[0] == atom && ! cmap->del ) {
							fprintf(file," %7d %7d %7d %7d %7d %7d %7d %7d\n",atom->atomid,
									cmap->atom[1]->atomid,cmap->atom[2]->atomid,
									cmap->atom[3]->atomid,cmap->atom[4]->atomid,
									cmap->atom[5]->atomid,cmap->atom[6]->atomid,
									cmap->atom[7]->atomid);
						}
					}
				}
			}
		}
		fprintf(file,"\n");
	}

	return 0;
}

static bool ends_with(const std::string& str, const std::string& end) {
    int pos = str.rfind(end);
    return  pos != std::string::npos && pos == (str.length() - end.length());
}


//  FILE *res_file;
//  const char *filename;
//  int charmmfmt, nocmap, i;
//  char msg[2048];
//  psfgen_data *psf = *(psfgen_data **)data;
//  PSFGEN_TEST_MOL(interp,psf);
//
//  if ( argc == 1 ) {
//    Tcl_SetResult(interp,"no psf file specified",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  if ( argc > 4 ) {
//    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  filename = argv[argc-1];
//
//
//  return TCL_OK;



/**
 * write an output file. The type of file is automatically determinied by the file name, or can be specified in
 * the optional file_type parameter.
 *
 * @param file_name: the file name, the type is guessed from the file extension, currently
 * supported files are .pdb and .psf
 * @param file_type: the file type, currently, either "pdb", "psf", or leave blank to have the file
 * type automatically determined from the file name.
 */
int topo_mol::write(const std::string& file_name, const std::string& file_type) {
	int result = -1;
    FILE *file = 0;
    const int npos = std::string::npos;
    std::string ufname  = file_name;
    std::string uftype = file_type;
    std::transform(ufname.begin(), ufname.end(),ufname.begin(), ::toupper);
    std::transform(uftype.begin(), uftype.end(),uftype.begin(), ::toupper);

    {
        std::stringstream ss;
        ss << "file_name: \"" << file_name << "\", file_type: \"" << file_type;
        PSFGEN_DEBUG(ss.str().c_str());
    }

	if(ends_with(ufname, ".PSF") || (uftype.find("PSF") != npos)) {

        int charmmfmt = 0;
        int nocmap = 0;
        
        if(uftype.find("CHARMM") != npos) charmmfmt = 1;
        if(uftype.find("X-PLOR") != npos) charmmfmt = 0;
        if(uftype.find("CMAP")   != npos) nocmap = 0;
        if(uftype.find("NOCMAP") != npos) nocmap = 1;

        {
            std::stringstream ss;
            ss << "writing psf file \"" << file_name << "\", ";;
            ss << (nocmap ? " without cross-terms, " : "with cross terms, ");
            ss << (charmmfmt ? " in CHARMM format" : "in x-plor format");
            PSFGEN_INFO(ss.str().c_str());
        }

        if ( ! ( file = fopen(file_name.c_str(),"w"))) {
            std::string msg = "Unable to open psf file ";
            msg += file_name; 
            msg += " for writing";
            PSFGEN_ERROR(msg.c_str());
        }
        if (result = topo_mol_write_psf(this,file,charmmfmt,nocmap)) {
            PSFGEN_ERROR("failed on writing structure to psf file");
        }
        fclose(file);
        PSFGEN_INFO("psf file complete.");
    } 
    else if(ends_with(ufname, ".PDB") || uftype == "PDB") {
        PSFGEN_INFO("writing pdb output");

        if ( ! ( file = fopen(file_name.c_str(),"w"))) {
            std::string msg = "Unable to open pdb file ";
            msg += file_name; 
            msg += " for writing";
            PSFGEN_ERROR(msg.c_str());
        }
        if (result = topo_mol_write_pdb(this,file)) {
            PSFGEN_ERROR("failed on writing coordinates to pdb file");
        }
        fclose(file);
        PSFGEN_INFO("pdb file complete.");
    }
    else {
        std::stringstream ss;
        ss << "Could not determine file type for file name \"" << file_name;
        ss << ", or given file_type \"" << file_type << ", currently only .psf and .pdb are supported";
        PSFGEN_ERROR(ss.str().c_str());
    }

    return result;

}


