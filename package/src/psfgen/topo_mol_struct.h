#ifndef TOPO_DEFS_MOL_H
#define TOPO_DEFS_MOL_H

#ifndef SWIG
#include <Python.h>
#include <vector>
#endif

#include "hasharray.h"
#include "memarena.h"
#include "topo_defs_struct.h"


#define NAMEMAXLEN 8

// TODO change to strnlen, 
// BSDs don't have this.
#define NAMETOOLONG(X) ( strlen(X) >= NAMEMAXLEN )

#ifdef SWIG
%feature("python:slot", "tp_str",  functype="reprfunc") topo_mol_bond::str;
%feature("python:slot", "tp_repr",  functype="reprfunc") topo_mol_bond::repr;
#endif

struct topo_mol_atom;

#ifndef SWIG
struct topo_mol_ident_t {
	const char *segid;
	const char *resid;
	const char *aname;
};
#endif

struct topo_mol_bond {
	struct topo_mol_bond *next[2];
	struct topo_mol_atom *atom[2];

	const char* __str__() {
		return "I'm a str";
	}

	const char* __repr__() {
		return "I'm a repr";
	}

#ifdef SWIG
	%rename(_del) del;
#endif

	int del;
};

struct topo_mol_angle {
	struct topo_mol_angle *next[3];
	struct topo_mol_atom *atom[3];

    //topo_mol_angle * next(topo_mol_atom *atom);


#ifdef SWIG
	%rename(_del) del;
#endif
	int del;
};

struct topo_mol_dihedral {
	struct topo_mol_dihedral *next[4];
	struct topo_mol_atom *atom[4];

#ifdef SWIG
	%rename(_del) del;
#endif
	int del;
};

struct topo_mol_improper {
	struct topo_mol_improper *next[4];
	struct topo_mol_atom *atom[4];

#ifdef SWIG
	%rename(_del) del;
#endif
	int del;
};

struct topo_mol_chargegroup {
	struct topo_mol_chargegroup *next;
	struct topo_mol_atom *atom;

#ifdef SWIG
	%rename(_del) del;
#endif
	int del;
};

struct topo_mol_cmap {
	struct topo_mol_cmap *next[8];
	struct topo_mol_atom *atom[8];

#ifdef SWIG
	%rename(_del) del;
#endif
	int del;
};

struct topo_mol_conformation {
	struct topo_mol_conformation *next[4];
	struct topo_mol_atom *atom[4];

#ifdef SWIG
	%rename(_del) del;
#endif
	int del;
	int improper;
	double dist12, angle123, dihedral, angle234, dist34;
};

#define TOPO_MOL_XYZ_VOID 0
#define TOPO_MOL_XYZ_SET 1
#define TOPO_MOL_XYZ_GUESS 2
#define TOPO_MOL_XYZ_BADGUESS 3

#ifdef SWIG
%feature("python:slot", "tp_str",  functype="reprfunc") topo_mol_atom::__str__;
%feature("python:slot", "tp_repr",  functype="reprfunc") topo_mol_atom::__str__;
#endif

struct topo_mol_atom {
	struct topo_mol_atom *next;
	struct topo_mol_atom *copy;
	topo_mol_bond *bonds;
	topo_mol_angle *angles;
	topo_mol_dihedral *dihedrals;
	topo_mol_improper *impropers;
	topo_mol_chargegroup *chargegroups;
	topo_mol_cmap *cmaps;
	topo_mol_conformation *conformations;
	char name[NAMEMAXLEN];
	char type[NAMEMAXLEN];
	char element[NAMEMAXLEN];
	double mass;
	double charge;
	double x,y,z;
	int xyz_state;
	int partition;
	int atomid;
    std::string __str__();
};


#ifdef SWIG
%feature("python:slot", "tp_str",  functype="reprfunc") topo_mol_residue_t::__str__;
%feature("python:slot", "tp_repr",  functype="reprfunc") topo_mol_residue_t::__str__;
#endif

struct topo_mol_residue_t {
	char resid[NAMEMAXLEN];
	char name[NAMEMAXLEN];
	char chain[NAMEMAXLEN];

    /**
     * get the list of atoms that belong to this residue.
     * @see topo_mol_atom
     * @returns: a list of atoms. 
     */
    std::vector<topo_mol_atom*> get_atoms();

    std::string __str__();

#ifndef SWIG
	topo_mol_atom *atoms;
#endif
    
};


#ifdef SWIG
%feature("python:slot", "tp_str",  functype="reprfunc") topo_mol_segment_t::__str__;
%feature("python:slot", "tp_repr",  functype="reprfunc") topo_mol_segment_t::__str__;
#endif

struct topo_mol_segment_t {

    /**
     * Segment ID, a numeric string, 1-4 characters.
     */
	char segid[NAMEMAXLEN];

    /**
     * Enable generation of angles from bonds.
     */
	int auto_angles;

    /**
     * Enable generation of dihedrals from angles.
     */
	int auto_dihedrals;

    /**
     * The name of the patch that was applied to the beginning of the  segment. 
     */
	char pfirst[NAMEMAXLEN];

    /**
     * The name of the patch that was applied to the end of the segment.
     */
	char plast[NAMEMAXLEN];

    /**
     * get a list of residues that belong to this segment.
     */
    std::vector<topo_mol_residue_t*> get_residues();
    std::string __str__();

#ifndef SWIG
	topo_mol_residue_t *residue_array;
    hasharray *residue_hash;
#endif
};

struct topo_mol_patchres_t {
	struct topo_mol_patchres_t *next;
	char segid[NAMEMAXLEN];
	char resid[NAMEMAXLEN];
};

struct topo_mol_patch_t {
	struct topo_mol_patch_t *next;
	char pname[NAMEMAXLEN];
	int npres;
	int deflt;
	topo_mol_patchres_t *patchresids;
};

class topo_mol {
public:

	topo_defs *defs;
	bool all_caps;

#ifndef SWIG
	int npatch;
	topo_mol_patch_t *patches;
	topo_mol_patch_t *curpatch;
	hasharray *segment_hash;
	memarena *arena;
	topo_mol_segment_t **segment_array;
#endif


	/**
	 * the segment currently being built.
	 *
	 * it's possible to modify fields of this variable.
	 *
	 * To override the default patch applied to first residue in segment, set the
	 * 'pfirst' field. The default is read from topology file and may be residue-specific.
	 *
	 * To override the default patch applied to last residue in segment, set the
	 * 'plast' filed.
	 */
	topo_mol_segment_t *buildseg;

	/**
	 * creates a new topo_mol object.
	 * does stuff.
	 */
	topo_mol();

	topo_mol(topo_defs *defs);
	~topo_mol();

	/**
	 * Get all the segments that belong to this molcule.
	 * Do not attempt to append to this list, to create a new segment,
	 * use the build_segment method.
	 * @return: a list of segments.
	 */
	std::vector<topo_mol_segment_t*> get_segments();

	/**
	 * Read a topology file.
	 * the patches of the topology file are added to the definitions and made available for use.
	 * This may be called multiple times, as the contents are simply added.
	 * @param filename: the name of the topo file
	 * @returns 0 on success, something else otherwise.
	 */
	int read_topology(const char* filename);

	/**
	 * read  a psf file.
	 * @param filename: the name of a file
	 * @returns: 0 on success, something else on failure.
	 */
	int read_psf(const char* filename);

    /**
     * Extract sequence information (residues) from PDB file when building segment. 
     * Residue IDs will be preserved, residue names must match entries in the 
     * topology file or should be aliased before pdb is called. 
     * @param filename: PDB file containing known or aliased residues. 
     */
	int read_residues(const char* filename);

    /**
     * Read coordinates from file (currently PDB), matching segment, residue and atom names. 
     * @param filename: PDB file containing known or aliased residues and atoms. 
     * @param segid: If specified override segment IDs in PDB file. 
     */
    int read_coordinates(const char* filename, const char* segid = 0);

    /**
	 * Currently not exactly sure what this does,
	 * something like counting the number of times a patch can be applied???
	 * @param presname: the name of a patch residue.
	 */
	int num_patch_targets(const char *presname);

	/**
	 * Apply a patch to the molecule.
	 * @param patch_res_name: the residiue name of the patch
	 * @param target: the segid/resid of the residue to apply the path to.
	 * @param prepend: optional, defaults to 0
	 * @param warn_angles: optional, defaults to 0
	 * @param warn_dihedrals: optional, defaults to 0
	 * @param deflt: optional, defaults to 0
	 */
	int patch(const char *patch_res_name, PyObject *target,
			int prepend = 0, int warn_angles = 0, int warn_dihedrals = 0, int deflt = 0);


	/**
	 * Select a set of atoms. This performs the same selections as delete_atom, except
	 * it returns the atoms instead of deletes them.
     *
     * @param target: a Dictionary with the keys {'segid', 'resid', 'aname'} where these
     * corespond to the segment id, residue id, and atom name. resid and aname may be null.
     *
     * @returns a list of topo_mol_atom objects.
     *
	 */
	std::vector<topo_mol_atom*> select_atom(const topo_mol_ident_t *target);

	/**
	 * write an output file. The type of file is automatically determinied by the file name, or can be specified in
	 * the optional file_type parameter.
	 *
	 * @param file_name: the file name, the type is guessed from the file extension, currently
	 * supported files are .pdb and .psf
	 * @param file_type: the file type, currently, either "pdb", "psf", or leave blank to have the file
	 * type automatically determined from the file name.
	 */
	int write(const std::string& file_name, const std::string& file_type="");

	/**
	 * delete one or more atoms based on the identifier dictionary.
	 * if the identifier only contains a segid, all atoms in the segment are deleted,
	 * if the identifier contains a segid and resid, all atoms in the resid are deleted, and
	 * if the identifier contans a segid, resid and atomid, a specific atom is deleted.
     *
     * @see topo_mol::select_atom
	 * @param mol: a topo_mol object
	 * @param target: a dictionary with the keys 'segid', 'resid', 'aname'";
	 */
	void delete_atom(topo_mol *mol, const topo_mol_ident_t *target);


	/**
	 *	delete one or more bonds between atoms based on a pair of identifier dictionaries.
	 *	if the identifier only contains a segid, all atoms in the segment are selected,
	 *	if the identifier contains a segid and resid, all atoms in the resid are selected, and
	 *	if the identifier contans a segid, resid and atomid, a specific atom is selected.
	 *
	 *	@param a: a dictionary identifier for the first group of atoms
	 *	@param b: a dictionary identifier for the second group of atoms.
	 */
	void delete_bonds(const topo_mol_ident_t *a, const topo_mol_ident_t *b);

	/**
	 * Build a segment.
	 * A new empty segment with the given segid is created, and buildseg now points
	 * to this new segment. If a new segment fails to be created for any reason,
	 * buildseg will be Null.
	 * @param segid: the segment id of the new segment
	 * @returns 0 on success, -2 if the segid name is too long, -3 if the segid
	 * already exists, -4 or -5 if the segid fails to be created.
	 */
	int build_segment(const char *segid);

	/**
	 * Build a new residue.
	 * A new empty residue is added to the segment currently being built.
	 * The residue type must be know by the defintions, i.e. it must exist in
	 * the definitions which may be loaded from a topology file.
	 * @param resid: the residue id, must be unique in a segment.
	 * @param rname: the residue name, must exist in the definitions.
	 * @param chain: ?
	 * @returns: 0 on success, negative otherwise.
	 */
	int build_residue(const char *resid, const char *rname, const char *chain);

	/**
	 * Finish off the molecule.
	 * Determines all bonds from the definitions, adds them to the molecule. Fills
	 * in missing atoms.
	 * Applies all queued patches.
	 * @returns 0 on success, negative otherwise.
	 */
	int end();

    /**
     * Change the type of a single residue in the current segment. 
     * @param resid: Unique name for residue, 1-5 characters, usually numeric.
     * @param resname: New residue type name from topology file.
     */
	int mutate(const char *resid, const char *resname);

    
    /**
     * Remove all angles and/or dihedrals and completely regenerate them using the segment 
     * automatic generation algorithms. This is only needed if patches were applied that 
     * do not correct angles and bonds. Segment and file defaults are ignored, and 
     * angles/dihedrals for the entire molecule are regenerated from scratch.
     * @param what: what to regenerate, currently either ANGLES or DIHEDRALS.
     */
	int regenerate(const char* what);

    /**
     * Create multiple images of a set of atoms for use in locally enhanced sampling. 
     * The beta column of the output pdb file is set to 1...ncopies for each image. 
     * Multiple copies of bonds, angles, etc. are created. Atom, residue or segment 
     * names are not altered; images are distinguished only by beta value. This is not a 
     * normal molecular structure and may confuse other tools.
     *
     * Should be called after one or more segments have been built, all patches applied, 
     * and coordinates guessed. The effects of this command may confuse other functions.
     * May be repeated as many times as necessary to include all atoms.
     *
     * @param target: a dictionary with keys {segid:,resid:,aname} identifying the  segment, 
     * residue, or atom to be multiplied. If :resid is omitted the entire
     * segment is multiplied; if :aname is omitted the entire residue is multiplied. 
     * @param ncopies: how many copies to make.
     */
    int multiply_atoms(const topo_mol_ident_t *target, int ncopies);

	int set_element(const topo_mol_ident_t *target,
			const char *element, int replace);

	int set_chain(const topo_mol_ident_t *target,
			const char *chain, int replace);

	int set_coordinates(const topo_mol_ident_t *target,
			double x, double y, double z);

    /**
     * Guesses coordinates of atoms for which they were not explicitly set. Calculation is
     * based on internal coordinate hints contained in toplogy definition files. When these
     * are insufficient, wild guesses are attempted based on bond lengths of 0.1 nm
     * and angles of 109 degrees.
     *
     * @returns 0 on success, negative on failure.
     */
	int guess_coordinates();


	/**
	 * DO NOT USE YET.
	 */
	int add_patch(const char *pname, int deflt);

	/**
	 * DO NOT USE YET.
	 */
	int add_patchres(const topo_mol_ident_t *targets);

	int validate_patchres(const char *pname, const char *segid, const char *resid);

    /**
     * Provide translations from residues found in PDB files to proper residue names read in from
     * topology definition files. Proper names from topology files will be used in generated PSF and PDB
     * files. 
     *
     * Should be called before reading sequence with read_residues. May be called multiple times.
     * @param altres: Residue name found in PDB file.
     * @param realres: Residue name found in topology file.
     */
    int alias_residue(const char* altres, const char* realres);


    /**
     * Provide translations from atom names found in PDB files to proper atom names read in
     * from topology definition files. Proper names from topology files will be used in generated PSF and
     * PDB files. 
     *
     * Should be called before reading coordinates with read_coordinates. May call multiple times.
     * @param resname: Proper or aliased residue name.
     * @param altname: Atom name found in PDB file.
     * @param realname: Atom name found in topology file.
     */
    int alias_atom(const char* resname, const char* altname, const char* realname);


private:
	/**
	 * private full ctor
	 */
	void _topo_mol(topo_defs *defs);

    struct stringhash *aliases;
};

/**
 * Set the psfgen logging level. This will eventually be replaced with the set_logger function.
 * @param lvl: A string which may be of ["DEBUG", "INFO", "WARN", "ERROR", "CRITICAL"].
 * The case is irrelevant.
 * @returns 0 on success, -1 if the string is not one of the above values.
 */
int logging_set_level(const std::string& lvl);

/**
 * Get the logging level.
 * @returns: A string of ["DEBUG", "INFO", "WARN", "ERROR", "CRITICAL"].
 */
std::string logging_get_level();

/**
 * Set the logger to a Python logging.Logger object.
 * NOT IMPLEMENTED YET!
 */
int set_logger(PyObject *logger);

topo_mol_bond * topo_mol_bond_next(
    topo_mol_bond *tuple, topo_mol_atom *atom);

topo_mol_angle * topo_mol_angle_next(
    topo_mol_angle *tuple, topo_mol_atom *atom);

topo_mol_dihedral * topo_mol_dihedral_next(
    topo_mol_dihedral *tuple, topo_mol_atom *atom);

topo_mol_improper * topo_mol_improper_next(
    topo_mol_improper *tuple, topo_mol_atom *atom);

topo_mol_cmap * topo_mol_cmap_next(
    topo_mol_cmap *tuple, topo_mol_atom *atom);

topo_mol_chargegroup * topo_mol_chargegroup_next(
    topo_mol_chargegroup *tuple, topo_mol_atom *atom);


#endif

