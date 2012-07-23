
// File: index.xml

// File: structtopo__defs.xml
%feature("docstring") topo_defs "";

%feature("docstring")  topo_defs::topo_defs "topo_defs::topo_defs()

creates a topo_defs object ";

%feature("docstring")  topo_defs::~topo_defs "topo_defs::~topo_defs()
";

%feature("docstring")  topo_defs::build_residue "int
topo_defs::build_residue(const char *rname, int patch)

checks if this defs has a residue with the given residue name. if one
exists, returns 0. If one does not exist, a new one is created and
inserted into the residue array and the buildres member points to this
new residue.

The new residue is intialized to all zeros except patch is set to the
given patch value.

buildres is set to 0 if either a residue already exists or fails to be
created.

0 on success, -2 if name too long, -4 if the name fails to insert into
hash array ";

%feature("docstring")  topo_defs::read_topology "int
topo_defs::read_topology(const char *fname, bool all_caps=false) ";

%feature("docstring")  topo_defs::get_residues "std::vector<topo_defs_residue*> topo_defs::get_residues() ";

%feature("docstring")  topo_defs::get_patches "std::vector<topo_defs_residue*> topo_defs::get_patches() ";

%feature("docstring")  topo_defs::get_topology_filenames "std::vector<std::string> topo_defs::get_topology_filenames() ";

%feature("docstring")  topo_defs::__str__ "std::string
topo_defs::__str__() ";


// File: structtopo__defs__angle__t.xml
%feature("docstring") topo_defs_angle_t "";


// File: structtopo__defs__atom__t.xml
%feature("docstring") topo_defs_atom_t "";


// File: structtopo__defs__bond__t.xml
%feature("docstring") topo_defs_bond_t "";


// File: structtopo__defs__cmap__t.xml
%feature("docstring") topo_defs_cmap_t "";


// File: structtopo__defs__conformation__t.xml
%feature("docstring") topo_defs_conformation_t "";


// File: structtopo__defs__dihedral__t.xml
%feature("docstring") topo_defs_dihedral_t "";


// File: structtopo__defs__improper__t.xml
%feature("docstring") topo_defs_improper_t "";


// File: structtopo__defs__residue.xml
%feature("docstring") topo_defs_residue "";

%feature("docstring")  topo_defs_residue::__str__ "std::string
topo_defs_residue::__str__() ";


// File: structtopo__defs__topofile__t.xml
%feature("docstring") topo_defs_topofile_t "";


// File: structtopo__defs__type__t.xml
%feature("docstring") topo_defs_type_t "";


// File: classtopo__mol.xml
%feature("docstring") topo_mol "";

%feature("docstring")  topo_mol::topo_mol "topo_mol::topo_mol()

creates a new topo_mol object. does stuff. ";

%feature("docstring")  topo_mol::topo_mol "topo_mol::topo_mol(topo_defs *defs) ";

%feature("docstring")  topo_mol::~topo_mol "topo_mol::~topo_mol() ";

%feature("docstring")  topo_mol::get_segments "std::vector<topo_mol_segment_t*> topo_mol::get_segments()

Get all the segments that belong to this molcule. Do not attempt to
append to this list, to create a new segment, use the build_segment
method. : a list of segments. ";

%feature("docstring")  topo_mol::read_topology "int
topo_mol::read_topology(const char *filename)

Read a topology file. the patches of the topology file are added to
the definitions and made available for use. This may be called
multiple times, as the contents are simply added.

Parameters:
-----------

filename:  ::  the name of the topo file

0 on success, something else otherwise. ";

%feature("docstring")  topo_mol::read_psf "int
topo_mol::read_psf(const char *filename)

read a psf file.

Parameters:
-----------

filename:  ::  the name of a file

: 0 on success, something else on failure. ";

%feature("docstring")  topo_mol::read_residues "int
topo_mol::read_residues(const char *filename)

Extract sequence information (residues) from PDB file when building
segment. Residue IDs will be preserved, residue names must match
entries in the topology file or should be aliased before pdb is
called.

Parameters:
-----------

filename:  ::  PDB file containing known or aliased residues. ";

%feature("docstring")  topo_mol::read_coordinates "int
topo_mol::read_coordinates(const char *filename, const char *segid=0)

Read coordinates from file (currently PDB), matching segment, residue
and atom names.

Parameters:
-----------

filename:  ::  PDB file containing known or aliased residues and
atoms.

segid:  ::  If specified override segment IDs in PDB file. ";

%feature("docstring")  topo_mol::num_patch_targets "int
topo_mol::num_patch_targets(const char *presname)

Currently not exactly sure what this does, something like counting the
number of times a patch can be applied???

Parameters:
-----------

presname:  ::  the name of a patch residue. ";

%feature("docstring")  topo_mol::patch "int topo_mol::patch(const
char *patch_res_name, PyObject *target, int prepend=0, int
warn_angles=0, int warn_dihedrals=0, int deflt=0)

Apply a patch to the molecule.

Parameters:
-----------

patch_res_name:  ::  the residiue name of the patch

target:  ::  the segid/resid of the residue to apply the path to.

prepend:  ::  optional, defaults to 0

warn_angles:  ::  optional, defaults to 0

warn_dihedrals:  ::  optional, defaults to 0

deflt:  ::  optional, defaults to 0 ";

%feature("docstring")  topo_mol::select_atom "std::vector<topo_mol_atom*> topo_mol::select_atom(const
topo_mol_ident_t *target)

Select a set of atoms. This performs the same selections as
delete_atom, except it returns the atoms instead of deletes them.

Parameters:
-----------

target:  ::  a Dictionary with the keys {'segid', 'resid', 'aname'}
where these corespond to the segment id, residue id, and atom name.
resid and aname may be null.

a list of topo_mol_atom objects. ";

%feature("docstring")  topo_mol::write "int topo_mol::write(const
std::string &file_name, const std::string &file_type=\"\")

write an output file. The type of file is automatically determinied by
the file name, or can be specified in the optional file_type
parameter.

Parameters:
-----------

file_name:  ::  the file name, the type is guessed from the file
extension, currently supported files are .pdb and .psf

file_type:  ::  the file type, currently, either \"pdb\", \"psf\", or
leave blank to have the file type automatically determined from the
file name. ";

%feature("docstring")  topo_mol::delete_atom "void
topo_mol::delete_atom(topo_mol *mol, const topo_mol_ident_t *target)

delete one or more atoms based on the identifier dictionary. if the
identifier only contains a segid, all atoms in the segment are
deleted, if the identifier contains a segid and resid, all atoms in
the resid are deleted, and if the identifier contans a segid, resid
and atomid, a specific atom is deleted.

See:  topo_mol::select_atom

Parameters:
-----------

mol:  ::  a topo_mol object

target:  ::  a dictionary with the keys 'segid', 'resid', 'aname'\";
";

%feature("docstring")  topo_mol::delete_bonds "void
topo_mol::delete_bonds(const topo_mol_ident_t *a, const
topo_mol_ident_t *b)

delete one or more bonds between atoms based on a pair of identifier
dictionaries. if the identifier only contains a segid, all atoms in
the segment are selected, if the identifier contains a segid and
resid, all atoms in the resid are selected, and if the identifier
contans a segid, resid and atomid, a specific atom is selected.

Parameters:
-----------

a:  ::  a dictionary identifier for the first group of atoms

b:  ::  a dictionary identifier for the second group of atoms. ";

%feature("docstring")  topo_mol::build_segment "int
topo_mol::build_segment(const char *segid)

Build a segment. A new empty segment with the given segid is created,
and buildseg now points to this new segment. If a new segment fails to
be created for any reason, buildseg will be Null.

Parameters:
-----------

segid:  ::  the segment id of the new segment

0 on success, -2 if the segid name is too long, -3 if the segid
already exists, -4 or -5 if the segid fails to be created. ";

%feature("docstring")  topo_mol::build_residue "int
topo_mol::build_residue(const char *resid, const char *rname, const
char *chain)

Build a new residue. A new empty residue is added to the segment
currently being built. The residue type must be know by the
defintions, i.e. it must exist in the definitions which may be loaded
from a topology file.

Parameters:
-----------

resid:  ::  the residue id, must be unique in a segment.

rname:  ::  the residue name, must exist in the definitions.

chain:  ::  ?

: 0 on success, negative otherwise. ";

%feature("docstring")  topo_mol::end "int topo_mol::end()

Finish off the molecule. Determines all bonds from the definitions,
adds them to the molecule. Fills in missing atoms. Applies all queued
patches. 0 on success, negative otherwise. ";

%feature("docstring")  topo_mol::mutate "int topo_mol::mutate(const
char *resid, const char *resname)

Change the type of a single residue in the current segment.

Parameters:
-----------

resid:  ::  Unique name for residue, 1-5 characters, usually numeric.

resname:  ::  New residue type name from topology file. ";

%feature("docstring")  topo_mol::regenerate "int
topo_mol::regenerate(RegenerateType what)

Remove all angles and/or dihedrals and completely regenerate them
using the segment automatic generation algorithms. This is only needed
if patches were applied that do not correct angles and bonds. Segment
and file defaults are ignored, and angles/dihedrals for the entire
molecule are regenerated from scratch.

Parameters:
-----------

what:  ::  what to regenerate, currently either ANGLES or DIHEDRALS.
";

%feature("docstring")  topo_mol::multiply_atoms "int
topo_mol::multiply_atoms(const topo_mol_ident_t *target, int ncopies)

Create multiple images of a set of atoms for use in locally enhanced
sampling. The beta column of the output pdb file is set to 1...ncopies
for each image. Multiple copies of bonds, angles, etc. are created.
Atom, residue or segment names are not altered; images are
distinguished only by beta value. This is not a normal molecular
structure and may confuse other tools.

Should be called after one or more segments have been built, all
patches applied, and coordinates guessed. The effects of this command
may confuse other functions. May be repeated as many times as
necessary to include all atoms.

Parameters:
-----------

target:  ::  a dictionary with keys {segid:,resid:,aname} identifying
the segment, residue, or atom to be multiplied. If :resid is omitted
the entire segment is multiplied; if :aname is omitted the entire
residue is multiplied.

ncopies:  ::  how many copies to make. ";

%feature("docstring")  topo_mol::set_element "int
topo_mol::set_element(const topo_mol_ident_t *target, const char
*element, int replace) ";

%feature("docstring")  topo_mol::set_chain "int
topo_mol::set_chain(const topo_mol_ident_t *target, const char *chain,
int replace) ";

%feature("docstring")  topo_mol::set_coordinates "int
topo_mol::set_coordinates(const topo_mol_ident_t *target, double x,
double y, double z) ";

%feature("docstring")  topo_mol::guess_coordinates "int
topo_mol::guess_coordinates()

Guesses coordinates of atoms for which they were not explicitly set.
Calculation is based on internal coordinate hints contained in toplogy
definition files. When these are insufficient, wild guesses are
attempted based on bond lengths of 0.1 nm and angles of 109 degrees.

0 on success, negative on failure. ";

%feature("docstring")  topo_mol::add_patch "int
topo_mol::add_patch(const char *pname, int deflt)

DO NOT USE YET. ";

%feature("docstring")  topo_mol::add_patchres "int
topo_mol::add_patchres(const topo_mol_ident_t *targets)

DO NOT USE YET. ";

%feature("docstring")  topo_mol::validate_patchres "int
topo_mol::validate_patchres(const char *pname, const char *segid,
const char *resid) ";

%feature("docstring")  topo_mol::alias_residue "int
topo_mol::alias_residue(const char *altres, const char *realres)

Provide translations from residues found in PDB files to proper
residue names read in from topology definition files. Proper names
from topology files will be used in generated PSF and PDB files.

Should be called before reading sequence with read_residues. May be
called multiple times.

Parameters:
-----------

altres:  ::  Residue name found in PDB file.

realres:  ::  Residue name found in topology file. ";

%feature("docstring")  topo_mol::alias_atom "int
topo_mol::alias_atom(const char *resname, const char *altname, const
char *realname)

Provide translations from atom names found in PDB files to proper atom
names read in from topology definition files. Proper names from
topology files will be used in generated PSF and PDB files.

Should be called before reading coordinates with read_coordinates. May
call multiple times.

Parameters:
-----------

resname:  ::  Proper or aliased residue name.

altname:  ::  Atom name found in PDB file.

realname:  ::  Atom name found in topology file. ";


// File: structtopo__mol__angle.xml
%feature("docstring") topo_mol_angle "";


// File: structtopo__mol__atom.xml
%feature("docstring") topo_mol_atom "";

%feature("docstring")  topo_mol_atom::__str__ "std::string
topo_mol_atom::__str__() ";


// File: structtopo__mol__bond.xml
%feature("docstring") topo_mol_bond "";

%feature("docstring")  topo_mol_bond::__str__ "const char*
topo_mol_bond::__str__() ";

%feature("docstring")  topo_mol_bond::__repr__ "const char*
topo_mol_bond::__repr__() ";


// File: structtopo__mol__chargegroup.xml
%feature("docstring") topo_mol_chargegroup "";


// File: structtopo__mol__cmap.xml
%feature("docstring") topo_mol_cmap "";


// File: structtopo__mol__conformation.xml
%feature("docstring") topo_mol_conformation "";


// File: structtopo__mol__dihedral.xml
%feature("docstring") topo_mol_dihedral "";


// File: structtopo__mol__ident__t.xml
%feature("docstring") topo_mol_ident_t "";


// File: structtopo__mol__improper.xml
%feature("docstring") topo_mol_improper "";


// File: structtopo__mol__patch__t.xml
%feature("docstring") topo_mol_patch_t "";


// File: structtopo__mol__patchres__t.xml
%feature("docstring") topo_mol_patchres_t "";


// File: structtopo__mol__residue__t.xml
%feature("docstring") topo_mol_residue_t "";

%feature("docstring")  topo_mol_residue_t::get_atoms "std::vector<topo_mol_atom*> topo_mol_residue_t::get_atoms()

get the list of atoms that belong to this residue. See:  topo_mol_atom

: a list of atoms. ";

%feature("docstring")  topo_mol_residue_t::__str__ "std::string
topo_mol_residue_t::__str__() ";


// File: structtopo__mol__segment__t.xml
%feature("docstring") topo_mol_segment_t "";

%feature("docstring")  topo_mol_segment_t::get_residues "std::vector<topo_mol_residue_t*> topo_mol_segment_t::get_residues()

get a list of residues that belong to this segment. ";

%feature("docstring")  topo_mol_segment_t::__str__ "std::string
topo_mol_segment_t::__str__() ";


// File: topo__defs__struct_8h.xml


// File: topo__mol__struct_8h.xml
%feature("docstring")  logging_set_level "int logging_set_level(const
std::string &lvl)

Set the psfgen logging level. This will eventually be replaced with
the set_logger function.

Parameters:
-----------

lvl:  ::  A string which may be of [\"DEBUG\", \"INFO\", \"WARN\",
\"ERROR\", \"CRITICAL\"]. The case is irrelevant.

0 on success, -1 if the string is not one of the above values. ";

%feature("docstring")  logging_get_level "std::string
logging_get_level()

Get the logging level. : A string of [\"DEBUG\", \"INFO\", \"WARN\",
\"ERROR\", \"CRITICAL\"]. ";

%feature("docstring")  set_logger "int set_logger(PyObject *logger)

Set the logger to a Python logging.Logger object. NOT IMPLEMENTED YET!
";

%feature("docstring")  topo_mol_bond_next "topo_mol_bond*
topo_mol_bond_next(topo_mol_bond *tuple, topo_mol_atom *atom) ";

%feature("docstring")  topo_mol_angle_next "topo_mol_angle*
topo_mol_angle_next(topo_mol_angle *tuple, topo_mol_atom *atom) ";

%feature("docstring")  topo_mol_dihedral_next "topo_mol_dihedral*
topo_mol_dihedral_next(topo_mol_dihedral *tuple, topo_mol_atom *atom)
";

%feature("docstring")  topo_mol_improper_next "topo_mol_improper*
topo_mol_improper_next(topo_mol_improper *tuple, topo_mol_atom *atom)
";

%feature("docstring")  topo_mol_cmap_next "topo_mol_cmap*
topo_mol_cmap_next(topo_mol_cmap *tuple, topo_mol_atom *atom) ";

%feature("docstring")  topo_mol_chargegroup_next "topo_mol_chargegroup* topo_mol_chargegroup_next(topo_mol_chargegroup
*tuple, topo_mol_atom *atom) ";

