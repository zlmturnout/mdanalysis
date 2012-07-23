%include std_string.i
%include doxydocs.i


%module (docstring="This is the example module's docstring") psfgen


/*
 number of segments
 nseg = hasharray_count(mol->segment_hash);

 seg = mol->segment_array[iseg];

 nres = hasharray_count(seg->residue_hash);
 res = &(seg->residue_array[ires]);
 */

// %ignore topo_mol_select_atom(topo_mol *, const topo_mol_ident_t *, int *, topo_mol_atom_t **);

//%feature("novaluewrapper") std::pair<int, int>;

%ignore topo_mol_ident_t;
%ignore hasharray;

%{
#include "pdb_file_extract.h"
#include "psf_file_extract.h"
#include "topo_defs.h"
#include "topo_mol.h"
#include "topo_mol_struct.h"
#include "charmm_parse_topo_defs.h"
#include "stringhash.h"
#include "extract_alias.h"
#include "psfgen_data.h"
#include "topo_defs_struct.h"
#include "string.h"
%}


//%include "std_vector.i"
//namespace std {
//   %template(topo_mol_atom_vector) vector<topo_mol_atom*>;
//}
//%feature("novaluewrapper") topo_mol_atom_vector;

%typemap(in) (const topo_mol_ident_t *targets, int ntargets) {
    int i;
    if (!PyList_Check($input)) {
        PyErr_SetString(PyExc_ValueError, "Expected a list of ident_t.");
        return NULL;
    }
    $2 = PyList_Size($input);
    $1 = (topo_mol_ident_t *)malloc($2*sizeof(topo_mol_ident_t));
    for (i=0; i<$2; i++) {
        PyObject *anameobj;
        PyObject *identobj = PyList_GET_ITEM($input, i);
        if (!PyDict_Check(identobj)) {
            free($1);
            PyErr_SetString(PyExc_ValueError, "Expected a list of dictionaries.");
            return NULL;
        }
        $1[i].segid = PyString_AsString(PyDict_GetItemString(identobj, "segid"));
        $1[i].resid = PyString_AsString(PyDict_GetItemString(identobj, "resid"));
        if ((anameobj = PyDict_GetItemString(identobj, "aname"))) {
            $1[i].aname = PyString_AsString(anameobj);
        } else {
			$1[i].aname = 0;
		}
    }
 }

%typemap(freearg) (const topo_mol_ident_t *targets, int ntargets) {
    free($1);
}

%typemap(in) const topo_mol_ident_t * {
    PyObject *tmpobj;
    if (!PyDict_Check($input)) {
        PyErr_SetString(PyExc_ValueError, "Expected a dictionary.");
        return NULL;
    }
    $1 = (topo_mol_ident_t *)malloc(sizeof(topo_mol_ident_t));
	tmpobj = PyDict_GetItemString($input, "segid");
	$1->segid = tmpobj ? PyString_AsString(tmpobj) : 0;
	tmpobj = PyDict_GetItemString($input, "resid");
	$1->resid = tmpobj ? PyString_AsString(tmpobj) : 0;
	tmpobj = PyDict_GetItemString($input, "aname");
	$1->aname = tmpobj ? PyString_AsString(tmpobj) : 0;
 }
%typemap(freearg) const topo_mol_ident_t * {
    free($1);
 }

%typemap(in) void (*print_msg)(void *, const char *) {
    $1=psfgen_msg;
 }


%typemap(in) struct topo_mol_atom *[ANY]() {
    
    topo_mol_atom *atmp[$1_dim0];
    
    if (PyTuple_Check($input)) {
        if(PyTuple_Size($input) == $1_dim0) {
            // get atom types from tuple
            int i;
            for(i = 0; i < $1_dim0; i++) {
                PyObject *obj = PyTuple_GetItem($input, i);
                if (SWIG_ConvertPtr(obj,(void **)&atmp[i], SWIGTYPE_p_topo_mol_atom,0) < 0 ) {
                    // bad pointer
                    char format[] = "could not convert tuple element %i to a topo_mol_atom *";
                    char *buffer = (char*)alloca(sizeof(format) + 10);
                    sprintf(buffer, format, i);
                    SWIG_exception_fail(SWIG_ValueError, buffer);
                }
            }
        } else {
            // bad size
            SWIG_exception_fail(SWIG_ValueError, "expected a tuple of size $1_dim0");
        }
    } else {
        // not a tuple
        SWIG_exception_fail(SWIG_ValueError, "argument is not a tuple");
    }
    
    // all good
    $1 = atmp;
 }

%typemap(out) struct topo_mol_atom *[ANY]() {
    
    // TODO add error handling
    PyObject *presult = PyTuple_New($1_dim0);
    
    // pack the tuple
    int i;
    for(i = 0; i < $1_dim0; i++) {
        PyObject *p = SWIG_NewPointerObj(SWIG_as_voidptr($1[i]), SWIGTYPE_p_topo_mol_atom, 0 |  0 );
        PyTuple_SetItem(presult, i, p);
    }
    
    // all good
    $result = presult;
 }


%typemap(in) struct topo_mol_bond *[ANY]() {
    
    topo_mol_bond *atmp[$1_dim0];
    
    if (PyTuple_Check($input)) {
        if(PyTuple_Size($input) == $1_dim0) {
            // get bond types from tuple
            int i;
            for(i = 0; i < $1_dim0; i++) {
                PyObject *obj = PyTuple_GetItem($input, i);
                if (SWIG_ConvertPtr(obj,(void **)&atmp[i], SWIGTYPE_p_topo_mol_bond,0) < 0 ) {
                    // bad pointer
                    char format[] = "could not convert tuple element %i to a topo_mol_bond *";
                    char *buffer = (char*)alloca(sizeof(format) + 10);
                    sprintf(buffer, format, i);
                    SWIG_exception_fail(SWIG_ValueError, buffer);
                }
            }
        } else {
            // bad size
            SWIG_exception_fail(SWIG_ValueError, "expected a tuple of size $1_dim0");
        }
    } else {
        // not a tuple
        SWIG_exception_fail(SWIG_ValueError, "argument is not a tuple");
    }
    
    // all good
    $1 = atmp;
}

%typemap(out) struct topo_mol_bond *[ANY]() {
    
    // TODO add error handling
    PyObject *presult = PyTuple_New($1_dim0);
    
    // pack the tuple
    int i;
    for(i = 0; i < $1_dim0; i++) {
        PyObject *p = SWIG_NewPointerObj(SWIG_as_voidptr($1[i]), SWIGTYPE_p_topo_mol_bond, 0 |  0 );
        PyTuple_SetItem(presult, i, p);
    }

    // all good
    $result = presult;
}

/**
 * make a typemap to only copy the chars up to null from a fixed array, 
 * othwise bunch of junk is returned.
 */
%typemap(out) char [ANY]() {
    
    // TODO add better error handling
    int len = ersatz_strnlen($1, $1_dim0);

    // all good
    $result = PyString_FromStringAndSize($1, len);
}

%typemap(out) std::vector<topo_mol_atom*> {
    // typemap std::vector<topo_mol_atom*> 

    const std::vector<topo_mol_atom*>& atoms = $1; 

    PyObject *list = PyList_New(atoms.size());
    
    for(int i = 0; i < atoms.size(); i++) {
        PyObject *p = SWIG_NewPointerObj(SWIG_as_voidptr(atoms[i]), SWIGTYPE_p_topo_mol_atom, 0 |  0 ); 
        PyList_SET_ITEM(list, i, p);
    }

    $result = list;
}

%typemap(out) std::vector<topo_mol_segment_t*> {
    // typemap std::vector<topo_mol_segment_t*> 

    const std::vector<topo_mol_segment_t*>& vec = $1; 

    PyObject *list = PyList_New(vec.size());
    
    for(int i = 0; i < vec.size(); i++) {
        PyObject *p = SWIG_NewPointerObj(SWIG_as_voidptr(vec[i]), SWIGTYPE_p_topo_mol_segment_t, 0 |  0 ); 
        PyList_SET_ITEM(list, i, p);
    }

    $result = list;
}


%typemap(out) std::vector<topo_defs_residue*> {
    // typemap std::vector<topo_defs_residue*> 

    const std::vector<topo_defs_residue*>& residues = $1; 

    PyObject *list = PyList_New(residues.size());
    
    for(int i = 0; i < residues.size(); i++) {
        PyObject *p = SWIG_NewPointerObj(SWIG_as_voidptr(residues[i]), SWIGTYPE_p_topo_defs_residue, 0 |  0 ); 
        PyList_SET_ITEM(list, i, p);
    }

    $result = list;
}


%typemap(out) std::vector<topo_mol_residue_t*> {
    // typemap std::vector<topo_mol_residue_t*> 

    const std::vector<topo_mol_residue_t*>& vec = $1; 

    PyObject *list = PyList_New(vec.size());
    
    for(int i = 0; i < vec.size(); i++) {
        PyObject *p = SWIG_NewPointerObj(SWIG_as_voidptr(vec[i]), SWIGTYPE_p_topo_mol_residue_t, 0 |  0 ); 
        PyList_SET_ITEM(list, i, p);
    }
    $result = list;
}


%typemap(out) std::vector<std::string> {
    // typemap std::vector<std::string>

    const std::vector<std::string>& vec = $1;

    PyObject *list = PyList_New(vec.size());
    
    for(int i = 0; i < vec.size(); i++) {
        PyObject *s = PyString_FromString(vec[i].c_str());
        PyList_SET_ITEM(list, i, s);
    }

    $result = list;
}



/**
 * input typemap for the natoms, atoms vars. This is needed
 * because the pointers that are passed to the function need to 
 * point to something, so they point to the tmp_* vars below.
 */
/* %typemap(in,numimputs=0) (int *natoms, topo_mol_atom_t ***atoms) { */
/* 	int tmp_natoms; */
/* 	topo_mol_atom_t **tmp_atoms; */
	
/* 	$1=&tmp_natoms; */
/* 	$2=&tmp_atoms; */
/* } */


/* %typemap(argout) (int *natoms, topo_mol_atom_t ***atoms) {   */
/*     printf("in the typemap, typemap(argout) (int *natoms, topo_mol_atom_t **atoms) {\n"); */

/*     // the func was called, output args now in $1 and $2. */
/*     int natoms = *$1; */
/*     topo_mol_atom_t **atoms = *$2; */
/*     PyObject *list = PyList_New(natoms); */
/*     int i; */
/*     for(i = 0; i < natoms; i++) { */
/*         PyObject *p = SWIG_NewPointerObj(SWIG_as_voidptr(atoms[i]), SWIGTYPE_p_topo_mol_atom_t, 0 |  0 ); */
/*         PyList_SET_ITEM(list, i, p); */
/*     } */
/*     free(atoms); */
/* 	$result = list; */
/* } */



%include pdb_file_extract.h
%include psf_file_extract.h
%include topo_defs.h
%include topo_mol.h
%include topo_mol_struct.h
%include charmm_parse_topo_defs.h
%include stringhash.h
%include extract_alias.h
%include psfgen_data.h
%include topo_defs_struct.h









