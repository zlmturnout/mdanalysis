#ifndef TOPO_MOL_H
#define TOPO_MOL_H
#ifndef SWIG

#pragma push_macro("HAVE_DIRENT_H")
#undef HAVE_DIRENT_H
#pragma push_macro( "HAVE_FSEEKO" )
#undef HAVE_FSEEKO
#pragma push_macro( "HAVE_FSYNC" )
#undef HAVE_FSYNC
#pragma push_macro( "HAVE_GETTIMEOFDAY" )
#undef HAVE_GETTIMEOFDAY
#pragma push_macro( "HAVE_LSTAT" )
#undef HAVE_LSTAT
#pragma push_macro( "HAVE_PTHREAD_H" )
#undef HAVE_PTHREAD_H
#pragma push_macro( "HAVE_SIGACTION" )
#undef HAVE_SIGACTION
#pragma push_macro( "HAVE_STRDUP" )
#undef HAVE_STRDUP
#pragma push_macro( "HAVE_SYSCONF" )
#undef HAVE_SYSCONF
#pragma push_macro( "HAVE_SYS_TIME_H" )
#undef HAVE_SYS_TIME_H
#pragma push_macro( "HAVE_SYS_TYPES_H" )
#undef HAVE_SYS_TYPES_H
#pragma push_macro( "HAVE_UNISTD_H" )
#undef HAVE_UNISTD_H
#pragma push_macro( "_POSIX_C_SOURCE" )
#undef _POSIX_C_SOURCE
#pragma push_macro( "_XOPEN_SOURCE" )
#undef _XOPEN_SOURCE

#include <Python.h>

#pragma pop_macro("HAVE_DIRENT_H")
#pragma pop_macro( "HAVE_FSEEKO" )
#pragma pop_macro( "HAVE_FSYNC" )
#pragma pop_macro( "HAVE_GETTIMEOFDAY" )
#pragma pop_macro( "HAVE_LSTAT" )
#pragma pop_macro( "HAVE_PTHREAD_H" )
#pragma pop_macro( "HAVE_SIGACTION" )
#pragma pop_macro( "HAVE_STRDUP" )
#pragma pop_macro( "HAVE_SYSCONF" )
#pragma pop_macro( "HAVE_SYS_TIME_H" )
#pragma pop_macro( "HAVE_SYS_TYPES_H" )
#pragma pop_macro( "HAVE_UNISTD_H" )
#pragma pop_macro( "_POSIX_C_SOURCE" )
#pragma pop_macro( "_XOPEN_SOURCE" )


#include "topo_defs.h"
#include "topo_mol_struct.h"
#include <string>

struct topo_mol;
typedef struct topo_mol topo_mol;

topo_mol * topo_mol_create(topo_defs *defs);
void topo_mol_destroy(topo_mol *mol);


int topo_mol_segment(topo_mol *mol, const char *segid);


int topo_mol_segment_auto_angles(topo_mol *mol, int autogen);
int topo_mol_segment_auto_dihedrals(topo_mol *mol, int autogen);

int topo_mol_residue(topo_mol *mol, const char *resid, const char *rname,
		const char *chain);
int topo_mol_mutate(topo_mol *mol, const char *resid, const char *rname);

int topo_mol_end(topo_mol *mol);


int topo_mol_regenerate_angles(topo_mol *mol);
int topo_mol_regenerate_dihedrals(topo_mol *mol);

void topo_mol_delete_atom(topo_mol *mol, const topo_mol_ident_t *target);


int topo_mol_multiply_atoms(topo_mol *mol, const topo_mol_ident_t *targets,
		int ntargets, int ncopies);

int topo_mol_set_element(topo_mol *mol, const topo_mol_ident_t *target,
		const char *element, int replace);

int topo_mol_set_chain(topo_mol *mol, const topo_mol_ident_t *target,
		const char *chain, int replace);

int topo_mol_set_xyz(topo_mol *mol, const topo_mol_ident_t *target,
		double x, double y, double z);

int topo_mol_guess_xyz(topo_mol *mol);

int topo_mol_add_patch(topo_mol *mol, const char *pname, int deflt);

int topo_mol_add_patchres(topo_mol *mol, const topo_mol_ident_t *targets);

int topo_mol_validate_patchres(topo_mol *mol, const char *pname, const char *segid, const char *resid);



enum psfgen_logging_level {LOG_DEBUG = 0, LOG_INFO = 1, LOG_WARN = 2, LOG_ERROR = 3, LOG_CRITICAL = 4};
/**
 * psfgen message handler
 * write stuff to a location specified by ...

 */
extern "C" void psfgen_msg(void* v, const char* s);

extern "C" void psfgen_log(psfgen_logging_level,const char* msg,const char* src=0);

#define PSFGEN_LOG(lvl,msg) psfgen_log(lvl,msg,__PRETTY_FUNCTION__)
#define PSFGEN_DEBUG(msg) psfgen_log(LOG_DEBUG,msg,__PRETTY_FUNCTION__)
#define PSFGEN_INFO(msg) psfgen_log(LOG_INFO,msg,__PRETTY_FUNCTION__)
#define PSFGEN_WARN(msg) psfgen_log(LOG_WARN,msg,__PRETTY_FUNCTION__)
#define PSFGEN_ERROR(msg) psfgen_log(LOG_ERROR,msg,__PRETTY_FUNCTION__)
#define PSFGEN_CRITICAL(msg) psfgen_log(LOG_CRITICAL,msg,__PRETTY_FUNCTION__)



/**
 * BSDs do not have strnlen, but Linux does.
 * need a macro around this to use the built in srtnlen on systems
 * that have it. 
 *
 *
 * Find the length of S, but scan at most MAXLEN characters.  If no '\0'
 * terminator is found within the first MAXLEN characters, return MAXLEN. 
 */
inline size_t ersatz_strnlen(const char* s, size_t maxlen)
{
    register const char *e;
    size_t n;

    for (e = s, n = 0; *e && n < maxlen; e++, n++)
        ;
    return n;
}


/**
 * create a string from a name beign safe as the name array *might* not have
 * a terminal null.
 */
inline std::string namestr(const char* name) {
    size_t len = ersatz_strnlen(name, NAMEMAXLEN);
    return std::string(name, len);
}


/**
 * populates an ident from a python dictionary.
 * returns 0 if OK, -1 on failure.
 */
int ident_from_dict(PyObject *dict, topo_mol_ident_t *ident);




#endif // TOPO_MOL_H
#endif // SWIG
