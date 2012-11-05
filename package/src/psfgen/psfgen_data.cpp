
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "psfgen_data.h"
#include "charmm_parse_topo_defs.h"
#include "topo_mol_output.h"
#include "pdb_file_extract.h"
#include "psf_file_extract.h"
#include "topo_defs_struct.h"
#include "topo_mol_struct.h"
#include "extract_alias.h"

#ifdef WIN32
#define strcasecmp  stricmp
#define strncasecmp strnicmp
#endif



/* 
 * Provide user feedback and warnings beyond result values.
 * If we are running interactively, Tcl_Main will take care of echoing results
 * to the console.  If we run a script, we need to output the results
 * ourselves.
 */
void newhandle_msg(void *v, const char *msg) {
    puts(msg);
}

//
//
//int psfgen_test_mol(Tcl_Interp *interp, psfgen_data *data) {
//  if (! data->mol) {
//        Tcl_AppendResult(interp,
//        "\nMOLECULE MISSING!  Use resetpsf to start over.",
//        NULL);
//    return -1;
//  }
//  return 0;
//}
//
//#define PSFGEN_TEST_MOL(INTERP,DATA) \
//  if ( psfgen_test_mol(INTERP,DATA) ) return TCL_ERROR
//


void psfgen_destroy(psfgen_data *data) {

    topo_mol_destroy(data->mol);
    topo_defs_destroy(data->defs);
    stringhash_destroy(data->aliases);
    free(data);

}

//
//void psfgen_data_delete_pointer(ClientData cd, Tcl_Interp *interp) {
//  psfgen_data **dataptr = (psfgen_data **)cd;
//  free(dataptr);
//}
//
//static void count_delete_proc(ClientData data, Tcl_Interp *interp) {
//  free(data);
//}
//
psfgen_data* psfgen_create() {
    char namebuf[128];
    int *countptr;
    int id = 0;
    psfgen_data *data;


    data = (psfgen_data *)malloc(sizeof(psfgen_data));
    data->defs = topo_defs_create();
    data->aliases = stringhash_create();
    data->mol = topo_mol_create(data->defs);
    data->id = id;
    data->in_use = 0;
    data->all_caps = 1;
    *countptr = id+1;

    return data;
}
//
//void psfgen_data_reset(Tcl_Interp *interp, psfgen_data *data) {
//  topo_mol_destroy(data->mol);
//  topo_defs_destroy(data->defs);
//  stringhash_destroy(data->aliases);
//  data->defs = topo_defs_create();
//  topo_defs_error_handler(data->defs,interp,newhandle_msg);
//  data->aliases = stringhash_create();
//  data->mol = topo_mol_create(data->defs);
//  topo_mol_error_handler(data->mol,interp,newhandle_msg);
//  data->all_caps = 1;
//}
//
//int tcl_psfcontext(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_topology(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_segment(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_residue(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_mutate(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_multiply(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_coord(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_auto(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_regenerate(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_alias(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_pdb(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_coordpdb(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_guesscoord(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_readpsf(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_readplugin(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_writepsf(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_writetop(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_writepdb(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_writetranslate(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_writeplugin(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_first(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_last(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_patch(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//int tcl_resetpsf(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]);
//
//
//#if defined(PSFGENTCLDLL_EXPORTS) && defined(_WIN32)
//#  undef TCL_STORAGE_CLASS
//#  define TCL_STORAGE_CLASS DLLEXPORT
//
//#define WIN32_LEAN_AND_MEAN /* Exclude rarely-used stuff from Windows headers */
//#include <windows.h>
//
//BOOL APIENTRY DllMain( HANDLE hModule,
//                       DWORD  ul_reason_for_call,
//                       LPVOID lpReserved
//                                         )
//{
//    return TRUE;
//}
//
//EXTERN int Psfgen_Init(Tcl_Interp *interp) {
//
//#else
//
//int Psfgen_Init(Tcl_Interp *interp) {
//
//#endif
//
//  /* Create psfgen data structures; keep in interp so that other libraries
//   * can access them.
//   */
//  psfgen_data **data;
//  Tcl_SetAssocData(interp, (char *)"Psfgen_count",0,(ClientData)0);
//  data = (psfgen_data **)malloc(sizeof(psfgen_data *));
//  Tcl_SetAssocData(interp, (char *)"Psfgen_pointer",
//		psfgen_data_delete_pointer,(ClientData)data);
//  *data = psfgen_data_create(interp);
//  (*data)->in_use++;
//
//  Tcl_CreateCommand(interp,"psfcontext",tcl_psfcontext,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"topology",tcl_topology,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"readpsf",tcl_readpsf,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"readmol",tcl_readplugin,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"segment",tcl_segment,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"residue",tcl_residue,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"mutate",tcl_mutate,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"multiply",tcl_multiply,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"coord",tcl_coord,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"auto",tcl_auto,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"regenerate",tcl_regenerate,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"alias",tcl_alias,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"pdbalias",tcl_alias,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"pdb",tcl_pdb,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"coordpdb",tcl_coordpdb,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"guesscoord",tcl_guesscoord,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"writepsf",tcl_writepsf,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"writetop",tcl_writetop,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"writepdb",tcl_writepdb,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"writetranslate",tcl_writetranslate,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"writemol",tcl_writeplugin,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"first",tcl_first,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"last",tcl_last,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"patch",tcl_patch,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"resetpsf", tcl_resetpsf,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//  Tcl_CreateCommand(interp,"delatom", tcl_delatom,
//	(ClientData)data, (Tcl_CmdDeleteProc*)NULL);
//
//  Tcl_PkgProvide(interp, "psfgen", "1.5.0");
//
//#ifdef NAMD_VERSION
//  {
//    char buf[1024];
//    sprintf(buf, "puts \"PSFGEN [package require psfgen] from NAMD %s for %s\"",
//                                             NAMD_VERSION, NAMD_PLATFORM);
//    Tcl_Eval(interp, buf);
//  }
//#endif
//
//  return TCL_OK;
//}
//
//
//char* splitcolon(char *s) {
//  if ( s ) {
//    while ( *s && *s != ':' ) { ++s; }
//    if ( *s ) *(s++) = 0; else s = 0;
//  }
//  return s;
//}
//
///*
//  Old-style calls:
//    set n [psfcontext new]   $n is old context (0)
//    psfcontext new delete    old context (1) is deleted
//    set m [psfcontext]       $m is current context (2)
//    set m [psfcontext $n]    $m is old context (2)
//    psfcontext $m delete     old context (0) is deleted
//
//  How they would have to be used:
//    set mycontext [psfcontext [psfcontext new]]
//    proc a { } {
//      global mycontext
//      set oldcontext [psfcontext $mycontext]
//      set retcode [catch {
//        ... error ...
//      } result] }
//      psfcontext $oldcontext
//      if { $retcode } { error $result; } else { return $result }
//    }
//    psfcontext [psfcontext $mycontext] delete
//
//  New-style calls and usage:
//    psfcontext reset    (clears all state from current context)
//
//    set mycontext [psfcontext create]
//    psfcontext eval $mycontext { ... }
//    psfcontext delete $mycontext
//
//    psfcontext stats    (returns numbers of contexts created and destroyed)
//*/
//
//int psfgen_psfcontext(ClientData data, Tcl_Interp *interp,
//					int argc, CONST84 char *argv[]) {
//
//  int oldid, newid;
//  int delold = 0;
//  psfgen_data **cur = (psfgen_data **)data;
//  char oldidstr[128];
//  oldid = (*cur)->id;
//  sprintf(oldidstr,"%d",oldid);
//
//  if ( argc == 1 ) {
//    Tcl_SetResult(interp,oldidstr,TCL_VOLATILE);
//    return TCL_OK;
//  }
//
//  if ( argc == 2 && ! strcmp(argv[1],"stats") ) {
//    char msg[128];
//    int nc, nd, *countptr;
//    nc = 0;  nd = 0;
//    countptr = Tcl_GetAssocData(interp, "Psfgen_count", 0);
//    if (countptr) {
//      nc = countptr[0];
//      nd = countptr[1];
//    }
//    sprintf(msg,"%d created %d destroyed",countptr[0],countptr[1]);
//    Tcl_SetResult(interp,msg,TCL_VOLATILE);
//    return TCL_OK;
//  }
//
//  if ( argc == 2 && ! strcmp(argv[1],"allcaps") ) {
//    newhandle_msg(interp,"mapping names to all caps on input");
//    (*cur)->all_caps = 1;
//    return TCL_OK;
//  }
//
//  if ( argc == 2 && ! strcmp(argv[1],"mixedcase") ) {
//    newhandle_msg(interp,"preserving case of names on input");
//    (*cur)->all_caps = 0;
//    return TCL_OK;
//  }
//
//  if ( argc == 2 && ! strcmp(argv[1],"reset") ) {
//    newhandle_msg(interp,"clearing structure, topology, and aliases");
//    if ( ! (*cur)->all_caps ) {
//      newhandle_msg(interp,"mapping names to all caps on input");
//    }
//    psfgen_data_reset(interp,*cur);
//    return TCL_OK;
//  }
//
//  if ( argc == 2 && ! strcmp(argv[1],"create") ) {
//    char msg[128];
//    psfgen_data *newdata = psfgen_data_create(interp);
//    sprintf(msg,"%d",newdata->id);
//    Tcl_SetResult(interp,msg,TCL_VOLATILE);
//    return TCL_OK;
//  }
//
//  if ( argc == 3 && ! strcmp(argv[1],"delete") ) {
//    if (Tcl_GetInt(interp,argv[2],&newid) == TCL_OK) {
//      char newkey[128];
//      psfgen_data *newdata;
//      sprintf(newkey,"Psfgen_%d",newid);
//      if ( newdata = Tcl_GetAssocData(interp,newkey,0) ) {
//        if ( newdata->in_use ) {
//          Tcl_SetResult(interp,"specified context in use",TCL_VOLATILE);
//          return TCL_ERROR;
//        }
//        Tcl_DeleteAssocData(interp,newkey);
//        sprintf(newkey,"deleted %d",newid);
//        Tcl_SetResult(interp,newkey,TCL_VOLATILE);
//        return TCL_OK;
//      }
//    }
//    Tcl_SetResult(interp,"specified context does not exist",TCL_VOLATILE);
//    return TCL_ERROR;
//  }
//
//  if ( argc > 1 && ! strcmp(argv[1],"eval") ) {
//    psfgen_data *newdata, *olddata;
//    char newkey[128];
//    int retval;
//    if ( argc != 4 ) {
//      Tcl_SetResult(interp,
//        "usage: psfcontext eval ?context? { ?commmands? }",TCL_VOLATILE);
//      return TCL_ERROR;
//    }
//    if (Tcl_GetInt(interp,argv[2],&newid) != TCL_OK) {
//      Tcl_SetResult(interp,"specified context does not exist",TCL_VOLATILE);
//      return TCL_ERROR;
//    }
//    sprintf(newkey,"Psfgen_%d",newid);
//    newdata = Tcl_GetAssocData(interp,newkey,0);
//    if ( ! newdata ) {
//      Tcl_SetResult(interp,"specified context does not exist",TCL_VOLATILE);
//      return TCL_ERROR;
//    }
//    olddata = *cur;
//    *cur = newdata;
//    (*cur)->in_use++;
//
//    newdata = 0;  /* Tcl_Eval might delete this context and change *cur */
//    retval = Tcl_Eval(interp,argv[3]);
//
//    (*cur)->in_use--;
//    *cur = olddata;
//    return retval;
//  }
//
//  if ( argc == 3 ) {
//    if ( strcmp(argv[2],"delete") == 0 ) {
//      delold = 1;
//    } else {
//      Tcl_SetResult(interp,"second argument must be delete",TCL_VOLATILE);
//      psfgen_kill_mol(interp,*cur);
//      return TCL_ERROR;
//    }
//  }
//
//  if ( delold && (*cur)->in_use > 1 ) {
//    Tcl_SetResult(interp,"current context in use",TCL_VOLATILE);
//    psfgen_kill_mol(interp,*cur);
//    return TCL_ERROR;
//  }
//
//  if ( argc > 3 ) {
//    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
//    psfgen_kill_mol(interp,*cur);
//    return TCL_ERROR;
//  }
//
//  if (strcmp(argv[1],"new") == 0) {
//    psfgen_data *newdata = psfgen_data_create(interp);
//    (*cur)->in_use--;
//    *cur = newdata;
//    (*cur)->in_use++;
//  } else if (Tcl_GetInt(interp,argv[1],&newid) == TCL_OK) {
//    psfgen_data *newdata;
//    char newkey[128];
//    if ( newid == oldid ) {
//      if ( delold ) {
//        Tcl_SetResult(interp,"specified context same as current",TCL_VOLATILE);
//        psfgen_kill_mol(interp,*cur);
//        return TCL_ERROR;
//      } else {
//        Tcl_SetResult(interp,oldidstr,TCL_VOLATILE);
//        return TCL_OK;
//      }
//    }
//    sprintf(newkey,"Psfgen_%d",newid);
//    if ( (newdata = Tcl_GetAssocData(interp,newkey,0)) ) {
//      (*cur)->in_use--;
//      *cur = newdata;
//      (*cur)->in_use++;
//    } else {
//      Tcl_SetResult(interp,"specified context does not exist",TCL_VOLATILE);
//      psfgen_kill_mol(interp,*cur);
//      return TCL_ERROR;
//    }
//  } else {
//    Tcl_SetResult(interp,"first argument must be existing context or new",TCL_VOLATILE);
//    psfgen_kill_mol(interp,*cur);
//    return TCL_ERROR;
//  }
//
//  if ( delold ) {
//    char oldkey[128];
//    sprintf(oldkey,"Psfgen_%d",oldid);
//    Tcl_DeleteAssocData(interp,oldkey);
//    sprintf(oldkey,"deleted %d",oldid);
//    Tcl_SetResult(interp,oldkey,TCL_VOLATILE);
//    return TCL_OK;
//  } else {
//    Tcl_SetResult(interp,oldidstr,TCL_VOLATILE);
//    return TCL_OK;
//  }
//
//}
//

int psfgen_readpsf (psfgen_data *psf, const char* filename) {
    FILE *psf_file;
    int retval;

    printf("%s, %s.\n", __PRETTY_FUNCTION__, filename);



    /* Open psf as a binary file because the reading code uses ftell and
     fseek which do not work properly if the file is opened as text
     on Windows.  fgetpos/fsetpos misbehave in the exact same way.    */
    if ( ! ( psf_file = fopen(filename,"rb") ) ) {
        printf("ERROR: Unable to open psf file %s\n",filename);
        return -1;

    } else {
        printf("reading structure from psf file %s\n",filename);

        retval = psf_file_extract(psf->mol, psf_file, 0, newhandle_msg);
        fclose(psf_file);
    }
    return retval;
}

//
////int tcl_readplugin(ClientData data, Tcl_Interp *interp,
////					int argc, CONST84 char *argv[]) {
////  const char *filename, *pluginname;
////  const char *coorpluginname=0;
////  const char *coorfilename=0;
////  char msg[2048];
////  psfgen_data *psf = *(psfgen_data **)data;
////  char *segid=NULL;
////  int curarg;
////  int coordinatesonly=0;
////  int residuesonly=0;
////  PSFGEN_TEST_MOL(interp,psf);
////
////  if ( argc < 3 ) {
////    Tcl_SetResult(interp,"missing file format and/or input filename",TCL_VOLATILE);
////    psfgen_kill_mol(interp,psf);
////    return TCL_ERROR;
////  }
////  pluginname = argv[1];
////  filename = argv[2];
////
////  sprintf(msg,"Info: reading file %s using plugin %s", filename, pluginname);
////  newhandle_msg(interp,msg);
////
////  for (curarg=3; curarg<argc; curarg++) {
////    if (!strcmp(argv[curarg], "segment")) {
////      curarg++;
////      if (curarg<argc) {
////        segid = strtoupper(argv[curarg], psf->all_caps);
////        sprintf(msg, "Info: read mode: coordinates for segment %s", segid);
////        newhandle_msg(interp,msg);
////      }
////    } else if (!strcmp(argv[curarg], "coordinatesonly")) {
////      coordinatesonly=1;
////      newhandle_msg(interp, "Info: read mode: coordinates only");
////    } else if (!strcmp(argv[curarg], "residuesonly")) {
////      residuesonly=1;
////      newhandle_msg(interp, "Info: read mode: residue sequence only");
////    } else { /* positional arguments for second coordinate file */
////      if ( curarg == 3 ) coorpluginname = argv[3];
////      if ( curarg == 4 ) coorfilename = argv[4];
////    }
////  }
////
////  if ( coorpluginname && coorpluginname ) {
////    sprintf(msg,"Info: reading coordinates from file %s using plugin %s",
////            coorfilename, coorpluginname);
////    newhandle_msg(interp,msg);
////  }
////
////  if ( topo_mol_read_plugin(psf->mol, pluginname, filename,
////                            coorpluginname, coorfilename,
////                            segid, psf->aliases, psf->all_caps,
////                            coordinatesonly, residuesonly,
////                            interp, newhandle_msg) ) {
////    if (segid != NULL)
////      free(segid);
////    Tcl_AppendResult(interp,"ERROR: failed reading file", NULL);
////    psfgen_kill_mol(interp,psf);
////    return TCL_ERROR;
////  }
////
////  if (segid != NULL)
////    free(segid);
////
////  return TCL_OK;
////}
//
//
//
//
//int tcl_segment(ClientData data, Tcl_Interp *interp,
//					int argc, CONST84 char *argv[]) {
//  char msg[2048];
//  char *seg;
//  psfgen_data *psf = *(psfgen_data **)data;
//  PSFGEN_TEST_MOL(interp,psf);
//
//  /*
//   * special case query commands: 'segment segids', 'segment <segid> first',
//   * 'segment <segid> last', 'segment <segid> resids',
//   * 'segment <segid> residue <resid>'
//   */
//  if (argc == 2 && !strcasecmp(argv[1], "segids")) {
//    topo_mol *mol = psf->mol;
//    if (mol) {
//      int i, n=hasharray_count(mol->segment_hash);
//      for (i=0; i<n; i++) {
//        Tcl_AppendElement(interp, mol->segment_array[i]->segid);
//      }
//      return TCL_OK;
//    }
//    /* Return nothing when there's no molecule */
//  } else if (argc == 3 && !strcasecmp(argv[1], "first")) {
//    topo_mol *mol = psf->mol;
//    int segindex = (mol ?
//        hasharray_index(mol->segment_hash, argv[2]) :
//        HASHARRAY_FAIL);
//    if (segindex != HASHARRAY_FAIL) {
//      topo_mol_segment_t *seg = mol->segment_array[segindex];
//      Tcl_SetResult(interp, seg->pfirst, TCL_VOLATILE);
//      return TCL_OK;
//    }
//    Tcl_AppendResult(interp, "Invalid segid: ", argv[2], NULL);
//    return TCL_ERROR;
//  } else if (argc == 3 && !strcasecmp(argv[1], "last")) {
//    topo_mol *mol = psf->mol;
//    int segindex = (mol ?
//        hasharray_index(mol->segment_hash, argv[2]) :
//        HASHARRAY_FAIL);
//    if (segindex != HASHARRAY_FAIL) {
//      topo_mol_segment_t *seg = mol->segment_array[segindex];
//      Tcl_SetResult(interp, seg->plast, TCL_VOLATILE);
//      return TCL_OK;
//    }
//    Tcl_AppendResult(interp, "Invalid segid: ", argv[2], NULL);
//    return TCL_ERROR;
//  } else if (argc == 3 && !strcasecmp(argv[1], "resids")) {
//    topo_mol *mol = psf->mol;
//    int segindex = (mol ?
//        hasharray_index(mol->segment_hash, argv[2]) :
//        HASHARRAY_FAIL);
//    if (segindex != HASHARRAY_FAIL) {
//      topo_mol_segment_t *seg = mol->segment_array[segindex];
//      int n = hasharray_count(seg->residue_hash);
//      int i;
//      for (i=0; i<n; i++) {
//        if (hasharray_index(seg->residue_hash, seg->residue_array[i].resid) != HASHARRAY_FAIL) {
//          Tcl_AppendElement(interp, seg->residue_array[i].resid);
//        }
//      }
//      return TCL_OK;
//    }
//    Tcl_AppendResult(interp, "Invalid segid: ", argv[2], NULL);
//    return TCL_ERROR;
//  } else if (argc == 4 && !strcasecmp(argv[1], "residue")) {
//    topo_mol *mol = psf->mol;
//    int segindex = (mol ?
//        hasharray_index(mol->segment_hash, argv[2]) :
//        HASHARRAY_FAIL);
//    if (segindex != HASHARRAY_FAIL) {
//      topo_mol_segment_t *seg = mol->segment_array[segindex];
//      int resindex = hasharray_index(seg->residue_hash, argv[3]);
//      if (resindex == HASHARRAY_FAIL) {
//        Tcl_AppendResult(interp, "Invalid resid '", argv[3], "' for segment '",
//            argv[1], "'.", NULL);
//        return TCL_ERROR;
//      }
//      Tcl_SetResult(interp, seg->residue_array[resindex].name, TCL_VOLATILE);
//      return TCL_OK;
//    }
//    Tcl_AppendResult(interp, "Invalid segid: ", argv[2], NULL);
//    return TCL_ERROR;
//  } else if (argc == 4 && !strcasecmp(argv[1], "atoms")) {
//    topo_mol *mol = psf->mol;
//    int segindex = (mol ?
//        hasharray_index(mol->segment_hash, argv[2]) :
//        HASHARRAY_FAIL);
//    if (segindex != HASHARRAY_FAIL) {
//      topo_mol_atom *atoms;
//      topo_mol_segment_t *seg = mol->segment_array[segindex];
//      int resindex = hasharray_index(seg->residue_hash, argv[3]);
//      if (resindex == HASHARRAY_FAIL) {
//        Tcl_AppendResult(interp, "Invalid resid '", argv[3], "' for segment '",
//            argv[1], "'.", NULL);
//        return TCL_ERROR;
//      }
//      atoms = seg->residue_array[resindex].atoms;
//      while (atoms) {
//        Tcl_AppendElement(interp, atoms->name);
//        atoms = atoms->next;
//      }
//      return TCL_OK;
//    }
//    Tcl_AppendResult(interp, "Invalid segid: ", argv[2], NULL);
//    return TCL_ERROR;
//  } else if (argc == 5 && !strcasecmp(argv[1], "coordinates")) {
//    topo_mol *mol = psf->mol;
//    int segindex = (mol ?
//        hasharray_index(mol->segment_hash, argv[2]) :
//        HASHARRAY_FAIL);
//    if (segindex != HASHARRAY_FAIL) {
//      topo_mol_atom *atoms;
//      topo_mol_segment_t *seg = mol->segment_array[segindex];
//      int resindex = hasharray_index(seg->residue_hash, argv[3]);
//      if (resindex == HASHARRAY_FAIL) {
//        Tcl_AppendResult(interp, "Invalid resid '", argv[3], "' for segment '",
//            argv[1], "'.", NULL);
//        return TCL_ERROR;
//      }
//      /*
//       * XXX Ouch, no hasharray for atom names
//       */
//      atoms = seg->residue_array[resindex].atoms;
//      while (atoms) {
//        if (!strcmp(atoms->name, argv[4])) {
//#if TCL_MINOR_VERSION >= 6
//          char buf[512];
//          sprintf(buf, "%f %f %f", atoms->x, atoms->y, atoms->z);
//          Tcl_AppendResult(interp, buf, NULL);
//#else
//          sprintf(interp->result, "%f %f %f", atoms->x, atoms->y, atoms->z);
//#endif
//          return TCL_OK;
//        }
//        atoms = atoms->next;
//      }
//      Tcl_AppendResult(interp, "Invalid atom name '", argv[4],
//          "' for segid '", argv[2], "', resid '", argv[3], "'.", NULL);
//      return TCL_ERROR;
//    }
//    Tcl_AppendResult(interp, "Invalid segid: ", argv[2], NULL);
//    return TCL_ERROR;
//  }
//
//  /*
//   * Fall through to segment-building commands
//   */
//
//  if ( argc < 3 ) {
//    Tcl_SetResult(interp,"arguments: segname { commmands }",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  if ( argc > 3 ) {
//    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  seg=strtoupper(argv[1], psf->all_caps);
//  if ( strlen(seg) > 7 ) {
//    Tcl_SetResult(interp,"segment name more than 7 characters",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//
//  sprintf(msg,"building segment %s",seg);
//  newhandle_msg(interp,msg);
//  if ( topo_mol_segment(psf->mol,seg) ) {
//    free(seg);
//    Tcl_AppendResult(interp,"ERROR: failed on segment",NULL);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  free(seg);
//
//  if ( Tcl_Eval(interp,argv[2]) != TCL_OK ) {
//    Tcl_AppendResult(interp,"\nERROR: failed while building segment",NULL);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//
//  newhandle_msg(interp,"Info: generating structure...");
//  if ( topo_mol_end(psf->mol) ) {
//    Tcl_AppendResult(interp,"ERROR: failed on end of segment",NULL);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  newhandle_msg(interp, "Info: segment complete.");
//  return TCL_OK;
//}
//
//int tcl_residue(ClientData data, Tcl_Interp *interp,
//					int argc, CONST84 char *argv[]) {
//  char *resid, *resname, *chain;
//  psfgen_data *psf = *(psfgen_data **)data;
//  PSFGEN_TEST_MOL(interp,psf);
//
//  if ( argc < 3 ) {
//    Tcl_SetResult(interp,"arguments: resid resname ?chain?",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  if ( argc > 4 ) {
//    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  resid=strtoupper(argv[1], psf->all_caps);
//  resname=strtoupper(argv[2], psf->all_caps);
//  chain=strtoupper(argc==4 ? argv[3] : "", psf->all_caps);
//
//  if ( topo_mol_residue(psf->mol,resid,resname,chain) ) {
//    free(resid);
//    free(resname);
//    Tcl_AppendResult(interp,"ERROR: failed on residue",NULL);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  free(resid);
//  free(resname);
//  free(chain);
//  return TCL_OK;
//}
//
//int tcl_mutate(ClientData data, Tcl_Interp *interp,
//					int argc, CONST84 char *argv[]) {
//  psfgen_data *psf = *(psfgen_data **)data;
//  char *resid, *resname;
//  PSFGEN_TEST_MOL(interp,psf);
//
//  if ( argc < 3 ) {
//    Tcl_SetResult(interp,"arguments: resid resname",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  if ( argc > 3 ) {
//    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  resid=strtoupper(argv[1], psf->all_caps);
//  resname=strtoupper(argv[2], psf->all_caps);
//
//  if ( topo_mol_mutate(psf->mol,resid, resname) ) {
//    free(resid);
//    free(resname);
//    Tcl_AppendResult(interp,"ERROR: failed on mutate",NULL);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  free(resid);
//  free(resname);
//
//  return TCL_OK;
//}
//
//int tcl_multiply(ClientData data, Tcl_Interp *interp,
//                                        int argc, CONST84 char *argv[]) {
//  int i, ncopies, ierr;
//  topo_mol_ident_t *targets;
//  char **tmp;
//  char msg[2048];
//  psfgen_data *psf = *(psfgen_data **)data;
//  PSFGEN_TEST_MOL(interp,psf);
//
//  if ( argc<3 || Tcl_GetInt(interp,argv[1],&ncopies) != TCL_OK || ncopies<2 ) {
//    Tcl_SetResult(interp,"arguments: ncopies segid?:resid?:atomname? ...",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//
//  targets = (topo_mol_ident_t *) Tcl_Alloc((argc-2)*sizeof(topo_mol_ident_t));
//  if ( ! targets ) {
//    Tcl_SetResult(interp,"memory allocation failed",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  tmp = (char **) Tcl_Alloc((argc-2)*sizeof(char *));
//  if (!tmp) {
//    Tcl_Free((char *)targets);
//    Tcl_SetResult(interp,"memory allocation failed",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  sprintf(msg,"generating %d copies of selected atoms",ncopies);
//  newhandle_msg(interp,msg);
//  for ( i=2; i<argc; ++i ) {
//    char *ctmp;
//    tmp[i-2] = strtoupper(argv[i], psf->all_caps);
//    targets[i-2].segid = ctmp = tmp[i-2];
//    targets[i-2].resid = ctmp = splitcolon(ctmp);
//    targets[i-2].aname = splitcolon(ctmp);
//  }
//  ierr = topo_mol_multiply_atoms(psf->mol,targets,(argc-2),ncopies);
//  for (i=2; i<argc; ++i) free(tmp[i-2]);
//  Tcl_Free((char *)tmp);
//  Tcl_Free((char *)targets);
//  if (ierr) {
//    sprintf(msg,"ERROR: failed to multiply atoms (error=%d)",ierr);
//    Tcl_SetResult(interp,msg,TCL_VOLATILE);
//    /* Tcl_AppendResult(interp,"ERROR: failed to multiply atoms",NULL); */
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//
//  return TCL_OK;
//}
//
////int tcl_coord(ClientData data, Tcl_Interp *interp,
////					int argc, CONST84 char *argv[]) {
////  double x,y,z;
////  topo_mol_ident_t target;
////  char *segid, *resid, *atomname;
////  int rc;
////  psfgen_data *psf = *(psfgen_data **)data;
////  PSFGEN_TEST_MOL(interp,psf);
////
////  if ( argc < 5 ) {
////    Tcl_SetResult(interp,"arguments: segid resid atomname { x y z }",TCL_VOLATILE);
////    psfgen_kill_mol(interp,psf);
////    return TCL_ERROR;
////  }
////  if ( argc > 5 ) {
////    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
////    psfgen_kill_mol(interp,psf);
////    return TCL_ERROR;
////  }
////  if ( sscanf(argv[4],"%lf %lf %lf",&x,&y,&z) != 3 ) {
////    Tcl_SetResult(interp,"arguments: segid resid atomname { x y z }",TCL_VOLATILE);
////    psfgen_kill_mol(interp,psf);
////    return TCL_ERROR;
////  }
////  segid=strtoupper(argv[1], psf->all_caps);
////  resid=strtoupper(argv[2], psf->all_caps);
////  atomname=strtoupper(argv[3], psf->all_caps);
////  target.segid = segid;
////  target.resid = resid;
////  target.aname = atomname;
////  rc = topo_mol_set_xyz(psf->mol,&target,x,y,z);
////  free(segid);
////  free(resid);
////  free(atomname);
////  if (rc) {
////    Tcl_AppendResult(interp,"ERROR: failed on coord",NULL);
////    psfgen_kill_mol(interp,psf);
////    return TCL_ERROR;
////  }
////
////  return TCL_OK;
////}
//
//
//int tcl_auto(ClientData data, Tcl_Interp *interp,
//					int argc, CONST84 char *argv[]) {
//  int i, angles, dihedrals;
//  psfgen_data *psf = *(psfgen_data **)data;
//  PSFGEN_TEST_MOL(interp,psf);
//
//  if ( argc < 2 ) {
//    Tcl_SetResult(interp,"arguments: ?angles? ?dihedrals? ?none?",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//
//  angles = 0;  dihedrals = 0;
//  for ( i = 1; i < argc; ++i ) {
//    if ( ! strcmp(argv[i],"angles") ) angles = 1;
//    else if ( ! strcmp(argv[i],"dihedrals") ) dihedrals = 1;
//    else if ( strcmp(argv[i],"none") ) {
//      Tcl_SetResult(interp,"arguments: ?angles? ?dihedrals? ?none?",TCL_VOLATILE);
//      psfgen_kill_mol(interp,psf);
//      return TCL_ERROR;
//    }
//  }
//
//  if ( angles ) newhandle_msg(interp,"enabling angle autogeneration");
//  else newhandle_msg(interp,"disabling angle autogeneration");
//  if ( topo_mol_segment_auto_angles(psf->mol,angles) ) {
//    Tcl_AppendResult(interp,"ERROR: failed setting angle autogen",NULL);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//
//  if ( dihedrals ) newhandle_msg(interp,"enabling dihedral autogeneration");
//  else newhandle_msg(interp,"disabling dihedral autogeneration");
//  if ( topo_mol_segment_auto_dihedrals(psf->mol,dihedrals) ) {
//    Tcl_AppendResult(interp,"ERROR: failed setting dihedral autogen",NULL);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//
//  return TCL_OK;
//}
//
//
//int tcl_regenerate(ClientData data, Tcl_Interp *interp,
//					int argc, CONST84 char *argv[]) {
//  int i, angles, dihedrals;
//  psfgen_data *psf = *(psfgen_data **)data;
//  PSFGEN_TEST_MOL(interp,psf);
//
//  if ( argc < 2 ) {
//    Tcl_SetResult(interp,"arguments: ?angles? ?dihedrals? ?none?",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//
//  angles = 0;  dihedrals = 0;
//  for ( i = 1; i < argc; ++i ) {
//    if ( ! strcmp(argv[i],"angles") ) angles = 1;
//    else if ( ! strcmp(argv[i],"dihedrals") ) dihedrals = 1;
//    else {
//      Tcl_SetResult(interp,"arguments: ?angles? ?dihedrals?",TCL_VOLATILE);
//      psfgen_kill_mol(interp,psf);
//      return TCL_ERROR;
//    }
//  }
//
//  if ( angles ) {
//    newhandle_msg(interp,"regenerating all angles");
//    if ( topo_mol_regenerate_angles(psf->mol) ) {
//      Tcl_AppendResult(interp,"ERROR: angle regeneration failed",NULL);
//      psfgen_kill_mol(interp,psf);
//      return TCL_ERROR;
//    }
//  }
//
//  if ( dihedrals ) {
//    newhandle_msg(interp,"regenerating all dihedrals");
//    if ( topo_mol_regenerate_dihedrals(psf->mol) ) {
//      Tcl_AppendResult(interp,"ERROR: dihedral regeneration failed",NULL);
//      psfgen_kill_mol(interp,psf);
//      return TCL_ERROR;
//    }
//  }
//
//  return TCL_OK;
//}
//
//int tcl_alias(ClientData data, Tcl_Interp *interp,
//					int argc, CONST84 char *argv[]) {
//  char msg[2048];
//  psfgen_data *psf = *(psfgen_data **)data;
//  int rc;
//  PSFGEN_TEST_MOL(interp,psf);
//
//  if ( argc < 2 ) {
//    Tcl_SetResult(interp,"arguments: atom | residue ...",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//
//  if ( ! strcmp(argv[1],"residue") ) {
//    char *altres, *realres;
//    if ( argc < 4 ) {
//      Tcl_SetResult(interp,"arguments: residue altres realres",TCL_VOLATILE);
//      psfgen_kill_mol(interp,psf);
//      return TCL_ERROR;
//    }
//    altres=strtoupper(argv[2], psf->all_caps);
//    realres=strtoupper(argv[3], psf->all_caps);
//    sprintf(msg,"aliasing residue %s to %s",argv[2],argv[3]);
//    newhandle_msg(interp,msg);
//    rc = extract_alias_residue_define(psf->aliases,altres, realres);
//    free(altres);
//    free(realres);
//    if (rc) {
//      Tcl_AppendResult(interp,"ERROR: failed on residue alias",NULL);
//      psfgen_kill_mol(interp,psf);
//      return TCL_ERROR;
//    }
//  } else if ( ! strcmp(argv[1],"atom") ) {
//    char *resname, *altatom, *realatom;
//    if ( argc < 5 ) {
//      Tcl_SetResult(interp,"arguments: atom resname altatom realatom",TCL_VOLATILE);
//      psfgen_kill_mol(interp,psf);
//      return TCL_ERROR;
//    }
//    resname=strtoupper(argv[2], psf->all_caps);
//    altatom=strtoupper(argv[3], psf->all_caps);
//    realatom=strtoupper(argv[4], psf->all_caps);
//    sprintf(msg,"aliasing residue %s atom %s to %s",argv[2],argv[3],argv[4]);
//    newhandle_msg(interp,msg);
//    rc=extract_alias_atom_define(psf->aliases,resname,altatom,realatom);
//    free(resname);
//    free(altatom);
//    free(realatom);
//    if (rc) {
//      Tcl_AppendResult(interp,"ERROR: failed on atom alias",NULL);
//      psfgen_kill_mol(interp,psf);
//      return TCL_ERROR;
//    }
//  }
//
//  return TCL_OK;
//}
//
//int tcl_pdb(ClientData data, Tcl_Interp *interp,
//					int argc, CONST84 char *argv[]) {
//  FILE *res_file;
//  const char *filename;
//  char msg[2048];
//  psfgen_data *psf = *(psfgen_data **)data;
//  PSFGEN_TEST_MOL(interp,psf);
//
//  if ( argc == 1 ) {
//    Tcl_SetResult(interp,"no pdb file specified",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  if ( argc > 2 ) {
//    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  filename = argv[1];
//  if ( ! ( res_file = fopen(filename,"r") ) ) {
//    sprintf(msg,"ERROR: Unable to open pdb file %s to read residues\n",filename);
//    Tcl_SetResult(interp,msg,TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  } else {
//    sprintf(msg,"reading residues from pdb file %s",filename);
//    newhandle_msg(interp,msg);
//    if ( pdb_file_extract_residues(psf->mol,res_file,psf->aliases,psf->all_caps,interp,newhandle_msg) ) {
//      Tcl_AppendResult(interp,"ERROR: failed on reading residues from pdb file",NULL);
//      fclose(res_file);
//      psfgen_kill_mol(interp,psf);
//      return TCL_ERROR;
//    }
//    fclose(res_file);
//  }
//
//  return TCL_OK;
//}
//
//int tcl_coordpdb(ClientData data, Tcl_Interp *interp,
//					int argc, CONST84 char *argv[]) {
//  FILE *res_file;
//  const char *filename;
//  char msg[2048];
//  int rc;
//  psfgen_data *psf = *(psfgen_data **)data;
//  PSFGEN_TEST_MOL(interp,psf);
//
//  if ( argc < 2 ) {
//    Tcl_SetResult(interp,"arguments: pdbfile ?segid?",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  if ( argc > 3 ) {
//    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  filename = argv[1];
//  if ( ! ( res_file = fopen(filename,"r") ) ) {
//    sprintf(msg,"ERROR: Unable to open pdb file %s to read coordinates\n",filename);
//    Tcl_SetResult(interp,msg,TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  } else {
//    char *segid;
//    if (argc == 3) {
//      /* Read only coordinates for given segid */
//      sprintf(msg,"reading coordinates from pdb file %s for segment %s",filename,argv[2]);
//      newhandle_msg(interp,msg);
//      segid = strtoupper(argv[2], psf->all_caps);
//    } else {
//      /* Read all segid's in pdb file */
//      sprintf(msg,"reading coordinates from pdb file %s",filename);
//      newhandle_msg(interp,msg);
//      segid = NULL;
//    }
//    rc=pdb_file_extract_coordinates(psf->mol,res_file,segid,psf->aliases,psf->all_caps,interp,newhandle_msg);
//    if (segid) free(segid);
//    if (rc) {
//      Tcl_AppendResult(interp,"ERROR: failed on reading coordinates from pdb file",NULL);
//      fclose(res_file);
//      psfgen_kill_mol(interp,psf);
//      return TCL_ERROR;
//    }
//    fclose(res_file);
//  }
//
//  return TCL_OK;
//
//}
//
//int tcl_guesscoord(ClientData data, Tcl_Interp *interp,
//					int argc, CONST84 char *argv[]) {
//  psfgen_data *psf = *(psfgen_data **)data;
//  PSFGEN_TEST_MOL(interp,psf);
//  if ( argc > 1 ) {
//    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  if ( topo_mol_guess_xyz(psf->mol) ) {
//    Tcl_AppendResult(interp,"ERROR: failed on guessing coordinates",NULL);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//
//  return TCL_OK;
//}
//
//int tcl_writepsf(ClientData data, Tcl_Interp *interp,
//					int argc, CONST84 char *argv[]) {
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
//  charmmfmt = 0;
//  nocmap = 0;
//  for ( i = 1; i < argc-1; ++i ) {
//    if ( strcmp(argv[i],"charmm") == 0 ) charmmfmt = 1;
//    else if ( strcmp(argv[i],"x-plor") == 0 ) charmmfmt = 0;
//    else if ( strcmp(argv[i],"cmap") == 0 ) nocmap = 0;
//    else if ( strcmp(argv[i],"nocmap") == 0 ) nocmap = 1;
//    else {
//      sprintf(msg,"ERROR: Unknown psf file format %s (not charmm or x-plor, cmap or nocmap).\n",argv[i]);
//      Tcl_SetResult(interp,msg,TCL_VOLATILE);
//      psfgen_kill_mol(interp,psf);
//      return TCL_ERROR;
//    }
//  }
//  filename = argv[argc-1];
//
//  if ( ! ( res_file = fopen(filename,"w") ) ) {
//    sprintf(msg,"ERROR: Unable to open psf file %s to write structure\n",filename);
//    Tcl_SetResult(interp,msg,TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  sprintf(msg,"Info: writing psf file %s%s%s",filename,
//                nocmap?" without cross-terms":"",
//                charmmfmt?" in CHARMM format":"");
//  newhandle_msg(interp,msg);
//  if ( topo_mol_write_psf(psf->mol,res_file,charmmfmt,nocmap,interp,newhandle_msg) ) {
//    Tcl_AppendResult(interp,"ERROR: failed on writing structure to psf file",NULL);
//    fclose(res_file);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  fclose(res_file);
//  newhandle_msg(interp, "Info: psf file complete.");
//
//  return TCL_OK;
//}
//
//int tcl_writetranslate(ClientData data, Tcl_Interp *interp,
//					int argc, CONST84 char *argv[]) {
//  FILE *res_file;
//  const char *filename;
//  int charmmfmt, nocmap, i;
//  char msg[2048];
//  psfgen_data *psf = *(psfgen_data **)data;
//  PSFGEN_TEST_MOL(interp,psf);
//
//  if ( argc == 1 ) {
//    Tcl_SetResult(interp,"no translate file specified",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  if ( argc > 2 ) {
//    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//
//  filename = argv[argc-1];
//
//  if ( ! ( res_file = fopen(filename,"w") ) ) {
//    sprintf(msg,"ERROR: Unable to open file %s to write translate table\n",filename);
//    Tcl_SetResult(interp,msg,TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//sprintf(msg,"Info: writing translate table to file %s",filename );
//  newhandle_msg(interp,msg);
//  if ( topo_mol_write_translate(psf->mol,res_file,interp,newhandle_msg) ) {
//    Tcl_AppendResult(interp,"ERROR: failed on writing translate table to file",NULL);
//    fclose(res_file);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  fclose(res_file);
//  newhandle_msg(interp, "Info: translate table complete.");
//
//  return TCL_OK;
//}
//
//
//int tcl_writetop(ClientData data, Tcl_Interp *interp,
//					int argc, CONST84 char *argv[]) {
//  FILE *res_file;
//  const char *filename;
//  char msg[2048];
//  psfgen_data *psf = *(psfgen_data **)data;
//  PSFGEN_TEST_MOL(interp,psf);
//
//  if ( argc == 1 ) {
//    Tcl_SetResult(interp,"no top file specified",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  if ( argc > 2 ) {
//    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//
//  filename = argv[argc-1];
//
//  if ( ! ( res_file = fopen(filename,"w") ) ) {
//    sprintf(msg,"ERROR: Unable to open top file %s to write structure\n",filename);
//    Tcl_SetResult(interp,msg,TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  sprintf(msg,"Info: writing top file %s",filename);
//  newhandle_msg(interp,msg);
//  if ( topo_mol_write_top(psf->mol,res_file,interp,newhandle_msg) ) {
//    Tcl_AppendResult(interp,"ERROR: failed on writing structure to top file",NULL);
//    fclose(res_file);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  fclose(res_file);
//  newhandle_msg(interp, "Info: top file complete.");
//
//  return TCL_OK;
//}
//
//

//
//
//int tcl_writeplugin(ClientData data, Tcl_Interp *interp,
//                    int argc, CONST84 char *argv[]) {
//  const char *filename, *pluginname;
//  char msg[2048];
//  struct image_spec images;
//  psfgen_data *psf = *(psfgen_data **)data;
//  PSFGEN_TEST_MOL(interp,psf);
//
//  images.na = 1; images.nb = 1; images.nc = 1;
//  images.ax = 0.; images.ay = 0.; images.az = 0.;
//  images.bx = 0.; images.by = 0.; images.bz = 0.;
//  images.cx = 0.; images.cy = 0.; images.cz = 0.;
//
//  if ( argc == 1 ) {
//    Tcl_SetResult(interp,"arguments: format filename ?na { x y z }? ?nb { x y z }? ?nc { x y z }?",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  if ( argc < 3 ) {
//    Tcl_SetResult(interp,"missing file format and/or output filename",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//
//  pluginname = argv[1];
//  filename = argv[2];
//
//  if ( argc > 3 ) {
//    if ( sscanf(argv[3],"%d",&images.na) != 1 || images.na < 1 ) {
//      Tcl_SetResult(interp,"image count not a positive integer",TCL_VOLATILE);
//      psfgen_kill_mol(interp,psf);
//      return TCL_ERROR;
//    }
//    if ( argc == 4 ) {
//      Tcl_SetResult(interp,"image count without offset vector",TCL_VOLATILE);
//      psfgen_kill_mol(interp,psf);
//      return TCL_ERROR;
//    }
//    if ( sscanf(argv[4],"%lf %lf %lf",&images.ax,&images.ay,&images.az) != 3 ) {
//      Tcl_SetResult(interp,"bad image offset vector format",TCL_VOLATILE);
//      psfgen_kill_mol(interp,psf);
//      return TCL_ERROR;
//    }
//  }
//
//  if ( argc > 5 ) {
//    if ( sscanf(argv[5],"%d",&images.nb) != 1 || images.nb < 1 ) {
//      Tcl_SetResult(interp,"image count not a positive integer",TCL_VOLATILE);
//      psfgen_kill_mol(interp,psf);
//      return TCL_ERROR;
//    }
//    if ( argc == 6 ) {
//      Tcl_SetResult(interp,"image count without offset vector",TCL_VOLATILE);
//      psfgen_kill_mol(interp,psf);
//      return TCL_ERROR;
//    }
//    if ( sscanf(argv[6],"%lf %lf %lf",&images.bx,&images.by,&images.bz) != 3 ) {
//      Tcl_SetResult(interp,"bad image offset vector format",TCL_VOLATILE);
//      psfgen_kill_mol(interp,psf);
//      return TCL_ERROR;
//    }
//  }
//
//  if ( argc > 7 ) {
//    if ( sscanf(argv[7],"%d",&images.nc) != 1 || images.nc < 1 ) {
//      Tcl_SetResult(interp,"image count not a positive integer",TCL_VOLATILE);
//      psfgen_kill_mol(interp,psf);
//      return TCL_ERROR;
//    }
//    if ( argc == 8 ) {
//      Tcl_SetResult(interp,"image count without offset vector",TCL_VOLATILE);
//      psfgen_kill_mol(interp,psf);
//      return TCL_ERROR;
//    }
//    if ( sscanf(argv[8],"%lf %lf %lf",&images.cx,&images.cy,&images.cz) != 3 ) {
//      Tcl_SetResult(interp,"bad image offset vector format",TCL_VOLATILE);
//      psfgen_kill_mol(interp,psf);
//      return TCL_ERROR;
//    }
//  }
//
//  if ( argc > 9 ) {
//    Tcl_SetResult(interp,"too many arguments specified",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//
//  sprintf(msg,"Info: writing file %s using plugin %s", filename, pluginname);
//  newhandle_msg(interp,msg);
//  if ( topo_mol_write_plugin(psf->mol, pluginname, filename, &images, interp, newhandle_msg) ) {
//    Tcl_AppendResult(interp,"ERROR: failed writing to file", NULL);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  newhandle_msg(interp, "Info: file complete.");
//
//  return TCL_OK;
//}
//
//
//
//
//int tcl_patch(ClientData data, Tcl_Interp *interp,
//					int argc, CONST84 char *argv[]) {
//  int i, j, rc, ipres, listall=0;
//  topo_mol_ident_t targets[10];
//  char *tmp[10];
//  char *pres;
//  char msg[2048];
//  topo_mol_patch_t *patch;
//  topo_mol_patchres_t *patchres;
//  Tcl_Obj *tcl_result;
//  psfgen_data *psf = *(psfgen_data **)data;
//  PSFGEN_TEST_MOL(interp,psf);
//
//  tcl_result = Tcl_NewListObj(0, NULL);
//
//  if (argc == 3 && !strcasecmp(argv[1], "targets")) {
//    return tcl_num_patch_targets(psf, interp, argv[2]);
//  }
//
//  if ( argc == 2 && (!strcasecmp(argv[1], "list") || !strcasecmp(argv[1], "listall"))) {
//    if (!strcasecmp(argv[1], "listall")) listall = 1;
//    for ( patch = psf->mol->patches; patch; patch = patch->next ) {
//      Tcl_Obj *patchlist = Tcl_NewListObj(0,NULL);
//      ipres = 0;
//      /* Only list all patches when 'patch listall' was invoked */
//      if (patch->deflt && !listall) continue;
//
//      for ( patchres = patch->patchresids; patchres; patchres = patchres->next ) {
//	/* Test the existence of segid:resid for the patch */
//	if (!topo_mol_validate_patchres(psf->mol,patch->pname,patchres->segid, patchres->resid)) {
//	  break;
//	};
//
//	if (ipres==0) {
//	  Tcl_ListObjAppendElement(interp, patchlist, Tcl_NewStringObj(patch->pname, -1));
//	}
//	Tcl_ListObjAppendElement(interp, patchlist, Tcl_NewStringObj(patchres->segid, -1));
//	Tcl_ListObjAppendElement(interp, patchlist, Tcl_NewStringObj(patchres->resid, -1));
//	ipres++;
//      }
//      Tcl_ListObjAppendElement(interp, tcl_result, patchlist);
//    }
//    Tcl_SetObjResult(interp, tcl_result);
//    return TCL_OK;
//  }
//
//  if ( argc < 2 ) {
//    Tcl_SetResult(interp,"arguments: list | presname segid:resid ...",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  if ( argc > 10 ) {
//    Tcl_SetResult(interp,"too many targets for patch",TCL_VOLATILE);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//  pres=strtoupper(argv[1], psf->all_caps);
//  sprintf(msg,"applying patch %s to %d residues",pres,(argc-2));
//  newhandle_msg(interp,msg);
//  for ( i=2; i<argc; ++i ) {
//    tmp[i-2]=strtoupper(argv[i], psf->all_caps);
//    targets[i-2].segid = tmp[i-2];
//    targets[i-2].resid = splitcolon(tmp[i-2]);
//    targets[i-2].aname = 0;
//    if ( ! targets[i-2].resid ) {
//      for (j=0; j<i-2; j++) free(tmp[j]);
//      sprintf(msg,"ERROR: resid missing from patch target %s",tmp[i-2]);
//      Tcl_SetResult(interp,msg,TCL_VOLATILE);
//      psfgen_kill_mol(interp,psf);
//      return TCL_ERROR;
//    }
//  }
//  rc=topo_mol_patch(psf->mol,targets,(argc-2),pres,0,0,0,0);
//  free(pres);
//  for (j=0; j<argc-2; j++) free(tmp[j]);
//  if (rc) {
//    Tcl_AppendResult(interp,"ERROR: failed to apply patch",NULL);
//    psfgen_kill_mol(interp,psf);
//    return TCL_ERROR;
//  }
//
//  return TCL_OK;
//}
//
//int tcl_resetpsf(ClientData data, Tcl_Interp *interp, int argc, CONST84 char *argv[]) {
//  psfgen_data *psf = *(psfgen_data **)data;
//
//  newhandle_msg(interp,"clearing structure, preserving topology and aliases");
//  topo_mol_destroy(psf->mol);
//  psf->mol = topo_mol_create(psf->defs);
//  topo_mol_error_handler(psf->mol,interp,newhandle_msg);
//
//  return TCL_OK;
//}
//



//
//#endif
//