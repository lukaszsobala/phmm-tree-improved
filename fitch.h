#ifndef _FITCH_H_
#define _FITCH_H_
int fitch_build_tree(const char *path_name_infile, const char *path_name_outfile, int tree_type, int num_threads);

// Fitch progress verbosity global (defined in fitch.c). 0=quiet,1=normal,2=verbose
extern int g_fitch_progress_level;

#endif /* _PHYLIP_H_ */
