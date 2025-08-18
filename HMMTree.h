#pragma once

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <stdio.h> /* fprintf, stderr, BUFSIZ */
#include <fcntl.h> /* open, flags */

//#include <hash_map>
#include <unordered_map>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <sstream>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <float.h>

//kitsch.h
extern "C"{
#include "kitsch.h"
#include "fitch.h"
#include "neighbor.h"
#include "upgma.h"
}

// Monotonic wall-clock time in milliseconds (for timing/logging only)
static inline long long now_millis() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return static_cast<long long>(tv.tv_sec) * 1000LL + static_cast<long long>(tv.tv_usec) / 1000LL;
}

const unsigned int MAX_BUF_LEN = 1024;             // buffer length for program existence checks (e.g., popen/which)
const unsigned int BUFFER_SIZE = 1024;             // buffer size for file copy
//const double CLUSTER_OK_THRESHOLD = 0.90;        // threshold to evaluate clustering quality
const double DBL_MAX_USER = 0.001;                 // sentinel value used as a "null" score
const unsigned int PTH_MATRIX_OUT = 0;             // output path for the final distance matrix (path_filename::path_get, parameter index_to_set)
const unsigned int PTH_UNALG_SEQ = 1;              // path for unaligned sequences (path_filename::path_get, parameter index_to_set)
const unsigned int PTH_F_ALG_CLU = 2;              // path for clusters and alignments (path_filename::path_get, parameter index_to_set)
const unsigned int PTH_F_HMM = 3;                  // path for HMM files (path_filename::path_get, parameter index_to_set)
const unsigned int PTH_F_EACH_2_RESULT = 4;        // path for pairwise HMM comparison result files (path_filename::path_get, parameter index_to_set)
const unsigned int PTH_F_1_LIARBRY_RESULT = 5;     // path for one-vs-library HMM comparison result files (path_filename::path_get, parameter index_to_set)

const std::string STR_ARGUMENTS_ERROR_MSG = "pHMM-Tree\n"
"(c) Yin Lab (NIU) and Zhang Lab (NKU) 2016\n"
"\n"
"Usage: phmm-tree <mode> <option> [-id <identity>] [-prc_hit <value>] [-acc] [analysis_methods] [file_path or file_path_name]\n"
"\n"
"USAGE (structured)\n"
"  phmm-tree <mode> <input_option> [flags] [analysis_methods] <path>\n"
"\n"
"MODES\n"
"  -prc      Compare profile HMMs in HMMER format using PRC.\n"
"  -hhsuite  Compare profile HMMs in HHM format using hhsuite.\n"
"\n"
"INPUT OPTIONS (choose one)\n"
"  -uals       Unaligned FASTA file with at least 9 sequences.\n"
"               - Requires -id <0.1..1.0> to set usearch identity for clustering.\n"
"               - Provide a file path/name.\n"
"               - prc_hit <value> sets PRC parameter (default: 10; 0 uses PRC default 100).\n"
"               - acc uses the 'ACC' key in matrix/tree if present (default is 'NAME').\n"
"               - lib runs PRC in library mode (default is pairwise mode).\n"
"\n"
"  -als        Directory containing at least three aligned FASTA files.\n"
"               - prc_hit and -acc are supported (same as in -uals).\n"
"               - Provide only a directory path.\n"
"\n"
"  -als_phmms  Two folders: one with alignments and one with HMMER2.x/3.x profile HMM files.\n"
"               - Options are the same as for -hmms.\n"
"\n"
"  -als_phhms  Two folders: one with alignments and one with HHsuite HHM files.\n"
"               - Options are the same as for -hhms.\n"
"\n"
"  -hmms       Directory with at least three HMMER2.x/3.x profile HMM files.\n"
"               - Options are the same as for -als.\n"
"\n"
"  -hhms       Directory with at least three HHsuite HHM files.\n"
"               - Options are the same as for -als.\n"
"\n"
"FLAGS\n"
"  -id <identity>         Usearch identity (0.1..1.0) for -uals.\n"
"  -prc_hit <value>       PRC parameter (default: 10; 0 uses PRC default 100).\n"
"  -acc                   Use ACC in matrix/tree if available (default: NAME).\n"
"  -lib                   PRC library mode (for -prc; default is pairwise).\n"
"  -prc_backend <auto|legacy|prcx>  Backend for PRC when HMMER3 HMMs are detected: \n"
"                               auto (default) uses prcX if available, else legacy convert+prc;\n"
"                               legacy forces convert-to-HMMER2 then prc;\n"
"                               prcx forces direct prcX without conversion.\n"
"\n"
"ANALYSIS METHODS (optional; if none specified, all run)\n"
"  -fitch    Run Fitch-Margoliash (both f-m and minimum evolution).\n"
"  -kitsch   Run Kitsch (contemporary tips; both f-m and minimum evolution).\n"
"  -upgma    Run UPGMA (Unweighted Pair Group Method with Arithmetic Mean).\n"
"  -nj       Run Neighbor-Joining.\n"
"  -fm       Run only f-m variants (affects Fitch and Kitsch).\n"
"  -min      Run only minimum evolution variants (affects Fitch and Kitsch).\n"
"\n"
"THREAD CONTROL (optional; 0 = auto-detect)\n"
"  -prc_threads <num>     Number of threads for PRC distance calculations (default: auto-detect).\n"
"  -phylo_threads <num>   Number of threads for phylogenetic algorithms (default: 1).\n"
"  -phylo_concurrent_threads <num>  Number of phylogenetic analyses to run concurrently (default: auto-detect).\n"

"\n"
"EXAMPLES\n"
"  PRC + unaligned FASTA:\n"
"    phmm-tree -prc -uals -id 0.9 ./some_path/a.fasta\n"
"    phmm-tree -prc -uals -id 0.9 -prc_hit 15 ./some_path/a.fasta\n"
"    phmm-tree -prc -uals -id 0.9 -prc_hit 15 -acc ./some_path/a.fasta\n"
"    phmm-tree -prc -uals -id 0.9 -fitch -fm ./some_path/a.fasta\n"
"    phmm-tree -prc -uals -id 0.9 -kitsch -min ./some_path/a.fasta\n"
"    phmm-tree -prc -uals -id 0.9 -fitch -upgma ./some_path/a.fasta\n"
"    phmm-tree -hhsuite -uals -id 0.9 -nj ./some_path/a.fasta\n"
"\n"
"  PRC + aligned FASTA folder:\n"
"    phmm-tree -prc -als ./some_path/\n"
"    phmm-tree -prc -als -prc_hit 15 ./some_path/\n"
"    phmm-tree -prc -als -prc_hit 15 -acc ./some_path/\n"
"    phmm-tree -prc -als -kitsch -nj ./some_path/\n"
"    phmm-tree -prc -als -fitch -fm ./some_path/\n"
"    phmm-tree -hhsuite -als ./some_path/\n"
"\n"
"  PRC + alignments with HMMER HMMs:\n"
"    phmm-tree -prc -als_phmms ./some_hmms_path/ ./some_als_path/\n"
"    phmm-tree -prc -als_phmms -prc_hit 15 ./some_hmms_path/ ./some_als_path/\n"
"    phmm-tree -prc -als_phmms -prc_hit 15 -acc ./some_hmms_path/ ./some_als_path/\n"
"    phmm-tree -prc -als_phmms -fitch -upgma ./some_hmms_path/ ./some_als_path/\n"
"    phmm-tree -prc -als_phmms -kitsch -min ./some_hmms_path/ ./some_als_path/\n"
"\n"
"  HHsuite + alignments with HHMs:\n"
"    phmm-tree -hhsuite -als_phhms ./some_hhms_path/ ./some_als_path/\n"
"    phmm-tree -hhsuite -als_phhms -nj -upgma ./some_hhms_path/ ./some_als_path/\n"
"\n"
"  PRC + HMMER HMM folder:\n"
"    phmm-tree -prc -hmms ./some_path/\n"
"    phmm-tree -prc -hmms -prc_hit 15 ./some_path/\n"
"    phmm-tree -prc -hmms -prc_hit 15 -acc ./some_path/\n"
"    phmm-tree -prc -hmms -fitch ./some_path/\n"
"    phmm-tree -prc -hmms -kitsch -fm ./some_path/\n"
"\n"
"  HHsuite + HHM folder:\n"
"    phmm-tree -hhsuite -hhms ./some_path/\n";
/* public_functions */
// Test whether the input string is numeric
int is_num_str(char *char_str_num);
// Split a string by any character in the delimiter list
std::vector<std::string> str_Split_by_char_list(std::string str, const char *pattern);
// Print the string if non-empty
void output_str(std::string  str);
// Convert unsigned int to string
std::string uint2str(unsigned int num);
// Check whether the directory exists
bool dir_exist_opendir(std::string path);
// Check whether the directory is non-empty
bool dir_noempty_opendir_readir(std::string path);
// Check that a file exists and is not empty
bool file_exists_and_empty_check(std::string str_path_name);
// Convert int to string
std::string int_2_string(int int_number);
// Convert double to string (6 decimal places)
std::string double_2_string(double double_x);
// List file names in a directory (optional extension filter)
int get_file_names(std::string path, std::vector<std::string> & vec_file_names, std::string extention);
// Delete all files in a directory
int delete_files(std::string path);
// Move files from path1 to path2
int mv_files(std::string path1, std::string path2, std::string extention);
// Copy one file from pathname1 to pathname2
int copy_one_file(std::string pathname1, std::string pathname2);
// Copy files from path1 to path2
int copy_files(std::string path1, std::string path2,std::string extention);
// Copy alignment-related files from path1 to path2
int copy_files_als_phmms_phhms_aligned(std::string path1, std::string path2,std::string extention);
// Handle return value from system()
int system_return(int status);
// Format a time duration as a human-readable string
std::string format_time_duration(long total_milliseconds);
// Print error and exit
void output_error_(std::string error_msg);
// Check if USEARCH exists (PATH or current directory)
int USEARCH_exist();
// Check if MAFFT exists (PATH or current directory)
int MAFFT_exist();
// Check if HMMER hmmbuild exists (PATH or current directory)
int HMMER_hmmbuild_exist();
// Check if HMMER hmmconvert exists (PATH or current directory)
int HMMER_hmmconvert_exist();
// Check if PRC exists (PATH or current directory)
int PRC_exist();
// Check if PRC-X (HMMER3 backend) exists (PATH or current directory)
int PRCX_exist();

// Check if HH-suite hhmake exists (PATH or current directory)
int hhsuite_hhmake_exist();
// Check if HH-suite hhalign exists (PATH or current directory)
int hhsuite_hhalign_exist();
// Replace all occurrences of s2 with s3 in s1, respecting protected prefixes in un_shorted
void string_replace(std::string &s1,const std::string&s2,const std::string&s3,unsigned int index_i, std::vector<std::string> un_shorted);

typedef struct stru_hmm_node {
	double M[20];
	double I[20];
	double D[20];
	double M_MID[3];
	double I_MID[3];
	double D_MID[3];
	stru_hmm_node() {
		for (size_t i = 0; i < 20; i++)
		{
			M[i] = 0.0;
			I[i] = 0.0;
			D[i] = 0.0;
		}

		for (size_t i = 0; i <3; i++)
		{
			M_MID[i] = 0.0;
			I_MID[i] = 0.0;
			D_MID[i] = 0.0;
		}
	}
}STRU_HMM_NODE;

typedef struct stru_hmm {
	std::string hmm_name;
	unsigned int length;
	std::vector<STRU_HMM_NODE> hmm_body;
}STRU_HMM;


// Structure to store per-alignment-line details
typedef struct struct_Result_hmm1_hmm2_each_line
{
    std::string hmm1_name;
    std::string hmm2_name;
	int begin_hmm1;         // begin position in hmm1
	int end_hmm1;           // end position in hmm1
	int begin_hmm2;         // begin position in hmm2
	int end_hmm2;           // end position in hmm2
	int hit_no;             // number of matching positions between the two HMMs
	double co_emis;         // co-emission score reported by PRC
	double simple;          // simple score reported by PRC
	double reverse;         // reverse score reported by PRC
	double E_value;         // E-value reported by PRC

	// Initialize members
	struct_Result_hmm1_hmm2_each_line()
	{
        hmm1_name="";
        hmm2_name="";
		begin_hmm1 = 0;
		end_hmm1 = 0;
		begin_hmm2 = 0;
		end_hmm2 = 0;
		hit_no = 0;
		co_emis = 0.0;
		simple = 0.0;
		reverse = 0.0;
	}
}STUC_RHHEL_NOTE;

// Structure that associates an HMM name with an ID
struct name_id
{
	std::string name;               // HMM name
	unsigned int identity_number;   // HMM ID
	name_id()
	{
		name = "";
		identity_number = 0;
	}
};

typedef struct struct_Result_hmm1_hmm2
{
	name_id hmm1;              // name and ID of HMM1
	name_id hmm2;              // name and ID of HMM2
	double distance;           // distance between the two HMMs
	struct_Result_hmm1_hmm2()
	{
		distance = 0.0;
	}
}STUC_RHH_NOTE;

// Structure to hold a path and its filenames
typedef struct strcut_path_filenames
{
	unsigned int count_of_filenames;
	std::string path;
	std::vector<std::string> filenames;
	strcut_path_filenames()
	{
		count_of_filenames = 0;
		path = "";
	}
}STUC_PFNS;

// Structure to hold HMM filename, NAME, and ACC fields
typedef struct hmm_name_acc {
    std::string filename;
	std::string NAME;
	std::string ACC;
	bool bool_flag_ACC;   // true if ACC is present/used
	hmm_name_acc(){
		NAME = "";
		ACC = "";
		filename="";
		bool_flag_ACC = false;
	}
}HMM_NAME_ACC;

// Structure to hold parsed command-line parameters
typedef struct cmd_params {
	bool prc_mode;
	bool prc_hmms;
	bool prc_prc_hit;
	bool prc_acc;
	bool prc_pair;
	// PRC backend selection: 0=auto, 1=legacy (convert+prc), 2=prcX direct
	int prc_backend;

	bool hhsuite_mode;
	bool hhsuite_hhms;

    bool uals;
    bool als;
	bool uals_id;
	bool als_phmms;
	bool als_phhms;

    // Phylogenetic analysis method selection
    bool run_fitch;
    bool run_kitsch;
    bool run_upgma;
    bool run_nj;

    // Fitch/Kitsch variant selection
    bool run_fm_only;       // run only f-m variants
    bool run_min_only;      // run only minimum evolution variants

    unsigned int prc_prc_hit_value;
    double uals_id_value;

    // Thread control parameters
    int prc_threads;        // number of threads for PRC distance calculations
	int phylo_threads;      // number of threads for phylogenetic algorithms (Fitch/Kitsch/NJ/UPGMA)
	int phylo_concurrent_threads; // number of concurrent phylogenetic analyses

	cmd_params(){
		prc_mode = false;
        prc_hmms = false;
        prc_prc_hit = false;
        prc_acc = false;
        prc_pair = false;
		prc_backend = 0; // auto by default

        hhsuite_mode = false;
        hhsuite_hhms = false;

        uals = false;
        als = false;
        uals_id = false;
        als_phmms = false;
        als_phhms = false;

        // By default, run all analyses if none are specified
        run_fitch = false;
        run_kitsch = false;
        run_upgma = false;
        run_nj = false;

        // By default, run both f-m and min variants
        run_fm_only = false;
        run_min_only = false;

        prc_prc_hit_value = 0;
        uals_id_value = 0.0;

	    // Default threads: PRC auto-detect (0); phylo default 1 (0 still means auto-detect)
        prc_threads = 0;
	    phylo_threads = 1;
		phylo_concurrent_threads = 0; // auto-detect concurrent processes
	}
}CMD_PARAMS;

// Main class
class HMMTree{
public:
	// Backend options for PRC processing
	enum PRCBackend { PRC_BACKEND_AUTO = 0, PRC_BACKEND_LEGACY = 1, PRC_BACKEND_PRCX = 2 };
	HMMTree() {
        hmm_a_mem_num=0;
        hmm_b_mem_num=0;
        bool_userach_in_folder = false;
        bool_MAFFT_in_folder = false;
        bool_HMMER_hmmbuild_in_folder = false;
        bool_PRC_in_folder = false;
		bool_PRCX_in_folder = false;
        bool_HMMER_hmmconvert_in_folder = false;
        NAME_OR_ACC = false;
        use_ACC = false;
        pairwise_mode = true;
        bool_hhsuite_hhmake_in_folder=false;
        bool_hhsuite_hhalign_in_folder=false;
        prc_hit_no = 10;
        double_user_id = 0.6;
        
        // Initialize analysis selection flags - default to all analyses
        run_fitch_analysis = true;
        run_kitsch_analysis = true;
        run_upgma_analysis = true;
        run_nj_analysis = true;
        
        // Initialize variant selection flags - default to both variants
        run_fm_only = false;
        run_min_only = false;
        
	    // Initialize thread control - PRC auto-detect (0), phylo default to 1
        prc_threads_count = 0;
	    phylo_threads_count = 1;
        phylo_concurrent_threads_count = 0; // auto-detect concurrent processes
        
        files_folder="";
        folder_hmms ="";
        folder_tree_files ="";
        folder_matrices ="";
        folder_clusters ="";
        folder_unalign_seqs ="";
        folder_invalid_clusters ="";
        folder_aligned ="";
        folder_hmmer2="";
        folder_hmmer3="";
        folder_hhms="";
        folder_prcfiles = "";
        folder_hmms_from_als="";
        folder_hhms_from_als="";
        prc_hhsuite = 0;
	prc_backend = PRC_BACKEND_AUTO;
	effective_prc_backend = PRC_BACKEND_LEGACY; // will be set during processing
	};
	~HMMTree() {};
    double CLUSTER_OK_THRESHOLD;                // threshold to evaluate clustering quality
    bool NAME_OR_ACC;                           // whether to use NAME or ACC in matrices when ACC is available
    bool use_ACC;                               // use the ACC key in PHYLIP/MEGA matrices if available
    bool pairwise_mode;                         // PRC pairwise mode (true) vs. library mode (false)
    int prc_hit_no;                             // PRC hit parameter (default 10; 0 uses PRC default 100)
    
    // Analysis selection flags
    bool run_fitch_analysis;
    bool run_kitsch_analysis;
    bool run_upgma_analysis;
    bool run_nj_analysis;
    
    // Variant selection flags for Fitch/Kitsch
    bool run_fm_only;        // run only f-m variants
    bool run_min_only;       // run only minimum evolution variants
    
    // Thread control parameters
    int prc_threads_count;    // number of threads for PRC analysis
    int phylo_threads_count;  // number of threads for phylogenetic analysis
	int phylo_concurrent_threads_count; // number of concurrent phylogenetic analyses (0 = auto)
	// PRC backend setting and effective backend selected for this run
	int prc_backend;                // user selection: auto/legacy/prcx
	int effective_prc_backend;      // resolved backend based on inputs and availability
    
    std::string files_folder;               // base output folder derived from the input path or filename
    std::string folder_hmms;
    std::string folder_tree_files;
    std::string folder_matrices;
    std::string folder_clusters;
    std::string folder_unalign_seqs;
    std::string folder_invalid_clusters;
    std::string folder_aligned;
    std::string folder_hmmer2;
    std::string folder_hmmer3;
    std::string folder_hhms;
    std::string folder_prcfiles;
    std::string folder_hmms_from_als;
    std::string folder_hhms_from_als;
    double double_user_id;
    int prc_hhsuite;
    static bool als_phmms_phhms;


	// Create an output folder named after the input file or path
	void create_files_folder(std::string input_file_or_folder_path, int prc_hhsuite, int int_run_mode);
	// Cluster sequences with USEARCH
    int usearch_cluster(std::string file_path_name, double identity);

	// Align all sequences in the current folder using MAFFT
	int align_do_mafft_all_from_file();

	// Build HMMs from all FASTA files in the folder using hmmbuild
	int hmm_do_hmmbuild_all_from_file();

	// Compute the distance between two HMMs
	int prc_each2();

	// Compute the distance between one HMM and a library of HMMs
	int prc_hmm1_library(std::string hmm1_name, std::string library, std::string output_path_filename);

	// Compute the distance between all HMMs and a library
	int prc_library();

	// Populate the matrix with distances between HMMs and a library
	int matrix_get_lib_hmms_result_2();

	// Populate the matrix with distances between each pair of HMMs
	int matrix_get_each2_hmms_result_2();

	// Populate the matrix with distances between each pair of HH-suite HMMs
	int matrix_get_each2_hhms_result_hhsuite();

	// Initialize the matrix, vector, and map for HMMs
	int matrix_init_matrix_vector_map();

	// Initialize the matrix, vector, and map for HH-suite HMMs
	int matrix_init_matrix_vector_map_hhsuite();

	// Process a FASTA sequence input for PRC
	void process_fasta_sequences(std::string file_path_name, double identity);

	// Process a FASTA sequence input for HH-suite
	void process_hhsuite_fasta_sequences(std::string file_path_name, double identity);


    // Process an alignments input with PRC
	void process_prc_alignments(std::string file_path_name);

	// Process alignments and profile HMMs (HMMER) with PRC
	void process_prc_alignments_phmms(std::string file_path_name, std::string file_path_name_2);

    // Process an alignments input with HH-suite
	void process_hhsuite_alignments(std::string file_path_name);

	// Process alignments and profile HHMs with HH-suite
	void process_hhsuite_alignments_phhms(std::string file_path_name, std::string file_path_name_2);

	// HMM 1000 occurrence test (?)
	int hmm_1000_A_B_occurence_test(std::string infile_path_and_name,int x);

	// Process a folder of HMMs (HMMER) with PRC
	void process_prc_HMMs(std::string infile_path_and_name);

	// Process a folder of HHMs (HH-suite) with HH-suite
	void process_hhsuite_HHMs(std::string infile_path_and_name);

	void draw_tree_test(){
	   draw_tree_selective(run_fitch_analysis, run_kitsch_analysis, run_upgma_analysis, run_nj_analysis, run_fm_only, run_min_only);
	}

	// Run selected phylogenetic analyses
	void draw_tree_selective(bool run_fitch, bool run_kitsch, bool run_upgma, bool run_nj, bool fm_only, bool min_only);
	// Backward-compatible entry point used elsewhere; delegates to draw_tree_selective
	void draw_tree_selective_concurrent(bool run_fitch, bool run_kitsch, bool run_upgma, bool run_nj, bool fm_only, bool min_only);
    // Replace shortened names in tree files with original names
	void trees_replace_shorted_names(std::string str_treefiles_path);



private:

    // USEARCH present in working directory?
    bool bool_userach_in_folder;
    // MAFFT present in working directory?
    bool bool_MAFFT_in_folder;
    // HMMER (hmmbuild) present in working directory?
    bool bool_HMMER_hmmbuild_in_folder;

    // HMMER (hmmconvert) present in working directory?
    bool bool_HMMER_hmmconvert_in_folder;
    // HH-suite (hhmake) present in working directory?
    bool bool_hhsuite_hhmake_in_folder;
    // HH-suite (hhalign) present in working directory?
    bool bool_hhsuite_hhalign_in_folder;
    // PRC present in working directory?
    bool bool_PRC_in_folder;
	// PRC-X present in working directory?
	bool bool_PRCX_in_folder;



    // Count of HMM-A members
    unsigned int hmm_a_mem_num;
    // Count of HMM-B members
    unsigned int hmm_b_mem_num;

    // Verify required programs exist
    void test_depend_programs(int mode_num);
    void test_depend_programs_2(int prc_hhsuite, int mode_num);
    // Cluster unaligned sequences with USEARCH
	int usearch_cluster_unalign_seqs(std::string file_path_name, double identity);
	// Rename USEARCH cluster result files
	int usearch_rename_result();
	// Ensure required directories exist; create if missing
	void dir_exist_or_create(int prc_hhsuite, int mode_num);

	// Verify FASTA file exists and format is valid
	bool align_fasta_file_exist_format_check(std::string infile_path_and_name);

	// Align a family with MAFFT L-INS-i approach
	int align_mafft_LINS_i_one(std::string family_name);

	// Move alignments to the aligned folder
	int align_move_aligned(std::string file_path_name);
	// Move alignments (for als_phmms/phhms modes) to the aligned folder
	int als_phmms_phhms_move_aligned(std::string file_path_name);

	// Verify HMM file exists and format is valid
	bool hmm_file_exist_format_check(std::string infile_path_and_name);

	// Split a multi-HMM file into single HMM files (each ends with "//")
	int hmm_divide_hmms_to_single_hmm(std::string infile_path_and_name, std::string outfiles_path, std::string str_key);

	// Build one HMM with hmmbuild
	int hmm_hmmbuild_one(std::string family_name);
	// Populate Pfam clan vector and HMM data vector/map
	int hmm_set_pfamclan_vector_hmmdat_vector_map();
	// Convert one HMMER3 HMM to HMMER2
	int hmm_hmmconvert_3_to_2_one(std::string hmmer3_name);
	// Convert all HMMER3 HMMs to HMMER2
	int hmm_hmmconvert_3_to_2();
	// HMMER hmmconvert dispatcher
	int hmm_hmmconvert();
    // Copy HMM files from user directory to hmms
    int hmm_copy_hmmfiles(std::string path);

    // Copy HMM files (als_phmms mode) from user directory to hmms
    int hmm_als_phmms_copy_hmmfiles(std::string path);

    // Build one HHM with hhmake
    int hhsuite_hhmake_one(std::string family_name);
    // Build all HHMs with hhmake
    int hhsuite_hhmake_all();
    // Pairwise hhalign over all HHMs
    int hhsuite_hhalign_each2();
    // Read hhalign results from files
    void hhsuite_read_result_from_file_hhalign(FILE * file_stream1, FILE * file_stream2);

    // Copy HHM files from user directory to hhms
    int hhsuite_copy_hhmfiles(std::string path);
    // Copy HHM files (als_phhms mode) from user directory to hhms
    int hhsuite_als_phhms_copy_hhmfiles(std::string path);

	// Clear PRC results
	void clear_Prc_Result()
	{
		res_msg.clear();
		for (int head_msg_i = 0; head_msg_i < 8; head_msg_i++)
			head_msg[head_msg_i] = "";
		hmm1 = "";
		hmm2 = "";
		length_hmm1 = 0;
		length_hmm2 = 0;
		//identity_number_hmm1 = -1;
		//identity_number_hmm2 = -1;
	}


    // Check whether profile HMM files are in HMMER2.0 format
    int prc_check_profile_HMM_format();

	// Read PRC pairwise result from file
	void prc_read_result_from_file(FILE * file_stream);
	// Thread-safe reader for PRC result
	void prc_read_result_from_file_threadsafe(FILE * file_stream, std::vector<STUC_RHHEL_NOTE>& local_res_msg);
	// Thread-safe distance computation from PRC lines
	double prc_hmms_dist_compute_threadsafe(const std::vector<STUC_RHHEL_NOTE>& local_res_msg);
	// Thread-safe setter for distance result
	void prc_set_STUC_RHH_NOTE_dist_threadsafe(STUC_RHH_NOTE* STUC_RHH_NOTE_res, const std::vector<STUC_RHHEL_NOTE>& local_res_msg, const std::string& local_hmm1, const std::string& local_hmm2);
	// Thread-safe aggregator into the pairwise matrix
	int matrix_get_each2_hmms_result_2_threadsafe(const std::vector<STUC_RHHEL_NOTE>& local_res_msg, const std::string& local_hmm1, const std::string& local_hmm2);
	
	// Helper to normalize PRC output identifiers to filenames
	std::string normalize_prc_identifier_to_filename(const std::string& prc_identifier);

	// Read PRC library-mode result file
	void prc_read_lib_result_from_file(std::string str_prc_lib_result_filename);

	// Compute distance from PRC lines
	double prc_hmms_dist_compute();

	// Set hmm1/hmm2 and distance into result struct
	void prc_set_STUC_RHH_NOTE_dist(STUC_RHH_NOTE* STUC_RHH_NOTE_res);

	// Output the distance matrix to a file in MEGA format
	int matrix_mega_output_dist_matrix_to_file();

	// Output the distance matrix to a file in PHYLIP format
	int matrix_phylip_output_dist_matrix_to_file();

	// Print the distance matrix to stdout
	int matrix_output_dist_matrix_to_window();

	// Validate matrix format
	bool matrix_check_matrix_format();

	// Compute average distance
	double matrix_compute_average_dist();

	// Clear matrix, maps, and vectors
	void all_clear_matrix_maps_vecotr(){
		//id_hmmname_hmmacc_list_vector.clear();
		hmm_NAME_id_list_unordered_map.clear();
		for (unsigned int i_matrix = 0; i_matrix < dist_matrix.size(); i_matrix++){
			dist_matrix[i_matrix].clear();
		}
		dist_matrix.clear();
	}

    // Output an HMM (debug)
    void output_hmm_test(STRU_HMM hmm);


    // Parse HMM data from file
    int hmm_get_data(std::string str_path_name, STRU_HMM *hmm) ;

    // Compute HMM node distance via KLD
    double hmm_compu_note(STRU_HMM_NODE note1, STRU_HMM_NODE note2);

    // Compute 2-HMM distance by longest subsequence
    double hmm_compu_dist(STRU_HMM hmm1, STRU_HMM hmm2, double EQUAL_LEVEL) ;

    // Compute pairwise distance between two HMM files
    double hmm_each2_dist(std::string str_path_name1 , std::string str_path_name2, double  EQUAL_LEVEL);

    // PHYLIP: draw tree
    int Phylip_draw_tree();

    // PHYLIP: draw tree (variant)
    int Phylip_draw_tree2();

	std::string PRC_CopyRight_license;
	std::string head_msg[8];           // header fields of PRC result
	std::string hmm1;                  // name of HMM1
	std::string hmm2;                  // name of HMM2
									   //int identity_number_hmm1;          // ID of HMM1
									   //int identity_number_hmm2;          // ID of HMM2
	int length_hmm1;                   // length of HMM1
	int length_hmm2;                   // length of HMM2

	std::vector<STUC_RHHEL_NOTE> res_msg;        // parsed PRC result lines

	//path_filename path_and_filenames;            // path and filename utilities

	//std::vector<HMM_NAME_ID_AC> vector_hmms_AB;  // selected HMM groups A and B

	//std::vector<PFAM_C> vector_pfam_clans;       // Pfam clan metadata
	//std::vector<HMM_DAT> vector_hmm_dat;         // Pfam hmm.dat metadata
	std::unordered_map <std::string, unsigned int> hmm_ACC_id_list_unordered_map;
	//std::vector<unsigned int> vector_index_right_accuracy; // accuracy indices
	//std::vector<PFAM_CLAN_AVG_DIS> vector_result_all;      // average distances per clan
    //std::vector<PFAM_CLAN> vector_rand_2clan_index;        // indices of random clan pairs
	std::vector<HMM_NAME_ACC> id_hmm_NAME_ACC_list_vector;  // list of HMM filename/NAME/ACC
    //std::vector<std::string> id_hmm_list_vector;           // list of HMM IDs
	std::unordered_map <std::string, unsigned int> hmm_NAME_id_list_unordered_map; // map NAME -> ID

	std::vector<std::vector<double> > dist_matrix;          // 2D array storing the final distance matrix

    // Names used in PHYLIP tree drawing before shortening
    std::vector <std::string> vec_unshorted_names;
    // Names used in PHYLIP tree drawing after shortening
    std::vector <std::string> vec_shorted_names;

	std::unordered_map <std::string, unsigned int> hmm_NAME_unordered_map;

	std::unordered_map <std::string, unsigned int> hmm_ACC_unordered_map;

};
