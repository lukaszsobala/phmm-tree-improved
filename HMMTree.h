#pragma once

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <stdio.h> /*fprintf(),stderr,BUFSIZ*/
#include <fcntl.h> /*open(),flag*/

//#include <hash_map>
#include <unordered_map>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <sstream>
#include <time.h>
#include <sys/timeb.h>
#include <sys/stat.h>
#include <float.h>

//kitsch.h
extern "C"{
#include "kitsch.h"
#include "fitch.h"
#include "neighbor.h"
#include "upgma.h"
}

const unsigned int MAX_BUF_LEN = 1024;             //program exists test
const unsigned int BUFFER_SIZE = 1024;             //file copy
//const double CLUSTER_OK_THRESHOLD = 0.90;                //a threshold to evaluate the cluster result
const double DBL_MAX_USER = 0.001;				//define the null score of the result
const unsigned int PTH_MATRIX_OUT = 0;				//set the last matrix out put path,used in function path_get of class path_filename ,the parametor of index_to_set
const unsigned int PTH_UNALG_SEQ = 1;				//set the unaligned sequences file path,used in function path_get of class path_filename ,the parametor of index_to_set
const unsigned int PTH_F_ALG_CLU = 2;				//set the clusters and aligenments path,used in function path_get of class path_filename ,the parametor of index_to_set
const unsigned int PTH_F_HMM = 3;					//set the hmm files path,used in function path_get of class path_filename ,the parametor of index_to_set
const unsigned int PTH_F_EACH_2_RESULT = 4;         //set the each two hmm compare result files path,used in function path_get of class path_filename ,the parametor of index_to_set
const unsigned int PTH_F_1_LIARBRY_RESULT = 5;      //set the one to a liarbry of hmms files path,used in function path_get of class path_filename ,the parametor of index_to_set

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
"               - -prc_hit <value> sets PRC parameter (default: 10; 0 uses PRC default 100).\n"
"               - -acc uses the 'ACC' key in matrix/tree if present (default is 'NAME').\n"
"               - -lib runs PRC in library mode (default is pairwise mode).\n"
"\n"
"  -als        Directory containing at least three aligned FASTA files.\n"
"               - -prc_hit and -acc are supported (same as in -uals).\n"
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
"  -prc_threads <num>     Number of threads for PRC distance calculations.\n"
"  -phylo_threads <num>   Number of threads for phylogenetic algorithms.\n"
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
//test the input string is num str or not
int is_num_str(char *char_str_num);
//split the string by char list
std::vector<std::string> str_Split_by_char_list(std::string str, const char *pattern);
//out put the string
void out_put_str(std::string  str);
//turn unsigned int to string
std::string uint2str(unsigned int num);
//function to test the dir exists or not
bool dir_exist_opendir(std::string path);
//function to test the folder is empty or not
bool dir_noempty_opendir_readir(std::string path);
//file exists and empty check
bool file_exists_and_empty_check(std::string str_path_name);
//convert int to string
std::string int_2_string(int int_number);
//convert double to string
std::string double_2_string(double double_x);
//get the file name in the folder
int get_file_names(std::string path, std::vector<std::string> & vec_file_names, std::string extention);
//delete the files in a path
int delete_files(std::string path);
//mv the files from path1 to path2
int mv_files(std::string path1, std::string path2, std::string extention);
//copy one file from pathname1 to pathname2
int copy_one_file(std::string pathname1, std::string pathname2);
//copy files from path1 to path2
int copy_files(std::string path1, std::string path2,std::string extention);
//copy files from path1 to path2
int copy_files_als_phmms_phhms_aligned(std::string path1, std::string path2,std::string extention);
//system return value deal
int system_return(int status);
//format time duration for human-readable output
std::string format_time_duration(long total_milliseconds);
//error deal
void output_error_(std::string error_msg);
//test the exists of USEARCH
int USEARCH_exist();
//test the exists of MAFFT
int MAFFT_exist();
//test the exists of HMMER hmmbuild
int HMMER_hmmbuild_exist();
//test the exists of HMMER hmmconvert
int HMMER_hmmconvert_exist();
//test the exists of PRC
int PRC_exist();

//test the exists of HMMER hmmbuild
int hhsuite_hhmake_exist();
//test the exists of HMMER hmmbuild
int hhsuite_hhalign_exist();
 //replace the string by another string in a string
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


//struct to save the the detail message of each line
typedef struct struct_Result_hmm1_hmm2_each_line
{
    std::string hmm1_name;
    std::string hmm2_name;
	int begin_hmm1;         //the bigin pot int the hmm1
	int end_hmm1;			//the end pot int the hmm2
	int begin_hmm2;			//the bigin pot int the hmm2
	int end_hmm2;			//the end pot int the hmm2
	int hit_no;				//the hit number of two hmms ,which means that there are hit_no of similar pot
	double co_emis;         //waitting for the comment
	double simple;		    //waitting for the comment
	double reverse;         //waitting for the comment
	double E_value;         //waitting for the comment

							//init the parameters
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

//struct to put the hmm name and the id number together
struct name_id
{
	std::string name;               //hmm name
	unsigned int identity_number;   // hmm id number
	name_id()
	{
		name = "";
		identity_number = 0;
	}
};

typedef struct struct_Result_hmm1_hmm2
{
	name_id hmm1;              //the first name and id of the hmm
	name_id hmm2;              //the second name and id of the hmm
	double distance;           //the diatance of two hmms
	struct_Result_hmm1_hmm2()
	{
		distance = 0.0;
	}
}STUC_RHH_NOTE;

//struct to save the path and the filenames
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

//struct to get the HMM name and acc structer
typedef struct hmm_name_acc {
    std::string filename;
	std::string NAME;
	std::string ACC;
	bool bool_flag_ACC;
	hmm_name_acc(){
		NAME = "";
		ACC = "";
		filename="";
		bool_flag_ACC = false;
	}
}HMM_NAME_ACC;

//struct to get the HMM name and acc structer
typedef struct cmd_params {
	bool prc_mode;
	bool prc_hmms;
	bool prc_prc_hit;
	bool prc_acc;
	bool prc_pair;

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
    bool run_fm_only;       // Run only f-m variants
    bool run_min_only;      // Run only minimum evolution variants

    unsigned int prc_prc_hit_value;
    double uals_id_value;

    // Thread control parameters
    int prc_threads;        // Number of threads for PRC distance calculations
    int phylo_threads;      // Number of threads for phylogenetic algorithms (Fitch/Kitsch/NJ/UPGMA)

	cmd_params(){
		prc_mode = false;
        prc_hmms = false;
        prc_prc_hit = false;
        prc_acc = false;
        prc_pair = false;

        hhsuite_mode = false;
        hhsuite_hhms = false;

        uals = false;
        als = false;
        uals_id = false;
        als_phmms = false;
        als_phhms = false;

        // By default, run all analyses if none specified
        run_fitch = false;
        run_kitsch = false;
        run_upgma = false;
        run_nj = false;

        // By default, run both f-m and min variants
        run_fm_only = false;
        run_min_only = false;

        prc_prc_hit_value = 0;
        uals_id_value = 0.0;

        // Default: auto-detect threads (0 means auto-detect)
        prc_threads = 0;
        phylo_threads = 0;

	}
}CMD_PARAMS;

//the main class
class HMMTree{
public:
	HMMTree() {
        hmm_a_mem_num=0;
        hmm_b_mem_num=0;
        bool_userach_in_folder = false;
        bool_MAFFT_in_folder = false;
        bool_HMMER_hmmbuild_in_folder = false;
        bool_PRC_in_folder = false;
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
        
        // Initialize thread control - default to auto-detect (0)
        prc_threads_count = 0;
        phylo_threads_count = 0;
        
        files_folder="";
        folder_hmms ="";
        folder_tree_files ="";
        folder_matrixs ="";
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
	};
	~HMMTree() {};
    double CLUSTER_OK_THRESHOLD;                //a threshold to evaluate the cluster result
    bool NAME_OR_ACC;                   //use 'NAME' or 'ACC' in the matrix when the hmm files have 'ACC' key
    bool use_ACC;                         //use the ACC key in the phylip matrix or not
    bool pairwise_mode;                    //use pairwise mode or not in PRC process
    int prc_hit_no;                         //use the ACC key in the phylip matrix or not
    
    // Analysis selection flags
    bool run_fitch_analysis;
    bool run_kitsch_analysis;
    bool run_upgma_analysis;
    bool run_nj_analysis;
    
    // Variant selection flags for Fitch/Kitsch
    bool run_fm_only;        // Run only f-m variants
    bool run_min_only;       // Run only minimum evolution variants
    
    // Thread control parameters
    int prc_threads_count;    // Number of threads for PRC analysis
    int phylo_threads_count;  // Number of threads for phylogenetic analysis
    
    std::string files_folder;               // the folder to save the files created by the input path or file name
    std::string folder_hmms;
    std::string folder_tree_files;
    std::string folder_matrixs;
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


	//create fils folder named after the input file or path
	void create_files_folder(std::string input_file_or_folder_path, int prc_hhsuite, int int_run_mode);
	//use usearch do clusters
    int usearch_cluster(std::string file_path_name, double identity);

	//Align all sequences in the current folder using MAFFT
	int align_do_mafft_all_from_file();

	//Build HMMs from all FASTA files in the folder using hmmbuild
	int hmm_do_hmmbuild_all_from_file();

	//Compute the distance between two HMMs
	int prc_each2();

	//Compute the distance between one HMM and a library of HMMs
	int prc_hmm1_library(std::string hmm1_name, std::string library, std::string output_path_filename);

	//Compute the distance between all HMMs and a library
	int prc_library();

	//Fill the matrix with distances between HMMs and a library
	int matrix_get_lib_hmms_result_2();

	//Fill the matrix with distances between each pair of HMMs
	int matrix_get_each2_hmms_result_2();

	//Fill the matrix with distances between each pair of HHsuite HMMs
	int matrix_get_each2_hhms_result_hhsuite();

	//Initialize the matrix, vector, and map for HMMs
	int matrix_init_matrix_vector_map();

	//Initialize the matrix, vector, and map for HHsuite HMMs
	int matrix_init_matrix_vector_map_hhsuite();

	//Process a FASTA sequence input for PRC
	void prc_fasta_seqs_deal(std::string file_path_name, double identity);

	//Process a FASTA sequence input for HHsuite
	void hhsuite_fasta_seqs_deal(std::string file_path_name, double identity);


    //function deal alignments input
	void prc_alignments_deal(std::string file_path_name);

	//function deal alignments and phmms input
	void prc_alignments_phmms_deal(std::string file_path_name, std::string file_path_name_2);

    //function deal alignments input
	void hhsuite_alignments_deal(std::string file_path_name);

	//function deal alignments and phhms input
	void hhsuite_alignments_phhms_deal(std::string file_path_name, std::string file_path_name_2);

	//hmm 1000 accurence test
	int hmm_1000_A_B_accurence_test(std::string infile_path_and_name,int x);

	//function deal HMMs input
	void prc_HMMs_deal(std::string infile_path_and_name);

	//function deal HMMs input
	void hhsuite_HHMs_deal(std::string infile_path_and_name);

	void draw_tree_test(){
	   draw_tree_selective(run_fitch_analysis, run_kitsch_analysis, run_upgma_analysis, run_nj_analysis, run_fm_only, run_min_only);
	}

	//function to run selected phylogenetic analyses
	void draw_tree_selective(bool run_fitch, bool run_kitsch, bool run_upgma, bool run_nj, bool fm_only, bool min_only);
    //function to replace the shortednames to the before-replaced names
	void trees_replace_shorted_names(std::string str_treefiles_path);



private:

    //usearch is in the first folder or not?
    bool bool_userach_in_folder;
    //MAFFT is in the first folder or not?
    bool bool_MAFFT_in_folder;
    //HMMER is in the first folder or not?
    bool bool_HMMER_hmmbuild_in_folder;

    //HMMER is in the first folder or not?
    bool bool_HMMER_hmmconvert_in_folder;
    //
    bool bool_hhsuite_hhmake_in_folder;
    //
    bool bool_hhsuite_hhalign_in_folder;
    //PRC is in the first folder or not?
    bool bool_PRC_in_folder;



    //count the hmmA hmm member' number
    unsigned int hmm_a_mem_num;
    //count the hmmB hmm member' number
    unsigned int hmm_b_mem_num;

    //test the depend programs exist or not
    void test_depend_programs(int mode_num);
    void test_depend_programs_2(int prc_hhsuite, int mode_num);
    //use usearch do clusters
	int usearch_cluster_unalign_seqs(std::string file_path_name, double identity);
	//rename the clusteer reslut files
	int usearch_rename_result();
	//test if the dir exists, or create them
	void dir_exist_or_create(int prc_hhsuite, int mode_num);

	//fasta file format and exists check
	bool align_fasta_file_exist_format_check(std::string infile_path_and_name);

	//seqs align by mafft LINS_i approch
	int align_mafft_LINS_i_one(std::string family_name);

	//move the alignments to the aligned folder
	int align_move_aligned(std::string file_path_name);
	//move the alignments to the aligned folder
	int als_phmms_phhms_move_aligned(std::string file_path_name);

	//HMM file format and exists check
	bool hmm_file_exist_format_check(std::string infile_path_and_name);

	//divide hmms into single hmm files and the hmms in the file must be ended by "//"
	int hmm_divide_hmms_to_single_hmm(std::string infile_path_and_name, std::string outfiles_path, std::string str_key);

	//hmmbuild a hmm model
	int hmm_hmmbuild_one(std::string family_name);
	//set the pfam clan vector and the hmms data vector and map
	int hmm_set_pfamclan_vector_hmmdat_vector_map();
	//hmmconvert one hmmer3 hmm model to hmmer2 hmm model
	int hmm_hmmconvert_3_to_2_one(std::string hmmer3_name);
	//hmmconvert all hmmer3 hmm models to hmmer2 hmm models
	int hmm_hmmconvert_3_to_2();
	//hmmconvert
	int hmm_hmmconvert();
    //copy the hmm files from user dir to hmms
    int hmm_copy_hmmfiles(std::string path);

    //copy the hmm files from user dir to hmms
    int hmm_als_phmms_copy_hmmfiles(std::string path);

    //hhmake build one HMM
    int hhsuite_hhmake_one(std::string family_name);
    //hhmake build all HMM
    int hhsuite_hhmake_all();
    //hhalign
    int hhsuite_hhalign_each2();
    //
    void hhsuite_read_result_from_file_hhalign(FILE * file_stream1, FILE * file_stream2);

    int hhsuite_copy_hhmfiles(std::string path);
    //copy the hhm files from user dir to hmms
    int hhsuite_als_phhms_copy_hhmfiles(std::string path);

	//clear the message
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


    //check the profile HMM files is HMMER2.O or not
    int prc_check_profile_HMM_format();

	//read into the  message from the prc_each2 file
	void prc_read_result_from_file(FILE * file_stream);
	void prc_read_result_from_file_threadsafe(FILE * file_stream, std::vector<STUC_RHHEL_NOTE>& local_res_msg);
	double prc_hmms_dist_compute_threadsafe(const std::vector<STUC_RHHEL_NOTE>& local_res_msg);
	void prc_set_STUC_RHH_NOTE_dist_threadsafe(STUC_RHH_NOTE* STUC_RHH_NOTE_res, const std::vector<STUC_RHHEL_NOTE>& local_res_msg, const std::string& local_hmm1, const std::string& local_hmm2);
	int matrix_get_each2_hmms_result_2_threadsafe(const std::vector<STUC_RHHEL_NOTE>& local_res_msg, const std::string& local_hmm1, const std::string& local_hmm2);
	
	// Helper function to normalize PRC output identifiers to filename
	std::string normalize_prc_identifier_to_filename(const std::string& prc_identifier);

	void prc_read_lib_result_from_file(std::string str_prc_lib_result_filename);

	//compute the distence
	double prc_hmms_dist_compute();

	//set the hmm1 and hmm2 message to the result matrix class
	void prc_set_STUC_RHH_NOTE_dist(STUC_RHH_NOTE* STUC_RHH_NOTE_res);

	//out put the distance matrix to a file in mega format
	int matrix_mega_out_put_dist_matrix_to_file();

	//out put the distance matrix to a file in phylip format
	int matrix_phylip_out_put_dist_matrix_to_file();

	//out put the distance matrix to the window
	int matrix_out_put_dist_matrix_to_window();

	//check the matrix is the right format or not
	bool matrix_check_matrix_format();

	//double compute average distance score
	double matrix_compute_average_dist();

	//clear the matrix, maps and the vector
	void all_clear_matrix_maps_vecotr(){
		//id_hmmname_hmmacc_list_vector.clear();
		hmm_NAME_id_list_unordered_map.clear();
		for (unsigned int i_matrix = 0; i_matrix < dist_matrix.size(); i_matrix++){
			dist_matrix[i_matrix].clear();
		}
		dist_matrix.clear();
	}

    //out put the hmm data to check
    void output_hmm_test(STRU_HMM hmm);


    //get data from the hmm file
    int hmm_get_data(std::string str_path_name, STRU_HMM *hmm) ;

    //compute the hmm note dist by KLD
    double hmm_compu_note(STRU_HMM_NODE note1, STRU_HMM_NODE note2);

    //compute the two hmms' distance by the longest subseq
    double hmm_compu_dist(STRU_HMM hmm1, STRU_HMM hmm2, double EQUAL_LEVEL) ;

    //compute the each two dist
    double hmm_each2_dist(std::string str_path_name1 , std::string str_path_name2, double  EQUAL_LEVEL);

    //Phylip draw tree
    int Phylip_draw_tree();

    //Phylip draw tree
    int Phylip_draw_tree2();

	std::string PRC_CopyRight_license;
	std::string head_msg[8];           //a string array to save the head message of the prc result
	std::string hmm1;                  //the first name of the hmm
	std::string hmm2;                  //the second name of the hmm
									   //int identity_number_hmm1;          //the identity number of hmm1
									   //int identity_number_hmm2;          //id number of hmm2
	int length_hmm1;                   //the length of the first hmm
	int length_hmm2;                   //the length of the second hmm

	std::vector<STUC_RHHEL_NOTE> res_msg;        //the last messages of the result

	//Prc_Result deal_prc_results;			//deal the prc reasult
	//path_filename path_and_filenames;		// deal with the path and the filenames

	//std::vector<HMM_NAME_ID_AC> vector_hmms_AB;      //VECTOR to save the shosen hmmAs and hmmBs

	//std::vector<PFAM_C> vector_pfam_clans;		//vector to save the clans' messages read from the Pfam-c file
	//std::vector<HMM_DAT> vector_hmm_dat;		//a vector to save the hmms data in Pfam.hmm.dat file
	std::unordered_map <std::string, unsigned int> hmm_ACC_id_list_unordered_map;
	//std::vector<unsigned int> vector_index_right_accuracy;		//a vector to save the right accuracy clans' indexes of compute times
	//std::vector<PFAM_CLAN_AVG_DIS> vector_result_all;		//a vector to save the average distance score of clan A, B and A_B
    //std::vector<PFAM_CLAN> vector_rand_2clan_index;			//A vector to save the indexs of the two rand clans
	std::vector<HMM_NAME_ACC> id_hmm_NAME_ACC_list_vector;			//a vector to save the pairs of hmm id and name
    //std::vector<std::string> id_hmm_list_vector;					//a vector to save the pairs of hmm id and name
	std::unordered_map <std::string, unsigned int> hmm_NAME_id_list_unordered_map;	//a unordered map to save the pairs of hmm name and id

	std::vector<std::vector<double> > dist_matrix;					 //two dimencions array two save the last result

    //a unordered map to save the names in phylip draw tree step before names short step
    std::vector <std::string> vec_unshorted_names;
        //a unordered map to save the names in phylip draw tree step after names short step
    std::vector <std::string> vec_shorted_names;

	std::unordered_map <std::string, unsigned int> hmm_NAME_unordered_map;


	std::unordered_map <std::string, unsigned int> hmm_ACC_unordered_map;

};
