#include "HMMTree.h"

/*
*-prc:
*    -uals -id <0.xxx>
*    -als
*    -hmms
*
*         -prc_hit <int>
*         -acc
*-hhsuite
*    -uals  -id <0.xxx>
*    -als
*    -hhms
*/
bool HMMTree::als_phmms_phhms = false;
int main(int argc, char * argv[]){

    // If no arguments are provided, print help/usage and exit successfully
    if (argc == 1) {
        std::cout << STR_ARGUMENTS_ERROR_MSG << std::endl;
        return 0;
    }

    // Hidden worker mode: run a single phylogenetic algorithm in isolation
    if (argc > 1 && strcmp(argv[1], "-phylo_worker") == 0) {
        if (argc != 6) {
            std::cerr << "Usage: phmm-tree -phylo_worker <algo> <matrix_file> <output_file> <threads>" << std::endl;
            return 1;
        }
        const char* algo = argv[2];
        const char* matrix = argv[3];
        const char* output = argv[4];
        int threads = atoi(argv[5]);

        int rc = 1;
        if (strcmp(algo, "kitsch_fm") == 0) {
            rc = kitsch_build_tree(matrix, output, 0, threads);
        } else if (strcmp(algo, "kitsch_min") == 0) {
            rc = kitsch_build_tree(matrix, output, 1, threads);
        } else if (strcmp(algo, "fitch_fm") == 0) {
            rc = fitch_build_tree(matrix, output, 0, threads);
        } else if (strcmp(algo, "fitch_min") == 0) {
            rc = fitch_build_tree(matrix, output, 1, threads);
        } else if (strcmp(algo, "nj") == 0) {
            rc = neighbor_build_tree(matrix, output, threads);
        } else if (strcmp(algo, "upgma") == 0) {
            rc = upgma_build_tree(matrix, output, threads);
        } else {
            std::cerr << "Unknown phylogenetic algorithm: " << algo << std::endl;
            return 2;
        }
        return rc;
    }

    if(argc == 2){
        const char* str_argv = argv[1];
        if(strcmp(str_argv, "-h") == 0){
            std::cout<< STR_ARGUMENTS_ERROR_MSG <<std::endl;
        }else{
            output_error_("Arguments number");
            exit(1);
        }
        exit(0);
    }
	if (argc < 3)
	{
		output_error_("Arguments number");
		return 0;
	}

	static const char* kArgs[] = {
        "-prc",       // 0
        "-uals",      // 1
        "-als",       // 2
        "-hmms",      // 3
        "-id",        // 4
        "-prc_hit",   // 5
        "-acc",       // 6
        "-lib",       // 7
        "-hhsuite",   // 8
        "-hhms",      // 9
        "-als_phmms", // 10
        "-als_phhms", // 11
        "-fitch",     // 12
        "-kitsch",    // 13
        "-upgma",     // 14
        "-nj",        // 15
        "-fm",        // 16
        "-min",       // 17
        "-prc_threads",   // 18
    "-phylo_threads",  // 19
    "-phylo_concurrent_threads", // 20
    "-prc_backend" // 21
    };
    const size_t kArgsCount = sizeof(kArgs) / sizeof(kArgs[0]);
	//string to save the compute model from the user input

    CMD_PARAMS strc_cmd;

    std::string str_file_path =argv[argc - 1];

    int arg_num = 1;
    while(arg_num < argc - 1 && (!strc_cmd.als_phmms || arg_num < argc - 2) && (!strc_cmd.als_phhms || arg_num < argc - 2)){
        const char* str_argu_temp = argv[arg_num];
        bool find_argu = false;
        size_t arg_vec_num = 0;
        while(arg_vec_num < kArgsCount){
            if(strcmp(str_argu_temp, kArgs[arg_vec_num]) == 0){
                switch(static_cast<int>(arg_vec_num)){

                    case 0:
                        if(!strc_cmd.prc_mode && !strc_cmd.hhsuite_mode){
                            strc_cmd.prc_mode = true;
                            find_argu = true;
                        }else{
                            output_error_("Mode '-prc' arguments( [-prc, -hhsuite] )");
                        }

                        break;
                    case 1:
                        if(!strc_cmd.uals && !strc_cmd.als && !strc_cmd.als_phmms && !strc_cmd.als_phhms && !strc_cmd.hhsuite_hhms  && !strc_cmd.prc_hmms){
                            strc_cmd.uals = true;
                            find_argu = true;
                        }else{
                            output_error_("Run mode '-uals' arguments( [-uals, -als, -als_phmms(prc), -als_phhms(hhsuite), -hmms(prc), -hhms(hhsuite)] )");
                        }

                        break;
                    case 2:
                        if(!strc_cmd.uals && !strc_cmd.als && !strc_cmd.als_phmms && !strc_cmd.als_phhms && !strc_cmd.hhsuite_hhms  && !strc_cmd.prc_hmms){
                            strc_cmd.als = true;
                            find_argu = true;
                        }else{
                            output_error_("Run mode '-als' arguments( [-uals, -als, -als_phmms(prc), -als_phhms(hhsuite), -hmms(prc), -hhms(hhsuite)] )");
                        }

                        break;
                    case 3:
                        if(!strc_cmd.uals && !strc_cmd.als && !strc_cmd.als_phmms && !strc_cmd.als_phhms && !strc_cmd.hhsuite_hhms  && !strc_cmd.prc_hmms){
                            strc_cmd.prc_hmms = true;
                            find_argu = true;
                        }else{
                            output_error_("Run mode '-hmms' arguments( [-uals, -als, -als_phmms(prc), -als_phhms(hhsuite), -hmms(prc), -hhms(hhsuite)] )");
                        }

                        break;
                    case 4:
                        arg_num++;
                        if(!strc_cmd.uals_id){
                            if(arg_num >= argc - 1 || !is_num_str(argv[arg_num])){
                                output_error_("'-id' value arguments");
                            }
                            strc_cmd.uals_id_value = atof(argv[arg_num]);
                            strc_cmd.uals_id = true;
                            find_argu = true;
                        }else{
                            output_error_("'-id' arguments");
                        }

                        break;
                    case 5:
                        arg_num++;
                        if(!strc_cmd.prc_prc_hit){
                            if(arg_num >= argc - 1 || !is_num_str(argv[arg_num])){
                                output_error_("'-prc_hit' value arguments");
                            }
                            strc_cmd.prc_prc_hit_value = atoi(argv[arg_num]);
                            strc_cmd.prc_prc_hit = true;
                            find_argu = true;
                        }else{
                            output_error_("'-prc_hit' arguments");
                        }

                        break;
                    case 6:
                        if(!strc_cmd.prc_acc){
                            strc_cmd.prc_acc = true;
                            find_argu = true;
                        }else{
                            output_error_("'-acc' arguments");
                        }

                        break;
                    case 7:
                        if(!strc_cmd.prc_pair){
                            strc_cmd.prc_pair = true;
                            find_argu = true;
                        }else{
                            output_error_("'-pair' arguments");
                        }
                        break;
                    case 8:
                        if(!strc_cmd.prc_mode && !strc_cmd.hhsuite_mode){
                            strc_cmd.hhsuite_mode = true;
                            find_argu = true;
                        }else{
                            output_error_("Mode '-hhsuite' arguments( [-prc, -hhsuite] )");
                        }
                        break;
                    case 9:
                        if(!strc_cmd.uals && !strc_cmd.als && !strc_cmd.als_phmms && !strc_cmd.als_phhms && !strc_cmd.hhsuite_hhms  && !strc_cmd.prc_hmms){
                            strc_cmd.hhsuite_hhms = true;
                            find_argu = true;
                        }else{
                            output_error_("Run mode '-hhms' arguments( [-uals, -als, -als_phmms(prc), -als_phhms(hhsuite), -hmms(prc), -hhms(hhsuite)] )");
                        }
                        break;
                    case 10:
                        if(!strc_cmd.als_phmms && !strc_cmd.uals && !strc_cmd.als && !strc_cmd.als_phhms && !strc_cmd.hhsuite_hhms  && !strc_cmd.prc_hmms){
                            strc_cmd.als_phmms = true;
                            find_argu = true;
                        }else{
                            output_error_("Run mode '-als_phmms' arguments( [-uals, -als, -als_phmms(prc), -als_phhms(hhsuite), -hmms(prc), -hhms(hhsuite)] )");
                        }
                        break;

                    case 11:
                        if(!strc_cmd.als_phhms && !strc_cmd.als_phmms && !strc_cmd.uals && !strc_cmd.als && !strc_cmd.hhsuite_hhms  && !strc_cmd.prc_hmms){
                            strc_cmd.als_phhms = true;
                            find_argu = true;
                        }else{
                            output_error_("Run mode '-als_phhms' arguments( [-uals, -als, -als_phmms(prc), -als_phhms(hhsuite), -hmms(prc), -hhms(hhsuite)] )");
                        }
                        break;
                    case 12:
                        strc_cmd.run_fitch = true;
                        find_argu = true;
                        break;
                    case 13:
                        strc_cmd.run_kitsch = true;
                        find_argu = true;
                        break;
                    case 14:
                        strc_cmd.run_upgma = true;
                        find_argu = true;
                        break;
                    case 15:
                        strc_cmd.run_nj = true;
                        find_argu = true;
                        break;
                    case 16:
                        if (!strc_cmd.run_min_only) {
                            strc_cmd.run_fm_only = true;
                            find_argu = true;
                        } else {
                            output_error_("Cannot specify both '-fm' and '-min' arguments");
                        }
                        break;
                    case 17:
                        if (!strc_cmd.run_fm_only) {
                            strc_cmd.run_min_only = true;
                            find_argu = true;
                        } else {
                            output_error_("Cannot specify both '-fm' and '-min' arguments");
                        }
                        break;
                    case 18:
                        // -prc_threads
                        arg_num++;
                        if(arg_num >= argc - 1 || !is_num_str(argv[arg_num])){
                            output_error_("'-prc_threads' argument must be followed by a positive integer");
                        }
                        strc_cmd.prc_threads = atoi(argv[arg_num]);
                        if(strc_cmd.prc_threads < 0){
                            output_error_("'-prc_threads' must be a positive integer (0 for auto-detect)");
                        }
                        find_argu = true;
                        break;
                    case 19:
                        // -phylo_threads
                        arg_num++;
                        if(arg_num >= argc - 1 || !is_num_str(argv[arg_num])){
                            output_error_("'-phylo_threads' argument must be followed by a positive integer");
                        }
                        strc_cmd.phylo_threads = atoi(argv[arg_num]);
                        if(strc_cmd.phylo_threads < 0){
                            output_error_("'-phylo_threads' must be a positive integer (0 for auto-detect)");
                        }
                        find_argu = true;
                        break;
                    case 20:
                        // -phylo_concurrent_threads
                        arg_num++;
                        if(arg_num >= argc - 1 || !is_num_str(argv[arg_num])){
                            output_error_("'-phylo_concurrent_threads' argument must be followed by a positive integer");
                        }
                        strc_cmd.phylo_concurrent_threads = atoi(argv[arg_num]);
                        if(strc_cmd.phylo_concurrent_threads < 0){
                            output_error_("'-phylo_concurrent_threads' must be a positive integer (0 for auto-detect)");
                        }
                        find_argu = true;
                        break;
                    case 21:
                        // -prc_backend <auto|legacy|prcx>
                        arg_num++;
                        if(arg_num >= argc - 1){
                            output_error_("'-prc_backend' requires a value: auto|legacy|prcx");
                        }
                        {
                            std::string v = argv[arg_num];
                            if (v == "auto") strc_cmd.prc_backend = 0;
                            else if (v == "legacy") strc_cmd.prc_backend = 1;
                            else if (v == "prcx") strc_cmd.prc_backend = 2;
                            else output_error_("'-prc_backend' value must be one of: auto|legacy|prcx");
                        }
                        find_argu = true;
                        break;
                    default:
                        output_error_("Arguments");
                        break;
                }
                // Found and processed a matching argument; stop scanning further flags for this token.
                break;
            }
            arg_vec_num++;
        }
        if(!find_argu){
            output_error_("Arguments");
            return 0;
        }
        arg_num++;
    }


    //test to compute
	HMMTree  test;
    if(strc_cmd.prc_mode){
        if(strc_cmd.hhsuite_hhms || strc_cmd.als_phhms){
            output_error_("'-prc' mode cannot run in '-hhms' or '-als_phhms' mode, arguments");
        }
        if(!strc_cmd.uals && !strc_cmd.als && !strc_cmd.prc_hmms && !strc_cmd.als_phmms){
            output_error_("Please at least one run mode('-uals', '-als', '-als_phmms' or '-hmms' ) to run in '-prc' mode ,arguments");
        }

        if(strc_cmd.uals){
            if(strc_cmd.uals_id ){
                if(strc_cmd.uals_id_value < 0.1 || strc_cmd.uals_id_value > 1.0 ){
                    std::cout<<"Error:  The user set '-id' must be a number between 0.1 to 1.0 (0.1<=x && x <= 1.0) , current input is '"<<strc_cmd.uals_id_value<<"' ."<<std::endl;
                    output_error_("'-prc' mode '-uals' argument '-id'");
                }
                test.double_user_id = strc_cmd.uals_id_value;
            }
        }else if(strc_cmd.uals_id ){
            output_error_("'-prc' mode argument '-id'");
        }

        if(strc_cmd.als || strc_cmd.prc_hmms || strc_cmd.als_phmms){
            if(strc_cmd.uals_id){
                output_error_("'-prc' mode '-als' arguments or '-als_phmms' arguments cannot set '-id', arguments");
            }
        }
        if(strc_cmd.prc_pair){
            test.pairwise_mode = false;
        }

        if(strc_cmd.prc_acc){
            test.NAME_OR_ACC = true;
        }

        if(strc_cmd.prc_prc_hit){
            if(test.prc_hit_no < 0 ){
                std::cout<<"Error:  The user set '-prc_hit' must be a non negative value , current input is '"<<test.prc_hit_no<<"' ."<<std::endl;
                output_error_("'-prc' mode '-prc_hit' arguments");
            }
            test.prc_hit_no = strc_cmd.prc_prc_hit_value;
        }


    }

    if(strc_cmd.hhsuite_mode){

        if(strc_cmd.prc_acc || strc_cmd.prc_hmms || strc_cmd.prc_pair || strc_cmd.prc_prc_hit || strc_cmd.als_phmms){
            output_error_("'-hhsuite' mode cannot run in '-acc', '-hmms', '-pair' or '-prc_hit' mode, arguments");
        }

        if(!strc_cmd.uals && !strc_cmd.als && !strc_cmd.hhsuite_hhms && !strc_cmd.als_phhms){
            output_error_("Please at least one run mode('-uals', '-als', '-als_phhms' or '-hhms' ) to run in '-hhsuite' mode, arguments");
        }

        if(strc_cmd.uals){
            if(strc_cmd.uals_id ){
                if(strc_cmd.uals_id_value < 0.1 || strc_cmd.uals_id_value > 1.0 ){
                    std::cout<<"Error:  The user set '-id' must be a number between 0.1 to 1.0 (0.1<=x && x <= 1.0) , current input is '"<<strc_cmd.uals_id_value<<"' ."<<std::endl;
                    output_error_("'-hhsuite' mode '-uals' argument '-id'");
                }
                test.double_user_id = strc_cmd.uals_id_value;
            }
        }else if(strc_cmd.uals_id){
            output_error_("'-hhsuite' mode argument '-id'");
        }

        test.prc_hhsuite = 1;
    }

    // Configure phylogenetic analysis selection
    // If specific analyses were requested, only run those
    bool any_analysis_specified = strc_cmd.run_fitch || strc_cmd.run_kitsch || strc_cmd.run_upgma || strc_cmd.run_nj;
    if (any_analysis_specified) {
        test.run_fitch_analysis = strc_cmd.run_fitch;
        test.run_kitsch_analysis = strc_cmd.run_kitsch;
        test.run_upgma_analysis = strc_cmd.run_upgma;
        test.run_nj_analysis = strc_cmd.run_nj;
        
        std::cout << "Selected phylogenetic analyses: ";
        if (strc_cmd.run_fitch) std::cout << "Fitch ";
        if (strc_cmd.run_kitsch) std::cout << "Kitsch ";
        if (strc_cmd.run_upgma) std::cout << "UPGMA ";
        if (strc_cmd.run_nj) std::cout << "Neighbor-Joining ";
        std::cout << std::endl;
    } else {
        // Default: run all analyses
        std::cout << "Running all phylogenetic analyses (default)" << std::endl;
    }

    // Configure Fitch/Kitsch variant selection
    test.run_fm_only = strc_cmd.run_fm_only;
    test.run_min_only = strc_cmd.run_min_only;
    if (strc_cmd.run_fm_only) {
        std::cout << "Variant selection: f-m only" << std::endl;
    } else if (strc_cmd.run_min_only) {
        std::cout << "Variant selection: minimum evolution only" << std::endl;
    } else if (any_analysis_specified && (strc_cmd.run_fitch || strc_cmd.run_kitsch)) {
        std::cout << "Variant selection: both f-m and minimum evolution (default)" << std::endl;
    }

    // Transfer thread control parameters to HMMTree instance
    test.prc_threads_count = strc_cmd.prc_threads;
    test.phylo_threads_count = strc_cmd.phylo_threads;
    test.phylo_concurrent_threads_count = strc_cmd.phylo_concurrent_threads;
    test.prc_backend = strc_cmd.prc_backend;
    std::cout << "PRC backend preference: " << (test.prc_backend == HMMTree::PRC_BACKEND_AUTO ? "auto" : (test.prc_backend == HMMTree::PRC_BACKEND_LEGACY ? "legacy" : "prcx")) << std::endl;
    
    std::cout << "Thread configuration - PRC analysis: " << 
        (test.prc_threads_count == 0 ? "auto-detect" : std::to_string(test.prc_threads_count)) << 
        ", Phylogenetic analysis: " << 
        (test.phylo_threads_count == 0 ? "auto-detect" : std::to_string(test.phylo_threads_count)) << 
        ", Concurrent phylogenetic executions: " <<
        (test.phylo_concurrent_threads_count == 0 ? "auto-detect" : std::to_string(test.phylo_concurrent_threads_count))
        << std::endl;

    std::string str_file_path_2="";
    if(strc_cmd.als_phmms || strc_cmd.als_phhms){
        str_file_path_2=argv[argc - 2];
    }

    //time count
    struct timeb startTime , endTime;

    if(test.prc_hhsuite == 0){
        if(strc_cmd.uals){
            ftime(&startTime);
            std::cout << " '-prc' mode, run in '-uals' mode !"<< std::endl;
            // Banner suppressed per request
            //test.CLUSTER_OK_THRESHOLD = double_user_id;
            test.process_fasta_sequences(str_file_path, test.double_user_id);
            ftime(&endTime);
            std::cout << "'-prc' mode in '-uals' mode runtime: " << format_time_duration((endTime.time-startTime.time)*1000 + (endTime.millitm - startTime.millitm)) << std::endl;
        }

        if(strc_cmd.als){
            ftime(&startTime);
            std::cout << " '-prc' mode, run in '-als' mode !" << std::endl;
            // Banner suppressed per request
            test.process_prc_alignments(str_file_path);
            ftime(&endTime);
            std::cout << "'-prc' mode in '-als' mode runtime: " << format_time_duration((endTime.time-startTime.time)*1000 + (endTime.millitm - startTime.millitm)) << std::endl;
		}

		if(strc_cmd.als_phmms){
            HMMTree::als_phmms_phhms = true;
            ftime(&startTime);
            std::cout << " '-prc' mode, run in '-als_phmms' mode !" << std::endl;
            // Banner suppressed per request
            test.process_prc_alignments_phmms(str_file_path, str_file_path_2);
            ftime(&endTime);
            std::cout << "'-prc' mode in '-als_phmms' mode runtime: " << format_time_duration((endTime.time-startTime.time)*1000 + (endTime.millitm - startTime.millitm)) << std::endl;
		}


		if(strc_cmd.prc_hmms){
           ftime(&startTime);
            std::cout << " '-prc' mode, run in '-hmms' mode !"<< std::endl;
            // Banner suppressed per request
            test.process_prc_HMMs(str_file_path);
            ftime(&endTime);
            std::cout << "'-prc' mode in '-hmms' mode runtime: " << format_time_duration((endTime.time-startTime.time)*1000 + (endTime.millitm - startTime.millitm)) << std::endl;

		}
    }else{
        if(strc_cmd.uals){
            ftime(&startTime);
            std::cout << " '-hhsuite' mode, run in '-uals' mode !"<< std::endl;
            std::cout << "usearch->mafft->hmmbuild->hhsuite->distance matrix" << std::endl;
            //test.CLUSTER_OK_THRESHOLD = double_user_id;
            test.process_hhsuite_fasta_sequences(str_file_path, test.double_user_id);
            ftime(&endTime);
            std::cout << "'-hhsuite' mode in '-uals' mode runtime: " << format_time_duration((endTime.time-startTime.time)*1000 + (endTime.millitm - startTime.millitm)) << std::endl;
        }

        if(strc_cmd.als){
            ftime(&startTime);
            std::cout <<" '-hhsuite' mode, run in '-als' mode !" << std::endl;
            std::cout << "hhmbuild->hhsuite->distance matrix" << std::endl;
            test.process_hhsuite_alignments(str_file_path);
            ftime(&endTime);
            std::cout << "'-hhsuite' mode in '-als' mode runtime: " << format_time_duration((endTime.time-startTime.time)*1000 + (endTime.millitm - startTime.millitm)) << std::endl;
		}


		if(strc_cmd.als_phhms){
		    HMMTree::als_phmms_phhms= true;
            ftime(&startTime);
            std::cout <<" '-hhsuite' mode, run in '-als_phhms' mode !" << std::endl;
            std::cout << "als->hhmbuild->hhm1s, (hhm1s + hhm2s)->hhsuite->distance matrix" << std::endl;
            test.process_hhsuite_alignments_phhms(str_file_path, str_file_path_2);
            ftime(&endTime);
            std::cout << "'-hhsuite' mode in '-als_phhms' mode runtime: " << format_time_duration((endTime.time-startTime.time)*1000 + (endTime.millitm - startTime.millitm)) << std::endl;
		}

		if(strc_cmd.hhsuite_hhms){
           ftime(&startTime);
            std::cout << " '-hhsuite' mode, run in '-hhms' mode !"<< std::endl;
            std::cout << "hhsuite->distance matrix" << std::endl;
            test.process_hhsuite_HHMs(str_file_path);
            ftime(&endTime);
            std::cout << "'-hhsuite' mode in '-hhms' mode runtime: " << format_time_duration((endTime.time-startTime.time)*1000 + (endTime.millitm - startTime.millitm)) << std::endl;

		}
    }



	return 0;
}
