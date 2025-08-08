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

    if(argc == 2){
        std::string str_argv = argv[1];
        if(str_argv == "-h"){
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

	std::vector<std::string> vec_arguments;
    vec_arguments.push_back("-prc"); //0
	vec_arguments.push_back("-uals"); //1
	vec_arguments.push_back("-als"); //2
	vec_arguments.push_back("-hmms"); //3
	vec_arguments.push_back("-id"); //4
	vec_arguments.push_back("-prc_hit"); //5
	vec_arguments.push_back("-acc"); //6
	vec_arguments.push_back("-lib"); //7

	vec_arguments.push_back("-hhsuite"); //8
	vec_arguments.push_back("-hhms"); //9aln_phmm

	vec_arguments.push_back("-als_phmms"); //10
	vec_arguments.push_back("-als_phhms"); //11
	//string to save the compute model from the user input

    CMD_PARAMS strc_cmd;

    std::string str_file_path =argv[argc - 1];

    int arg_num = 1;
    while(arg_num < argc - 1 && (!strc_cmd.als_phmms || arg_num < argc - 2) && (!strc_cmd.als_phhms || arg_num < argc - 2)){
        std::string str_argu_temp =argv[arg_num];
        bool find_argu = false;
        int arg_vec_num = 0;
        while(arg_vec_num <= vec_arguments.size()-1){
            if(strcmp(str_argu_temp.c_str(),vec_arguments[arg_vec_num].c_str()) == 0){
                switch(arg_vec_num){

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
                            output_error_("Run style '-uals' arguments( [-uals, -als, -als_phmms(prc), -als_phhms(hhsuite), -hmms(prc), -hhms(hhsuite)] )");
                        }

                        break;
                    case 2:
                        if(!strc_cmd.uals && !strc_cmd.als && !strc_cmd.als_phmms && !strc_cmd.als_phhms && !strc_cmd.hhsuite_hhms  && !strc_cmd.prc_hmms){
                            strc_cmd.als = true;
                            find_argu = true;
                        }else{
                            output_error_("Run style '-als' arguments( [-uals, -als, -als_phmms(prc), -als_phhms(hhsuite), -hmms(prc), -hhms(hhsuite)] )");
                        }

                        break;
                    case 3:
                        if(!strc_cmd.uals && !strc_cmd.als && !strc_cmd.als_phmms && !strc_cmd.als_phhms && !strc_cmd.hhsuite_hhms  && !strc_cmd.prc_hmms){
                            strc_cmd.prc_hmms = true;
                            find_argu = true;
                        }else{
                            output_error_("Run style '-hmms' arguments( [-uals, -als, -als_phmms(prc), -als_phhms(hhsuite), -hmms(prc), -hhms(hhsuite)] )");
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
                            output_error_("Run style '-hhms' arguments( [-uals, -als, -als_phmms(prc), -als_phhms(hhsuite), -hmms(prc), -hhms(hhsuite)] )");
                        }
                        break;
                    case 10:
                        if(!strc_cmd.als_phmms && !strc_cmd.uals && !strc_cmd.als && !strc_cmd.als_phhms && !strc_cmd.hhsuite_hhms  && !strc_cmd.prc_hmms){
                            strc_cmd.als_phmms = true;
                            find_argu = true;
                        }else{
                            output_error_("Run style '-als_phmms' arguments( [-uals, -als, -als_phmms(prc), -als_phhms(hhsuite), -hmms(prc), -hhms(hhsuite)] )");
                        }
                        break;

                    case 11:
                        if(!strc_cmd.als_phhms && !strc_cmd.als_phmms && !strc_cmd.uals && !strc_cmd.als && !strc_cmd.hhsuite_hhms  && !strc_cmd.prc_hmms){
                            strc_cmd.als_phhms = true;
                            find_argu = true;
                        }else{
                            output_error_("Run style '-als_phhms' arguments( [-uals, -als, -als_phmms(prc), -als_phhms(hhsuite), -hmms(prc), -hhms(hhsuite)] )");
                        }
                        break;
                    default:
                        output_error_("Arguments");
                        break;
                }
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
            output_error_("'-prc' mode cannot run in '-hhms' or '-als_phhms' style, arguments");
        }
        if(!strc_cmd.uals && !strc_cmd.als && !strc_cmd.prc_hmms && !strc_cmd.als_phmms){
            output_error_("Please at least one run style('-uals', '-als', '-als_phmms' or '-hmms' ) to run in '-prc' mode ,arguments");
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
            output_error_("'-hhsuite' mode cannot run in '-acc', '-hmms', '-pair' or '-prc_hit' style, arguments");
        }

        if(!strc_cmd.uals && !strc_cmd.als && !strc_cmd.hhsuite_hhms && !strc_cmd.als_phhms){
            output_error_("Please at least one run style('-uals', '-als', '-als_phhms' or '-hhms' ) to run in '-hhsuite' mode, arguments");
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

    std::string str_file_path_2="";
    if(strc_cmd.als_phmms || strc_cmd.als_phhms){
        str_file_path_2=argv[argc - 2];
    }

    //time count
    struct timeb startTime , endTime;

    if(test.prc_hhsuite == 0){
        if(strc_cmd.uals){
            ftime(&startTime);
            std::cout << " '-prc' mode, run in '-uals' style !"<< std::endl;
            std::cout << "usearch->mafft->hmmbuild->prc->distance matrix" << std::endl;
            //test.CLUSTER_OK_THRESHOLD = double_user_id;
            test.prc_fasta_seqs_deal(str_file_path, test.double_user_id);
            ftime(&endTime);
            std::cout << "'-prc' mode in '-uals' style run time: " << format_time_duration((endTime.time-startTime.time)*1000 + (endTime.millitm - startTime.millitm)) << std::endl;
        }

        if(strc_cmd.als){
            ftime(&startTime);
            std::cout << " '-prc' mode, run in '-als' style !" << std::endl;
            std::cout << "hhmbuild->prc->distance matrix" << std::endl;
            test.prc_alignments_deal(str_file_path);
            ftime(&endTime);
            std::cout << "'-prc' mode in '-als' style run time: " << format_time_duration((endTime.time-startTime.time)*1000 + (endTime.millitm - startTime.millitm)) << std::endl;
		}

		if(strc_cmd.als_phmms){
            HMMTree::als_phmms_phhms = true;
            ftime(&startTime);
            std::cout << " '-prc' mode, run in '-als_phmms' style !" << std::endl;
            std::cout << "als->hmmbuild->hmm1s, (hmm1s+hmm2s)->prc->distance matrix " << std::endl;
            test.prc_alignments_phmms_deal(str_file_path, str_file_path_2);
            ftime(&endTime);
            std::cout << "'-prc' mode in '-als_phmms' style run time: " << format_time_duration((endTime.time-startTime.time)*1000 + (endTime.millitm - startTime.millitm)) << std::endl;
		}


		if(strc_cmd.prc_hmms){
           ftime(&startTime);
            std::cout << " '-prc' mode, run in '-hmms' style !"<< std::endl;
            std::cout << "prc->distance matrix" << std::endl;
            test.prc_HMMs_deal(str_file_path);
            ftime(&endTime);
            std::cout << "'-prc' mode in '-hmms' style run time: " << format_time_duration((endTime.time-startTime.time)*1000 + (endTime.millitm - startTime.millitm)) << std::endl;

		}
    }else{
        if(strc_cmd.uals){
            ftime(&startTime);
            std::cout << " '-hhsutie' mode, run in '-uals' style !"<< std::endl;
            std::cout << "usearch->mafft->hmmbuild->hhsutie->distance matrix" << std::endl;
            //test.CLUSTER_OK_THRESHOLD = double_user_id;
            test.hhsuite_fasta_seqs_deal(str_file_path, test.double_user_id);
            ftime(&endTime);
            std::cout << "'-hhsutie' mode in '-uals' style run time: " << format_time_duration((endTime.time-startTime.time)*1000 + (endTime.millitm - startTime.millitm)) << std::endl;
        }

        if(strc_cmd.als){
            ftime(&startTime);
            std::cout <<" '-hhsutie' mode, run in '-als' style !" << std::endl;
            std::cout << "hhmbuild->hhsuite->distance matrix" << std::endl;
            test.hhsuite_alignments_deal(str_file_path);
            ftime(&endTime);
            std::cout << "'-hhsutie' mode in '-als' style run time: " << format_time_duration((endTime.time-startTime.time)*1000 + (endTime.millitm - startTime.millitm)) << std::endl;
		}


		if(strc_cmd.als_phhms){
		    HMMTree::als_phmms_phhms= true;
            ftime(&startTime);
            std::cout <<" '-hhsutie' mode, run in '-als_phhms' style !" << std::endl;
            std::cout << "als->hhmbuild->hhm1s, (hhm1s + hhm2s)->hhsuite->distance matrix" << std::endl;
            test.hhsuite_alignments_phhms_deal(str_file_path, str_file_path_2);
            ftime(&endTime);
            std::cout << "'-hhsutie' mode in '-als_phhms' style run time: " << format_time_duration((endTime.time-startTime.time)*1000 + (endTime.millitm - startTime.millitm)) << std::endl;
		}

		if(strc_cmd.hhsuite_hhms){
           ftime(&startTime);
            std::cout << " '-hhsutie' mode, run in '-hhms' style !"<< std::endl;
            std::cout << "hhsuite->distance matrix" << std::endl;
            test.hhsuite_HHMs_deal(str_file_path);
            ftime(&endTime);
            std::cout << "'-hhsutie' mode in '-hhms' style run time: " << format_time_duration((endTime.time-startTime.time)*1000 + (endTime.millitm - startTime.millitm)) << std::endl;

		}
    }



	return 0;
}
