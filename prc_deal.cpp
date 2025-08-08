#include "HMMTree.h"
#include <omp.h>
#include <sstream>

//compute the distance of two hmms
int HMMTree::prc_each2(){
	std::vector<std::string> vec_hmms_filenames;
	if(!get_file_names(folder_hmms,vec_hmms_filenames,"")){
        return 0;
	}
	if(vec_hmms_filenames.size() < 3){
        return 0;
	}
	size_t unint_hmms = 0;
	std::vector <std::string>  hmm_names;
	while (unint_hmms < vec_hmms_filenames.size())
	{
		if(!hmm_file_exist_format_check(folder_hmms + vec_hmms_filenames[unint_hmms]))
		{
			vec_hmms_filenames.clear();
			return 0;
		}
		hmm_names.push_back(folder_hmms+ vec_hmms_filenames[unint_hmms]);
		unint_hmms++;
	}
	vec_hmms_filenames.clear();

	//init the matrix, vector and map;
    if(!matrix_init_matrix_vecotr_map()){
        std::cout<<"Fatal error: matrix_init_matrix_vecotr_map() failed !"<<std::endl;
        exit(1);
	}

    // Create list of all pairs to process
    std::vector<std::pair<unsigned int, unsigned int>> pairs_to_process;
    for (unsigned int i = 0; i < hmm_names.size() - 1; i++) {
        for (unsigned int j = i + 1; j < hmm_names.size(); j++) {
            pairs_to_process.push_back(std::make_pair(i, j));
        }
    }

    std::cout << "PRC deal: Processing " << pairs_to_process.size() << " pairwise comparisons..." << std::endl;
    
    // Track progress with thread-safe counter
    int completed_count = 0;
    int total_pairs = pairs_to_process.size();

    // Parallelize the pairwise comparisons
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t pair_idx = 0; pair_idx < pairs_to_process.size(); pair_idx++) {
        unsigned int i_hmm_names1 = pairs_to_process[pair_idx].first;
        unsigned int i_hmm_names2 = pairs_to_process[pair_idx].second;
        
        std::string str_name_hmm1 = hmm_names[i_hmm_names1];
        std::string str_name_hmm2 = hmm_names[i_hmm_names2];

        // Get thread ID for unique temporary handling
        int thread_id = omp_get_thread_num();
        
        // Thread-local variables to avoid race conditions
        std::vector<STUC_RHHEL_NOTE> local_res_msg;
        // Extract just the filename from the full path for matrix lookup
        std::string local_hmm1 = str_name_hmm1.substr(str_name_hmm1.find_last_of('/') + 1);
        std::string local_hmm2 = str_name_hmm2.substr(str_name_hmm2.find_last_of('/') + 1);
        
        std::string str_cmd = "";
        if(prc_hit_no != 0){
            std::string str_prc_hit_num = int_2_string(prc_hit_no);
            if(!bool_PRC_in_folder){
                str_cmd = "OMP_NUM_THREADS=1 prc -hits " + str_prc_hit_num + " " + str_name_hmm1 + " " + str_name_hmm2;
            }else{
                str_cmd = "OMP_NUM_THREADS=1 ./prc -hits " + str_prc_hit_num + " " + str_name_hmm1 + " " + str_name_hmm2;
            }
        }else{
            if(!bool_PRC_in_folder){
                str_cmd = "OMP_NUM_THREADS=1 prc " + str_name_hmm1 + " " + str_name_hmm2;
            }else{
                str_cmd = "OMP_NUM_THREADS=1 ./prc " + str_name_hmm1 + " " + str_name_hmm2;
            }
        }

        FILE *stream = popen(str_cmd.c_str(), "r");
        if (stream != nullptr) {
            prc_read_result_from_file_threadsafe(stream, local_res_msg);
            pclose(stream);
            
            // Thread-safe matrix update
            #pragma omp critical
            {
                matrix_get_each2_hmms_result_2_threadsafe(local_res_msg, local_hmm1, local_hmm2);
                completed_count++;
                
                // Thread-safe progress output
                std::string hmm1_name = str_name_hmm1.substr(str_name_hmm1.find_last_of('/') + 1);
                std::string hmm2_name = str_name_hmm2.substr(str_name_hmm2.find_last_of('/') + 1);
                
                std::cout << "[" << completed_count << "/" << total_pairs << "] " 
                         << hmm1_name << " vs " << hmm2_name 
                         << " (Thread " << thread_id << ")" << std::endl;
            }
        }
    }
    
    std::cout << "PRC deal completed: " << completed_count << " comparisons processed." << std::endl;
	return 1;
}

//prc hmm1 to a library
int HMMTree::prc_hmm1_library(std::string str_name_hmm1, std::string str_library,std::string prc_out_path_filename)
{
	std::string str_cmd = "";
    if(prc_hit_no != 0){
        std::string str_prc_hit_num =int_2_string(prc_hit_no);
        if(!bool_PRC_in_folder){
            str_cmd = "OMP_NUM_THREADS=1 prc -hits "+str_prc_hit_num+" " + str_name_hmm1 + " " + str_library + "  " + prc_out_path_filename;
        }else{
            str_cmd = "OMP_NUM_THREADS=1 ./prc -hits "+str_prc_hit_num+" " + str_name_hmm1 + " " + str_library + "  "+ prc_out_path_filename;
        }
    }else{
        if(!bool_PRC_in_folder){
            str_cmd = "OMP_NUM_THREADS=1 prc " + str_name_hmm1 + " " + str_library + "  "+ prc_out_path_filename;
        }else{
            str_cmd = "OMP_NUM_THREADS=1 ./prc " + str_name_hmm1 + " " + str_library + "  "+ prc_out_path_filename;
        }
    }

	return system(str_cmd.c_str());
}

//compute the distance of hmm1 to a library
int HMMTree::prc_library()
{
	std::vector<std::string> vec_hmms_filenames;
	if(!get_file_names(folder_hmms,vec_hmms_filenames,"")){
        return 0;
	}
	if(vec_hmms_filenames.size() < 3){
        return 0;
	}
	size_t unint_hmms = 0;
	std::vector <std::string>  hmm_names;
	while (unint_hmms < vec_hmms_filenames.size())
	{
		if(!hmm_file_exist_format_check(folder_hmms + vec_hmms_filenames[unint_hmms]))
		{
			vec_hmms_filenames.clear();
			return 0;
		}
		hmm_names.push_back(folder_hmms+ vec_hmms_filenames[unint_hmms]);
		unint_hmms++;
	}
	vec_hmms_filenames.clear();

	// check the hmm files version is HMMER2 or not
	std::ifstream file_test_hmm_version;
	std::string str_oneline_file="";
	file_test_hmm_version.open(hmm_names[0].c_str());
	if (!file_test_hmm_version.is_open())
	{
		std::cout << "Failed to open hmm file !" << std::endl;
		return 0;
	}
	if(file_test_hmm_version.eof())
	{
		std::cout << "The hmm file is empty !" << std::endl;
		return 0;
	}
	std::getline(file_test_hmm_version, str_oneline_file);
	if(str_oneline_file.find("HMMER2") == -1)
	{
		std::cout<<"HMMER version error, PRC can only deal hmm files of HMMER2 version!"<<std::endl;
		file_test_hmm_version.close();
		return 2;
	}
	file_test_hmm_version.close();

	//init the matrix, vector and map;
    if(!matrix_init_matrix_vecotr_map()){
        std::cout<<"Fatal error: matrix_init_matrix_vecotr_map() failed !"<<std::endl;
        exit(1);
	}

	//if the prcfiles folder is not empty then empty it
	if(dir_noempty_opendir_readir(folder_prcfiles))
	{
		//system_return(system("rm ./prcfiles/*.txt"));
		if(!delete_files(folder_prcfiles)){
            return 0;
		}
	}
    unsigned int uint_computed_num = 0;
    std::string str_lib_list=files_folder+"list.lib";
    std::cout<<"PRC deal: "<<std::endl;
	for (unsigned int i_hmm_names1 = 0; i_hmm_names1 < hmm_names.size() - 1; i_hmm_names1++)
	{
        std::string str_hmm_names = "";
        str_hmm_names = hmm_names[i_hmm_names1].substr(hmm_names[i_hmm_names1].find_last_of('/')+1);
        std::cout<<std::endl<<std::setiosflags(std::ios::left)<<std::setw(10)<<str_hmm_names<<" :  "<<i_hmm_names1<<std::endl;

        std::string str_prc_num =uint2str(i_hmm_names1);
        std::string str_lib_list =folder_prcfiles + "list.lib";

        std::ofstream if_library;
        if_library.open(str_lib_list.c_str());
        if(!if_library.is_open()){
            std::cout<<"prc_library(): Failed to oprn list.lib"<<std::endl;
            hmm_names.clear();
            return 0;
        }
		for (unsigned int i_hmm_names2 = i_hmm_names1 + 1; i_hmm_names2 < hmm_names.size(); i_hmm_names2++)
		{
            if_library<<hmm_names[i_hmm_names2]<<std::endl;
		}
		if_library.close();
        std::string str_lib_result_file =folder_prcfiles + str_prc_num+".scores";
        bool error_flag = false;
        int return_flag = 0;
        return_flag = prc_hmm1_library(hmm_names[i_hmm_names1], str_lib_list,folder_prcfiles + str_prc_num);
        if(system_return(return_flag)){
            std::cout<<std::endl;
            return 0;
        }
        prc_read_lib_result_from_file(str_lib_result_file);

        matrix_get_lib_hmms_result_2();

        system_return(system(("rm -f "+str_lib_result_file).c_str()));
        system_return(system(("rm -f "+str_lib_list).c_str()));
	}
	return 1;
}

void HMMTree::prc_read_lib_result_from_file(std::string str_prc_lib_result_filename){
    std::ifstream i_result_file;
	i_result_file.open(str_prc_lib_result_filename.c_str());
	//std::cout << str_result_file_path.c_str() << std::endl;
	if (!i_result_file.is_open())
	{
		std::cout << std::endl << "Failed to reading data from result file, please be sure that the file path is right!" << std::endl;
		exit(0);
	}

	std::string str_msg_temp;                 //temp variable to save the msg readed from the file
	int line_num_temp = 0;                        //temp variable to count the line of the data
	bool bool_flag = false;
	while (std::getline(i_result_file, str_msg_temp))
	{
		line_num_temp++;
        if(str_msg_temp.find("# hmm1") != -1){
            bool_flag = true;
            continue;
        }
		if (bool_flag)
		{
			if (str_msg_temp.find("# END") != -1)
			{
				break;
			}
			//temp variable to save the line messages of the result file
			STUC_RHHEL_NOTE struct_line_temp;

			//push the varibale to the calss variable to save it
			res_msg.push_back(struct_line_temp);

			//split the string to get messages
			std::vector<std::string> result_line_split;
			result_line_split = str_Split_by_char_list(str_msg_temp, " \t");

			// Skip malformed lines that don't have enough fields
			if(result_line_split.size() < 12){
			    res_msg.pop_back();
			    continue;
			}

			//put the message into the class result array
			res_msg.back().hmm1_name = normalize_prc_identifier_to_filename(result_line_split[0]);
			res_msg.back().begin_hmm1 = atoi(result_line_split[1].c_str());
			res_msg.back().end_hmm1 = atoi(result_line_split[2].c_str());
			res_msg.back().hit_no = atoi(result_line_split[4].c_str());
			res_msg.back().hmm2_name = normalize_prc_identifier_to_filename(result_line_split[5]);
			res_msg.back().begin_hmm2 = atoi(result_line_split[6].c_str());
			res_msg.back().end_hmm2 = atoi(result_line_split[7].c_str());
			res_msg.back().co_emis = atof(result_line_split[9].c_str());
			res_msg.back().simple = atof(result_line_split[10].c_str());
			res_msg.back().reverse = atof(result_line_split[11].c_str());

			// Skip entries with empty HMM names
			if(res_msg.back().hmm1_name.empty() || res_msg.back().hmm2_name.empty()){
			    res_msg.pop_back();
			}

		}
		if (line_num_temp > 4 && line_num_temp < 13)
		{
			head_msg[line_num_temp - 5] = str_msg_temp;
		}
		if (line_num_temp < 4)
		{
			PRC_CopyRight_license += str_msg_temp;
		}

		//std::cout << line_num_temp << "  " << str_msg_temp << std::endl;
	}

	i_result_file.close();
	return ;
}

//read into the  message from the file
void HMMTree::prc_read_result_from_file(FILE * file_stream)
{
    char   buf[1024];
    bool bool_flag = false;
    bool get_names = false;
    int line_num_temp = 0;
    while(fgets(buf,1024,file_stream)){
        std::string str_msg_temp = buf;
        line_num_temp++;

        if(str_msg_temp.find("# hmm1") != -1){
            bool_flag = true;
            continue;
        }
		if (bool_flag)
		{
			if (str_msg_temp.find("# END") != -1)
			{
				break;
			}
			//temp variable to save the line messages of the result file
			STUC_RHHEL_NOTE struct_line_temp;

			//push the varibale to the calss variable to save it
			res_msg.push_back(struct_line_temp);

			//split the string to get messages
			std::vector<std::string> result_line_split;
			result_line_split = str_Split_by_char_list(str_msg_temp, " \t");

			// Skip malformed lines that don't have enough fields
			if(result_line_split.size() < 12){
			    res_msg.pop_back();
			    continue;
			}

			//put the message into the class result array
			res_msg.back().hmm1_name = normalize_prc_identifier_to_filename(result_line_split[0]);
			res_msg.back().begin_hmm1 = atoi(result_line_split[1].c_str());
			res_msg.back().end_hmm1 = atoi(result_line_split[2].c_str());
			res_msg.back().hit_no = atoi(result_line_split[4].c_str());
			res_msg.back().hmm2_name = normalize_prc_identifier_to_filename(result_line_split[5]);
			res_msg.back().begin_hmm2 = atoi(result_line_split[6].c_str());
			res_msg.back().end_hmm2 = atoi(result_line_split[7].c_str());
			res_msg.back().co_emis = atof(result_line_split[9].c_str());
			res_msg.back().simple = atof(result_line_split[10].c_str());
			res_msg.back().reverse = atof(result_line_split[11].c_str());

			if(!get_names){
                hmm1 = res_msg.back().hmm1_name;
                hmm2 = res_msg.back().hmm2_name;
                get_names = true;
			}
			
			// Skip entries with empty HMM names
			if(res_msg.back().hmm1_name.empty() || res_msg.back().hmm2_name.empty()){
			    res_msg.pop_back();
			}
		}
		if (line_num_temp > 4 && line_num_temp < 13)
		{
			head_msg[line_num_temp - 5] = str_msg_temp;
		}
		if (line_num_temp < 4)
		{
			PRC_CopyRight_license += str_msg_temp;
		}
    }
	return ;
}

//compute the distence
double HMMTree::prc_hmms_dist_compute()
{
    double answer=0;
    answer=res_msg[0].simple;

    //std::cout<<answer<<"   "<<std::endl;

    if(answer <= 0){
        answer = DBL_MAX_USER;
    }
	return 1.0/answer;
}

//Thread-safe version that uses local res_msg instead of class member
double HMMTree::prc_hmms_dist_compute_threadsafe(const std::vector<STUC_RHHEL_NOTE>& local_res_msg)
{
    double answer=0;
    if(!local_res_msg.empty()){
        answer=local_res_msg[0].simple;
    }

    //std::cout<<answer<<"   "<<std::endl;

    if(answer <= 0){
        answer = DBL_MAX_USER;
    }
	return 1.0/answer;
}

//set the hmm1 and hmm2 message to the result matrix class
void HMMTree::prc_set_STUC_RHH_NOTE_dist(STUC_RHH_NOTE* STUC_RHH_NOTE_res)
{
	STUC_RHH_NOTE_res->hmm1.name = hmm1;
	STUC_RHH_NOTE_res->hmm2.name = hmm2;
	STUC_RHH_NOTE_res->distance = prc_hmms_dist_compute();
}

//Thread-safe version that takes local parameters instead of using class members
void HMMTree::prc_set_STUC_RHH_NOTE_dist_threadsafe(STUC_RHH_NOTE* STUC_RHH_NOTE_res, const std::vector<STUC_RHHEL_NOTE>& local_res_msg, const std::string& local_hmm1, const std::string& local_hmm2)
{
	STUC_RHH_NOTE_res->hmm1.name = local_hmm1;
	STUC_RHH_NOTE_res->hmm2.name = local_hmm2;
	STUC_RHH_NOTE_res->distance = prc_hmms_dist_compute_threadsafe(local_res_msg);
}

//check the profile HMM files is HMMER2.O or not
int HMMTree::prc_check_profile_HMM_format(){
    std::vector<std::string> vec_hmm3_2_filenames;
	if(!get_file_names(folder_hmms,vec_hmm3_2_filenames,"")){
	    std::cout<<"Get the file name error, in prc_check_profile_HMM_format()"<<std::endl;
        return 0;
	}
	if(vec_hmm3_2_filenames.size() < 3){
	    std::cout<<"Please be sure there are at least three Profile HMM files !"<<std::endl;
        return 0;
	}

	for(unsigned int uint_hmms_i=0; uint_hmms_i < vec_hmm3_2_filenames.size();uint_hmms_i++){

        std::string str_path_name = folder_hmms+vec_hmm3_2_filenames[uint_hmms_i];
        // check the hmm files version is HMMER2 or not
        std::ifstream file_test_hmm_version;
        std::string str_oneline_file="";
        file_test_hmm_version.open(str_path_name.c_str());
        if (!file_test_hmm_version.is_open())
        {
            std::cout << "Failed to open hmm file, in prc_check_profile_HMM_format()" << std::endl;
            return 0;
        }
        if(file_test_hmm_version.eof())
        {
            std::cout << "The hmm file is empty, in prc_check_profile_HMM_format()" << std::endl;
            return 0;
        }
        std::getline(file_test_hmm_version, str_oneline_file);
        if(str_oneline_file.find("HMMER2") == -1)
        {
            std::cout<<"HMMER version error, PRC can only deal hmm files of HMMER2 verson !"<<std::endl;
            file_test_hmm_version.close();
            return 2;
        }
        file_test_hmm_version.close();

	}

	return 1;
}

//Thread-safe version that uses local variables instead of class members
void HMMTree::prc_read_result_from_file_threadsafe(FILE * file_stream, std::vector<STUC_RHHEL_NOTE>& local_res_msg)
{
    char   buf[1024];
    bool bool_flag = false;
    int line_num_temp = 0;
    while(fgets(buf,1024,file_stream)){
        std::string str_msg_temp = buf;
        line_num_temp++;

        if(str_msg_temp.find("# hmm1") != -1){
            bool_flag = true;
            continue;
        }
		if (bool_flag)
		{
			if (str_msg_temp.find("# END") != -1)
			{
				break;
			}
			//temp variable to save the line messages of the result file
			STUC_RHHEL_NOTE struct_line_temp;

			//push the variable to the local vector instead of class member
			local_res_msg.push_back(struct_line_temp);

			//split the string to get messages
			std::vector<std::string> result_line_split;
			result_line_split = str_Split_by_char_list(str_msg_temp, " \t");

			// Skip malformed lines that don't have enough fields
			if(result_line_split.size() < 12){
			    local_res_msg.pop_back();
			    continue;
			}

			//put the message into the local result array
			local_res_msg.back().hmm1_name = normalize_prc_identifier_to_filename(result_line_split[0]);
			local_res_msg.back().begin_hmm1 = atoi(result_line_split[1].c_str());
			local_res_msg.back().end_hmm1 = atoi(result_line_split[2].c_str());
			local_res_msg.back().hit_no = atoi(result_line_split[4].c_str());
			local_res_msg.back().hmm2_name = normalize_prc_identifier_to_filename(result_line_split[5]);
			local_res_msg.back().begin_hmm2 = atoi(result_line_split[6].c_str());
			local_res_msg.back().end_hmm2 = atoi(result_line_split[7].c_str());
			local_res_msg.back().co_emis = atof(result_line_split[9].c_str());
			local_res_msg.back().simple = atof(result_line_split[10].c_str());
			local_res_msg.back().reverse = atof(result_line_split[11].c_str());
			
			// Skip entries with empty HMM names
			if(local_res_msg.back().hmm1_name.empty() || local_res_msg.back().hmm2_name.empty()){
			    local_res_msg.pop_back();
			}
		}
    }
	return ;
}

