#include "HMMTree.h"

//hhmake build one hhm
int HMMTree::hhsuite_hhmake_one(std::string family_name){

    std::string str_cmd = "";
	std::string str_pickname_outof_ext = "";
	int begin = family_name.find_last_of('/') + 1;
	int end = family_name.find_last_of('_');
	str_pickname_outof_ext = family_name.substr(begin, end - begin);
	//str_pickname_outof_ext = family_name.substr(family_name.find_last_of('/') + 1, family_name.find_last_of('_'));
	if(!bool_hhsuite_hhmake_in_folder){
        str_cmd = "hhmake  -i " + family_name + " -name " + str_pickname_outof_ext + " -o " +folder_hhms+ str_pickname_outof_ext + ".hhm ";
	}else{
        str_cmd = "./hhmake  -i " + family_name + " -name " + str_pickname_outof_ext + " -o " +folder_hhms+ str_pickname_outof_ext + ".hhm ";
	}

	return system(str_cmd.c_str());
}
//hhmake build all hhm
int HMMTree::hhsuite_hhmake_all(){

    if(dir_noempty_opendir_readir(folder_hhms))
	{
		//system_return(system("rm ./hhms/*.hhm"));
		if(!delete_files(folder_hhms)){
            return 0;
		}
	}
	std::vector<std::string> vec_aligned_filenames;
	if(!get_file_names(folder_aligned,vec_aligned_filenames,"")){
        return 0;
	}

	if(vec_aligned_filenames.size() < 3 && (!HMMTree::als_phmms_phhms)){
        return 0;
	}

	size_t unint_aligned = 0;
	bool error_flag = false;
	while (unint_aligned < vec_aligned_filenames.size())
	{
		if(!align_fasta_file_exist_format_check(folder_aligned+vec_aligned_filenames[unint_aligned]))
		{
			vec_aligned_filenames.clear();
			error_flag = true;
			break;
		}
		int return_flag = 0;
		return_flag = hhsuite_hhmake_one(folder_aligned+vec_aligned_filenames[unint_aligned]);
		if(system_return(return_flag)){
            error_flag = true;
		}
		if (error_flag)
		{
			vec_aligned_filenames.clear();
			break;
		}
		unint_aligned++;
	}
	if (!error_flag)
	{
		vec_aligned_filenames.clear();
	}
	else {
		return 0;
	}
	return 1;
}

//hhalign
int HMMTree::hhsuite_hhalign_each2(){

    std::vector<std::string> vec_hhms_filenames;
	if(!get_file_names(folder_hhms,vec_hhms_filenames,"")){
        return 0;
	}

	if(vec_hhms_filenames.size() < 3){
        return 0;
	}

	size_t unint_hhms = 0;
	std::vector <std::string>  hhm_names;
	while (unint_hhms < vec_hhms_filenames.size())
	{
		if(!file_exists_and_empty_check(folder_hhms + vec_hhms_filenames[unint_hhms]))
		{
			vec_hhms_filenames.clear();
			return 0;
		}
		hhm_names.push_back(folder_hhms+ vec_hhms_filenames[unint_hhms]);
		unint_hhms++;
	}
	vec_hhms_filenames.clear();


	//init the matrix, vector and map;
	if(!matrix_init_matrix_vector_map_hhsuite()){
		std::cout<<"Fatal error: matrix_init_matrix_vector_map_hhsuite() failed !"<<std::endl;
        exit(1);
	}

	/* removed unused counter uint_computed_num */
    std::cout<<"hhsuite processing: ";
	for (unsigned int i_hhm_names1 = 0; i_hhm_names1 < hhm_names.size() - 1; i_hhm_names1++)
	{
        std::string str_hhm_names = "";
        str_hhm_names = hhm_names[i_hhm_names1].substr(hhm_names[i_hhm_names1].find_last_of('/')+1);
        std::cout<<std::endl<<std::setiosflags(std::ios::left)<<std::setw(10)<<str_hhm_names<<" "<<std::endl;
		for (unsigned int i_hhm_names2 = i_hhm_names1 + 1; i_hhm_names2 < hhm_names.size(); i_hhm_names2++)
		{
            std::string str_name_hhm1=hhm_names[i_hhm_names1];
            std::string str_name_hhm2=hhm_names[i_hhm_names2];

            std::string str_cmd1 = "";
            std::string str_cmd2 = "";
            if(!bool_hhsuite_hhalign_in_folder){
                str_cmd1 = "hhalign -i " + str_name_hhm1 + " -t " + str_name_hhm2;
                str_cmd2 = "hhalign -i " + str_name_hhm2 + " -t " + str_name_hhm1;
            }else{
                str_cmd2 = "./hhalign -i " + str_name_hhm1 + " -t " + str_name_hhm2;
                str_cmd2 = "./hhalign -i " + str_name_hhm2 + " -t " + str_name_hhm1;
            }

            FILE   *stream1, *stream2;

            stream1 = popen(str_cmd1.c_str(), "r" ); //
            stream2 = popen(str_cmd2.c_str(), "r" ); //

            hhsuite_read_result_from_file_hhalign(stream1,stream2);
            pclose( stream1 );
            pclose( stream2 );
            matrix_get_each2_hhms_result_hhsuite();
			std::cout<<". "<<std::flush;
		}
		std::cout<<std::endl;
	}
	return 1;
}

void HMMTree::hhsuite_read_result_from_file_hhalign(FILE * file_stream1, FILE * file_stream2){
    char   buf[1024];
    bool get_names = false;
    double stream1_score = 0.0;
    double stream1_evalue = 0.0;
    while(fgets(buf,1024,file_stream1)){

        std::string str_msg_temp = buf;

        if(str_msg_temp.find("Score = ") != -1){
            //temp variable to save the line messages of the result file
			STUC_RHHEL_NOTE struct_line_temp;

			//push the varibale to the calss variable to save it
			res_msg.push_back(struct_line_temp);

			//split the string to get messages
			std::vector<std::string> result_line_split;
			result_line_split = str_Split_by_char_list(str_msg_temp, " \t");

			//put the message into the class result array
			res_msg.back().hmm1_name = result_line_split[1];
			//std::cout<<"hmm1_name: "<<res_msg.back().hmm1_name<<std::endl;
			result_line_split[3] = result_line_split[3].substr(0,result_line_split[3].length()-1);
			res_msg.back().hmm2_name = result_line_split[3];
			//std::cout<<"hmm2_name: "<<res_msg.back().hmm2_name<<std::endl;

			stream1_score= atof(result_line_split[6].c_str());
			//std::cout<<"simple: "<<res_msg.back().simple<<std::endl;
			stream1_evalue = atof(result_line_split[9].c_str());
			//std::cout<<"simple: "<<res_msg.back().E_value<<std::endl;


			if(!get_names){
                hmm1 = res_msg.back().hmm1_name;
                hmm2 = res_msg.back().hmm2_name;
                get_names = true;
			}
            break;
        }
    }

     double stream2_score = 0.0;
     double stream2_evalue = 0.0;
     while(fgets(buf,1024,file_stream2)){

        std::string str_msg_temp = buf;

        if(str_msg_temp.find("Score = ") != -1){


			//split the string to get messages
			std::vector<std::string> result_line_split;
			result_line_split = str_Split_by_char_list(str_msg_temp, " \t");

			//res_msg.back().simple = atof(result_line_split[6].c_str());  
			stream2_score = atof(result_line_split[6].c_str());
			//std::cout<<"simple: "<<  res_msg.back().simple<<std::endl;
			stream2_evalue = atof(result_line_split[9].c_str());
			//std::cout<<"simple: "<<res_msg.back().E_value<<std::endl;


			if(!get_names){
                hmm1 = res_msg.back().hmm1_name;
                hmm2 = res_msg.back().hmm2_name;
                get_names = true;
			}
            break;
        }
    }
    res_msg.back().simple = (stream1_score+stream2_score)/2.0;
    res_msg.back().E_value = (stream1_evalue+stream2_evalue)/2.0;
	return ;
}


//copy the hhm files from user dir to hmms
int HMMTree::hhsuite_copy_hhmfiles(std::string path){

    if (!dir_exist_opendir(path))
	{
        std::cout<<"Failed to open the file directory !"<<std::endl;
		return 0;
	}
    create_files_folder(path,prc_hhsuite,2);

    //test if the dir exists, or create them
    dir_exist_or_create(prc_hhsuite, 2);

    if(dir_noempty_opendir_readir(folder_hhms))
	{
		//system_return(system("rm ./hmms/*.hmm"));
		if(!delete_files(folder_hhms)){
            return 0;
		}
	}


	std::string  path_temp=path;
    if(path_temp[path_temp.length()-1] != '/'){
        path_temp = path_temp+"/";
    }
    copy_files(path_temp,folder_hhms,".hhm");

    return 1;
}


//copy the hhm files from user dir to hmms
int HMMTree::hhsuite_als_phhms_copy_hhmfiles(std::string path){

    if (!dir_exist_opendir(path))
	{
        std::cout<<"Failed to open the file directory !"<<std::endl;
		return 0;
	}
    //create_files_folder(path,prc_hhsuite,2);

    //test if the dir exists, or create them
    //dir_exist_or_create(prc_hhsuite, 2);
/*
    if(dir_noempty_opendir_readir(folder_hhms))
	{
		//system_return(system("rm ./hmms/*.hmm"));
		if(!delete_files(folder_hhms)){
            return 0;
		}
	}

*/

        //test the hmms folder
    if(!dir_exist_opendir(folder_hhms_from_als))
    {
        system_return(system(("mkdir  -m 777 "+folder_hhms_from_als).c_str()));
    }

    copy_files(folder_hhms,folder_hhms_from_als,".hhm");

	std::string  path_temp=path;
    if(path_temp[path_temp.length()-1] != '/'){
        path_temp = path_temp+"/";
    }
    copy_files(path_temp,folder_hhms,".hhm");

    return 1;
}

