#include "HMMTree.h"

//init the matrix 、 vector  and map
int HMMTree::matrix_init_matrix_vecotr_map()
{
	std::vector<std::string> vec_matrix_hmm_filenames;
	if(!get_file_names(folder_hmms,vec_matrix_hmm_filenames,".hmm")){
        return 0;
	}
	if(vec_matrix_hmm_filenames.size() < 3){
        return 0;
	}

	size_t unint_matrix_hmm = 0;
	size_t unint_name_num = 0;
	size_t unint_acc_num = 0;

	while (unint_matrix_hmm < vec_matrix_hmm_filenames.size())
	{
		std::ifstream  if_name_acc_temp;
		if_name_acc_temp.open(folder_hmms+vec_matrix_hmm_filenames[unint_matrix_hmm]);
		if(!if_name_acc_temp.is_open()){
            std::cout<<"matrix_init_matrix_vecotr_map(): Open "<<vec_matrix_hmm_filenames[unint_matrix_hmm]<<"failed !"<<std::endl;
            return 0;
		}
        HMM_NAME_ACC hmm_name_acc_temp;
        id_hmm_NAME_ACC_list_vector.push_back(hmm_name_acc_temp);
        id_hmm_NAME_ACC_list_vector.back().filename=folder_hmms+vec_matrix_hmm_filenames[unint_matrix_hmm];
		std::string str_online_temp = "";
		size_t unint_temp = 0;
		while(unint_temp < 5 && std::getline(if_name_acc_temp,str_online_temp) ){

            if(str_online_temp.find("NAME") != -1){
                std::string str_temp_0 = str_online_temp;
                str_temp_0 = str_temp_0.substr(str_temp_0.find_first_not_of(' '));
                str_temp_0 = str_temp_0.substr(str_temp_0.find_first_of(' '));
                str_temp_0 = str_temp_0.substr(str_temp_0.find_first_not_of(' '));
                std::vector<std::string> result_split =str_Split_by_char_list(str_temp_0,"  ");
                //str_temp_0 = str_temp_0.substr(str_temp_0.find_first_of(' '));
                //str_temp_0 = str_temp_0.substr(str_temp_0.find_first_not_of(' '), str_temp_0.find_first_of(' ') - str_temp_0.find_first_not_of(' '));

                str_temp_0 =result_split[0];
                if(hmm_NAME_unordered_map.find(str_temp_0) != hmm_NAME_unordered_map.end()){
                    std::cout<<"File '"<<id_hmm_NAME_ACC_list_vector.back().filename <<"' has the same 'NAME' with file '"<<id_hmm_NAME_ACC_list_vector[hmm_NAME_unordered_map[str_temp_0]].filename <<"' !"<<std::endl;
                    id_hmm_NAME_ACC_list_vector.clear();
                    return 0;
                }
                id_hmm_NAME_ACC_list_vector.back().NAME = str_temp_0;
                hmm_NAME_unordered_map[str_temp_0]=unint_name_num;
                unint_name_num++;
            }

            if(str_online_temp.find("ACC") != -1){
                std::string str_temp_0 = str_online_temp;
                str_temp_0 = str_temp_0.substr(str_temp_0.find_first_not_of(' '));
                str_temp_0 = str_temp_0.substr(str_temp_0.find_first_of(' '));
                str_temp_0 = str_temp_0.substr(str_temp_0.find_first_not_of(' '));
                std::vector<std::string> result_split =str_Split_by_char_list(str_temp_0,"  ");
                //str_temp_0 = str_temp_0.substr(str_temp_0.find_first_of(' '));
                //str_temp_0 = str_temp_0.substr(str_temp_0.find_first_not_of(' '), str_temp_0.find_first_of(' ') - str_temp_0.find_first_not_of(' '));

                str_temp_0 =result_split[0];
                if(hmm_ACC_unordered_map.find(str_temp_0) != hmm_ACC_unordered_map.end()){
                    std::cout<<"File '"<<id_hmm_NAME_ACC_list_vector.back().filename <<"' has the same 'ACC' with file '"<<id_hmm_NAME_ACC_list_vector[hmm_ACC_unordered_map[str_temp_0]].filename <<"' !"<<std::endl;
                    id_hmm_NAME_ACC_list_vector.clear();
                    return 0;
                }
                id_hmm_NAME_ACC_list_vector.back().ACC = str_temp_0;
                hmm_ACC_unordered_map[str_temp_0] =unint_acc_num;
                id_hmm_NAME_ACC_list_vector.back().bool_flag_ACC = true;
                unint_acc_num++;
            }
            unint_temp++;
		}
		if_name_acc_temp.close();
		unint_matrix_hmm++;
	}
	vec_matrix_hmm_filenames.clear();

    if(use_ACC){
        if(unint_name_num != unint_acc_num){
            NAME_OR_ACC = false;
            std::cout<<"Ignore the '-acc' option !"<<std::endl;
        }else {
            if(use_ACC){
                NAME_OR_ACC = true;
            }
        }
    }

    //update the unordered map
	for (unsigned int i_vector = 0; i_vector < id_hmm_NAME_ACC_list_vector.size(); i_vector++)
	{

		if(id_hmm_NAME_ACC_list_vector[i_vector].bool_flag_ACC){

            hmm_ACC_id_list_unordered_map[id_hmm_NAME_ACC_list_vector[i_vector].ACC] = i_vector;
		}else{

            hmm_NAME_id_list_unordered_map[id_hmm_NAME_ACC_list_vector[i_vector].NAME] = i_vector;
		}
		//std::cout << hmm_NAME_id_list_unordered_map[id_hmm_NAME_ACC_list_vector[i_vector]] << std::endl;
	}

	//update the matrix

	for (unsigned int i_matrix_row = 0; i_matrix_row < id_hmm_NAME_ACC_list_vector.size(); i_matrix_row++)
	{
		std::vector<double> vector_distance_temp;
		dist_matrix.push_back(vector_distance_temp);
		for (unsigned int i_matrix_col = 0; i_matrix_col < id_hmm_NAME_ACC_list_vector.size(); i_matrix_col++)
		{
			dist_matrix.back().push_back(0.0);
		}
	}

	return 1;
}

//init the matrix 、 vector  and map
int HMMTree::matrix_init_matrix_vecotr_map_hhsuite()
{
	std::vector<std::string> vec_matrix_hmm_filenames;

	if(!get_file_names(folder_hhms,vec_matrix_hmm_filenames,".hhm")){
        return 0;
	}
	if(vec_matrix_hmm_filenames.size() < 3){
        return 0;
	}

	size_t unint_matrix_hmm = 0;
	size_t unint_name_num = 0;
	size_t unint_acc_num = 0;

	while (unint_matrix_hmm < vec_matrix_hmm_filenames.size())
	{
		std::ifstream  if_name_acc_temp;
		if_name_acc_temp.open(folder_hhms+vec_matrix_hmm_filenames[unint_matrix_hmm]);

		if(!if_name_acc_temp.is_open()){
            std::cout<<"matrix_init_matrix_vecotr_map(): Open "<<vec_matrix_hmm_filenames[unint_matrix_hmm]<<"failed !"<<std::endl;
            return 0;
		}
        HMM_NAME_ACC hmm_name_acc_temp;
        id_hmm_NAME_ACC_list_vector.push_back(hmm_name_acc_temp);
        id_hmm_NAME_ACC_list_vector.back().filename=folder_hhms+vec_matrix_hmm_filenames[unint_matrix_hmm];
		std::string str_online_temp = "";
		size_t unint_temp = 0;
		while(unint_temp < 5 && std::getline(if_name_acc_temp,str_online_temp) ){

            if(str_online_temp.find("NAME") != -1){
                std::string str_temp_0 = str_online_temp;
                str_temp_0 = str_temp_0.substr(str_temp_0.find_first_not_of(' '));
                str_temp_0 = str_temp_0.substr(str_temp_0.find_first_of(' '));
                str_temp_0 = str_temp_0.substr(str_temp_0.find_first_not_of(' '));
                std::vector<std::string> result_split =str_Split_by_char_list(str_temp_0,"  ");
                //str_temp_0 = str_temp_0.substr(str_temp_0.find_first_of(' '));
                //str_temp_0 = str_temp_0.substr(str_temp_0.find_first_not_of(' '), str_temp_0.find_first_of(' ') - str_temp_0.find_first_not_of(' '));

                str_temp_0 =result_split[0];
                //char * result_split = strtok(str_temp_0,"   ");
                //std::cout<<result_split[0]<<"-------------------"<<std::endl;
                if(hmm_NAME_unordered_map.find(str_temp_0) != hmm_NAME_unordered_map.end()){
                    std::cout<<"File '"<<id_hmm_NAME_ACC_list_vector.back().filename <<"' has the same 'NAME' with file '"<<id_hmm_NAME_ACC_list_vector[hmm_NAME_unordered_map[str_temp_0]].filename <<"' !"<<std::endl;
                    id_hmm_NAME_ACC_list_vector.clear();
                    return 0;
                }
                id_hmm_NAME_ACC_list_vector.back().NAME = str_temp_0;
                hmm_NAME_unordered_map[str_temp_0]=unint_name_num;
                unint_name_num++;
            }

            if(str_online_temp.find("ACC") != -1){
                std::string str_temp_0 = str_online_temp;
                str_temp_0 = str_temp_0.substr(str_temp_0.find_first_not_of(' '));
                str_temp_0 = str_temp_0.substr(str_temp_0.find_first_of(' '));
                str_temp_0 = str_temp_0.substr(str_temp_0.find_first_not_of(' '));
                std::vector<std::string> result_split =str_Split_by_char_list(str_temp_0,"  ");
                //str_temp_0 = str_temp_0.substr(str_temp_0.find_first_of(' '));
                //str_temp_0 = str_temp_0.substr(str_temp_0.find_first_not_of(' '), str_temp_0.find_first_of(' ') - str_temp_0.find_first_not_of(' '));

                str_temp_0 =result_split[0];
                if(hmm_ACC_unordered_map.find(str_temp_0) != hmm_ACC_unordered_map.end()){
                    std::cout<<"File '"<<id_hmm_NAME_ACC_list_vector.back().filename <<"' has the same 'ACC' with file '"<<id_hmm_NAME_ACC_list_vector[hmm_ACC_unordered_map[str_temp_0]].filename <<"' !"<<std::endl;
                    id_hmm_NAME_ACC_list_vector.clear();
                    return 0;
                }
                id_hmm_NAME_ACC_list_vector.back().ACC = str_temp_0;
                 hmm_ACC_unordered_map[str_temp_0] =unint_acc_num;
                id_hmm_NAME_ACC_list_vector.back().bool_flag_ACC = true;
                unint_acc_num++;
            }
            unint_temp++;
		}
		if_name_acc_temp.close();
		unint_matrix_hmm++;
	}
	vec_matrix_hmm_filenames.clear();

    if(use_ACC){
        if(unint_name_num != unint_acc_num){
            NAME_OR_ACC = false;
            std::cout<<"Ignore the '-acc' option !"<<std::endl;
        }else {
            if(use_ACC){
                NAME_OR_ACC = true;
            }
        }
    }

    //update the unordered map
	for (unsigned int i_vector = 0; i_vector < id_hmm_NAME_ACC_list_vector.size(); i_vector++)
	{

		if(id_hmm_NAME_ACC_list_vector[i_vector].bool_flag_ACC){

            hmm_ACC_id_list_unordered_map[id_hmm_NAME_ACC_list_vector[i_vector].ACC] = i_vector;
		}else{

            hmm_NAME_id_list_unordered_map[id_hmm_NAME_ACC_list_vector[i_vector].NAME] = i_vector;
		}
		//std::cout << hmm_NAME_id_list_unordered_map[id_hmm_NAME_ACC_list_vector[i_vector]] << std::endl;
	}

	//update the matrix

	for (unsigned int i_matrix_row = 0; i_matrix_row < id_hmm_NAME_ACC_list_vector.size(); i_matrix_row++)
	{
		std::vector<double> vector_distance_temp;
		dist_matrix.push_back(vector_distance_temp);
		for (unsigned int i_matrix_col = 0; i_matrix_col < id_hmm_NAME_ACC_list_vector.size(); i_matrix_col++)
		{
			dist_matrix.back().push_back(0.0);
		}
	}

	return 1;
}


//temple test function to compare each 2 hmms
int  HMMTree::matrix_get_each2_hmms_result_2(){

    //temp variable to save the distance
    STUC_RHH_NOTE STUC_RHH_NOTE_distanc_temp;
    //compute the distance  and get the distance
    prc_set_STUC_RHH_NOTE_dist(&STUC_RHH_NOTE_distanc_temp);
    //clear the result class
    clear_Prc_Result();

    int row_number_matrix = 0;
    int col_number_matrix = 0;
    if(hmm_NAME_id_list_unordered_map.find(STUC_RHH_NOTE_distanc_temp.hmm1.name) == hmm_NAME_id_list_unordered_map.end()){
        row_number_matrix = hmm_ACC_id_list_unordered_map[STUC_RHH_NOTE_distanc_temp.hmm1.name];
    }else{
        row_number_matrix = hmm_NAME_id_list_unordered_map[STUC_RHH_NOTE_distanc_temp.hmm1.name];
    }
    if(hmm_NAME_id_list_unordered_map.find(STUC_RHH_NOTE_distanc_temp.hmm2.name) == hmm_NAME_id_list_unordered_map.end()){
        col_number_matrix = hmm_ACC_id_list_unordered_map[STUC_RHH_NOTE_distanc_temp.hmm2.name];
    }else{
        col_number_matrix = hmm_NAME_id_list_unordered_map[STUC_RHH_NOTE_distanc_temp.hmm2.name];
    }
    //put the distance into the matrix
    dist_matrix[row_number_matrix][col_number_matrix] = STUC_RHH_NOTE_distanc_temp.distance;
    dist_matrix[col_number_matrix][row_number_matrix] = STUC_RHH_NOTE_distanc_temp.distance;


	return 1;
}

//
int HMMTree::matrix_get_lib_hmms_result_2(){

    std::vector<int> vec_int;
    for(unsigned int i = 0; i< dist_matrix.size(); i++)
    {
        vec_int.push_back(0);
    }
    unsigned int uint_res_msg = 0;
    while(uint_res_msg < res_msg.size()){
        std::string str_name1 = res_msg[uint_res_msg].hmm1_name;
        std::string str_name2 = res_msg[uint_res_msg].hmm2_name;
        double distance = res_msg[uint_res_msg].simple;
        double answer=0;
        if(distance <= 0){
            distance = DBL_MAX_USER;
        }
        answer = 1.0/distance;

        int row_number_matrix = 0;
        int col_number_matrix = 0;
        if(hmm_NAME_id_list_unordered_map.find(str_name1) == hmm_NAME_id_list_unordered_map.end()){
            row_number_matrix = hmm_ACC_id_list_unordered_map[str_name1];
        }else{
            row_number_matrix = hmm_NAME_id_list_unordered_map[str_name1];
        }
        if(hmm_NAME_id_list_unordered_map.find(str_name2) == hmm_NAME_id_list_unordered_map.end()){
            col_number_matrix = hmm_ACC_id_list_unordered_map[str_name2];
        }else{
            col_number_matrix = hmm_NAME_id_list_unordered_map[str_name2];
        }
        if(vec_int[col_number_matrix] == 0 ){
            //put the distance into the matrix
            dist_matrix[row_number_matrix][col_number_matrix] = answer;
            dist_matrix[col_number_matrix][row_number_matrix] = answer;
            vec_int[col_number_matrix] = 1;
        }else{
            if(dist_matrix[row_number_matrix][col_number_matrix] > answer){
                //put the distance into the matrix
                dist_matrix[row_number_matrix][col_number_matrix] = answer;
                dist_matrix[col_number_matrix][row_number_matrix] = answer;
            }
        }
        uint_res_msg++;
    }
    vec_int.clear();
    //clear the result class
    clear_Prc_Result();
	return 1;
}

//out put the distance matrix
int HMMTree::matrix_mega_out_put_dist_matrix_to_file()
{
	//matrix_out_put_dist_matrix_to_window();
	//temp file stream to output the matrix to a phiply file
	std::ofstream file_dist_matrix_out;

	/*
	char c_int2str[255];                   //temp varibale to save the result of turn unsigned int to str

	//turn unsigned int to str
	_itoa_s(dist_matrix.size(), c_int2str, 10);
	append
	*/

	std::string str_name_temp = folder_matrixs+"file_dist_matrix_out_mega.meg";                 //temp string variable to name the out distance file
	file_dist_matrix_out.open(str_name_temp.c_str(),std::ios_base::out|std::ios_base::trunc);

	if (!file_dist_matrix_out.is_open())
	{
		std::cout << "There is an error when open the matrix output file!" << std::endl;
		return 0;
	}

	//out put the count of the hmms to the file
	//file_dist_matrix_out << std::left << std::setw(8) << dist_matrix.size() << std::endl;
	/*
	#mega
	!Title: xxx--x;
	!Format DataType = Distance DataFormat = LowerLeft NTaxa = 6;

	[1] #Rodent
	[2] #Primate
	[3] #Lagomorpha
	[4] #Artiodactyla
	[5] #Carnivora
	[6] #Perissodactyla

	[  1  2  3  4  5  6 ]
	[1]
	[2]  0.514
	[3]  0.535 0.436
	[4]  0.530 0.388 0.418
	[5]  0.521 0.353 0.417 0.345
	[6]  0.500 0.331 0.402 0.327 0.349
	*/
	/*
	#mega
	!Title: xxx--x;
	!Format DataType = Distance DataFormat = LowerLeft NTaxa = 6;

	*/
	file_dist_matrix_out << std::left << "#mega"<<std::endl;
	file_dist_matrix_out << std::left << "!Title: ;"<< std::endl;
	file_dist_matrix_out << std::left << "!Format DataType = Distance DataFormat = LowerLeft NTaxa = "<< dist_matrix.size()<<";" << std::endl;
	file_dist_matrix_out << std::left << std::endl;

	/*
	[1] #Rodent
	[2] #Primate
	[3] #Lagomorpha
	[4] #Artiodactyla
	[5] #Carnivora
	[6] #Perissodactyla

	*/
	//vector to save the number string "[x]"
    std::vector<std::string> vec_str_tags;
	if(NAME_OR_ACC){
        for (unsigned int i_vec_row_dist_matrix = 0; i_vec_row_dist_matrix < dist_matrix.size(); i_vec_row_dist_matrix++)
        {
            std::string  str_temp = "";
            str_temp = "[" + uint2str(i_vec_row_dist_matrix + 1) + "]";
            vec_str_tags.push_back(str_temp);
            file_dist_matrix_out << std::left << std::setw(6) << str_temp;
            file_dist_matrix_out << std::left <<"#" << std::left << std::setw(16) << id_hmm_NAME_ACC_list_vector[i_vec_row_dist_matrix].ACC << std::endl;
        }
	}else{
        for (unsigned int i_vec_row_dist_matrix = 0; i_vec_row_dist_matrix < dist_matrix.size(); i_vec_row_dist_matrix++)
        {
            std::string  str_temp = "";
            str_temp = "[" + uint2str(i_vec_row_dist_matrix + 1) + "]";
            vec_str_tags.push_back(str_temp);
            file_dist_matrix_out << std::left << std::setw(6) << str_temp;
            file_dist_matrix_out << std::left <<"#" << std::left << std::setw(16) << id_hmm_NAME_ACC_list_vector[i_vec_row_dist_matrix].NAME << std::endl;
        }
	}

	file_dist_matrix_out << std::left << std::endl;
	/*
	[  1  2  3  4  5  6 ]
	*/
	file_dist_matrix_out << std::left << std::setw(6) << "[" ;
	for (unsigned int i_vec_row_dist_matrix = 1; i_vec_row_dist_matrix <= dist_matrix.size(); i_vec_row_dist_matrix++)
	{
		file_dist_matrix_out << std::left << std::setw(18) <<i_vec_row_dist_matrix ;
	}
	file_dist_matrix_out << std::left << std::setw(6) << "]"<< std::endl;
	/*
	[1]
	[2]  0.514
	[3]  0.535 0.436
	[4]  0.530 0.388 0.418
	[5]  0.521 0.353 0.417 0.345
	[6]  0.500 0.331 0.402 0.327 0.349
	*/
	for (unsigned int i_vec_row_dist_matrix = 0; i_vec_row_dist_matrix < dist_matrix.size(); i_vec_row_dist_matrix++)
	{
		file_dist_matrix_out << std::left << std::setw(6) << vec_str_tags[i_vec_row_dist_matrix];
		for (unsigned int i_vec_coll_dist_matrix = 0; i_vec_coll_dist_matrix < i_vec_row_dist_matrix; i_vec_coll_dist_matrix++)
		{
		    //<< std::fixed<<std::setprecision(10)
			file_dist_matrix_out << std::left << std::setw(18) << (dist_matrix[i_vec_row_dist_matrix][i_vec_coll_dist_matrix]+ dist_matrix[i_vec_coll_dist_matrix][i_vec_row_dist_matrix]) /2.0;
		}
		file_dist_matrix_out << std::endl;
	}
	vec_str_tags.clear();
	file_dist_matrix_out.close();
	chmod((folder_matrixs+"file_dist_matrix_out_mega.meg").c_str(), S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IWOTH | S_IROTH);
	return 1;
}

//out put the distance matrix
int HMMTree::matrix_phylip_out_put_dist_matrix_to_file()
{
	if(!matrix_check_matrix_format()){
		std::cout<<"Format error!!!"<<std::endl;
	}
	//matrix_out_put_dist_matrix_to_window();
	//temp file stream to output the matrix to a phiply file
	std::ofstream file_dist_matrix_out;

	/*
	char c_int2str[255];                   //temp varibale to save the result of turn unsigned int to str

	//turn unsigned int to str
	_itoa_s(dist_matrix.size(), c_int2str, 10);
	append
	*/

	std::string str_name_temp = folder_matrixs+"file_dist_matrix_out_phylip.txt";                 //temp string variable to name the out distance file
	file_dist_matrix_out.open(str_name_temp.c_str(),std::ios_base::out|std::ios_base::trunc);
	if (!file_dist_matrix_out.is_open())
	{
		std::cout << "There is an error when open the matrix output file!" << std::endl;
		return 0;
	}

	//out put the count of the hmms to the file
	file_dist_matrix_out << std::left << std::setw(8) << dist_matrix.size() << std::endl;

    std::ofstream of_shorted_names;
    of_shorted_names.open(folder_matrixs+"shorted_names.txt",std::ios_base::out|std::ios_base::trunc);
   // of_shorted_names.
   // chmod("./matrixs/shorted_names.txt",S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP);

    if (!of_shorted_names.is_open())
	{
		std::cout << "There is an error when open the 'shorted_names.txt' file !" << std::endl;
		return 0;
	}
	bool shorted_name_flag=false;
    //a unordered map to save the shorted names in phylip draw tree step
    std::unordered_map <std::string, int> shorted_names_map;
    int all_names_num = 0;
    int shorted_num = 0;
    if(NAME_OR_ACC){
        for (unsigned int i_vec_row_dist_matrix = 0; i_vec_row_dist_matrix < dist_matrix.size(); i_vec_row_dist_matrix++)
        {
            std::string str_hmms_names=id_hmm_NAME_ACC_list_vector[i_vec_row_dist_matrix].ACC;
            if(str_hmms_names.length() > 10){
                vec_unshorted_names.push_back(str_hmms_names);
                shorted_name_flag = true;
                of_shorted_names<<str_hmms_names<<"==";
                
                // Create base shortened name: first 3 + last 7 characters
                std::string base_short_name = str_hmms_names.substr(0,3) + str_hmms_names.substr(str_hmms_names.length()-7,7);
                str_hmms_names = base_short_name;
                
                // Handle collisions by appending a unique suffix
                int collision_counter = 1;
                while(shorted_names_map.find(str_hmms_names) != shorted_names_map.end() ){
                    std::string suffix = "_" + int_2_string(collision_counter);
                    
                    // If the suffix would make the name too long, truncate the base name
                    if(base_short_name.length() + suffix.length() > 10) {
                        int available_length = 10 - suffix.length();
                        if(available_length < 3) {
                            // If suffix too long, use shorter suffix format
                            suffix = int_2_string(collision_counter);
                            available_length = 10 - suffix.length();
                        }
                        str_hmms_names = base_short_name.substr(0, available_length) + suffix;
                    } else {
                        str_hmms_names = base_short_name + suffix;
                    }
                    
                    collision_counter++;
                    std::cout<<"Collision resolved: " << str_hmms_names << std::endl;
                }
                of_shorted_names<<str_hmms_names<<std::endl;
                vec_shorted_names.push_back(str_hmms_names);
                shorted_num++;
            }
            shorted_names_map[str_hmms_names] = all_names_num;
            all_names_num++;
            file_dist_matrix_out << std::left << std::setw(16) << str_hmms_names.c_str() << "    ";

            for (unsigned int i_vec_coll_dist_matrix = 0; i_vec_coll_dist_matrix < dist_matrix.size(); i_vec_coll_dist_matrix++)
            {

                file_dist_matrix_out << std::left << std::setw(18) <<dist_matrix[i_vec_row_dist_matrix][i_vec_coll_dist_matrix];

            }

            file_dist_matrix_out << std::endl;
        }
    }else{
        for (unsigned int i_vec_row_dist_matrix = 0; i_vec_row_dist_matrix < dist_matrix.size(); i_vec_row_dist_matrix++)
        {
            std::string str_hmms_names=id_hmm_NAME_ACC_list_vector[i_vec_row_dist_matrix].NAME;
            if(str_hmms_names.length() > 10){
                vec_unshorted_names.push_back(str_hmms_names);
                shorted_name_flag = true;
                of_shorted_names<<str_hmms_names<<"==";
                
                // Create base shortened name: first 3 + last 7 characters
                std::string base_short_name = str_hmms_names.substr(0,3) + str_hmms_names.substr(str_hmms_names.length()-7,7);
                str_hmms_names = base_short_name;
                
                // Handle collisions by appending a unique suffix
                int collision_counter = 1;
                while(shorted_names_map.find(str_hmms_names) != shorted_names_map.end() ){
                    std::string suffix = "_" + int_2_string(collision_counter);
                    
                    // If the suffix would make the name too long, truncate the base name
                    if(base_short_name.length() + suffix.length() > 10) {
                        int available_length = 10 - suffix.length();
                        if(available_length < 3) {
                            // If suffix too long, use shorter suffix format
                            suffix = int_2_string(collision_counter);
                            available_length = 10 - suffix.length();
                        }
                        str_hmms_names = base_short_name.substr(0, available_length) + suffix;
                    } else {
                        str_hmms_names = base_short_name + suffix;
                    }
                    
                    collision_counter++;
                    std::cout<<"Collision resolved: " << str_hmms_names << std::endl;
                }
                of_shorted_names<<str_hmms_names<<std::endl;
                vec_shorted_names.push_back(str_hmms_names);
                shorted_num++;
            }
            shorted_names_map[str_hmms_names] = all_names_num;
            all_names_num++;
            file_dist_matrix_out << std::left << std::setw(16) << str_hmms_names.c_str() << "    ";

            for (unsigned int i_vec_coll_dist_matrix = 0; i_vec_coll_dist_matrix < dist_matrix.size(); i_vec_coll_dist_matrix++)
            {

                file_dist_matrix_out << std::left << std::setw(18) <<dist_matrix[i_vec_row_dist_matrix][i_vec_coll_dist_matrix];

            }

            file_dist_matrix_out << std::endl;
        }
    }
    of_shorted_names.close();
	file_dist_matrix_out.close();
	
	// Set file permissions after closing
	chmod((folder_matrixs+"file_dist_matrix_out_phylip.txt").c_str(), S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IWOTH | S_IROTH);
	if(shorted_name_flag) {
		chmod((folder_matrixs+"shorted_names.txt").c_str(), S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IWOTH | S_IROTH);
	}
	
	if(!shorted_name_flag){
        system_return(system(("rm -f "+folder_matrixs+"shorted_names.txt").c_str()));
	}
	return 1;
}

//matrix_out_put_dist_matrix_to_window
int HMMTree::matrix_out_put_dist_matrix_to_window()
{

	for (unsigned int i_vec_row_dist_matrix = 0; i_vec_row_dist_matrix < dist_matrix.size(); i_vec_row_dist_matrix++)
	{

		for (unsigned int i_vec_coll_dist_matrix = 0; i_vec_coll_dist_matrix < dist_matrix[i_vec_row_dist_matrix].size(); i_vec_coll_dist_matrix++)
		{

			std::cout << dist_matrix[i_vec_row_dist_matrix][i_vec_coll_dist_matrix] << "  ";
		}

		std::cout << std::endl;
	}


	return 1;

}


bool HMMTree::matrix_check_matrix_format(){
	for (unsigned int i_vec_row_dist_matrix = 0; i_vec_row_dist_matrix < dist_matrix.size(); i_vec_row_dist_matrix++)
	{
		for (unsigned int i_vec_coll_dist_matrix = 0; i_vec_coll_dist_matrix <= i_vec_row_dist_matrix; i_vec_coll_dist_matrix++)
		{
			if(dist_matrix[i_vec_row_dist_matrix][i_vec_coll_dist_matrix] != dist_matrix[i_vec_coll_dist_matrix][i_vec_row_dist_matrix])
			{
				return false;
			}

		}
	}

	return true;
}



//double compute average distance score
double HMMTree::matrix_compute_average_dist(){
	double result=0;
	for (unsigned int i_vec_row_dist_matrix = 0; i_vec_row_dist_matrix < dist_matrix.size(); i_vec_row_dist_matrix++)
	{
		for (unsigned int i_vec_coll_dist_matrix = 0; i_vec_coll_dist_matrix < i_vec_row_dist_matrix; i_vec_coll_dist_matrix++)
		{
			result += dist_matrix[i_vec_row_dist_matrix][i_vec_coll_dist_matrix];

		}
	}

	return result*2/(dist_matrix.size()*dist_matrix.size());
}

//
int HMMTree::matrix_get_each2_hhms_result_hhsuite(){
    //temp variable to save the distance
    STUC_RHH_NOTE STUC_RHH_NOTE_distanc_temp;
    //compute the distance  and get the distance
    prc_set_STUC_RHH_NOTE_dist(&STUC_RHH_NOTE_distanc_temp);
    //std::cout<<"555555555555555555555555555555"<<std::endl;
    //clear the result class
    clear_Prc_Result();
    //std::cout<<STUC_RHH_NOTE_distanc_temp.distance<<"   "<<std::endl;
    //std::cout<<STUC_RHH_NOTE_distanc_temp.hmm1.name<<"   "<<STUC_RHH_NOTE_distanc_temp.hmm2.name<<"  :   "<<STUC_RHH_NOTE_distanc_temp.distance<<std::endl;

    int row_number_matrix = 0;
    int col_number_matrix = 0;
    if(hmm_NAME_id_list_unordered_map.find(STUC_RHH_NOTE_distanc_temp.hmm1.name) == hmm_NAME_id_list_unordered_map.end()){
        row_number_matrix = hmm_ACC_id_list_unordered_map[STUC_RHH_NOTE_distanc_temp.hmm1.name];
    }else{
        row_number_matrix = hmm_NAME_id_list_unordered_map[STUC_RHH_NOTE_distanc_temp.hmm1.name];
    }
    if(hmm_NAME_id_list_unordered_map.find(STUC_RHH_NOTE_distanc_temp.hmm2.name) == hmm_NAME_id_list_unordered_map.end()){
        col_number_matrix = hmm_ACC_id_list_unordered_map[STUC_RHH_NOTE_distanc_temp.hmm2.name];
    }else{
        col_number_matrix = hmm_NAME_id_list_unordered_map[STUC_RHH_NOTE_distanc_temp.hmm2.name];
    }
    //std::cout<<"7777777777777777777777777777777"<<std::endl;
    //put the distance into the matrix
    dist_matrix[row_number_matrix][col_number_matrix] = STUC_RHH_NOTE_distanc_temp.distance;
    dist_matrix[col_number_matrix][row_number_matrix] = STUC_RHH_NOTE_distanc_temp.distance;
    //std::cout<<"88888888888888888888888888888888"<<std::endl;

	return 1;
}
