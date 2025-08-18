#include "HMMTree.h"
#include <string.h>
#include <sstream>
#include <iomanip>
#include <cerrno>
#include <sys/stat.h>
#include <unistd.h>

// Test whether the input string is numeric.
int is_num_str(char *char_str_num){

    std::string str_temp = char_str_num;
    std::cout<<str_temp.length()<<std::endl;
    if(str_temp.length() <= 0){
        return 0;
    }
    if(str_temp[0] < '0' || str_temp[0] > '9'){
        return 0;
    }
    unsigned int uint_i = 1;
    bool point_flag = false;
    while(uint_i < str_temp.length()-1){
        if(str_temp[uint_i] > '0' && str_temp[uint_i] < '9'){
            uint_i++;
        }else if(str_temp[uint_i] == '.' && !point_flag){
            uint_i++;
            point_flag = true;
        }else{
            return 0;
        }
    }
    if(uint_i == str_temp.length()-1 && (str_temp[uint_i] < '0' || str_temp[uint_i] > '9')){
        return 0;
    }
    return 1;
}

// Split a string using any character from the given list.
std::vector<std::string> str_Split_by_char_list(std::string str, const char *pattern)
{
	// Build a string containing the delimiter characters.
	std::string str_pat = "";

	for (unsigned int i_pat = 0; i_pat < strlen(pattern); i_pat++)
	{
		str_pat += pattern[i_pat];
	}

	// Result vector.
	std::vector<std::string> result;

	// Temporary buffer: original string plus pattern to simplify loop termination.
	std::string temp = str + pattern; // buffer combining original string with pattern to simplify splitting logic

	// Iterator over the buffer string.
	std::string::iterator itr_str = temp.begin();

	// Current token buffer.
	std::string str_temp = "";

	// Loop to split the string.
    while (itr_str != temp.end())
	{
		// If current char is a delimiter, end the current token.
        if (str_pat.find(*itr_str) != std::string::npos)
		{
			// If the current token is non-empty, push it to the result vector.
			if (str_temp.length() > 0)
			{
				result.push_back(str_temp);
				str_temp = "";
			}
		}
		else
		{// If not a delimiter, append it to the current token.
			str_temp += *itr_str;
		}
		// Advance iterator.
		itr_str++;
	}

	return result;
}

// Print the string if non-empty.
void output_str(std::string  str)
{
	if ((int)str.length() > 0)
	{
		std::cout << str << std::endl;
	}
}

// Convert unsigned int to a decimal string.
std::string uint2str(unsigned int num)
{
	std::string result = "";
	if (num == 0)
	{
		result.insert(0, 1, '0');
		return result;
	}

	unsigned int num_temp = num;
	while (num_temp > 0)
	{
		int x = 0;
		x = num_temp % 10;
		result.insert(0, 1, (x + '0'));
		num_temp = num_temp / 10;
	}
	return result;
}

// Check whether a directory exists.
bool dir_exist_opendir(std::string path)
{
	DIR *dirptr = NULL;
	if(path.length() == 0)
	{
		return false;
	}

	if((dirptr=opendir(path.c_str())) == NULL)
	{
		return false;
	}

	closedir(dirptr);
	return true;
}

// Check whether a directory is non-empty.
bool dir_noempty_opendir_readir(std::string path)
{
    if (path.empty()) {
        return false;
    }
    DIR* dirptr = opendir(path.c_str());
    if (!dirptr) {
        return false;
    }
    struct dirent* p_dirent;
    int file_num = 0;
    while ((p_dirent = readdir(dirptr)) != NULL) {
        file_num++;
    }
    closedir(dirptr);
    return file_num > 2; // '.' and '..' only -> empty
}

// Check that a file exists and is not empty.
bool file_exists_and_empty_check(std::string str_path_name)
{
	std::ifstream ifstream_empty_test;
	ifstream_empty_test.open(str_path_name.c_str());
	if(!ifstream_empty_test.is_open())
	{
		std::cout<<"Failed to open file: "<<str_path_name<<std::endl;
		return false;
	}
    if (ifstream_empty_test.peek() == EOF)
	{
		ifstream_empty_test.close();
		std::cout<<"File '"<<str_path_name<<"' is empty!"<<std::endl;
		return false;
	}
	ifstream_empty_test.close();
	return true;
}

// Check whether a file exists.
bool file_exists_als_phmms_phhms(std::string str_path_name)
{
	std::ifstream ifstream_empty_test;
	ifstream_empty_test.open(str_path_name.c_str());
	if(!ifstream_empty_test.is_open())
	{
		return false;
	}
	ifstream_empty_test.close();
	return true;
}

// Convert int to string.
std::string int_2_string(int int_number) {
	char a[10];
	std::string str;

	sprintf(a,"%d",int_number);
	str = a;
	return str;
}


std::string double_2_string(double double_x){
    char a[40];
	std::string str;

	sprintf(a,"%.6f",double_x);
	str = a;
	return str;
}

int system_return(int status){
    int result = 0;
    if (-1 == status)
    {
        std::cerr << "system() failed: errno=" << errno << std::endl;
        result = 1;
    }
    else
    {
        // printf("exit status value = [0x%x]\n", status);

        if (WIFEXITED(status))  // Process exited normally.
        {
            if (0 == WEXITSTATUS(status)) // Success.
            {
                // printf("run shell script successfully.\n");
            }
            else  // Non-zero exit status.
            {
                std::cout << "run shell command failed, exit code: " << WEXITSTATUS(status) << std::endl;
                result = 1;
            }
        }
        else  // Process did not exit normally.
        {
            std::cout << "process did not exit normally (status=" << status << ")" << std::endl;
            result = 1;
        }
    }
    return result;
}

// List file names in a directory, optionally filtering by extension.
int get_file_names(std::string path, std::vector<std::string> & vec_file_names, std::string extention){
    if(!dir_exist_opendir(path)){
        return 0;
    }
    bool ext_filter = false;
    if(extention.length() > 0){
        ext_filter = true;
    }
    DIR * dir;
    struct dirent * ptr;
    // removed unused variable 'i'
    dir = opendir(path.c_str()); //open a dir
    while((ptr = readdir(dir)) != NULL) // Iterate over directory entries.
    {
        if(ptr->d_type == 8){
            if(!ext_filter){
                vec_file_names.push_back(ptr->d_name);
            }else{
                std::string str_temp = ptr->d_name;
                if(!strcmp((str_temp.substr(str_temp.length()-extention.length(),extention.length())).c_str(), extention.c_str())){
                    vec_file_names.push_back(str_temp);
                }
            }
        }
    }
    closedir(dir);//
    return 1;
}

// Delete all files in a directory.
int delete_files(std::string path){
    std::vector<std::string> vec_file_names;
    get_file_names(path,vec_file_names,"");
    size_t uint_files_num = 0;
    while(uint_files_num < vec_file_names.size()){
        std::string path_temp=path + vec_file_names[uint_files_num];
        unlink(path_temp.c_str());
        uint_files_num++;
    }
    vec_file_names.clear();
    return 1;
}

// Move files from path1 to path2.
int mv_files(std::string path1, std::string path2, std::string extention){
    std::vector <std::string> vec_filenames_path1;
    get_file_names(path1,vec_filenames_path1,extention);
    size_t uint_file_num = 0;
    while(uint_file_num < vec_filenames_path1.size()){
        std::string path1name = path1+vec_filenames_path1[uint_file_num];
        std::string path2name = path2+vec_filenames_path1[uint_file_num];
        if(rename(path1name.c_str(),path2name.c_str())){
            std::cout<<"mv_files(): error!"<<std::endl;
            return 0;
        }
        uint_file_num++;
    }
    return 1;
}

// Copy files from path1 to path2.
int copy_files(std::string path1, std::string path2,std::string extention){

    std::vector <std::string> vec_names_path1;
    get_file_names(path1,vec_names_path1,extention);
    if(vec_names_path1.size() < 3 && (!HMMTree::als_phmms_phhms)){
        vec_names_path1.clear();
        output_error_("'copy_files()' ");
        return 0;
    }
    size_t uint_file_num = 0;
    while(uint_file_num < vec_names_path1.size()){
        std::string path1name = path1 + vec_names_path1[uint_file_num];
        std::string path2name = path2 + vec_names_path1[uint_file_num];
        if(file_exists_als_phmms_phhms(path2name)){
            std::string path_name2=path2name.substr(0,path2name.find_last_of("."));
            std::string path_extention=path2name.substr(path2name.find_last_of("."),path2name.length()-path2name.find_last_of("."));
            path2name=path_name2+"_"+path_extention;
        }
        copy_one_file(path1name,path2name);
        uint_file_num++;
    }
    vec_names_path1.clear();
    return 1;
}

// Copy a single file from pathname1 to pathname2.
int copy_one_file(std::string pathname1, std::string pathname2){
    int from_fd,to_fd;
    int bytes_read,bytes_write;
    char buffer[BUFFER_SIZE];
    char *ptr;

    /* Open the source file. */
    if((from_fd=open(pathname1.c_str(),O_RDONLY))==-1) /* open file read-only; return -1 on error */
    {
        fprintf(stderr,"Open %s Error:%s\n",pathname1.c_str(),strerror(errno));
        exit(1);
    }

    /* Create the target file. */
    /* When using O_CREAT, open() needs a third argument;
       mode S_IRUSR|S_IWUSR means user read/write. */
    if((to_fd=open(pathname2.c_str(),O_WRONLY|O_CREAT,S_IRUSR|S_IWUSR))==-1)
    {
        fprintf(stderr,"Open %s Error:%s\n",pathname2.c_str(),strerror(errno));
        exit(1);
    }

    /* Copy the file. */
    while ((bytes_read = read(from_fd, buffer, BUFFER_SIZE)))
    {
        /* Handle fatal error. */
        if((bytes_read==-1)&&(errno!=EINTR))
            break;
        else if(bytes_read>0)
        {
            ptr=buffer;
            while ((bytes_write = write(to_fd, ptr, bytes_read)))
            {
                /* Handle fatal error. */
                if((bytes_write==-1)&&(errno!=EINTR))
                    break;
                /* All bytes copied. */
                else if(bytes_write==bytes_read)
                    break;
                /* If partially copied, continue. */
                else if(bytes_write>0)
                {
                    ptr+=bytes_write;
                    bytes_read-=bytes_write;
                }
            }
            /* Handle fatal error while writing. */
            if(bytes_write==-1)
            break;
        }
    }
    close(from_fd);
    close(to_fd);
    return 1;
}

void output_error_(std::string error_msg){
    std::cout << error_msg << " error!" << std::endl;
    // Also print the full usage to help the user immediately
    std::cout << STR_ARGUMENTS_ERROR_MSG << std::endl;
    exit(1);
}

// Check if USEARCH exists.
int USEARCH_exist(){
    bool usr_bin_path = true;
    char *buf = NULL;
    buf = (char*)malloc(sizeof(char) * MAX_BUF_LEN);
    FILE *file = popen("which usearch", "r");
    memset(buf, 0, MAX_BUF_LEN);
    if(fgets(buf, MAX_BUF_LEN, file) == NULL){
        usr_bin_path = false;
    }
    pclose(file);
    free(buf);
    if(!usr_bin_path){
        if(access("usearch",F_OK) == 0){
            if(access("usearch",X_OK) == 0){
                return 1;
            }
        }
    }else{
        return -1;
    }

    return 0;
}

// Check if PRC exists.
int PRC_exist(){
    bool usr_bin_path = true;
    char *buf = NULL;
    buf = (char*)malloc(sizeof(char) * MAX_BUF_LEN);
    FILE *file = popen("which prc", "r");
    memset(buf, 0, MAX_BUF_LEN);
    if(fgets(buf, MAX_BUF_LEN, file) == NULL){
        usr_bin_path = false;
    }
    pclose(file);
    free(buf);
    if(!usr_bin_path){
        if(access("prc",F_OK) == 0){
            if(access("prc",X_OK) == 0){
                return 1;
            }
        }
    }else{
        return -1;
    }

    return 0;
}

// Check if PRC-X exists.
int PRCX_exist(){
    bool usr_bin_path = true;
    char *buf = NULL;
    buf = (char*)malloc(sizeof(char) * MAX_BUF_LEN);
    FILE *file = popen("which prcX", "r");
    memset(buf, 0, MAX_BUF_LEN);
    if(fgets(buf, MAX_BUF_LEN, file) == NULL){
        usr_bin_path = false;
    }
    pclose(file);
    free(buf);
    if(!usr_bin_path){
        if(access("prcX",F_OK) == 0){
            if(access("prcX",X_OK) == 0){
                return 1;
            }
        }
    }else{
        return -1;
    }

    return 0;
}

// Check if MAFFT exists.
int MAFFT_exist(){
    bool usr_bin_path = true;
    char *buf = NULL;
    buf = (char*)malloc(sizeof(char) * MAX_BUF_LEN);
    FILE *file = popen("which mafft", "r");
    memset(buf, 0, MAX_BUF_LEN);
    if(fgets(buf, MAX_BUF_LEN, file) == NULL){
        usr_bin_path = false;
    }
    pclose(file);
    free(buf);
    if(!usr_bin_path){
        if(access("mafft",F_OK) == 0){
            if(access("mafft",X_OK) == 0){
                return 1;
            }
        }
    }else{
        return -1;
    }

    return 0;
}

// Check if HMMER hmmbuild exists.
int HMMER_hmmbuild_exist(){
    bool usr_bin_path = true;
    char *buf = NULL;
    buf = (char*)malloc(sizeof(char) * MAX_BUF_LEN);
    FILE *file = popen("which hmmbuild", "r");
    memset(buf, 0, MAX_BUF_LEN);
    if(fgets(buf, MAX_BUF_LEN, file) == NULL){
        usr_bin_path = false;
    }
    pclose(file);
    free(buf);
    if(!usr_bin_path){
        if(access("hmmbuild",F_OK) == 0){
            if(access("hmmbuild",X_OK) == 0){
                return 1;
            }
        }
    }else{
        return -1;
    }

    return 0;

}

// Check if HMMER hmmconvert exists.
int HMMER_hmmconvert_exist(){
    bool usr_bin_path = true;
    char *buf = NULL;
    buf = (char*)malloc(sizeof(char) * MAX_BUF_LEN);
    FILE *file = popen("which hmmconvert", "r");
    memset(buf, 0, MAX_BUF_LEN);
    if(fgets(buf, MAX_BUF_LEN, file) == NULL){
        usr_bin_path = false;
    }
    pclose(file);
    free(buf);
    if(!usr_bin_path){
        if(access("hmmconvert",F_OK) == 0){
            if(access("hmmconvert",X_OK) == 0){
                return 1;
            }
        }
    }else{
        return -1;
    }
    return 0;
}

// Check if HH-suite hhmake exists.
int hhsuite_hhmake_exist(){
    bool usr_bin_path = true;
    char *buf = NULL;
    buf = (char*)malloc(sizeof(char) * MAX_BUF_LEN);
    FILE *file = popen("which hhmake", "r");
    memset(buf, 0, MAX_BUF_LEN);
    if(fgets(buf, MAX_BUF_LEN, file) == NULL){
        usr_bin_path = false;
    }
    pclose(file);
    free(buf);
    if(!usr_bin_path){
        if(access("hhmake",F_OK) == 0){
            if(access("hhmake",X_OK) == 0){
                return 1;
            }
        }
    }else{
        return -1;
    }

    return 0;

}

// Check if HH-suite hhalign exists.
int hhsuite_hhalign_exist(){
    bool usr_bin_path = true;
    char *buf = NULL;
    buf = (char*)malloc(sizeof(char) * MAX_BUF_LEN);
    FILE *file = popen("which hhalign", "r");
    memset(buf, 0, MAX_BUF_LEN);
    if(fgets(buf, MAX_BUF_LEN, file) == NULL){
        usr_bin_path = false;
    }
    pclose(file);
    free(buf);
    if(!usr_bin_path){
        if(access("hhalign",F_OK) == 0){
            if(access("hhalign",X_OK) == 0){
                return 1;
            }
        }
    }else{
        return -1;
    }

    return 0;
}

// Replace occurrences of s2 with s3 in s1, considering un_shorted as protected prefixes.
void string_replace(std::string &s1,const std::string&s2,const std::string&s3,unsigned int index_i, std::vector<std::string> un_shorted)
{	/*
	std::string::size_type pos=0;
	std::string::size_type a=s2.size();
	std::string::size_type b=s3.size();
	while((pos=s1.find(s2,pos))!=std::string::npos)
	{
		s1.replace(pos,a,s3);
		pos+=b;
	}
	return;
*/
	std::string::size_type pos=0;
    // removed unused variable 'pos1'
	std::string::size_type a=s2.size();
	std::string::size_type b=s3.size();

	while((pos=s1.find(s2,pos))!=std::string::npos)
	{
         bool flag = false;
         unsigned int max_len = 0;
         unsigned int max_index= 0;
	    std::string s4 = s1.substr(pos,s1.length()-pos);
	    unsigned int index_num = 0;
        while(index_num< un_shorted.size()){
            if(index_num != index_i){
                if(s4.find(un_shorted[index_num]) == 0){
                    flag = true;
                    if(un_shorted[index_num].length() > max_len){
                       max_len = un_shorted[index_num].length();
                       max_index = index_num;
                    }
                }
            }
            index_num++;
        }
        if(flag){
            pos+=un_shorted[max_index].size();
        }else{
             s1.replace(pos,a,s3);
            pos+=b;
        }

	}
	return;

}

// Format a time duration as a human-readable string.
std::string format_time_duration(long total_milliseconds) {
    // Convert milliseconds to seconds
    double total_seconds = total_milliseconds / 1000.0;
    
    // If less than 60 seconds, show seconds with decimal places
    if (total_seconds < 60.0) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(3) << total_seconds << " seconds";
        return oss.str();
    }
    
    // Convert to integer seconds for longer durations
    long seconds = static_cast<long>(total_seconds);
    
    // Calculate days, hours, minutes, seconds
    long days = seconds / (24 * 3600);
    seconds %= (24 * 3600);
    long hours = seconds / 3600;
    seconds %= 3600;
    long minutes = seconds / 60;
    seconds %= 60;
    
    std::ostringstream oss;
    bool first = true;
    
    if (days > 0) {
        oss << days << " day" << (days > 1 ? "s" : "");
        first = false;
    }
    
    if (hours > 0) {
        if (!first) oss << ", ";
        oss << hours << " hour" << (hours > 1 ? "s" : "");
        first = false;
    }
    
    if (minutes > 0) {
        if (!first) oss << ", ";
        oss << minutes << " minute" << (minutes > 1 ? "s" : "");
        first = false;
    }
    
    if (seconds > 0 || first) {  // Always show seconds if it's the only unit
        if (!first) oss << ", ";
        oss << seconds << " second" << (seconds > 1 ? "s" : "");
    }
    
    return oss.str();
}
