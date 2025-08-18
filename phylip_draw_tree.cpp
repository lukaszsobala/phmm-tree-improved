#include "HMMTree.h"
#include <omp.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <algorithm>

int HMMTree::Phylip_draw_tree2(){
     kitsch_build_tree((folder_matrices+"file_dist_matrix_out_phylip.txt").c_str(),(folder_tree_files+"f-m").c_str(), 0, phylo_threads_count);
     kitsch_build_tree((folder_matrices+"file_dist_matrix_out_phylip.txt").c_str(),(folder_tree_files+"min").c_str(), 1, phylo_threads_count);
     fitch_build_tree((folder_matrices+"file_dist_matrix_out_phylip.txt").c_str(),(folder_tree_files+"f-m").c_str(), 0, phylo_threads_count);
     fitch_build_tree((folder_matrices+"file_dist_matrix_out_phylip.txt").c_str(),(folder_tree_files+"min").c_str(), 1, phylo_threads_count);
     neighbor_build_tree((folder_matrices+"file_dist_matrix_out_phylip.txt").c_str(),(folder_tree_files+"neighbor").c_str(), phylo_threads_count);
     upgma_build_tree((folder_matrices+"file_dist_matrix_out_phylip.txt").c_str(),(folder_tree_files+"upgma").c_str(), phylo_threads_count);
     //UPGMA_build_tree("ssssss");
return 1;
}

//function to run selected phylogenetic analyses
void HMMTree::draw_tree_selective(bool run_fitch, bool run_kitsch, bool run_upgma, bool run_nj, bool fm_only, bool min_only) {
    std::string matrix_file = folder_matrices + "file_dist_matrix_out_phylip.txt";
    
    // Determine concurrency settings
    int concurrent = phylo_concurrent_threads_count;
    if (concurrent <= 0) {
        // Auto-detect similar to PRC: use available cores, but cap at number of tasks later
        #ifdef OPENMP_ENABLED
        concurrent = omp_get_max_threads();
        #else
        concurrent = 1;
        #endif
    }

    // Build task list (each as an external worker call) to avoid static-state races in PHYLIP
    struct Task { std::string algo; std::string matrix; std::string output; int threads; };
    std::vector<Task> tasks;

    auto queue_task = [&](const std::string& algo, const std::string& suffix){
        tasks.push_back(Task{algo, matrix_file, folder_tree_files + suffix, phylo_threads_count <= 0 ? 1 : phylo_threads_count});
    };

    if (run_kitsch) {
        std::cout << "Running Kitsch analysis (contemporary tips method)..." << std::endl;
        // Timing moved to worker completion reporting
        
        // Determine which variants to run
        bool run_fm = !min_only;    // Run f-m unless min_only is specified
        bool run_min = !fm_only;    // Run min unless fm_only is specified
        
    if (run_fm) queue_task("kitsch_fm", "kitsch_f-m");
    if (run_min) queue_task("kitsch_min", "kitsch_min");
        // Completion time is reported per worker when it finishes
    }
    
    if (run_fitch) {
        std::cout << "Running Fitch-Margoliash analysis..." << std::endl;
        // Timing moved to worker completion reporting
        
        // Determine which variants to run
        bool run_fm = !min_only;    // Run f-m unless min_only is specified
        bool run_min = !fm_only;    // Run min unless fm_only is specified
        
    if (run_fm) queue_task("fitch_fm", "fitch_f-m");
    if (run_min) queue_task("fitch_min", "fitch_min");
        // Completion time is reported per worker when it finishes
    }
    
    if (run_nj) {
        std::cout << "Running Neighbor-Joining analysis..." << std::endl;
        // Timing moved to worker completion reporting
        queue_task("nj", "neighbor");
        // Completion time is reported per worker when it finishes
    }
    
    if (run_upgma) {
        std::cout << "Running UPGMA analysis..." << std::endl;
        // Timing moved to worker completion reporting
        queue_task("upgma", "upgma");
        // Completion time is reported per worker when it finishes
    }

    // Nothing to do?
    if (tasks.empty()) return;

    // Cap concurrency to number of tasks
    if (concurrent > (int)tasks.size()) concurrent = (int)tasks.size();

    std::cout << "Executing " << tasks.size() << " phylogenetic analyses with up to " << concurrent << " concurrent worker(s)..." << std::endl;

    // Run tasks using a simple worker-pool with system() calls (process isolation)
    size_t next = 0; size_t done = 0; int active = 0;
    std::vector<pid_t> pids(tasks.size(), -1);
    std::vector<long long> start_times(tasks.size());
    std::vector<int> results(tasks.size(), -1);

    auto pretty_name = [&](const std::string& algo) -> std::string {
        if (algo == "nj") return "Neighbor-joining";
        if (algo == "upgma") return "UPGMA";
        if (algo == "fitch_fm") return "Fitch-Margoliash (f-m)";
        if (algo == "fitch_min") return "Fitch-Margoliash (min)";
        if (algo == "kitsch_fm") return "Kitsch (f-m)";
        if (algo == "kitsch_min") return "Kitsch (min)";
        return algo;
    };

    auto spawn = [&](size_t idx){
        const auto& t = tasks[idx];
        pid_t pid = fork();
        if (pid == 0) {
            // Child: exec self with -phylo_worker
            std::string threads_str = std::to_string(std::max(1, t.threads));
            execlp("./phmm-tree", "phmm-tree", "-phylo_worker", t.algo.c_str(), t.matrix.c_str(), t.output.c_str(), threads_str.c_str(), (char*)NULL);
            // If exec fails
            _exit(127);
        } else if (pid > 0) {
            pids[idx] = pid;
            start_times[idx] = now_millis();
            active++;
        } else {
            std::cerr << "Failed to fork for task " << idx << std::endl;
            results[idx] = 127;
        }
    };

    // Start initial batch
    while (active < concurrent && next < tasks.size()) spawn(next++);

    // Reap and keep spawning
    while (done < tasks.size()) {
        int status = 0; pid_t pid = wait(&status);
        if (pid > 0) {
            // Find which task
            for (size_t i = 0; i < pids.size(); ++i) if (pids[i] == pid) {
                results[i] = WIFEXITED(status) ? WEXITSTATUS(status) : 128;
                active--; done++;
                // Report accurate per-task runtime upon completion
                long diff_ms = static_cast<long>(now_millis() - start_times[i]);
                std::cout << '\r' << pretty_name(tasks[i].algo) << " completed in: " << format_time_duration(diff_ms) << std::endl;
                break;
            }
        }
        while (active < concurrent && next < tasks.size()) spawn(next++);
    }

    // Summary
    int ok = 0; for (int r : results) if (r == 0) ok++;
    std::cout << "Parallel execution completed: " << ok << "/" << tasks.size() << " analyses successful." << std::endl;
}

// Backward-compatible entry point; currently uses the same scheduling
void HMMTree::draw_tree_selective_concurrent(bool run_fitch, bool run_kitsch, bool run_upgma, bool run_nj, bool fm_only, bool min_only) {
    draw_tree_selective(run_fitch, run_kitsch, run_upgma, run_nj, fm_only, min_only);
}

//function to replace the shortednames to the before-replaced names
void HMMTree::trees_replace_shorted_names(std::string str_treefiles_path){
    if(!dir_exist_opendir(str_treefiles_path)){
        return;
    }

    std::vector<std::string>  vec_treefiles_names;
    if(!get_file_names(str_treefiles_path,vec_treefiles_names,"")){
        std::cout<<"Failed to get the tree files names in shorted names replace step !"<<std::endl;
        return;
	}

	for(unsigned int int_shortednames_i = 0; int_shortednames_i<vec_treefiles_names.size();int_shortednames_i++){
        std::ifstream input_treefile;
        std::string infile_path_name = folder_tree_files+vec_treefiles_names[int_shortednames_i];
        input_treefile.open(infile_path_name.c_str());
        if(!input_treefile.is_open()){
            std::cout<<"Failed to open the file: '"<<infile_path_name<<"' in the shorted names replace step !"<<std::endl;
            return;
        }



        std::string str_treefile_name_after_replacing_shorted_names = vec_treefiles_names[int_shortednames_i]+"_renamed";
        std::ofstream output_treefile;
        std::string outfile_path_name = folder_tree_files +str_treefile_name_after_replacing_shorted_names;
        output_treefile.open(outfile_path_name.c_str());
        if(!output_treefile.is_open()){
            std::cout<<"Failed to open the file: '"<<outfile_path_name<<"' in the shorted names replace step !"<<std::endl;
            return;
        }

        std::string str_one_line="";
        while(std::getline(input_treefile,str_one_line)){
            for(size_t uint_shorted_names_i = 0; uint_shorted_names_i < vec_unshorted_names.size(); ++uint_shorted_names_i){
                std::string str_unshorted_name = vec_unshorted_names[uint_shorted_names_i];
                std::string str_shorted_name = vec_shorted_names[uint_shorted_names_i];

                string_replace(str_one_line,str_shorted_name,str_unshorted_name,uint_shorted_names_i,vec_unshorted_names);
            }
            output_treefile<<str_one_line<<std::endl;
        }
        output_treefile.close();
        input_treefile.close();
        unlink(vec_treefiles_names[int_shortednames_i].c_str());
        if(rename(outfile_path_name.c_str(),infile_path_name.c_str())){
            std::cout<<"Failed to rename the file '"<<outfile_path_name <<"' to '"<< infile_path_name<<"' in the step of replacing the shorted names in tree files !"<<std::endl;
            return;
        }
	}
	return;
 }



