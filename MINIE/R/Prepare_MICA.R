#Create shell script to run scMINER
#sh file

generate_MICA_cmd<-function(save_sh_at,
                            input_file,
                            project_name,
                            num_cluster,
                            output_path,
                            host="lsf",
                            queue=NULL,
                            memory=NULL,
                            threads=10,
                            bootstrap=10,
                            dim_reduction_method="MDS",
                            visualization="tsne",
                            perplexity=30,
                            min_dist=0.01,
                            slice_size=1000){
                            
  file.sh<-file.path(save_sh_at,paste0("01_run_MICA_",project_name,'.sh'))
  if(file.exists(file.sh)) stop("File already existed!")
  
  sink(file.sh)
  
  if (tolower(host)=="lsf"){
    project<-ifelse(is.null(project_name),'',paste('#BSUB -P ',project_name,' \n'))
    queue<-ifelse(is.null(queue),'',paste('#BSUB -q ', q.job,' \n'))
    
    sh.scminer<-paste0(
      '#!/bin/env bash\n',
      project,
      '#BSUB -oo ',project_name,'.sh.out \n',
      '#BSUB -eo ',project_name,'.sh.err \n',
      '#BSUB -R \"rusage[mem=2000]\" \n',
      queue,
      "mica lsf ",
      ifelse(is.null(memory),"",paste0("-r ",paste(memory, collapse = " "), " ")))
  } else if (tolower(host)=="local") {
    sh.scminer<-paste0(
      '#!/bin/env bash\n' ,
      "mica local ")
  }
  
  #add generic attributes
  sh.scminer<-paste0(sh.scminer,
    paste0("-i ", normalizePath(input_file), " "),
    paste0("-p ", project_name, " "),
    paste0("-k ",  paste0(num_cluster, collapse=" ", " ")),
    paste0("-o ",  normalizePath(output_path), " "),
    ifelse(is.null(bootstrap),"",paste0("-b ",bootstrap," ")),
    ifelse(is.null(dim_reduction_method),"",paste0("-dr ",dim_reduction_method," ")),
    ifelse(is.null(visualization),"",paste0("-v ",visualization, " ")),
    ifelse(is.null(perplexity),"",paste0("-pp ",perplexity," ")),
    ifelse(is.null(min_dist),"",paste0("-d ",min_dist," ")),
    ifelse(is.null(slice_size),"",paste0("-sn ", slice_size," ")),
    ifelse(is.null(threads),"",paste0("-t ", threads," "))
  )
  cat(sh.scminer)
  sink()
  
  cat(basename(file.sh),'is generated!\n')
  
}
