pipelineFolder="/users/fischerd/git/Pipeline-GBS"
projectFolder="/scratch/project_2001746/TestProject"

module load bioconda/3

# Create the rulegraph
snakemake -s $pipelineFolder/GBS-pipeline.smk \
          --configfile $projectFolder/GBS-pipeline_config.yaml \
          --rulegraph | dot -T png > $projectFolder/workflow.png

snakemake -s $pipelineFolder/GBS-pipeline.smk \
          -j 150 \
          --use-singularity \
          --singularity-args "-B /scratch,/projappl,/dev/shm:/tmp" \
          --configfile $projectFolder/GBS-pipeline_config.yaml \
          --latency-wait 60 \
          --cluster-config $pipelineFolder/GBS-pipeline_server-config.yaml \
          --cluster "sbatch -t {cluster.time} --account={cluster.account} --gres=nvme:{cluster.nvme} --job-name={cluster.job-name} --tasks-per-node={cluster.ntasks} --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu} --mail-user={cluster.mail-user} --mail-type={cluster.mail-type} -p {cluster.partition} -D {cluster.working-directory}" \
          --scheduler greedy \
          $@

