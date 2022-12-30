pipelineFolder="/users/fischerd/git/Pipeline-GBS"
projectFolder="/scratch/project_2001746/TestProject"

module load snakemake

export SINGULARITY_TMPDIR="/scratch/project/tmp"
export SINGULARITY_CACHEDIR="/scratch/project/tmp"

# Create the rulegraph
snakemake -s $pipelineFolder/GBS-pipeline.smk \
          --configfile $projectFolder/GBS-pipeline_config.yaml \
          --rulegraph | dot -T png > $projectFolder/workflow.png

snakemake -s $pipelineFolder/GBS-pipeline.smk \
          -j 150 \
          --use-singularity \
          --singularity-args "-B /scratch,/projappl,/users,/dev/shm:/tmp,/run:/run" \
          --configfile $projectFolder/GBS-pipeline_config.yaml \
          --latency-wait 60 \
          --cluster-config $pipelineFolder/GBS-pipeline_server-config.yaml \
          --cluster "sbatch -t {cluster.time} --account={cluster.account} --gres=nvme:{cluster.nvme} --job-name={cluster.job-name} --tasks-per-node={cluster.ntasks} --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu} --mail-user={cluster.mail-user} --mail-type={cluster.mail-type} -p {cluster.partition} -D {cluster.working-directory} --parsable" \
          --scheduler greedy \
          --cluster-cancel scancel \
          $@

