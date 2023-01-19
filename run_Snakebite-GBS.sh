pipelineFolder="~/git/Snakebite-GBS"
projectFolder="/scratch/project_2001746/GBS_Example"

module load snakemake

export APPTAINER_TMPDIR="/scratch/project_2001746/tmp"
export APPTAINER_CACHEDIR="/scratch/project_2001746/tmp"

mkdir -p $APPTAINER_TMPDIR

# Create the rulegraph
snakemake -s $pipelineFolder/Snakebite-GBS.smk \
          --configfile $projectFolder/Snakebite-GBS_config.yaml \
          --rulegraph | dot -T png > $projectFolder/workflow.png

snakemake -s $pipelineFolder/Snakebite-GBS.smk \
          -j 150 \
          --use-singularity \
          --singularity-args "-B /scratch,/projappl,/users,/dev/shm:/tmp,/run:/run" \
          --configfile $projectFolder/Snakebite-GBS_config.yaml \
          --latency-wait 60 \
          --cluster-config $projectFolder/Snakebite-GBS_server-config.yaml \
          --cluster "sbatch -t {cluster.time} --account={cluster.account} --gres=nvme:{cluster.nvme} --job-name={cluster.job-name} --tasks-per-node={cluster.ntasks} --cpus-per-task={cluster.cpus-per-task} --mem-per-cpu={cluster.mem-per-cpu} --mail-user={cluster.mail-user} --mail-type={cluster.mail-type} -p {cluster.partition} -D {cluster.working-directory --parsable}" \
          --scheduler greedy \
	  --cluster-cancel scancel \
          $@

