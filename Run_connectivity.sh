#!/bin/bash
#SBATCH --account=OD-232538
#SBATCH --time=1-00:00
#SBATCH --mem=40g
#SBATCH --cpus-per-task=10
#SBATCH --job-name='Coral - connectivity'
#SBATCH --array=0-21

module load python



dates=(
    "2015-10-29"
    "2015-11-28"
    "2015-12-27"
    "2016-10-18"
    "2016-11-17"
    "2016-12-16"
    "2017-10-08"
    "2017-11-06"
    "2017-12-06"
    "2018-10-27"
    "2018-11-25"
    "2018-12-25"
    "2019-10-16"
    "2019-11-15"
    "2019-12-14"
    "2020-10-04"
    "2020-11-03"
    "2020-12-02"
    "2021-01-01"
    "2021-10-23"
    "2021-11-21"
    "2021-12-21"

)

(( release = SLURM_ARRAY_TASK_ID ))

python matrix_calculations.py ${dates[release]}
