import subprocess
import os, sys
from simulation_params import numbers_of_cohorts, filter_sizes, primary_decisions, secondary_decisions, primary_threshold_values, secondary_threshold_values, populations, samples, num_jobs
import itertools
import numpy as np
job_number = int(os.environ['SGE_TASK_ID'])
decode_string = f'decode_{job_number}.csv'
params_string = f'params_{job_number}.csv'
pop_params_string = f'pop_params_{job_number}.csv'

def delete_map():
    try:
        os.remove('map_{job_number}.Rdata')
    except FileNotFoundError: 
        print("Error: couldn't remove the map file, sim probabily failed", file=sys.stderr)


def write_params(params, pop_params):
    # Write params to file
    with open('params.csv', 'r') as file:
        data = file.readlines()

    data[1] = ','.join([str(x) for x in params])
    with open(f'{params_string}', 'w') as file:
        file.writelines(data)
        file.write('\n')

    with open('pop_params.csv', 'r') as file:
        data = file.readlines()

    data[1] = ','.join([str(x) for x in pop_params])
    with open(f'{pop_params_string}', 'w') as file:
        file.writelines(data)
        file.write('\n')
   
def write_decode_params(decode_params):
    with open('decode.csv', 'r') as file:
        data = file.readlines()

    decode_params = [str(x) for x in decode_params]
    data[1] = ','.join(decode_params)
    with open(f'{decode_string}', 'w') as file:
        file.writelines(data)
        file.write('\n')

        
def run_sim(decode_params):

    new_env = os.environ.copy()
    new_env["RAPPOR_REPO"] = "../../"
    bashCommand = f"Rscript cli.R map{job_number}.Rdata {params_string} {pop_params_string} {decode_string}"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, env=new_env)
    output, error = process.communicate()

    results = output.decode('utf').split('\n')
    results = [x.split(',')[1:] for x in results]
    results_check = results[6:8]
    results_check = sum([int(x[1]) for x in results_check])
    if results_check == 0:
        return "Signal too poor"
    results = [results[12],results[14]]
    return (', '.join([': '.join([y.strip() for y in x[1:]]) for x in results]))


# "k","h","m","p","q","f"
# 256,2,8,0.5,0.75,0
sim_results = []

# Initial default params which we will change and pass to the sim
params = [256,2,8,0.5,0.75,0]
pop_params = [300, 0.7, 0.3, 0,"Zipf", 10, 0.05, 300000, 0]
decode_params = [0.001, 0.001, "majority", "any"]

# Obviously this next bit is completely horrific but it was quick to do
params_set = itertools.product(filter_sizes, numbers_of_cohorts, populations, samples, primary_threshold_values,secondary_threshold_values, primary_decisions, secondary_decisions)
params_list = list(params_set)
job_params_list = np.array_split(params_list, num_jobs)

job_number = int(os.environ['SGE_TASK_ID'])

for p in job_params_list[job_number-1]:
    params[0] = p[0]
    params[2] = p[1]
    pop_params[0] = p[2]
    pop_params[7] = p[3]
    decode_params[0] = p[4]
    decode_params[1] = p[5]
    decode_params[2] = p[6]
    decode_params[3] = p[7]
    print(decode_params)
    # Deletion of the map: things that change the map go above this line
    delete_map()
    write_params(params,pop_params)
    write_decode_params(decode_params) 

    #decode_params_set = itertools.product(primary_threshold_values, secondary_threshold_values, primary_decisions, secondary_decisions)
    sim_results.append((run_sim(decode_params), params.copy(), list(decode_params).copy(), pop_params.copy()))

print("disasters, success, filter_size, cohort_size, primary_threshold, secondary_threshold, primary dec, secondary dec, population size, sample size")
for x in sim_results:
    print(x[0], end=',')
    print(str(x[1][0]) + ', ' + str(x[1][2]) + ", " + str(', '.join([str(s) for s in x[2]])) + ', ' + str(x[3][0]) + ', ' + str(x[3][7]))


