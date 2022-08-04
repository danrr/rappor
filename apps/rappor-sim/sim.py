import subprocess
import os
from simulation_params import numbers_of_cohorts, filter_sizes, primary_decisions, secondary_decisions, threshold_values
import itertools
def delete_map():
    try:
        os.remove('map.Rdata')
    except FileNotFoundError: 
        print("Error: couldn't remove the map file, sim probabily failed")


def write_params(params):
    # Write params to file
    with open('params.csv', 'r') as file:
        data = file.readlines()

    data[1] = ','.join([str(x) for x in params])
    with open('params.csv', 'w') as file:
        file.writelines(data)
        file.write('\n')

def write_decode_params(decode_params):
    with open('decode.csv', 'r') as file:
        data = file.readlines()

    decode_params = [str(x) for x in decode_params]
    data[1] = ','.join(decode_params)
    with open('decode.csv', 'w') as file:
        file.writelines(data)
        file.write('\n')

def run_sim(decode_params):

    write_decode_params(decode_params) 
    new_env = os.environ.copy()
    new_env["RAPPOR_REPO"] = "../../"
    bashCommand = "Rscript cli.R map.Rdata params.csv pop_params.csv decode.csv"
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

params = [2048,2,8,0.5,0.75,0]
decode_params = [0.001, "majority", "any"]

# Obviously this next bit is completely horrific but it was quick to do
for filter_size in filter_sizes:
    params[0]=filter_size
    for size in numbers_of_cohorts:
        params[2] = size
        # Deletion of the map: things that change the map go above this line
        delete_map()
        write_params(params)
        # Create all possible permutations of decode params
        decode_params_set = itertools.product(threshold_values, primary_decisions, secondary_decisions)
        for d in decode_params_set:
            sim_results.append((run_sim(d), params.copy(), list(d).copy()))

print("disasters, success, filter_size, cohort_size, threshold, primary dec, secondary dec")
for x in sim_results:
    print(x[0], end=',')
    print(str(x[1][0]) + ', ' + str(x[1][2]) + str(','.join([str(s) for s in x[2]])))


