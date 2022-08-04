import subprocess
import os

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

def write_decode_params(decode):
    with open('decode.csv', 'r') as file:
        data = file.readlines()

    data[1] = str(decode) + ",majority, majority"
    with open('decode.csv', 'w') as file:
        file.writelines(data)
        file.write('\n')

def run_sim(params, threshold):

    write_params(params)
    write_decode_params(threshold) 
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
thresholds = [0.002, 0.005, 0.001]
filter_sizes = range(0,2)
filter_sizes = [2**x for x in filter_sizes]
cohort_sizes = [256]#,64,128,256,512,1024,2048,4096]


for filter_size in filter_sizes:
    params[0]=filter_size
    for size in cohort_sizes:
        delete_map()
        for t in thresholds:
            params[2] = size
            sim_results.append((run_sim(params, t), params.copy()))

print("disasters, success, filter_size, cohort_size")
for x in sim_results:
    print(x[0], end=',')
    print(str(x[1][0]) + ', ' + str(x[1][2]))


