filter_sizes = [2**x for x in range(8,9)]
numbers_of_cohorts = [2]
#threshold_values = [0.00005,0.0001,0.0005,0.001,0.005, 0.01, 0.05, 0.1,0.5]
primary_threshold_values = [0.0001]
secondary_threshold_values = [0.0001]
primary_decisions = ['majority', 'any', 'all']
secondary_decisions = ['any','all']

populations = [300]
samples = [30000]

num_jobs = 1
