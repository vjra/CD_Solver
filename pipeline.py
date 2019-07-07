import subprocess

program_list = ['19_4_19_experiments_EV.py', '20_4_19_experiments.py']

for program in program_list:
    subprocess.call(['python3', program])
    print("Finished:" + program)
