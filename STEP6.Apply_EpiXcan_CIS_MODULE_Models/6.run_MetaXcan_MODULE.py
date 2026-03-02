
import os
import subprocess
import time
import sys

# Check if Joblib is installed, and install it if not
try:
    import joblib
except ImportError:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "joblib"])

# import it
import joblib
from joblib import Parallel, delayed

# Function to execute shell commands
def execute_shell_command(command):
    subprocess.run(command, shell=True)


start_all = time.time()

# Set the brain regions
regions = ["dlpfc"] ## your regions here
dataset = "your_testing_genotype_name" ## genotype to make predictions 

for region in regions:
        # Filter sh files based on brain region
        directory = "./STEP6.Apply_EpiXcan_CIS_MODULE_Models/MetaXcan_sh_scripts"
        files_sh = [file for file in os.listdir(directory) if "MODULE_%s_MetaXcan_script_%s" % (region,dataset) in file]

        files_sh = [os.path.join(directory, file_sh) for file_sh in files_sh]
        region_files_sh = [file_sh for file_sh in files_sh if region in file_sh]

        start = time.time()
        # Set the number of CPU cores
        num_cores = 22 ### Your number of cores


        Parallel(n_jobs=num_cores)(delayed(execute_shell_command)(script) for script in region_files_sh)

        end = time.time()
        print("Total time region %s : %s" % (region, end - start))

end_all = time.time()

print("Total time all regions: %s" % (end_all - start_all))
