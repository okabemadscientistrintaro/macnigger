import os
import sys
import glob
import time
import tqdm


def set_abcluster_files(params: dict):
	natoms = 0
	elements = []
	for key, value in params.items():
		if key.startswith("NUMELEM"):
			natoms += value
		if key.startswith("ELEM"):
			elements.append(value)

	for i, elem in enumerate(elements):
		fname = elem+str(natoms)
		composition = elem+' '+str(natoms)
		tmpfolder = params["MOD1"]["TMP_FOLDER"]

		with open("./core/abcluster.inp" ,'r+') as f:   
			lines = f.readlines()
		f.close()

		lines = [w.replace('FNAME', fname) for w in lines]
		lines = [w.replace('COMPOSITION', composition) for w in lines]
		lines = [w.replace('MNUMCALC', str(params["MOD1"]["NUMCALC"])) for w in lines]
		lines = [w.replace('PATH_XTB', str(params["PATH_XTB"])) for w in lines]


		with open("./core/job.sh" ,'r+') as fi:   
			lin = fi.readlines()
		fi.close()
		lin = [w.replace('FNAME', fname) for w in lin]
		lin = [w.replace('PATH_ABCLUSTER', str(params["PATH_ABCLUSTER"])) for w in lin]

		################### CUBE
		print(f"\t\t\t- ABCluster Structures Generation: CUBE ({fname})")
		finp = [w.replace('STRUCTYPE', 'cube '+str(params["MOD1"]["CUBE"])) for w in lines]
		for k in range(0,params["MOD1"]["NUMGEN"]):
			folder = f"{tmpfolder}/RC_{fname}_{k:04d}"
			os.system("mkdir "+folder)
			fil = open(folder+'/'+fname+'.inp', 'w+')
			for l in finp:
				fil.write(l) 
			fil.close()

			jobfile = open(folder+'/job.sh', 'w+')
			for l in lin:
				jobfile.write(l) 
			jobfile.close()


		###################

		################### SPHERE
		print(f"\t\t\t- ABCluster Structures Generation: SPHERE ({fname})")
		finp = [w.replace('STRUCTYPE', 'sphere ') for w in lines]
		for k in range(0,params["MOD1"]["NUMGEN"]):
			folder = f"{tmpfolder}/RS_{fname}_{k:04d}"
			os.system("mkdir "+folder)
			fil = open(folder+'/'+fname+'.inp', 'w+')
			for l in finp:
				fil.write(l) 
			fil.close()

			jobfile = open(folder+'/job.sh', 'w+')
			for l in lin:
				jobfile.write(l) 
			jobfile.close()


		###################

def abcluster_run(params: dict):
	for folder in sorted(glob.glob("varmod1/"+"R*/")):
		print(folder)
		os.system("cd "+folder+"&& sh ./job.sh")
		subfolder = glob.glob(os.path.join(folder, '*/'))[0]
		basefolder = os.path.basename(os.path.normpath(folder))
		print(subfolder)
		print(basefolder)
		os.system("cd "+subfolder+" && for file in *.xyz; do cp $file ../../unfiltered/"+basefolder+"_$file; done > /dev/null")



# import subprocess

# command = "your_command_here"

# # Run the command and capture the output
# result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

# # Get the standard output and error
# output = result.stdout
# error = result.stderr

# # Print or process the output and error
# print("Output:", output)
# print("Error:", error)


####Run in cluster
#for folder in glob.glob("core/*/"):
#	os.system("cd "+folder+"&& qsub job.pbs")
#time.sleep(3600)  # tune this value so the cluster waits the computations finish 
#i=1
#for folder in glob.glob("core/*/"):
#	os.system("cd "+folder+cont+"_LM/ && for file in *.xyz; do cp $file ../../../all_xyz/"+str(i)+"$file; done")
#	i=i+1

