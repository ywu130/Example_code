#Submit multiple files to XSEDE
#Author: Gaurav MOghe
import sys, os

def help1():
    print ("Submit each line in the commands file as a separate job")
    print ("REQ:")
    print ("-commands:  File containing list of commands (see below)")
    print ("-jobBase:   Base name for the job")
    print ("-mail-user: Email address")
    print ("-ntasks:    Always 1")
    print ("-cpus-per-task: Upto 60")
    print ("-mem:   Base name for the job")
    print ("-time: Time in HH:MM:SS format")
    print ("-outputFile:    Name of log file")
    print ("##########")
    print ("NOTE: The commands file should have all export paths in the beginning" \
          "starting with the # sign")
    print ("##########")

    sys.exit()

fn=jobBase=mailUser=ntasks=cpus=mem=time=log=""
for i in range(1,len(sys.argv)):
    if sys.argv[i]=='-commands':
        fn=sys.argv[i+1]
    if sys.argv[i]=='-jobBase':
        jobBase=(sys.argv[i+1])
    if sys.argv[i]=='-mail-user':
        mailUser=sys.argv[i+1]
    if sys.argv[i]=='-ntasks':
        ntasks=sys.argv[i+1]
    if sys.argv[i]=='-cpus-per-task':
        cpus=sys.argv[i+1]
    if sys.argv[i]=='-mem':
        mem=sys.argv[i+1]
    if sys.argv[i]=='-time':
        time=sys.argv[i+1]
        time=sys.argv[i+1]
#out1=open(fn+".log",'w')

if fn=="" or jobBase=="" or mailUser=="" or cpus=="" or mem=="" or time=="":
    help1()
else:
    file1=open(fn,'r')
    line1=file1.readline()
    plist=[]; i=0
    while line1:
        if line1.startswith('#'):
            plist.append(line1[1:])
        else:
            #Build job
            i+=1
            tmp1=open('{}{}.sh'.format(jobBase,i),'w')
            tmp1.write('#!/bin/sh\n\n')
            tmp1.write('#SBATCH --job-name={}{}.sh # Job name\n'.format(jobBase,i))
            tmp1.write('#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)\n')
            tmp1.write('#SBATCH --mail-user={} # Where to send mail\n'.format(mailUser))
            tmp1.write('#SBATCH --ntasks=1 # Run a single task\n')
            tmp1.write('#SBATCH --cpus-per-task={} # Number of CPU cores per task\n'.format(cpus))
            tmp1.write('#SBATCH --mem={}gb # Job memory request\n'.format(mem))
            tmp1.write('#SBATCH --time={} # Time limit hrs:min:sec\n'.format(time))
            tmp1.write('#SBATCH --output={}{}.log # Standard output and error log\n'.format(jobBase,i))
            tmp1.write('pwd; hostname; date\n')
            for path in plist:
                tmp1.write(path)
            tmp1.write(line1)
            tmp1.close()

            #Submit job
            os.system('sbatch {}{}.sh'.format(jobBase,i))
            #sys.exit()
        line1=file1.readline()
    file1.close()
print ("Jobs submitted: ", i)
print ("Done!")
