#!/usr/bin/ruby

=begin
This is a sequential queue submit script.
=end

#Select your system.
System = "xd"
#System = "opt"

#If Restart is "No", this script submit a job with "subfind.param" 
#If Restart is "Yes", this script submit a job with "subfind2.param"
#Case (1): RestartFromQueue = "No", simply submit a job
#Case (2): RestartFromQueue = "Yes", submit a job with -W depend=afterany:RestartID.  
Restart = "Yes"
#Restart = "No"
#RestartFromQueue = "No"
RestartFromQueue = "Yes"
RestartID = "5943"

NumberofQueue=10
BinaryName1 = "./ramses_dice/bin/ramses3d"
BinaryName2 = "./ramses/bin/ramses3d"
NumNodes=4
NumThredsPerProc=1
NumProc = NumNodes * 112 / NumThredsPerProc
NumTasksPerNode=112/NumThredsPerProc
NumTasksPerSocket = 56 / NumThredsPerProc

RunName = "rt_new"
StartUp = "srun"
#QName="M-large-b"
#QName="M-test-b"
#QName="P-test-b"
QName="P-large-b"
Lcommand='nodes=' + NumNodes.to_s

def write_script(str, snapnum, new)

    outputfile = open(str,"w")
    outputfile.print "#!/bin/bash\n"     
    outputfile.print '#SBATCH --job-name=', RunName, "\n"
    outputfile.print '#SBATCH --partition=', QName, "\n"
    outputfile.print '#SBATCH --nodes=', NumNodes.to_s, "\n"
    outputfile.print '#SBATCH --ntasks=', NumProc.to_s, "\n"
    outputfile.print '#SBATCH --ntasks-per-node=', NumTasksPerNode.to_s, "\n"
    outputfile.print '#SBATCH --ntasks-per-socket=', NumTasksPerSocket.to_s, "\n"
    outputfile.print "#SBATCH -o flat-%J.out\n" 
    outputfile.print "#SBATCH --mem=461G\n"    # ノード当たりのメモリ容量
    #outputfile.print "#SBATCH --mem=115G\n"    # ノード当たりのメモリ容量
    outputfile.print "#SBATCH --time=24:00:00\n"   
    #outputfile.print "#SBATCH --time=00:30:00\n"   
    outputfile.print ". /work/opt/local/bin/enable-oneapi.sh \n"

    outputfile.print "\n"
    outputfile.print "cd ${SLURM_SUBMIT_DIR}\n"
    if new == 1
        outputfile.print StartUp, " ", BinaryName1, " " 
        outputfile.print "xd_rt.nml > log.", snapnum.to_s, " 2>&1" "\n"
    else
        outputfile.print StartUp, " ", BinaryName2, " " 
        outputfile.print "xd_rt_rs.nml > log.", snapnum.to_s, " 2>&1" "\n"
    end

    outputfile.flush
    outputfile.close
end

def Attention 
  if Restart == "Yes" 
    if RestartFromQueue == "YES" 
      if RestartID == "" 
        print "Incorrect RestartID\n" 
        exit(1) 
      end 
    end
  end
end

Attention()

result = 0
(NumberofQueue).times{ |i| 
  if i == 0 
    if Restart == "No" 
      write_script("tmp", i, 1) 
      `chmod 744 ./tmp` 
      result = `sbatch ./tmp` 
      print result
      print "Submit first script:",result.chomp,"\n" 
      `rm -rf ./tmp` 
    else 
      if RestartFromQueue == "No" 
        write_script("tmp_", i, 2) 
        `chmod 744 ./tmp_` 
        result = `sbatch ./tmp_` 
        print "Submit first script:",result.chomp,"\n" 
        `rm -rf ./tmp_` 
      else 
        write_script("tmp_", i, 2) 
        `chmod 744 ./tmp_` 
        command = "sbatch -d afterany:" + RestartID + " ./tmp_" 
        result = `#{command}` 
        print "Submit first script:",result.chomp,"\n" 
        `rm -rf ./tmp_` 
      end 
    end 
  else 
    print "hoge\n"
    print result, "\n" 
    print "moge\n"
    str = result.chomp!
    qid = str[/\d+$/]
    print qid, "\n"
    write_script("tmp_", i,2) 
    `chmod 744 ./tmp_` 
    command = "sbatch -d afterany:" + qid.to_s + " ./tmp_" 
    print "hoge\n"
    print command, "\n"
    print "moge\n"
    result = `#{command}` 
    print "Submit ",i+1,"-th script:",result.chomp,"\n" 
    `rm -rf ./tmp_` 
  end
}

