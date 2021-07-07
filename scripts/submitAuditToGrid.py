#!/bin/sh
import argparse
import ifdh
from functionsForSubmissionScripts import *

## README
# The 'only' assumptions this script makes about your personal set up is that you have a test release of the Minerva software framework, which corresponds in this code to `testReleaseName` and lives in an area pointed at by `baseDirRoot`. The only required packages are `Tools/ProductionScriptsLite` and `Ana/PlotUtils`. The user is invited to customize as suits their needs.
# The merging is presumed to be at the playlist level for data, and at the run level for MC. That is to say, all data tuples in each playlist get merged together and all MC tuple in each run of each playlist get merged together.
######################

def getMergedTuples( mergedDir, toolName, dataSwitch, template=None ):
  if template:
    tupleBaseName=template
  else:
    tupleBaseName = '{tool}_{dataTag}'.format( tool=toolName, dataTag=dataSwitch )
  tupleList = glob.glob('{0}/{1}*.root'.format( mergedDir, tupleBaseName ))

  tupleList = makeXROOTD( tupleList )

  return tupleList

def getUnmergedDirectoryName( unmergedDir, dataSwitch, runNum, releaseName, process ):
  if True: #dataSwitch != 'data':
    runNum       = "{:0>8}".format(runNum)
    runNumString = '/'.join([runNum[i:i+2] for i in range(0, len(runNum), 2)])

    version = re.search("v[0-9]+r[0-9]+p[0-9]+",releaseName)
    if not version:
      print "Version number {0} not in standard format".format(releaseName)
      return
  
  if dataSwitch == 'data':
    baseDir = '{directory}/grid/minerva/ana/*/{release}/'.format(directory=unmergedDir, release=version.group(0))
    unmergedDirectoryPath = baseDir + runNumString
  elif dataSwitch == 'mc':
    baseDir = '{directory}/grid/{process}/minerva/ana/{release}/'.format(directory=unmergedDir,release=version.group(0),  process=process)
    unmergedDirectoryPath = baseDir + runNumString

  #This is terrible, but merging non numibeam files seems very remote.  A
  fillWildcard = glob.glob(unmergedDirectoryPath)

  return fillWildcard[0]

def submitJob( outDir , tmpDir , outDirLogs , mergedFilePath , subPlaylist , tarballName , treeName , releaseName , release_dir , memory , disk , bTest ):
  # Probably specific to name construction of merged tuples
  shortName = mergedFilePath.split('/')[-1].split('.')[0]

  ifdh_handle = ifdh.ifdh()
  # Create wrapper
  wrapper_name = "{0}/wrapper_{1}.sh".format(outDirLogs,shortName) 
  tmp_wrapper  = "{0}/wrapper_{1}.sh".format(tmpDir,mergedFileShortName)  
  logname = "$CONDOR_DIR_LOGS/{0}.log".format(mergedFileShortName) 
  
  my_wrapper = open(tmp_wrapper,"w")
  my_wrapper.write("#!/bin/sh\n")
  addBashLine(my_wrapper,'env',logname)
  addBashLine(my_wrapper,'pwd',logname)
  addBashLine(my_wrapper,'date',logname)
  unpackTarball(my_wrapper,tarballName,release_dir)
  addBashLine(my_wrapper,'pwd',logname)
  addBashLine(my_wrapper,'ls',logname)
  addBashLine(my_wrapper,"echo $HOSTNAME",logname)

  # This is the bash line that will be executed on the grid
  my_wrapper.write( "python $CONDOR_DIR_INPUT/Minerva_{release}/Ana/PlotUtils/scripts/AuditMergedFiles.py -m '{merge}' -p '{playlist}' -t '{tool}' --grid &>> $CONDOR_DIR_LOGS/{name}.log\n".format(release=release_dir,merge=mergedFilePath,playlist=subPlaylist,tool=treeName,name=shortName) )
  
  addBashLine(my_wrapper,'date >> $CONDOR_DIR_LOGS/{0}.log'.format(shortName))
  my_wrapper.close()
  
  ifdh_handle.cp([tmp_wrapper,wrapper_name])
  os.remove(tmp_wrapper)

  os.system( "chmod 777 %s" % wrapper_name )
  
  version = re.search("v[0-9]+r[0-9]+p[0-9]+",releaseName)
  if not version:
    print "Version number {0} not in standard format".format(releaseName)
    return

  #Executables compiled on a specific OS are only guaranteed to work on that OS on the grid.  Infer the OS from CMT.
  cmtConfig = os.environ["CMTCONFIG"]
  whichOS = cmtConfig.split('-')[1]
  whichOS = re.sub("slc", "sl", whichOS)

  #Don't let the user submit jobs with v21r1p1 to SL7 Singularity containers.
  if(version.group(0) == "v21r1p1" and whichOS != "sl6"):
    print "MINERvA framework version {version} has not been compiled for OS {whichOS}, so I will not submit this job!".format(version=version.group(0), whichOS=whichOS)
    sys.exit(0) #TODO: We should really be exiting with something nonzero everywhere that there is a problem.  Keeping with convention for now.

  #cmd = "jobsub_submit --group=minerva --cmtconfig=x86_64-slc6-gcc49-opt --OS=SL6 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --role=Analysis --memory {memory}MB --disk {disk}GB -f {outDirLogs}/{sub} -f /pnfs/minerva/resilient/tarballs/{tarball} -d LOGS {outDir} -r {release} -i /cvmfs/minerva.opensciencegrid.org/minerva/software_releases/{release}/ file://{wrapper}".format(memory=memory,disk=disk,outDir=outDir,outDirLogs=outDirLogs,sub=subPlaylist,tarball=tarballName,wrapper=wrapper_name,release = version.group(0))
  cmd = "jobsub_submit --group=minerva --cmtconfig={cmtConfig} --append_condor_requirements='(TARGET.HAS_SINGULARITY=?=true)' --lines '+SingularityImage=\\\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-{whichOS}:latest\\\"' --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --role=Analysis --memory {memory}MB --disk {disk}GB -f {outDirLogs}/{sub} -f /pnfs/minerva/resilient/tarballs/{tarball} -d LOGS {outDir} -r {release} -i /cvmfs/minerva.opensciencegrid.org/minerva/software_releases/{release}/ file://{wrapper}".format(cmtConfig=cmtConfig,memory=memory,disk=disk,outDir=outDir,outDirLogs=outDirLogs,sub=subPlaylist,tarball=tarballName,wrapper=wrapper_name,release = version.group(0), whichOS=whichOS)
 
  print cmd
  if not bTest:
    os.system(cmd)
 
  # Don't pound the grid; spread out submissions slightly 
  sleepTime = 2 
  print "Sleeping for %i seconds\n" % sleepTime
  time.sleep(sleepTime)

##############################################################
#submitAuditToGrid.py
##############################################################
if __name__ == '__main__': 
  userEnv = os.environ["USER"]

  ##############################################################
  #Get all the user inputs
  ##############################################################
  parser = argparse.ArgumentParser(description="Submit directories of merged files and the unmerged files used to create them to the grid in order to audit.  Currently, this only works with files on persistent")
  #Set things up so required arguments first
  optional_args = parser._action_groups.pop()
  required_args = parser.add_argument_group('required arguments')
  parser._action_groups.append(optional_args)

  required_args.add_argument("--tool",required=True,
                             help="The name of the analysis tool")
  required_args.add_argument("--release",required=True,
                             help ="Minerva release (vXrYpZ<additional identifier>)")
  required_args.add_argument("--unmerged_dir",required=True,
                             help="Directory in /pnfs/minerva/persistent/users/$USER where the unmerged files live")
  required_args.add_argument("--merged_dir",required=True,
                             help="Directory in /pnfs/minerva/persistent/users/$USER where the merged files live")

  #Flags for data switches
  optional_args.add_argument("--data",action="store_true",
                             help="Unmerged directory are from data")
  optional_args.add_argument("--mc",action="store_true",
                             help="Unmerged directory are from mc")
                              
  optional_args.add_argument("--scratch",action="store_true",
                             help="Change the base directory from /pnfs/minerva/persistent/users/$USER to /pnfs/minerva/scratch/users/$USER")
  optional_args.add_argument("--tree",default=None,
                             help="Name of the analysis tree.  (Default: <name of analysis tool>)")
                              
  required_args.add_argument("--unmerged_release",default="",
                             help ="Minerva release for the release files (vXrYpZ<additional identifier>)")
  optional_args.add_argument("--template",default=None,
                             help="Allows for different naming schemes for the merged files.  You will need to provide a run number if mc (Default: None)")
  optional_args.add_argument("--run",nargs="+",default=-1,
                             help="Run numbers if you want to look at specific runs. Can input multiple runs (ex. --run 1000 1250)  Required when using template")
  optional_args.add_argument("--processing",default="central_value",
                             help="Allows for merging when look at specifically processings (ex. CCCohSignal) (default: central_value")
  optional_args.add_argument("--base_dir",default="/minerva/app/users/$USER/cmtuser",
                             help="This points to base directory for your Minerva release ( default: /minerva/app/users/$USER/cmtuser )")
  optional_args.add_argument("--release_dir",default="",
                             help="This points to base directory for your Minerva release ( default :/minerva/app/users/$USER/cmtuser )")
  optional_args.add_argument("--user",default =userEnv,
                             help="The name of the user who owns the unmerged files (Default: $USER)")
  optional_args.add_argument("--memory",type=int,default=500,
                             help="This is the amount of memory to use for your job ( default: 500 MB )") 
  optional_args.add_argument("--disk",type=int,default=10,
                             help="This is the amount of disk limit to use for your job in GB ( default: 10 GB )") 
  optional_args.add_argument("--test",action="store_true",
                             help="Run this without submitting this to the grid")

  args = parser.parse_args()
  if not (args.data or args.mc):
    print "submitAuditToGrid.py: error: --data or --mc is required"
    sys.exit(0)

  if args.data and args.mc:
    print "submitAuditToGrid.py: error: cannot have both --data or --mc"
    sys.exit(0)

  if args.data:
    dataSwitch="data"
  elif args.mc:
    dataSwitch="mc"

  if not args.tree:
    anaTree = args.tool
  else:
    anaTree = args.tree

  if args.template and args.run==-1:
    print "Need a run number if using template"
    sys.exit(0)

  if args.template and len(args.run)>1:
    print "Cannot currently audit templated files with multiple runs"
    sys.exit(0)

  mergeDir = "/pnfs/minerva/{dcache}/users/{user}/{outdir}".format(dcache="scratch" if args.scratch else "persistent" , user=args.user, outdir=args.merged_dir )

  if not os.path.isdir(mergeDir):
    print "{0} does not exist!  Exiting".format(mergeDir)
    sys.exit(0)

  ##############################################################
  #Get the merged tuples
  ##############################################################
  mergedTuples = getMergedTuples( mergeDir, args.tool, dataSwitch, args.template )
  if len(mergedTuples) == 0:
    print "No files in {0}.  Exiting before creating tarball or output directory".format(mergeDir)
    sys.exit(0)

  print 'mergedTuples: ' 
  for _ in mergedTuples:
    print _

  outDirRuns = ""
  if args.run != -1:
    outDirRuns = "_"+"_".join(args.run[:3])#Giving a sensible title
    print "Looking at runs: ", args.run
  ##############################################################
  #Make the audit output directory
  ##############################################################
  outDir = "{mergeDir}/audit{runs}".format(mergeDir=mergeDir, runs=outDirRuns )
  if not os.path.isdir(outDir):
    print "Making output directory for audit {0}".format(outDir)
    os.system( "mkdir -p %s" % outDir )
  else:
    print "{0} already exists! Exiting before repopulating".format(outDir)
    sys.exit(0)

  #Need a place to make the wrappers and playlists before copying them over
  #Putting them into a directory in the users data area
  blueArcDir = "/minerva/data/users/{user}/grid_merge".format(user=userEnv)
  if not os.path.isdir(blueArcDir):
    print "Making BlueArc directory {0}".format(blueArcDir)
    os.system("mkdir -p %s" % blueArcDir)

  tag = 'audit_{0}-Tuples'.format( args.tool )
  processingID = '%s_%s-%s00' % (tag , dt.date.today() , dt.datetime.today().strftime("%H") )
  outDirLogs = "{0}/{1}_wrapper_playlist".format(outDir,processingID)
  os.system( "mkdir -p %s" % outDirLogs )

  ##############################################################
  #Create the tarball
  ##############################################################
  if len(args.release_dir)==0:
    args.release_dir = args.release
  if len(args.unmerged_release)==0:
    args.unmerged_release = args.release
  print args.release_dir
  tarballName = createTarball( args.base_dir, args.release_dir )

  ##############################################################
  #Get a list of unmerged runs so we can check whether it is there later
  ##############################################################
  unmergedRuns = []
  tupleDirsToMerge,baseDir = gatherUnmergedDirs( dataSwitch, "/pnfs/minerva/{0}/users/{1}/{2}".format("scratch" if args.scratch else "persistent",args.user,args.unmerged_dir), args.processing, args.run ) 
  for unmergedDir in tupleDirsToMerge:
    unmergedDirComponents = unmergedDir.split('/')
    runNum = ''.join(unmergedDirComponents[-4:])
    unmergedRuns.append(runNum)
  unmergedRuns.sort()

  for mergedFilePath in mergedTuples:
    mergedFileShortName = mergedFilePath.split('/')[-1].split('.')[0]
    runNum = -1
    #Using template
    # with run: grab run number from the user input
    if args.template:
      runNum = args.run[0] 
    else:
      #If not using template:
      # no run: use the run number on the merged Tuples
      runMatch = re.search("\d{8}",mergedFilePath[mergedFilePath.find("run")+3:-5])
      runNum = runMatch.group(0)
      # with run: see if the merged run num is in the runs given
      if type(args.run)==list: 
        formatRuns = ["{:0>8}".format(_) for _ in args.run]
        if runNum not in formatRuns:
          continue

    if runNum in unmergedRuns:
      unmergedRuns.remove(runNum)
    else:
      print "Somehow, you have a merged run",runNum,"without a unmerged directory"

    
    unmergedDir = getUnmergedDirectoryName( "/pnfs/minerva/{0}/users/{1}/{2}".format("scratch" if args.scratch else "persistent", args.user,args.unmerged_dir), dataSwitch, runNum, args.unmerged_release, args.processing )

    if not os.path.isdir(unmergedDir):
      print "{0} is not a directory!  Skipping".format(unmergedDir)
      continue
    subPlaylist = mapUnmergedDirToSubPlaylist( unmergedDir, mergedFileShortName, blueArcDir, outDirLogs )
    submitJob( outDir , blueArcDir , outDirLogs , mergedFilePath , subPlaylist , tarballName , anaTree , args.release , args.release_dir , args.memory , args.disk , args.test )

  if len(unmergedRuns) > 0:
    miss_log ="{0}/missing_merged_runs.log".format(outDir) 
    print "\n#########################################################################################"
    print "You have unmerged runs without merged files! Look at {0} for a list".format(miss_log)
    print "#########################################################################################"
    with open(miss_log,"w") as f:
      for runNum in unmergedRuns:
        f.write("{0}\n".format(runNum))
