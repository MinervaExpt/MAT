#!/usr/bin/env python
import ifdh
import datetime as dt
import subprocess, os, time, glob, itertools, sys, re, shutil


## Functions that are common to submitMergeToGrid.py and submitAuditToGrid.py

def createTarball( baseDirRoot , testReleaseName ):
  print "I'm inside createTarball()"
  print "Using Minerva_{0} in {1}".format(testReleaseName,baseDirRoot)
  os.environ["IFDH_DEBUG"] = "0"
  os.environ["IFDH_CP_MAXRETRIES"] = "0"

  ifdh_handle = ifdh.ifdh()
  
  user = os.environ["USER"]
  resilient_tardir = "/pnfs/minerva/resilient/tarballs"
  tarballName = '{0}_gridMerge-{1}.tgz'.format(user,dt.datetime.now().strftime("%Y%m%d_%H%M%S"))
  found = os.path.isfile("{0}/{1}".format(resilient_tardir,tarballName))
  if(not found):
    #Can't get this to work
    #cmd = "tar -czf /pnfs/minerva/resilient/tarballs/{0} -C {1} Minerva_{2}/cmt Minerva_{2}/Ana/PlotUtils Minerva_{2}/Tools/ProductionScriptsLite".format(tarballName,baseDirRoot,testReleaseName)
    #cmd = "tar -czf /pnfs/minerva/resilient/tarballs/{0} -C {1} Minerva_{2}".format(tarballName,baseDirRoot,testReleaseName)

    tmpTarball = "/minerva/data/users/{0}/{1}".format(user,tarballName) 
    cmd = "tar -czf {0} -C {1} Minerva_{2}".format(tmpTarball,baseDirRoot,testReleaseName)
    print cmd
    if (os.system(cmd)!=0):
      print "Failed to create the tarball."
      sys.exit(1)
    ifdh_handle.mv([tmpTarball,resilient_tardir+"/"+tarballName])

  print "I'm done creating the tarball {0}/{1}".format(resilient_tardir,tarballName)
  return tarballName

def unpackTarball( mywrapper , tarballName , testReleaseName ):
  # Add lines to wrapper that wil unpack tarball; add additional setup steps here if necessary  
  mywrapper.write("cd $CONDOR_DIR_INPUT\n")
  mywrapper.write("tar -xvzf {0}\n".format(tarballName))
  # Set up test release
  mywrapper.write("cd Minerva_{0}/Tools/ProductionScriptsLite/cmt/\n".format(testReleaseName))
  mywrapper.write("cmt config\n") #This step is necessary so that environment variables are constructed relative to $CONDOR_INPUT_DIR, as opposed to /minerva/app/<...>
  mywrapper.write("source setup.sh\n") #This step is necessary so that the environment variables get defined
  mywrapper.write("cd ../../../Ana/PlotUtils/cmt/\n")
  mywrapper.write("cmt config\n")
  mywrapper.write("source setup.sh\n")
  # Go back to the top-level working directory
  mywrapper.write("cd $CONDOR_DIR_INPUT\n")

def addBashLine( wrapper , command , logname ="" ):
  logger = ""
  if logname != "":
    logger = "&>> %s" % logname
  wrapper.write("echo '---------------' %s\n"% logger)
  wrapper.write("echo '-------%s'%s\n" % (command,logger) )
  wrapper.write("%s %s\n" % (command,logger) )
  wrapper.write("echo '---------------' %s\n"% logger)

def makeXROOTD( fileList ):
  for i,fileName in enumerate(fileList): fileList[i] = 'root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/' + fileName.lstrip('/pnfs/')

  return fileList

def gatherUnmergedDirs( dataSwitch , inputDir, processing, runs ):
  if dataSwitch == 'data':
    switchDir = '{indir}/grid/minerva/ana/*/*/{run}'
  elif dataSwitch == 'mc':
    switchDir = '{indir}/grid/{processing}/minerva/ana/*/{run}'

  runStr = "*/*/*/*"
  if runs ==-1:
    baseDir = switchDir.format(indir=inputDir,processing=processing,run=runStr)
    runList = glob.glob(baseDir) 
  elif type(runs)==list:
    runList=[]
    for run in runs:
      run = "{:0>8}".format(run)
      runStr = '/'.join([run[i:i+2] for i in range(0, len(run), 2)])
      baseDir = switchDir.format(indir=inputDir,processing=processing,run=runStr)
      tmpList = glob.glob(baseDir)
      if len(tmpList) == 0:
        print "Run {0} is not in {1}".format(run,switchDir.format(indir=inputDir,processing=processing,run=""))
        continue
      for _ in tmpList:
        runList.append(_)
        
  return runList,baseDir

def mapUnmergedDirToSubPlaylist( unmergedDir , mergedTupleShortName , tmpDir, outDirLogs, bDeleteFile = True ):
  ifdh_handle = ifdh.ifdh()

  unmergedFiles = os.popen("find %s -type f "%unmergedDir).readlines()
  unmergedFiles = makeXROOTD(unmergedFiles)  

  # I think it would be better (less error-prone) to not write directly to PNFS, but rather to write the file locally, then mv it to the PNFS location. However, every which way I tried to do the move, every file consistently got truncated to 404 lines...I tried various combination of 'cp', 'mv', 'ifdh cp' using subprocess and shutil
  outFileName = "playlist_{0}.txt".format(mergedTupleShortName)
  outFilePath = "{0}/{1}".format(outDirLogs,outFileName)
  tmpFilePath = "{0}/{1}".format(tmpDir,outFileName)

  outFile = open(tmpFilePath,'w')
  for unmergedFile in sorted(unmergedFiles):
    outFile.write(unmergedFile)

  outFile.close()

  ifdh_handle.cp([tmpFilePath,outFilePath])
  if bDeleteFile:
    os.remove(tmpFilePath)

  return outFileName

## Functions that are unique to submitMergeToGrid.py

def constructMergedFileShortName( toolName , dataSwitch , unmergedDir ):

  shortNameBase = '{0}_{1}_AnaTuple'.format(toolName,dataSwitch)

  unmergedDirComponents = unmergedDir.split('/')
  runNum = ''.join(unmergedDirComponents[-4:])
  shortName = '{0}_run{1}'.format(shortNameBase,runNum)

  ##XXX DATA!
  #if dataSwitch == 'data':
  #  shortName = shortNameBase
  #else:
  #  unmergedDirComponents = unmergedDir.split('/')
  #  runNum = ''.join(unmergedDirComponents[-4:])
  #
  #  shortName = '{0}_run{1}'.format(shortNameBase,runNum)

  return shortName

## Functions that are unique to submitAuditToGrid.py
def copyFile(path,outDirLogs,name) :
  if path:
    try:
      shutil.copyfile(path,"{dir}/{name}".format(dir=outDirLogs,name=name))
    except IOError:
      print "Failed to copy input file."
      sys.exit(1)
    return name
  else:
    return None

