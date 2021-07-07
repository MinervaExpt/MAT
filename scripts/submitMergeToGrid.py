#!/usr/bin/env python
import argparse
import ifdh
from functionsForSubmissionScripts import *

## README
# The 'only' assumptions this script makes about your personal set up is that
# you have a test release of the Minerva software framework, which corresponds
# in this code to `releaseName` and lives in an area pointed at by
# `baseDirRoot`. The only required packages are `Tools/ProductionScriptsLite`
# and `Ana/PlotUtils`. The user is invited to customize as suits their needs.

# The merging is presumed to be at the playlist level for data, and at the run
# level for MC. That is to say, all data tuples in each playlist get merged
# together and all MC tuple in each run of each playlist get merged together.

################################################################################
# Create grid script/wrapper, formulate jobsub_submit command, submit
################################################################################
def submitJob(
    outDir,
    stagingArea,
    outDirPlaylist,      # playlist subdir within the outdir
    outDirLogs,          # logfile subdir within the outdir
    mergedFileShortName, # "<anatool>_<data/mc>_AnaTuple_runXXXXXXXX"
    subPlaylist,         # "playlist_<anatool>_<data/mc>_AnaTuple_runXXXXXXXX.txt"
    dataSwitch,          # <data/mc>
    args,
    tarballName,
    variables,
):
    ifdh_handle = ifdh.ifdh()

    ############################################################################
    # Process args
    ############################################################################
    toolName = args.tool

    if not args.tree:
        treeName = args.tool
    else:
        treeName = args.tree

    releaseName = args.release
    release_dir = args.release_dir

    memory = args.memory
    disk = args.disk

    with_var = args.with_var

    bTest = args.test
    bDebug = args.debug
    bInteractive = args.interactive

    version = re.search("v[0-9]+r[0-9]+p[0-9]+", releaseName)
    if not version:
        print "Version number {0} not in standard format".format(releaseName)
        return

    if bDebug:
        os.environ["IFDH_CP_MAXRETRIES"] = "0"

    if bInteractive:
        subPlaylist = "{0}/{1}".format(stagingArea, subPlaylist)
        activeOutDir = outDir
        logname = "{0}/{1}.log".format(stagingArea, mergedFileShortName)
    else:
        activeOutDir = "$CONDOR_DIR_MERGE"
        logname = "$CONDOR_DIR_LOGS/{0}.log".format(mergedFileShortName)

    ############################################################################
    # Create grid script/wrapper
    ############################################################################
    # Create it in the staging area and then copy it to the outdir.
    # TODO: the wrapper can come directly from BlueArc/local and jobsub will
    # automatically copy it to your outdir.
    wrapper_name = "{0}/wrapper_{1}.sh".format(outDirPlaylist, mergedFileShortName)
    tmp_wrapper = "{0}/wrapper_{1}.sh".format(stagingArea, mergedFileShortName)

    my_wrapper = open(tmp_wrapper, "w")
    my_wrapper.write("#!/bin/sh\n")
    addBashLine(my_wrapper, "env", logname)
    addBashLine(my_wrapper, "pwd", logname)
    addBashLine(my_wrapper, "date", logname)
    if not bInteractive:
        unpackTarball(my_wrapper, tarballName, release_dir)
    addBashLine(my_wrapper, "pwd", logname)
    addBashLine(my_wrapper, "ls", logname)
    addBashLine(my_wrapper, "which MergeFiles.exe", logname)
    addBashLine(my_wrapper, "which MergeFilesWithVar.exe", logname)
    addBashLine(my_wrapper, "echo $HOSTNAME", logname)

    # This is the bash line that will be executed on the grid
    if with_var:
        my_wrapper.write(
            "MergeFilesWithVar.exe -g -o {active}/{name}.root -p '{sub}' {var} &>> {logname}\n".format(
                sub=subPlaylist,
                active=activeOutDir,
                name=mergedFileShortName,
                logname=logname,
                var=(variables if variables else ""),
            )
        )
    else:
        if dataSwitch == "data":
            my_wrapper.write(
                "MergeFiles.exe -data -e '{sub}' -o '{active}' -s '{name}' -a '{tool}' -t '{tree}' &>> {logname}\n".format(
                    sub=subPlaylist,
                    active=activeOutDir,
                    name=mergedFileShortName,
                    logname=logname,
                    tool=toolName,
                    tree=treeName,
                )
            )
        else:
            my_wrapper.write(
                "MergeFiles.exe -e '{sub}' -o '{active}' -s '{name}' -a '{tool}' -t '{tree}' {only_tracker_truth} &>> {logname}\n".format(
                    sub=subPlaylist,
                    active=activeOutDir,
                    name=mergedFileShortName,
                    logname=logname,
                    tool=toolName,
                    tree=treeName,
                    only_tracker_truth="-r" if only_tracker_truth else "",
                )
            )

    my_wrapper.write("exit_code=$?\n")
    my_wrapper.write(
        'echo "Exit Code of Merge is $exit_code." &>> {logname}\n'.format(
            logname=logname
        )
    )
    addBashLine(my_wrapper, "date", logname)
    my_wrapper.write("exit $exit_code")
    my_wrapper.close()

    # copy from staging area to outdir
    ifdh_handle.cp([tmp_wrapper, wrapper_name])

    ############################################################################
    # Submit
    ############################################################################
    if bInteractive:
        os.system("chmod 777 %s" % tmp_wrapper)
        cmd = "bash {0}".format(tmp_wrapper)
        print "Running interactively"
    else:
        os.system("chmod 777 %s" % wrapper_name)

        # Executables compiled on a specific OS are only guaranteed to work on
        # that OS on the grid.  Infer the OS from CMT.
        cmtConfig = os.environ["CMTCONFIG"]
        whichOS = cmtConfig.split("-")[1]
        whichOS = re.sub("slc", "sl", whichOS)

        # Don't let the user submit jobs with v21r1p1 to SL7 Singularity
        # containers.
        if version.group(0) == "v21r1p1" and whichOS != "sl6":
            print "MINERvA framework version {version} has not been compiled for OS {whichOS}, so I will not submit this job!".format(
                version=version.group(0), whichOS=whichOS
            )
            sys.exit(
                0
            )  # TODO: We should really be exiting with something nonzero everywhere that there is a problem.  Keeping with convention for now.

        cmd = ("jobsub_submit -e USER --group=minerva "
              "--cmtconfig={cmtConfig} "
              "--append_condor_requirements='(TARGET.HAS_SINGULARITY=?=true)' "
              "-l '+SingularityImage=\\\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-{whichOS}:latest\\\"' "
              "--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC "
              "--role=Analysis "
              "--memory {memory}MB --disk {disk}GB "
              "-f {outDirPlaylist}/{sub} "
              "-f /pnfs/minerva/resilient/tarballs/{tarball} "
              "-d MERGE {outDir} -d LOGS {outDirLogs} "
              "-r {release} {transfer_var_file} "
              "-i /cvmfs/minerva.opensciencegrid.org/minerva/software_releases/{release}/ "
              "{debug} "
              "file://{wrapper}".format(
                  cmtConfig=cmtConfig,
                  whichOS=whichOS,
                  memory=memory,
                  disk=disk,
                  outDir=outDir,
                  outDirLogs=outDirLogs,
                  outDirPlaylist=outDirPlaylist,
                  sub=subPlaylist,
                  tarball=tarballName,
                  wrapper=wrapper_name,
                  release=version.group(0),
                  transfer_var_file=(
                      "-f {0}/Variables.txt".format(outDirPlaylist) if variables else ""
                  ),
                  debug="-e IFDH_CP_MAXRETRIES" if bDebug else "",
              )
        )

    print cmd
    if not bTest:
        os.system(cmd)

    # Cleanup
    os.remove(tmp_wrapper)
    if bInteractive:
        ifdh_handle.cp([logname, "{0}/{1}.log".format(outDirLogs, mergedFileShortName)])
        os.remove(logname)

    # Don't pound the grid; spread out submissions slightly
    sleepTime = 2
    print "Sleeping for %i seconds\n" % sleepTime
    time.sleep(sleepTime)

##############################################################
# Get all the user inputs
##############################################################
def getArgs():
    parser = argparse.ArgumentParser(
        description="Submit directories of unmerged files to the grid in order to merge.  Currently, this only works with files on persistent"
    )
    # Set things up so required arguments first
    optional_args = parser._action_groups.pop()
    required_args = parser.add_argument_group("required arguments")
    parser._action_groups.append(optional_args)

    required_args.add_argument(
        "--tool", required=True, help="The name of the analysis tool"
    )
    required_args.add_argument(
        "--release",
        required=True,
        help="Minerva release (vXrYpZ<additional identifier>)",
    )
    required_args.add_argument(
        "--inputdir",
        required=True,
        help="Directory in /pnfs/minerva/persistent/users/$USER where the unmerged files live",
    )
    required_args.add_argument(
        "--outputdir",
        required=True,
        help="Directory to place merged files in /pnfs/minerva/persistent/users/$USER",
    )

    # Flags for data switches
    optional_args.add_argument(
        "--data", action="store_true", help="Unmerged files are from data"
    )
    optional_args.add_argument(
        "--mc", action="store_true", help="Unmerged files are from mc"
    )

    optional_args.add_argument(
        "--scratch",
        action="store_true",
        help="Change the base directory from /pnfs/minerva/persistent/users/$USER to /pnfs/minerva/scratch/users/$USER",
    )
    optional_args.add_argument(
        "--tree",
        default=None,
        help="Name of the analysis tree.  (Default: <name of analysis tool>)",
    )

    optional_args.add_argument(
        "--run",
        nargs="+",
        default=-1,
        help="Run numbers if you want to look at specific runs. Can input multiple runs (ex. --run 1000 1250)",
    )
    optional_args.add_argument(
        "--processing",
        default="central_value",
        help="Allows for merging when look at specific MC processings (ex. CCCohSignal) (default: central_value)",
    )
    optional_args.add_argument(


        "--base_dir",
        default="/minerva/app/users/$USER/cmtuser",
        help="This points to base directory for your Minerva release ( default: /minerva/app/users/$USER/cmtuser )",
    )
    optional_args.add_argument(
        "--release_dir",
        default="",
        help="Alternative method to set your release directory for your Minerva release",
    )
    optional_args.add_argument(
        "--user",
        default=userEnv,
        help="The name of the user who owns the input files (Default: $USER)",
    )
    optional_args.add_argument(
        "--memory",
        type=int,
        default=500,
        help="This is the amount of memory to use for your job in MB ( default: 500 MB )",
    )
    optional_args.add_argument(
        "--disk",
        type=int,
        default=10,
        help="This is the amount of disk limit to use for your job in GB ( default: 10 GB )",
    )
    optional_args.add_argument(
        "--test",
        action="store_true",
        help="Run this without submitting this to the grid",
    )
    optional_args.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="This raises some debug flags",
    )
    optional_args.add_argument(
        "--interactive",
        action="store_true",
        default=False,
        help="Run the job interactively",
    )
    optional_args.add_argument(
        "--tracker_only",
        action="store_true",
        default=False,
        help="Only save truth infomation of events with vertex in tracker(i.e. mc_vtx[2]>5891 && mc_vtx[2]<8439). "
        "Doesn't work with --with_var",
    )

    # Flags for using MergeToolWithVar
    optional_args.add_argument(
        "--with_var",
        action="store_true",
        default=False,
        help="Run MergeToolWithVar that support variables stripping."
        "Warning: Experimental feature. Please use 'another_merge_tool' branch of MergeToolWithVar",
    )

    optional_args.add_argument(
        "--variables",
        default=None,
        help="Path to Variable file for using in MergeToolWithVar",
    )

    args = parser.parse_args()
    return args

################################################################################
# submitMergeToGrid.py
################################################################################
if __name__ == "__main__":
    userEnv = os.environ["USER"]

    ############################################################################
    # Get args and do some arg validation/processing
    ############################################################################
    args = getArgs()
    if not (args.data or args.mc):
        print "submitMergeToGrid.py: error: --data or --mc is required"
        sys.exit(0)

    if args.data and args.mc:
        print "submitMergeToGrid.py: error: cannot have both --data or --mc"
        sys.exit(0)

    if args.data:
        dataSwitch = "data"
    elif args.mc:
        dataSwitch = "mc"

    if not args.tree:
        anaTree = args.tool
    else:
        anaTree = args.tree

    only_tracker_truth = args.tracker_only

    ############################################################################
    # Create a list of dirs that contain tuples - one dir for each run. The
    # tuples in each dir will be merged together, run-by-run. One job per
    # dir/run.
    # gatherUnmergedDirs starts in your inputTopDir and assumes data/mc
    # ProcessAna dir structure to find run dirs in the form XX/XX/XX/XX.  The
    # baseDir returned is like (e.g. for data):
    # <inputTopDir>/grid/minerva/ana/*/*/*/*/*/*
    ############################################################################
    inputTopDir = "/pnfs/minerva/{0}/users/{1}/{2}".format(
            "scratch" if args.scratch else "persistent", args.user, args.inputdir
        )
    tupleDirsToMerge, baseDir = gatherUnmergedDirs(
        dataSwitch,
        inputTopDir,
        args.processing,
        args.run,
    )
    if len(tupleDirsToMerge) == 0:
        print "No directories in {0}.  Exiting before creating tarball or output directory".format(
            baseDir
        )
        sys.exit(0)

    print "tupleDirsToMerge: "
    for _ in tupleDirsToMerge:
        print _

    ############################################################################
    # Create the output directories
    ############################################################################
    outDirRuns = ""
    if args.run != -1:
        outDirRuns = "/" + "_".join(args.run[:3])  # Giving a sensible title
        print "Looking at runs: ", args.run
    outDir = "/pnfs/minerva/{dcache}/users/{user}/{outdir}{runs}".format(
        dcache="scratch" if args.scratch else "persistent",
        user=userEnv,
        outdir=args.outputdir,
        runs=outDirRuns,
    )
    if not os.path.isdir(outDir):
        print "Making output directory {0}".format(outDir)
        os.system("mkdir -p %s" % outDir)
    else:
        print "{0} already exists! Exiting before repopulating".format(outDir)
        sys.exit(0)

    # Create output subdirectories
    tag = "merge_wrapper-{0}".format(args.tool)
    processingID = "%s_%s-%s00" % (
        tag,
        dt.date.today(),
        dt.datetime.today().strftime("%H"),
    )
    outDirPlaylist = "{0}/{1}_wrapper_playlist".format(outDir, processingID)
    os.system("mkdir -p %s" % outDirPlaylist)
    outDirLogs = "{0}/{1}_logs".format(outDir, processingID)
    os.system("mkdir -p %s" % outDirLogs)

    # Copy input files to output dir
    varfile = copyFile(args.variables, outDirPlaylist, "Variables.txt")

    ############################################################################
    # Create BlueArc staging area to which we can write the playlists and grid
    # scripts before we copy them to the grid.
    # TODO this strategy is no longer best practice - should create them
    # locally, then playlist is jammed in the tarball and grid script can be
    # passed to jobsub_submit directly.
    # Also: not cool side effect to write to user's BlueArc.
    ############################################################################
    blueArcDir = "/minerva/data/users/{user}/grid_merge".format(user=userEnv)
    if not os.path.isdir(blueArcDir):
        print "Making BlueArc directory {0}".format(blueArcDir)
        os.system("mkdir -p %s" % blueArcDir)


    ############################################################################
    # Create the tarball in the correct place: /pnfs/minerva/resilient/tarballs
    ############################################################################
    if len(args.release_dir) == 0:
        args.release_dir = args.release
    if args.interactive:
        tarballName = "interactive_notarball.tar.gz"
    else:
        tarballName = createTarball(args.base_dir, args.release_dir)

    ##############################################################
    # Loop each run's directory
    ##############################################################
    for unmergedDir in tupleDirsToMerge:
        mergedFileShortName = constructMergedFileShortName(
            args.tool, dataSwitch, unmergedDir
        ) # "<anatool>_<data/mc>_AnaTuple_runXXXXXXXX"

        # Search the dir for tuples, create an xrootd playlist of them, and
        # copy it to the outdir. 
        subPlaylist = mapUnmergedDirToSubPlaylist(
            unmergedDir, mergedFileShortName, blueArcDir, outDirPlaylist, False
        ) # "playlist_<mergedFileShortName>.txt"

        submitJob(
            outDir,
            blueArcDir,
            outDirPlaylist,
            outDirLogs,
            mergedFileShortName,
            subPlaylist,
            dataSwitch,
            args,
            tarballName,
            varfile,
        )
        os.remove("{0}/{1}".format(blueArcDir, subPlaylist))
