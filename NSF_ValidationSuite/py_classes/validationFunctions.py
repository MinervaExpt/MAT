import ROOT,os

def passCuts( chain, treeName , entry ):

  if treeName == 'CCQENu': return passCuts_CCQENu(chain,entry)
  else: print "I'm inside passCuts() and I don't recognize the toolName you passed in!"

def passCuts_CCQENu( chain , entry ):

  if not chain.GetValue("has_interaction_vertex",entry) == 1:                               return False 
  if not chain.GetValue("phys_n_dead_discr_pair_upstream_prim_track_proj",entry) <= 1:      return False
  if not chain.GetValue("muon_theta",entry) <= 0.349:                                       return False 
  if ROOT.TMath.IsNaN(chain.GetValue("muon_theta",entry)):                                  return False
  if not chain.GetValue("CCQENu_nuHelicity",entry) == 1:                                    return False 
  if not chain.GetValue("mc_current",entry) == 1:                                           return False
  if not chain.GetValue("mc_incoming",entry) == 14:                                         return False

  return True

def passCuts_NukeCC( chain , entry ):

  if not chain.GetValue("muon_theta",entry) <= 0.296706:                                    return False
  if not chain.GetValue("NukeCC_nuHelicity",entry) == 1:                                    return False 
  #if not chain.GetValue("mc_current",entry) == 1:                                           return False
  #if not chain.GetValue("mc_incoming",entry) == 14:                                         return False

  return True

def passCuts( chain , toolName , entry ):
  
  if toolName == "CCQENu": return passCuts_CCQENu(chain,entry)
  elif toolName == "NukeCC": return passCuts_NukeCC(chain,entry)
  else:
      print "I don't know how to apply cuts for your chosen toolName, so I'm going to exit gracefully"
      sys.exit()

def copyMetaTreeToOutput( OPTS_VEC , tuplePath , outputFile ):

  inputFile = ROOT.TFile( tuplePath )
 
  try:
    metaTree = inputFile.Get( "Meta" ).CloneTree()
  except:
    print "I couldn't find a Meta tree in this file. What's up with that? Not going to try to write it to output..."
    return

  outputFile.cd()
  metaTree.Write()
 
def constructOutputFilePath( tuplePath , outDir , toolName ):

  fileName = tuplePath.split('/')[-1]
  fileNameComponents = fileName.split('_')
  dataSwitch = fileNameComponents[1]
  playlistName = fileNameComponents[2]

  print 'fileNameComponents: ' , fileNameComponents
  
  if toolName == 'NukeCC' or toolName == 'MECAna':
    subPlaylistIt = fileNameComponents[-1].split('.')[0] + '_' + fileNameComponents[-3] + '_' + fileNameComponents[-2]
  elif toolName == 'CCQENu':
    subPlaylistIt = fileNameComponents[-1].split('.')[0] 
 
  outPath = "{0}/validationHists_{1}_{2}_{3}.root".format(outDir,dataSwitch,playlistName,subPlaylistIt)

  # Create outDir if it doesn't exist
  if not os.path.isdir(outDir):
    print "Making plot directory {0}".format(outDir)
    os.system( "mkdir %s" % outDir )
  
 
  return outPath


