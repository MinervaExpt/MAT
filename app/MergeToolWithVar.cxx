#include "TChain.h"
#include "TCollection.h"
#include "TFile.h"
#include "TFileMerger.h"
#include "TKey.h"
#include "TSystem.h"
#include "TStopwatch.h"

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

//#include <cstdlib>
#include <ctime>      //provide time stamp for default file name
#include <sys/stat.h> //make directory
#include <unistd.h>   // getopts

int MAX_OPENED_INPUTFILE = 5;

// Copied from MergeTool.cxx
bool isGoodFile(const char *const file) {
  TFile *f = TFile::Open(file);
  if (!f) {
    return false;
  } else if (f->IsZombie()) {
    delete f;
    return false;
  } else {
    bool good = true;
    TTree *meta = (TTree *)f->Get("Meta");
    good = good && meta;
    good = good && meta->GetEntriesFast() != 0;
    good = good && meta->GetBranch("POT_Total");
    good = good && meta->GetBranch("POT_Used");
    delete f;
    return good;
  }
}

/*
bool quickCheckConsistency(TTree *const new_tree, TTree *const old_tree) {
  //If there is no branch need for a tree, don't check.
  if (new_tree->GetEntries() != 0) {
    std::cout << "Checking Number of Entries." << std::endl;
    long int Nold = old_tree->GetEntries();
    delete old_tree; // Disconnect the new tree to old tree.
    long int Nnew = new_tree->GetEntries();
    std::cout << "new: " << Nnew << ", and old: " << Nold << std::endl;
    if (Nnew != Nold) {
      return false;
    }
  }
  // std::cout << "Checking Address of first entry" <<std::endl;
  // std::cout << "new: " << new->getEntries() << ", and old: " <<
  // old->getEntries() <<std::endl;
  return true;
}
*/

// TODO Print a better formated help info
void printHelp() {
  std::cout
      << "Usage: $MergeToolWithVar [OPTIONS]... Playlist.txt Variables.txt "
      << std::endl;
  std::cout << "Playlist.txt: list of anatuples files you want to merge, one "
               "file per line."
            << std::endl;
  std::cout << "Variables.txt: list of anatuple variables you want to save."
            << std::endl;
  std::cout << "Format of Variables.txt: specify name of TTree by ending with "
               "a ':', followed by names of branch variables, seperated by "
               "space or newline ."
            << std::endl;
  std::cout << "Here is an example:" << std::endl;
  std::cout << "\tTTree_1: TBranch_1_of_TTree_1 TBranch_2_of_TTree_1 ... "
            << std::endl;
  std::cout << "\tTTree_2: TBranch_1_of_TTree_2 " << std::endl;
  std::cout << "\tTBranch_2_of_TTree_2 " << std::endl;
  std::cout << "Check Variables.txt for a template." << std::endl;
  std::cout << "[OPTIONS]: " << std::endl;
  std::cout << "-h          : Print this message." << std::endl;
  std::cout << "-o [Path]   : Specify ABSOLUTE path of output directory or "
               "output file name, or both."
            << std::endl;

  std::cout << "-g          : Run this tool on grid. " << std::endl;
  std::cout << "Warning: You shouldn't use this option directly. Please use "
    "python script to submit merge job to grid." <<std::endl;
  std::cout << "-p [Path]   : Specify path of Playlist.txt file. " <<std::endl;
  std::cout << " Please omit the Playlist.txt argument if you used this option." <<std::endl;
  std::cout << "-v [Path]   : Specify path of Variable.txt file. " <<std::endl;
  std::cout << " Please omit the Variable.txt argument if you used this option." <<std::endl;
  std::cout << "-n [Number] : Specify number of input files to be merged at the same time." << std::endl;
  std::cout << " Default is " << MAX_OPENED_INPUTFILE << ". Set to 1 will merge files one by one" <<std::endl;
}

/*
bool copyOneFile(const std::map<std::string,TTree*>& trees, const char* const
rootfile ,TFile * fout) { TFile* f = TFile::Open(rootfile); if (!f ||
f->IsZombie()) { std::cerr << "Failed to open file: " << rootfile << std::endl;
    return false;
  }
  for (auto it = trees.begin(); it != trees.end(); ++it) {
    TTree *toBeCopied = (TTree*) f->Get(it->first.c_str());
    if (!toBeCopied) {
      std::cerr << "TTree " << it->first << "not found in this file." <<
std::endl; return false; } else {
      //Connect branchaddress before copy entries, which is assumed for
CopyEntries(). it->second->CopyAddresses(toBeCopied); if
(it->second->CopyEntries(toBeCopied)>0) { toBeCopied->ResetBranchAddresses();
//Disconnect two tree. delete toBeCopied; } else { std::cerr << "TTree not
copied!" << std::endl; return false;
      }
    }
    fout=it->second->GetCurrentFile();
  }
  f->Close();
  delete f;
  return true;
}
*/
int copyStruct(std::map<std::string, TTree *> &ttrees, TTree *const tmp) {
  if (!tmp)
    return 1;
  ttrees.insert(std::make_pair(std::string(tmp->GetName()), tmp->CloneTree(0)));
  delete tmp;
  return 0;
}

int initializeTTrees(std::map<std::string, TTree *> &ttrees,
                     const char *const rootfile, const char *const vars,
                     TFile *fout) {
  std::ifstream branch_names(vars);
  if (!branch_names.is_open()) {
    std::cerr << "Variables file not found: " << vars << std::endl;
    exit(1);
  }
  TFile *tree_struct_file = TFile::Open(rootfile);
  if (!tree_struct_file || tree_struct_file->IsZombie()) {
    std::cerr << "Failed to read first file, exiting." << std::endl;
    exit(1);
  }
  fout->cd();
  TTree *tmp = NULL;
  std::string name;

  while (1) {
    branch_names >> name;
    if (branch_names.eof()) {
      copyStruct(ttrees, tmp);
      break;
    }
    if (name.find(":") != std::string::npos) {
      std::string tree_name = name.substr(0, name.find(":"));
      // Moving from one TTree to another. Assuming tree status is set, make a
      // copy of tree structure and insert to map.
      if (tmp) {
        copyStruct(ttrees, tmp);
      }
      tmp = (TTree *)tree_struct_file->Get(tree_name.c_str());
      if (!tmp) {
        std::cerr << "TTree: " << tree_name << " not found." << std::endl;
        exit(1);
      }
      tmp->SetBranchStatus("*", 0);
    } else if (tmp) {
      // Do we want to check if branch exist or not?
      // It will prevent user from wildcarding.
      // if (!it->second->GetBranch(name.c_str())) {
      //  std::cerr << "Rrequested TBranch not found: " << name << " in TTree: "
      //  << it->first << std::endl; exit(1);
      //}
      tmp->SetBranchStatus(name.c_str(), 1);
    } else {
      std::cerr << "Please specify TTree name by ending with a ':' before "
                   "entering TBranch name: "
                << name << std::endl;
      exit(1);
    }
  }
  tree_struct_file->Close();
  return 0;
}

int merge(const std::string &playlist, const std::string &vars,
          const std::string &merged_file) {
  std::cout << "Enter merge()" << std::endl;
  std::ifstream playliststream(playlist.c_str());
  if (!playliststream.is_open()) {
    std::cerr << "Playlist file not found: " << playlist << std::endl;
    exit(1);
  }
  TStopwatch TS;
  bool first = true;
  TFileMerger TFM(false);
  TFM.SetMaxOpenedFiles(2);
  TFM.SetPrintLevel(1);
  TFM.OutputFile(merged_file.c_str());
  int EPartialMergeOption = TFileMerger::kAll | TFileMerger::kIncremental;
  auto fout = TFM.GetOutputFile();
  while (1) {
    std::string rootfile;
    playliststream >> rootfile;
    if (playliststream.eof())
      break;
    // Do we want to make sure file is good before adding?
    if (!isGoodFile(rootfile.c_str())) {
      std::cerr << "File: " << rootfile << " is not a good file, skipping."
                << std::endl;
      continue;
    }
    if (first && vars.size() != 0) {
      std::map<std::string, TTree *> TTrees;
      initializeTTrees(TTrees, rootfile.c_str(), vars.c_str(), fout);
      for (auto it : TTrees) {
        TFM.AddObjectNames(it.first.c_str());
      }
      EPartialMergeOption = EPartialMergeOption | TFileMerger::kOnlyListed;
    }
    TFM.AddFile(rootfile.c_str());
    first = first && false;
  }
  TS.Stop();
  std::cout << "Time spent on checking input files and initialize output file." <<std::endl;
  TS.Print();
  TS.Continue();
  TFM.SetMaxOpenedFiles(1+MAX_OPENED_INPUTFILE);
  std::cout << "Start merging." << std::endl;
  bool r = TFM.PartialMerge(EPartialMergeOption);
  TS.Stop();
  std::cout << "The merged anatuple is located at: " << merged_file
            << std::endl;
  std::cout << "Total time spent:" <<std::endl;
  TS.Print();
  return r;
}

std::string verifyOutputPath(const char *const output, bool skip_checking) {
  if (skip_checking) {
    std::cout <<output <<std::endl;
    return std::string(output);
  }
  // Attach a time stamp to default file name in order to prevent accidental
  // overwrite.
  // TODO Better ideas on naming?
  const int MAX_FILE_NAME = 100;
  char default_merged_file_name[MAX_FILE_NAME];
  sprintf(default_merged_file_name, "%s%u.root", "merged-", time(NULL));

  // const char * const default_merged_file_name = "merged.root";
  const char *const blue_arc_dir = "/minerva/data/users/";
  const char *const dCache_persistent_dir = "/pnfs/minerva/persistent/users/";
  const char *const default_folder_name = "/merged_files";
  const char *const default_dir = blue_arc_dir;
  const char *const dCache_xroot_head =
      "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr";
  bool xrootd = false;

  // Extract output directory and file name from user input.
  // TODO is there better way to do this?
  std::stringstream ss;
  if (!output) {
    // return default output path
    ss << default_dir << std::getenv("USER") << default_folder_name << "/"
       << default_merged_file_name;
  } else {
    std::string tmpout(output);
    if (tmpout[0] == '/') {
      // check if the output path given by user is in persistent or data
      if (tmpout.find(blue_arc_dir) == 0) {
        std::cout << "Saving to blueArc. Warning: it won't work on grid. "
                  << std::endl;
      } else if (tmpout.find(dCache_persistent_dir) == 0) {
        std::cout << "Saving to dCache. Maybe change path to url to use xrootd?"
                  << std::endl;
        xrootd = true;
      } else {
        std::cerr << "Warning: Please save merged file in dCache persistent or "
                     "BlueArc data. "
                  << std::endl;
        exit(1);
      }
      if (tmpout.back() == '/') {
        ss << tmpout << default_merged_file_name;
      } else {
        if (tmpout.rfind(".root") == std::string::npos) {
          ss << tmpout << "/" << default_merged_file_name;
        } else {
          ss << tmpout;
        }
      }
    } else {
      // User only give me file name, use default path:
      ss << default_dir << std::getenv("USER") << default_folder_name << "/"
         << tmpout;
    }
  }

  // Try creating directory if not already exists.
  // TODO This only create one subdirectory, not a tree of directories. may want
  // to mkdir recursively?
  std::string final_output = ss.str();
  const std::string directory = final_output.substr(0, final_output.rfind('/'));
  struct stat info;
  if (stat(directory.c_str(), &info) != 0 &&
      mkdir(directory.c_str(),
            S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH)) {
    std::cerr << "The output directory: " << directory
              << " doesn't exist and can't be created. " << std::endl;
    exit(1);
  }
  if (stat(final_output.c_str(), &info) == 0) {
    std::cout << "Warning: The output file: " << final_output
              << " Already exists, will be overwrited." << std::endl;
  }

  return xrootd ? final_output.replace(0, 5, dCache_xroot_head) : final_output;
}

int main(int argc, char **argv) {
  // TODO Add more options, like MergeTool.cxx provides
  if (argc < 3) {
    printHelp();
    return 1;
  }
  // Added this because root complains about missing library. But this only
  // shows up after I switched branch. Don't know why.
  gSystem->Load("libTree");
  std::string playlist_file;
  std::string variable_file;

  // Parsing options
  char opt;
  const char *output = NULL;
  bool grid = false;
  while ((opt = getopt(argc, argv, "gho:v:p:n:")) != -1) {
    switch (opt) {
    case 'h':
      printHelp();
      exit(0);
      break;
    case 'o':
      output = optarg;
      break;
    case 'g':
      grid = true;
      break;
    case 'v':
      variable_file = optarg;
      break;
    case 'p':
      playlist_file = optarg;
      break;
    case 'n':
      MAX_OPENED_INPUTFILE = atoi(optarg);
      break;
    default:
      std::cerr << "Unknown option: " << opt << std::endl;
      exit(1);
    }
  }

  while (argc - optind != 0) {
    if (playlist_file.size() == 0) {
      playlist_file = argv[optind];
      ++optind;
    } else if (variable_file.size() == 0) {
      variable_file = argv[optind];
      ++optind;
    } else {
      std::cerr << "Undefined behavior for non-option argument: "
                << argv[optind] << std::endl;
      std::cerr << "May be you have tell me input file by option?" << std::endl;
      ++optind;
    }
  }

  if (playlist_file.size() == 0) {
    std::cerr << "Please tell me your playlist." << std::endl;
  }

  if (variable_file.size() == 0) {
    std::cout << "No variable file, will merge all branches." << std::endl;
  }

  auto output_file = verifyOutputPath(output, grid);
  return merge(playlist_file.c_str(), variable_file.c_str(),
               output_file.c_str())? 0 : 1;
}
