#include <iostream>
#include <fstream>
#include <cstdlib>
#include "glob.h"

#include "ChainWrapper.h"
#include "ROOTglob.h"

#include "TSystem.h"

//Stuff to let me handle errors differently based on compiler flags

#ifdef PLOTUTILS_THROW_EXCEPTIONS
namespace MAT
{
  ChainWrapper::BadFile::BadFile(const std::string& fileName): std::runtime_error("No such file or directory: " + fileName), file(fileName)
  {
  }
}

#define FATAL(fileName)\
  throw MAT::ChainWrapper::BadFile(fileName);

#else //if !PLOTUTILS_THROW_EXCEPTIONS
#define FATAL(fileName)\
  std::cerr << "File not found: " << fileName << std::endl;\
  exit(1); //TODO: Don't exit() ever!

#endif //PLOTUTILS_THROW_EXCEPTIONS

MAT::ChainWrapper::ChainWrapper(const char* name)
  : TreeWrapper(new TChain(name))
{
  wrappingChain=true;
}

//===========================================================================

int MAT::ChainWrapper::AddFiles(const char* name)
{
  // Urgh, downcast
  //TChain* ch=(TChain*)tree;
  TChain* ch= dynamic_cast<TChain*>(tree);

  //If name seems to be an xrootd URL
  if(std::string(name).find("root:") != std::string::npos)
  {
    const std::vector<std::string> fileNames = MAT::glob(name, *gSystem);
    if(fileNames.empty())
    {
      FATAL(name);
    }

    for(std::vector<std::string>::const_iterator file = fileNames.begin();
        file != fileNames.end(); ++file)
    {
      ch->Add(file->c_str());
      fListOfFiles.push_back(*file);
    }

    return fileNames.size();
  }

  //else if not xrootd URL:
  glob_t g;
  glob(name, 0, 0, &g);
  //TODO: Check errno()?

  if ((int)g.gl_pathc == 0){
    std::cout << " ChainWrapper is unhappy with "<< name << std::endl;
   // FATAL(name)
  }

  for(int i=0; i<(int)g.gl_pathc; ++i){
    const char* filename=g.gl_pathv[i];
    ch->Add(filename);
    fListOfFiles.push_back(filename);
  }
  int ret=g.gl_pathc;

  globfree(&g);

  return ret;
}

//===========================================================================

int MAT::ChainWrapper::AddPlaylist(const char* name)
{
  // Urgh, downcast
  //TChain* ch=(TChain*)tree;
  TChain* ch= dynamic_cast<TChain*>(tree);

  std::ifstream input_stream(name);

  if (!input_stream.is_open())
  {
    FATAL(name);
  }
  
  int nfiles = 0;

  while (!input_stream.eof()) {

      std::string item;
      input_stream >> item;

      if (input_stream.eof()) break;

      std::cout << "\tadding file: " << item << std::endl;
      fListOfFiles.push_back(item);
      nfiles += ch->Add(item.c_str());
  }

  return nfiles;
}

//===========================================================================

int MAT::ChainWrapper::Add(std::string name){
  if (name.find(".root") != std::string::npos) 
    return MAT::ChainWrapper::AddFiles(name.c_str());
  else return MAT::ChainWrapper::AddPlaylist(name.c_str());
}
