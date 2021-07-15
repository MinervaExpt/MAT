#ifndef CHAINWRAPPER_H
#define CHAINWRAPPER_H

#include "TreeWrapper.h"
#include "TChain.h"

#ifdef PLOTUTILS_THROW_EXCEPTIONS
#include <stdexcept>
#endif

namespace MAT {

    class ChainWrapper : public TreeWrapper
    {
      public:
        
        ChainWrapper(const char* name);
        
        virtual ~ChainWrapper() { delete tree; }
        
        int AddPlaylist(const char* file_path);

        int AddFiles(const char* file_path);

        int Add(const std::string file_path);

        TChain* GetChain() { return (TChain*)tree; }

        const std::vector<std::string>& GetListOfFiles() const {return fListOfFiles; }

        #ifdef PLOTUTILS_THROW_EXCEPTIONS
          class BadFile: public std::runtime_error
          {
            public:
              BadFile(const std::string& fileName);
              virtual ~BadFile() throw() {}

              const std::string file;
          };
        #endif

        ClassDef(ChainWrapper, 0);

      private:
        std::vector<std::string> fListOfFiles;

        ChainWrapper() {}
        
    };
    
}
#endif
