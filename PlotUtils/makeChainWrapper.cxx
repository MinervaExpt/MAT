#include <iostream>

#include "makeChainWrapper.h"


MAT::ChainWrapper* makeChainWrapperPtr(const std::string& playlist,
                                             const std::string& name)
{
    MAT::ChainWrapper* chw = new MAT::ChainWrapper(name.c_str());
    int nfiles = chw->Add(playlist.c_str());
    
    std::cout << "Added " << nfiles << " files from playlist " << playlist << std::endl;
    
    return chw;
}


MAT::ChainWrapper& makeChainWrapper(const std::string& playlist,
                                          const std::string& name)
{
    MAT::ChainWrapper* chw = makeChainWrapperPtr(playlist, name); 
    
    return *chw;
}
