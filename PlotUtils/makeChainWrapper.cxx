#include <iostream>

#include "PlotUtils/makeChainWrapper.h"


PlotUtils::ChainWrapper* makeChainWrapperPtr(const std::string& playlist,
                                             const std::string& name)
{
    PlotUtils::ChainWrapper* chw = new PlotUtils::ChainWrapper(name.c_str());
    int nfiles = chw->Add(playlist.c_str());
    
    std::cout << "Added " << nfiles << " files from playlist " << playlist << std::endl;
    
    return chw;
}


PlotUtils::ChainWrapper& makeChainWrapper(const std::string& playlist,
                                          const std::string& name)
{
    PlotUtils::ChainWrapper* chw = makeChainWrapperPtr(playlist, name); 
    
    return *chw;
}
