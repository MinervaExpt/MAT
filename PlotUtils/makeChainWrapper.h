#ifndef make_chain_wrapper_h
#define make_chain_wrapper_h

#include <string>

#include "PlotUtils/ChainWrapper.h"

PlotUtils::ChainWrapper* makeChainWrapperPtr(const std::string& playlist,
                                             const std::string& name);

PlotUtils::ChainWrapper& makeChainWrapper(const std::string& playlist,
                                          const std::string& name);

#endif
