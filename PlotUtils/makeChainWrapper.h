#ifndef make_chain_wrapper_h
#define make_chain_wrapper_h

#include <string>

#include "ChainWrapper.h"

MAT::ChainWrapper* makeChainWrapperPtr(const std::string& playlist,
                                             const std::string& name);

MAT::ChainWrapper& makeChainWrapper(const std::string& playlist,
                                          const std::string& name);

#endif
