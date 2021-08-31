#ifndef MacroUtil_cxx
#define MacroUtil_cxx

//PlotUtils includes
#include "PlotUtils/MacroUtil.h"
#include "PlotUtils/makeChainWrapper.h"

//ROOT includes
#include "TH1.h"

//c++ includes
#include <iostream>
#include <cassert>

using namespace PlotUtils;

// Accumulate all POT in a "playlist" file
double CountPOT(const std::string& fileName)
{
  PlotUtils::ChainWrapper meta("Meta");
  meta.Add(fileName);

  const int nEntries = meta.GetEntries();
  assert((nEntries >= meta.GetChain()->GetListOfFiles()->GetEntriesFast())
         && "Each AnaTuple file should have exactly at least 1 entry for POT!");

  double pot_used = 0;
  for(int entry = 0; entry < nEntries; ++entry) pot_used += meta.GetValue("POT_Used", entry);

  return pot_used;
}

// Data tree only
MacroUtil::MacroUtil(const std::string& reco_tree_name, const std::string& data_file_list,
                     const std::string& plist_name, const bool is_grid): m_data(makeChainWrapperPtr(data_file_list, reco_tree_name)),
                                                                         m_mc(new ChainWrapper(reco_tree_name.c_str())),
                                                                         m_truth(new ChainWrapper("Truth")),
                                                                         m_data_pot(CountPOT(data_file_list)),
                                                                         m_mc_pot(0)
{
  CommonInitialization(plist_name, is_grid);
}

// MC reco tree.  Choose whether the Truth tree also is loaded
MacroUtil::MacroUtil(const std::string& reco_tree_name, const std::string& mc_file_list,
                     const std::string& plist_name, const bool wantsTruth, const bool is_grid)
                    : m_data(new ChainWrapper(reco_tree_name.c_str())),
                      m_mc(makeChainWrapperPtr(mc_file_list, reco_tree_name)),
                      m_truth(wantsTruth?makeChainWrapperPtr(mc_file_list, "Truth"):(new ChainWrapper("Truth"))),
                      m_data_pot(0), m_mc_pot(CountPOT(mc_file_list))
{
  CommonInitialization(plist_name, is_grid);
}

//Data, MC reco, and Truth trees
MacroUtil::MacroUtil(const std::string& reco_tree_name, const std::string& mc_file_list,
                     const std::string& data_file_list, const std::string& plist_name,
                     const bool wantsTruth, const bool is_grid): m_data(makeChainWrapperPtr(data_file_list, reco_tree_name)),
                                                                 m_mc(makeChainWrapperPtr(mc_file_list, reco_tree_name)),
                                                                 m_truth(wantsTruth?makeChainWrapperPtr(mc_file_list, "Truth"):(new ChainWrapper("Truth"))),
                                                                 m_data_pot(CountPOT(data_file_list)),
                                                                 m_mc_pot(CountPOT(mc_file_list))
{
  CommonInitialization(plist_name, is_grid);
}

// Delegating constructor to reduce code duplication
void MacroUtil::CommonInitialization(const std::string& plist_name, const bool is_grid)
{
  m_plist_string = plist_name;
  m_is_grid = is_grid;
  TH1::AddDirectory(false);
}

// Print this macro's configuration
void MacroUtil::PrintMacroConfiguration(std::string macro_name) {
  std::cout << "\nMacroUtil configuration of this macro " << macro_name
            << "\n** Number of data files is "
            << m_data->GetChain()->GetListOfFiles()->GetEntriesFast()
            << "\n** Number of MC reco files is "
            << m_mc->GetChain()->GetListOfFiles()->GetEntriesFast()
            << "\n** Number of Truth files is "
            << m_truth->GetChain()->GetListOfFiles()->GetEntriesFast()
            << "\n** Data POT extracted from anatuples is " << m_data_pot
            << "\n** MC POT extracted from anatuples is " << m_mc_pot
            << "\n** Playlist string is " << m_plist_string
            << "\n** Grid is "
            << std::boolalpha << m_is_grid << "\n\n";
}

#endif
