#ifndef MacroUtil_h
#define MacroUtil_h

//PlotUtils includes
#include "PlotUtils/ChainWrapper.h"

//c++ includes
#include <string>

namespace PlotUtils
{
  class MacroUtil {
    public:
      // Data tree only
      MacroUtil(const std::string& reco_tree_name, const std::string& data_file_list,
                const std::string& plist_name, const bool is_grid);
  
      // MC reco tree.  Choose whether the Truth tree also is loaded
      MacroUtil(const std::string& reco_tree_name, const std::string& mc_file_list,
                const std::string& plist_name, const bool wantsTruth, const bool is_grid);
  
      //Data, MC reco, and Truth trees
      MacroUtil(const std::string& reco_tree_name, const std::string& mc_file_list,
                const std::string& data_file_list, const std::string& plist_name,
                const bool wantsTruth, const bool is_grid);
  
      // All of these ChainWrapper pointers will be non-NULL.  Some of them
      // just might be empty.
      PlotUtils::ChainWrapper* m_data;
      PlotUtils::ChainWrapper* m_mc;
      PlotUtils::ChainWrapper* m_truth;
    
      // Program conditions
      std::string m_plist_string;
      bool m_is_grid;
    
      // POT counting
      double m_data_pot;
      double m_mc_pot;
    
      // Getters
      long int GetDataEntries() const { return m_data->GetEntries(); }
      long int GetMCEntries() const { return m_mc->GetEntries(); }
      long int GetTruthEntries() const { return m_truth->GetEntries(); }
    
      // Print this macro's configuration
      virtual void PrintMacroConfiguration(std::string macro_name = "");
    
    private:
      MacroUtil();
  
      // Delegating some constructor work to reduce code duplication
      void CommonInitialization(const std::string& plist_name, const bool is_grid);
  };
} // namespace PlotUtils

#endif  // MacroUtil_h
