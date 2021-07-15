
#ifndef TREEWRAPPER_H
#define TREEWRAPPER_H

#include "TTree.h"
#include "TLeaf.h"
#include "TBranch.h"
#include <map>
#include <string>
#include <vector>
#include <cassert>


struct LeafAndBranch
{
  TLeaf* leaf;
  TBranch* branch;
};

/**
   Typical use:

   TTree *t = (TTree*)somefile.Get("sometree");
   TreeWrapper w(t);
   double energy=w.GetValue("evt.energy", 1);

*/

namespace MAT {

  //Gory details: type traits to map c++ type back to TLeaf type names
  namespace detail
  {
    //Base behavior: don't compile because this type doesn't have a name!
    template <class T>
    struct typeName
    {
    };

    //Specializations for types both MinervaAnalysisTool and ROOT know about.
    //Pre-c++11 (i.e. interpreted ROOT 5 scripts), I have to define these at the top of the .cxx because of https://stackoverflow.com/questions/1639154/how-to-declare-a-static-const-char-in-your-header-file
    template <>
    struct typeName<int>
    {
      static const char* name;
    };

    template <>
    struct typeName<double>
    {
      static const char* name;
    };

    template <>
    struct typeName<bool>
    {
      static const char* name;
    };

    template <>
    struct typeName<std::vector<std::vector<double> > > // Need this for 2D vector getter.
    {
      static const char* name;
    };
  }

    class TreeWrapper : public TObject {
      public:
            /**
               Create a TreeWrapper for TTree t. TChain derives from TTree so
               can pass a TChain here too. On second thoughts that might not
               work. Who knows?
            */
        TreeWrapper(TTree* t);
        
        virtual ~TreeWrapper() {}
            /**
               Register a branch that you want to use, a la
               TTree::SetBranchAddress. All you need here is the branch name
               though - no pointers to pointers etc
               
               It is now unneccessary to call this function. It'll happen automagically
               on first call of GetValue.
            */
        virtual bool AddBranch(const std::string& branchName);
            /**
               Get the value of branch \a branchName in tree entry \a ientry.
               Use non-zero \a leafVal for leaves that are arrays
            */

        virtual double GetValue(const std::string& branchName, Long64_t ientry, int leafVal=0)
        {
          return (double) GetLeaf(branchName,ientry)->GetValue(leafVal);
        }

        virtual int GetInt(const std::string& branchName, Long64_t ientry, int leafVal=0)
        {
          return (int) GetLeaf(branchName,ientry)->GetValue(leafVal); 
        }

        virtual bool GetBool(const std::string& branchName, Long64_t ientry, int leafVal=0)
        {
          return (bool) GetLeaf(branchName,ientry)->GetValue(leafVal);
        }
            /**
               Number of entries in the tree
            */
        virtual Long64_t GetEntries() const { return tree->GetEntries(); }
        
        virtual TTree* GetTree() const { return tree; }
        
        

        template<class T>
        T* GetValueArray(const std::string& branchName, Long64_t ientry)
        {
          return GetValueVector<T>(branchName, ientry).data();
        }

        template<class T>
        std::vector<T> GetValueVector(const std::string& branchName, Long64_t ientry)
        {
          TLeaf* thisleaf = GetLeaf(branchName,ientry);
          //TODO: Sort branches by type instead so I can afford to throw an exception?  Users already pay for the lookup, and I'd probably be making it faster.
          assert(!strcmp(detail::typeName<T>::name, thisleaf->GetTypeName()) && "Stored vector branch type does not match type user requested!  "
                                                                               "You were getting junk here before I saved you.");
          int N=thisleaf->GetLen();
          T* d=(T*)thisleaf->GetValuePointer();
          std::vector<T> ret(d, d+N);
          return ret;
        }

        template<class T>
        T GetValueNDVector(const std::string& branchName, Long64_t ientry)
        {
          TLeaf* thisleaf = GetLeaf(branchName,ientry);
          //TODO: Sort branches by type instead so I can afford to throw an exception?  Users already pay for the lookup, and I'd probably be making it faster.
          assert(!strcmp(detail::typeName<T>::name, thisleaf->GetTypeName()) && "Stored vector branch type does not match type user requested!  "
                                                                               "You were getting junk here before I saved you.");
          return T(*((T*)thisleaf->GetValuePointer()));
        }

        virtual void PrintKnownBranchNames();

          /**
              Expert Interface. Don't use this unless you really need the TLeaf rather than its value.
          */
        virtual TLeaf* GetLeaf(const std::string& branchName, Long64_t ientry);

      protected:
        TTree* tree;
        
        std::map<std::string, LeafAndBranch> leavesAndBranches;

        typedef std::map<std::string, LeafAndBranch>::iterator itLaB;

            // Are we wrapping a TChain?
        bool wrappingChain;
            // If we're on a TChain, the tree number of the last entry we read
        int lastReadTree;
            // The entry offset for the current tree
        int currentOffset;
            /**
               Reset the branch addresses because we're on a new TTree (if wrapping a chain)
            */
        virtual Bool_t Notify();

            // Find the tree number containing entry
        int GetTree(Long64_t entry, int guess = -1);
            // Set the branch and leaf addresses to the ones appropriate for tree number treeNum
        bool SetBranchAddresses();
            // Get the entry offset for tree treeNum
        int GetOffset(int treeNum);
        
        ClassDef(TreeWrapper, 0);

            // Stupid default constructor forced by ROOT.
        TreeWrapper() {}


    };

}

#endif
