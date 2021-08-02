//macro to add histogram files
//NOTE: This macro is kept for back compatibility only.
//Use instead the executable $ROOTSYS/bin/hadd
//
//This macro will add histograms from a list of root files and write them
//to a target root file. The target file is newly created and must not be
//identical to one of the source files.
//
//Author: Sven A. Schmidt, sven.schmidt@cern.ch
//Date:   13.2.2001

//This code is based on the hadd.C example by Rene Brun and Dirk Geppert,
//which had a problem with directories more than one level deep.
//(see macro hadd_old.C for this previous implementation).
//
//The macro from Sven has been enhanced by
//   Anne-Sylvie Nicollerat <Anne-Sylvie.Nicollerat@cern.ch>
// to automatically add Trees (via a chain of trees).
//
//To use this macro, modify the file names in function hadd.
//
//NB: This macro is provided as a tutorial.
//    Use $ROOTSYS/bin/hadd to merge many histogram files



#include <string.h>
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"

#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif

#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>

TList *FileList;
TFile *Target;

void MergeRootfile( TDirectory *target, TList *sourcelist );


void mnvadd() {
   // in an interactive ROOT session, edit the file names
   // Target and FileList, then
   // root > .L mnvadd.cxx
   // root > mnvadd()

   Target = TFile::Open( "result.root", "RECREATE" );

   FileList = new TList();
   FileList->Add( TFile::Open("hsimple1.root") );
   FileList->Add( TFile::Open("hsimple2.root") );

   MergeRootfile( Target, FileList );

}

void MergeRootfile( TDirectory *target, TList *sourcelist ) {

   //  cout << "Target path: " << target->GetPath() << endl;
   TString path( (char*)strstr( target->GetPath(), ":" ) );
   path.Remove( 0, 2 );

   TFile *first_source = (TFile*)sourcelist->First();
   first_source->cd( path );
   TDirectory *current_sourcedir = gDirectory;
   //gain time, do not add the objects in the list in memory
   Bool_t status = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   // loop over all keys in this directory
   TChain *globChain = 0;
   TIter nextkey( current_sourcedir->GetListOfKeys() );
   TKey *key, *oldkey=0;
   while ( (key = (TKey*)nextkey())) {

      //keep only the highest cycle number for each key
      if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;

      // read object from first source file
      first_source->cd( path );
      TObject *obj = key->ReadObj();

      PlotUtils::MnvH1D* mnvh1d = dynamic_cast<PlotUtils::MnvH1D*>(obj);
      if (mnvh1d) printf("MnvH1D object\n");
      PlotUtils::MnvH2D* mnvh2d = dynamic_cast<PlotUtils::MnvH2D*>(obj);
      if (mnvh2d) printf("MnvH2D object\n");

      if (mnvh1d) {
              //if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {
         // descendant of TH1 -> merge it

         //      cout << "Merging histogram " << obj->GetName() << endl;
              //TH1 *h1 = (TH1*)obj;

         
         // loop over all source files and add the content of the
         // correspondant histogram to the one pointed to by "h1"
         TFile *nextsource = (TFile*)sourcelist->After( first_source );
         while ( nextsource ) {

            // make sure we are at the correct directory level by cd'ing to path
            nextsource->cd( path );
            TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(key->GetName());
            if (key2) {
                PlotUtils::MnvH1D *h2 = dynamic_cast<PlotUtils::MnvH1D*>(key2->ReadObj());
               mnvh1d->Add( h2 );
               delete h2;
            }

            nextsource = (TFile*)sourcelist->After( nextsource );
         }
      }
      else if (mnvh2d) {
              //if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {
         // descendant of TH1 -> merge it

         //      cout << "Merging histogram " << obj->GetName() << endl;
              //TH1 *h1 = (TH1*)obj;

         
         // loop over all source files and add the content of the
         // correspondant histogram to the one pointed to by "h1"
         TFile *nextsource = (TFile*)sourcelist->After( first_source );
         while ( nextsource ) {

            // make sure we are at the correct directory level by cd'ing to path
            nextsource->cd( path );
            TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(key->GetName());
            if (key2) {
                PlotUtils::MnvH2D *h2 = dynamic_cast<PlotUtils::MnvH2D*>(key2->ReadObj());
               mnvh2d->Add( h2 );
               delete h2;
            }

            nextsource = (TFile*)sourcelist->After( nextsource );
         }
      }

      else if ( obj->IsA()->InheritsFrom( TTree::Class() ) ) {

         // loop over all source files create a chain of Trees "globChain"
         const char* obj_name= obj->GetName();

         globChain = new TChain(obj_name);
         globChain->Add(first_source->GetName());
         TFile *nextsource = (TFile*)sourcelist->After( first_source );
         //      const char* file_name = nextsource->GetName();
         // cout << "file name  " << file_name << endl;
         while ( nextsource ) {

            globChain->Add(nextsource->GetName());
            nextsource = (TFile*)sourcelist->After( nextsource );
         }

      } else if ( obj->IsA()->InheritsFrom( TDirectory::Class() ) ) {
         // it's a subdirectory

         std::cout << "Found subdirectory " << obj->GetName() << std::endl;

         // create a new subdir of same name and title in the target file
         target->cd();
         TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );

         // newdir is now the starting point of another round of merging
         // newdir still knows its depth within the target file via
         // GetPath(), so we can still figure out where we are in the recursion
         MergeRootfile( newdir, sourcelist );

      } else {

         // object is of no type that we know or can handle
         std::cout << "Unknown object type, name: "
         << obj->GetName() << " title: " << obj->GetTitle() << std::endl;
      }

      // now write the merged histogram (which is "in" obj) to the target file
      // note that this will just store obj in the current directory level,
      // which is not persistent until the complete directory itself is stored
      // by "target->Write()" below
      if ( obj ) {
         target->cd();

         //!!if the object is a tree, it is stored in globChain...
         if(obj->IsA()->InheritsFrom( TTree::Class() ))
            globChain->Merge(target->GetFile(),0,"keep");
         else
            obj->Write( key->GetName() );
      }

   } // while ( ( TKey *key = (TKey*)nextkey() ) )

   // save modifications to target file
   target->SaveSelf(kTRUE);
   TH1::AddDirectory(status);
}

int main( int argc, char *argv[] )
{
  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable();
  #endif

  if( argc < 2 )
  {
    std::cout << std::endl;
    std::cout << "Usage: " << std::endl;
    std::cout << "    " << argv[0] << " targetfile source1 [source2 ... sourceN]" << std::endl;
    std::cout << std::endl;
    std::cout << "This is hadd with support for PlotUtils classes." << std::endl;
    std::cout << "This program will add histograms from a list of root files and write them" << std::endl;
    std::cout << "to a target root file. The target file is newly created and must not exist." << std::endl;
    std::cout << "Supply at least two source files for this to make sense... ;-)" << std::endl;
    std::cout << std::endl;
    std::cout << "Note: the hadd flags are not implemented yet." << std::endl;
    return 1;
  }

  Target = TFile::Open( argv[1], "RECREATE" );

  FileList = new TList();
  for( int i = 2; i != argc; ++i )
    FileList->Add( TFile::Open( argv[i] ) );

  MergeRootfile( Target, FileList );
  Target->Close();
  delete Target;
  return 0;
}
