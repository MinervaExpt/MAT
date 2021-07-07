#ifndef ARACHNEUTILS_cxx
#define ARACHNEUTILS_cxx

#include "ArachneUtils.h"

#include <TError.h>

#include <iostream>
#include <sstream>

using namespace PlotUtils;

const std::string MnvEVD::data_defaultdir = "/minerva/data/data_processing/minerva/dst/numibeam/v10r6";
const std::string MnvEVD::mc_defaultdir = "/minerva/data/mc_production/central_value/minerva/dst/v10r6";
bool MnvEVD::mc_useNthEventMethod = true;

std::ostream& operator<<( std::ostream& os, const PlotUtils::MnvEVD::Event& evt )
{
  os
    << "Event = { "
    << "Is MC: "    << (evt.isMC?"T":"F")
    << "; Run: "    << evt.run
    << "; Subrun: " << evt.subrun
    << "; Gate: "   << evt.gate
    << "; Slice: "  << evt.slice
    << ";Info: \"" << evt.info << "\""
    << "}";

  return os;
}

string MnvEVD::FormArachneLink( const MnvEVD::Event& evt, const string& basedir ){
  return FormArachneLink( evt.run, evt.subrun, evt.gate, evt.slice, evt.isMC, basedir );
}

string MnvEVD::FormArachneLink( int run, int subrun, int gate, int slice,  bool isMC, const string& basedir, int datarun, int datasubrun ){

  //-------------------------------
  // add on run number part of path
  //-------------------------------
  string dir = basedir;
  if( ! dir.size() )
    dir = isMC ? mc_defaultdir : data_defaultdir;

  const char *runstr = Form("%08d", run);
  for( unsigned int i = 0; i < 8; i += 2)
    dir += Form("/%c%c/", runstr[i], runstr[i+1] );

  //--------------------------------
  // get subrun string
  //--------------------------------
  TString subrun_str1 = Form("_%04d_", subrun );
  TString subrun_str2 = Form("-%04d", subrun );
  TString subrun_str3 = Form("_%04d-", subrun );

  //--------------------------------
  // Loop over files and find the
  // one we want
  //--------------------------------
  bool foundDST = false;
  DIR *dp;
  struct dirent *ep;
  TString file;
  string filename;

  dp = opendir (dir.c_str());
  if (dp != NULL)
  {
    while ( ( ep = readdir (dp) ) ){
      file = (ep->d_name);
      if( file.Contains("bad") ) continue;
      if( ( file.Contains(subrun_str1.Data()) || file.Contains(subrun_str2.Data()) || file.Contains(subrun_str3.Data()) ) && file.Contains("_DST_") && file.Contains(".root") )
        break;
    }
    if( !ep )
    {
      if( basedir.empty() && isMC )
      {
        //TODO: this is a temprorary fix because MC DSTs are in two locations
        return FormArachneLink(run, subrun, gate, slice, isMC, "/minerva/data/users/minervapro/mc_production_v10r4p4_runbase_00010200/grid/central_value/minerva/dst/v10r4p4" );
        //return FormArachneLink(run, subrun, gate, slice, isMC, "/minerva/data/users/minervapro/mc_production_v10r4p4_runbase_00050200/grid/central_value/minerva/dst/v10r4p4" );
      }

      Warning( "FormArachneLink", Form("Did not find DST for run %d, subrun %d", run, subrun ) );
    }
    else
      foundDST = true;

    (void) closedir (dp);
  }
  else if( basedir.empty() && isMC ) 
  {
    //TODO: this is a temprorary fix because MC DSTs are in two locations
    return FormArachneLink(run, subrun, gate, slice, isMC, "/minerva/data/users/minervapro/mc_production_v10r4p4_runbase_00010200/grid/central_value/minerva/dst/v10r4p4" );
    //return FormArachneLink(run, subrun, gate, slice, isMC, "/minerva/data/users/minervapro/mc_production_v10r4p4_runbase_00050200/grid/central_value/minerva/dst/v10r4p4" );
  }
  else
    perror (Form("Couldn't open the directory %s", dir.c_str()));

  if( ! foundDST )
  {
    return "";
  }

  string dstfile = string(dir) + string(file);

  int entry = -1;
  if( isMC && mc_useNthEventMethod ) 
  {
    entry = gate; //gate for MC means nthEventInFile
  }
  else
  {
    TFile *f = new TFile(dstfile.c_str());

    Int_t ev_run, ev_subrun, ev_gate, mc_run, mc_subrun, mc_nthEvtInFile;
    TBranch *b_ev_run, *b_ev_subrun, *b_ev_gate, *b_mc_run, *b_mc_subrun, *b_mc_nthEvtInFile;

    if( f != NULL && f->IsOpen() ){
      TTree *tree = (TTree*)f->Get("minerva");
      tree->SetBranchStatus("*",0); //disable all branches
      tree->SetBranchStatus("ev_gate",1);

      tree->SetBranchAddress("ev_gate", &ev_gate, &b_ev_gate);
      
      if (isMC) {
        tree->SetBranchStatus("mc_run",1);
        tree->SetBranchStatus("mc_subrun",1);
	tree->SetBranchStatus("mc_int_nevFile",1);
      
        tree->SetBranchAddress("mc_run", &mc_run, &b_mc_run);
	tree->SetBranchAddress("mc_subrun", &mc_subrun, &b_mc_subrun);
	tree->SetBranchAddress("mc_int_nevFile", &mc_nthEvtInFile, &b_mc_nthEvtInFile);
      }

      if( datarun && datasubrun ){
        tree->SetBranchStatus("ev_run",1);
        tree->SetBranchStatus("ev_sub_run",1);

        tree->SetBranchAddress("ev_run", &ev_run, &b_ev_run);
        tree->SetBranchAddress("ev_sub_run", &ev_subrun, &b_ev_subrun);
      }

      for( int i = 0; i < tree->GetEntries(); i++ ){
        tree->GetEntry(i);
	
	if (isMC && (mc_run!=run || mc_subrun!=subrun) ) continue;

        if( datarun && datasubrun ){
          if( ev_run!=datarun || ev_subrun!=datasubrun ) continue;
        }

        if( !isMC && ev_gate == gate ){
          entry = i;
          break;
        }
	if (isMC && mc_nthEvtInFile == gate ) {
	  entry = i;
	  break;
	}
      }
    }
    f->Close();
    delete f;

    if( entry < 0 )
      Warning( "FormArachneLink", Form("Could not find ev_gate = %d in DST file %s", gate, dstfile.c_str() ) );
  }

  string link = Form("http://minerva05.fnal.gov/Arachne/arachne.html?filename=%s&entry=%i&slice=%i",
      dstfile.c_str(),entry,slice);

  return link;
}

string MnvEVD::FormArachneLinks( const MnvEVD::Events& events, const bool useInfo, const std::string& basedir ) {
  string rval;

  for( Events::const_iterator i = events.begin(); i != events.end(); ++i )
  {
    rval += FormArachneLink( *i, basedir );
    if( useInfo )
      rval += Form( " %s", i->info.c_str() );
    rval += "\n";
  }

  return rval;
}

string MnvEVD::GetLinkAnchor( const MnvEVD::Event& evt, bool fullLink, const std::string& basedir ) {
  return GetLinkAnchor( evt.run, evt.subrun, evt.gate, evt.slice, evt.isMC, fullLink, basedir );
}


string MnvEVD::GetLinkAnchor( int run, int subrun, int gate, int slice, bool isMC, bool fullLink, const std::string& basedir ) {
  string link = FormArachneLink( run, subrun, gate, slice, isMC, basedir );
  string name;
  if( fullLink )
    name  = link;
  else
    name = Form("Run %d, Subrun %d, Gate %d, Slice %d", run, subrun, gate, slice );

  string rval = Form("<a href=\"%s\">%s</a>", link.c_str(), name.c_str() );
  return rval;
}

string MnvEVD::MakeLinksTable( const std::vector< MnvEVD::Event >& events, const string& basedir ) {
  string rval = "<table border=\"1\">\n";
  rval += "<tr><th>Event Idx</th><th>Event</th><th>Info</th></tr>\n";

  int idx = 1;
  for( MnvEVD::Events::const_iterator i = events.begin(); i != events.end(); ++i ) {
    rval += Form("<tr><td>%d</td>"   , idx++ );
    rval += Form("<td>%s</td>"       , GetLinkAnchor( *i, false, basedir ).c_str() );
    rval += Form("<td>%s</td></tr>\n", i->info.c_str() );
  }

  rval += "</table>\n";
  return rval;
}

string MnvEVD::MakeLinksTableCompact( const std::vector< MnvEVD::Event >& events, const string& basedir ) {
  string rval = "<html>\n";
  rval += "<head>\n";
  rval += "<title>Arachne Links</title>\n";
  rval += "</head>\n";
  rval += "<body bgcolor='#2B547E' vlink='#FF0000' link='#FFFFFF' text='#C0C0C0'>\n";
  rval += "1)";

  int row = 2;
  int idx = 0;
  for( MnvEVD::Events::const_iterator i = events.begin(); i != events.end(); ++i ) {

    rval += Form("<a href=\"%s\">%i %i %i %i<\a>\n", GetLinkAnchor( *i, false, basedir).c_str(),i->run,i->subrun,i->gate,i->slice);
    if( idx % 10 == 0 ){
      rval += Form("<br><br>\n %i)",row);
      row++;
    }
    idx++;
  }

  rval += "<br><br>";
  rval +=  Form("total entries = %i",idx);
  rval += "</body>\n";
  rval += "</html>\n";
  rval += "</table>\n";
  return rval;
}

void MnvEVD::str_replace_all( std::string& s, const std::string& takeOut, const std::string& putIn /*="_"*/ ) {
  size_t nOut = takeOut.size();
  size_t nIn  = putIn.size();

  size_t pos = s.find( takeOut, 0 );
  while( pos != std::string::npos ) {
    s.replace( pos, nOut, putIn );
    pos = s.find( takeOut, pos + nIn );
  }
}


std::string MnvEVD::MakeLinksPage( const MnvEVD::EventGroupMap& groupMap, const std::string& basedir /*= ""*/) {
  std::stringstream rval;

  //add links to the top of the page
  for( EventGroupMap::const_iterator it = groupMap.begin(); it != groupMap.end(); ++it ) {
    std::string tableName = it->first;
    str_replace_all( tableName, " ", "_" );
    rval << Form("<a href='#%s'>%s</a> <br />", tableName.c_str(), it->first.c_str() ) << std::endl;
  }

  //add the events tables
  for( EventGroupMap::const_iterator it = groupMap.begin(); it != groupMap.end(); ++it )
  {
    std::string tableName = it->first;
    str_replace_all( tableName, " ", "_" );

    //!add a header
    rval << Form("<h2 id='%s'>%s</h2>", tableName.c_str(), it->first.c_str() ) << std::endl;

    //!add the table
    rval << MnvEVD::MakeLinksTable( it->second, basedir ) << std::endl;
  }

  return rval.str();
}

std::string MnvEVD::MakeScanSheet( const MnvEVD::Events& events )
{
  std::vector<std::string> xData;
  return MakeScanSheet(events, xData);
}

std::string MnvEVD::MakeScanSheet( const MnvEVD::Events& events, const std::vector<std::string>& xData, const std::string& basedir /* = "" */){
  std::stringstream rval;

  typedef std::vector<std::string>::const_iterator sVecItr;
  //! Make a header that contains the extra data names
  for( sVecItr i = xData.begin(); i != xData.end(); ++i )
    rval << *i << ", ";

  //! Also add a header for a link to the event
  rval << "Arachne Link" << std::endl;

  //! Loop over all events
  for( Events::const_iterator event = events.begin(); event != events.end(); ++event )
  {
    const std::string link = FormArachneLink( *event, basedir );
    if( "" == link )
    {
      Warning( "MnvEVD::MakeScanSheet", "Arachne link could not be formed.  Skip this entry." );
      continue;
    }

    //! Add the xData in the order supplied by the caller (check both int and double extra data)
    for( sVecItr i = xData.begin(); i != xData.end(); ++i )
    {
      std::map<std::string, int>::const_iterator itInt = event->extraIntData.find(*i);
      if( itInt != event->extraIntData.end() )
      {
        rval << itInt->second << ", ";
        continue;
      }

      std::map<std::string, double>::const_iterator itDouble = event->extraDoubleData.find(*i);
      if( itDouble != event->extraDoubleData.end() )
      {
        rval << itDouble->second << ", ";
        continue;
      }      

      //extra data not found, don't write anything.
      rval << ", ";
    }

    //! Add the Arachne link to the last column
    rval << link << std::endl;
  }

  return rval.str();
}

#endif

//#############################################################################
//
// END 
//
//#############################################################################
