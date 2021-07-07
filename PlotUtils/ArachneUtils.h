#ifndef ARACHNEUTILS_H
#define ARACHNEUTILS_H

#include "TString.h"
#include "TFile.h"
#include "TTree.h"

#include <string>
#include <vector>
#include <sys/types.h>
#include <dirent.h>

using std::string;

namespace PlotUtils {

  class MnvEVD {
    public:
    
    static const std::string data_defaultdir;
    static const std::string mc_defaultdir;
    static bool mc_useNthEventMethod;

    struct Event {

      //! Full constructor
      Event( bool i_isMC, int i_run, int i_subrun, int i_gate, int i_slice, std::string i_info ) :
        isMC( i_isMC ),
        run(i_run),
        subrun(i_subrun),
        gate(i_gate),
        slice(i_slice),
        info(i_info)
      {};


      //! Constructor with no info
      Event( bool i_isMC, int i_run, int i_subrun, int i_gate, int i_slice ) :
        isMC( i_isMC ),
        run(i_run),
        subrun(i_subrun),
        gate(i_gate),
        slice(i_slice),
        info()
      {};

      //! Constructor for all slices with no info
      Event( bool i_isMC, int i_run, int i_subrun, int i_gate) :
        isMC( i_isMC ),
        run(i_run),
        subrun(i_subrun),
        gate(i_gate),
        slice(-1),
        info()
      {};

      //! Default constructor
      Event():
        isMC(false),
        run(-1),
        subrun(-1),
        gate(-1),
        slice(-1),
        info()
      {};


      bool isMC;
      int run;
      int subrun;
      int gate;
      int slice;
      std::string info;
      std::map<std::string, int> extraIntData;
      std::map<std::string, double> extraDoubleData;
    };

    //! Vector of Events
    typedef std::vector<Event> Events;

    //! Pair of string description of event group and vector of events
    typedef std::pair<std::string, Events> EventGroup;
    //! Map from string description of event group to vector of events
    typedef std::map<std::string, Events> EventGroupMap;

    static string FormArachneLink( const Event& evt, 
        const string& basedir = "" );

    static string FormArachneLink( int run, 
        int subrun, 
        int gate, 
        int slice = -1, 
        bool isMC = false,
        const string& basedir = "",
        int datarun = 0,
        int datasubrun = 0 ); 

    /*! Create a list of links to be used for Aranchne/roundup.
      @param[in] events the vector of events
      @param useInfo[in] Do you want to use info as a short name for the event (default=false= show run/sub/gate/slice)
      @param[in] basedir Top level DST directory
     */
    static string FormArachneLinks( const Events& events,
        const bool useInfo = false,
        const string& basedir = "" );



    //! Get the html code to create a link to an event
    static string GetLinkAnchor( const Event& evt,
        bool fulllink = false,
        const string& basedir = "");

    //! Get the html code to create a link to an event
    static string GetLinkAnchor( int run,
        int subrun,
        int gate,
        int slice = -1,
        bool isMC = false,
        bool fulllink = false,
        const string& basedir = "");

    //! Get the html code for a table of event links with descriptions for each event
    //! @return A complete html table object
    static string MakeLinksTable( const Events& events,
        const string& basedir = "" );

    //! Get the html code for a compact table of event links
    //! @return html table object
    static string MakeLinksTableCompact( const Events& events,
        const string& basedir = "" );

    //! String replace utility which is helpful for making links websafe
    static void str_replace_all( std::string& s, const std::string& takeOut, const std::string& putIn = "_" );

    //! Create a bunch of arachne links tables and a header with links to those tables
    //! @return html of a bunch of h2's that link to tables
    static string MakeLinksPage( const EventGroupMap& groupMap,
        const string& basedir = "");

    //! Create a csv page suitable for performing an arachne scan adding selected extra data to the table
    static string MakeScanSheet( const Events& events, const std::vector<std::string>& xData,
        const string& basedir = "");
    static string MakeScanSheet( const Events& events );

  }; // end of MnvEVD

} //end of PlotUtils

std::ostream& operator<<( std::ostream& os, const PlotUtils::MnvEVD::Event& evt );

#endif


