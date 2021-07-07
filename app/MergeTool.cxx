#ifndef _MERGETOOL_
#define _MERGETOOL_

#include "TString.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TList.h"
#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include <stdlib.h>// or cstdlib is c++
#include <unistd.h>

#ifndef __CINT__
#include "glob.h"
#endif

using namespace std;

class MergeTool {
    public:
        MergeTool() : m_realdata(false), m_fullpath(false), m_checkmeta(true) {};
        ~MergeTool(){};
        
        void CombineMergedFiles(const char* outDir, const char* tag, const char* save_name = "");
        void SingleMergeFromExternalPlaylist(const char* extPlaylistPath, const char* outDir, const char* save_name = "");
        void SingleMergeRuns(const char* inDirBase, const char* tag, int first_run, int last_run, const char* outDir, const char* save_name = "");
        void MergeEachRun(const char* inDirBase, const char* tag, int first_run, int last_run, const char* outDir, const char* save_name = "");

        void AddListOfFiles( const char* merge_name, const vector< string > &listOfFiles ){ m_listOfAllFiles[string(merge_name)] = listOfFiles; }

        void SetMinervaRelease(char const * var){ minerva_release = string(var); }
        void SetCheckMetaData(bool var){ m_checkmeta = var; }
        void SetIsRealData(bool var){ m_realdata = var; }
        void SetFullPath(bool var){ m_fullpath = var; }
        void SetAnaTree(string anatree){ m_anatree = anatree; }
        void SetAnaTool(string anatool){ m_anatool = anatool; }
        void SetTrackerCut(bool var){ m_nontracker_truthcut = var; }
        
    private:
        void VerifyFilesandGetMeta();        
        void Merge( const char* outDir, const char* save_name );

        string minerva_release;
        string m_anatree;
        string m_anatool;

        bool m_realdata;
        bool m_fullpath;
        bool m_checkmeta;
        bool m_nontracker_truthcut;
      
        map< string, vector< string > > m_listOfAllFiles;
        map< string, vector< string > > m_listOfGoodFiles;

        map< string, map< string, double > > m_metaInfo;

        const string pnfs = "/pnfs";
        const int n_pnfs = pnfs.size();
        const string xrootd = "root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr";
        const string tracker_cut = "mc_vtx[2]>5891 && mc_vtx[2]<8439";
};

#endif

void MergeTool::CombineMergedFiles(const char* outDir, const char* tag, const char* save_name)
{    
    //Grab the files from the following wildcard
    TString inGlob(TString::Format("%s/merged_%s_%s_run*.root", outDir, tag, save_name));
    cout << "glob :" << inGlob.Data() << endl;
    glob_t g;
    glob(inGlob.Data(), 0, 0, &g);

    vector< string > tmpFiles;
    for(int i=0; i<(int)g.gl_pathc; ++i){
        string s(g.gl_pathv[i]);
        s.replace(0, n_pnfs, xrootd);
        tmpFiles.push_back(s);
    }

    AddListOfFiles("Full",tmpFiles);
    VerifyFilesandGetMeta();        
    Merge( outDir, save_name );
}

void MergeTool::MergeEachRun(const char* inDirBase, const char* tag, int first_run, int last_run, const char* outDir, const char* save_name )
{
    //Loop over the runs
    for( int run = first_run; run < last_run + 1; ++run ){
        //Seperate string into run directories
        string runStr(TString::Format("%08d", run));
        string runStrParts[4];
        for(int i=0; i<4; ++i) runStrParts[i]=runStr.substr(i*2, 2);

        //Set the subrun path
        string subpath;
        if(!m_fullpath) subpath = m_realdata? "grid/minerva/ana/numibeam/":"grid/central_value/minerva/ana/";
        else subpath = "";

        TString inGlob(TString::Format("%s/%s%s/%s/%s/%s/%s/%s_*%s_*_%s*.root",
                                       inDirBase,
                                       subpath.c_str(),
                                       minerva_release.c_str(),
                                       runStrParts[0].c_str(),
                                       runStrParts[1].c_str(),
                                       runStrParts[2].c_str(),
                                       runStrParts[3].c_str(),
                                       m_realdata ? "MV" : "SIM",
                                       runStr.c_str(),
                                       tag));
        cout << "glob :" << inGlob.Data() << endl;
        glob_t g;
        glob(inGlob.Data(), 0, 0, &g);

        vector< string > tmpFiles;
        for(int i=0; i<(int)g.gl_pathc; ++i){
            string s(g.gl_pathv[i]);
            s.replace(0, n_pnfs, xrootd);
            tmpFiles.push_back(s);
        }

        AddListOfFiles(runStr.c_str(),tmpFiles);
    }

    VerifyFilesandGetMeta();        
    Merge( outDir, save_name );
}

void MergeTool::SingleMergeFromExternalPlaylist(const char* extPlaylistPath, const char* outDir, const char* save_name)
{
    ifstream fin(extPlaylistPath);
    assert(fin.is_open());
    std::string filename;

    //Grab the files from the playlist
    vector< string > tmpFiles;
    while (getline(fin, filename))
    {
        tmpFiles.push_back(filename);
    }
    fin.close();

    AddListOfFiles("Playlist",tmpFiles);
    VerifyFilesandGetMeta();        
    Merge( outDir, save_name );
}

void MergeTool::SingleMergeRuns(const char* inDirBase,  const char* tag, int first_run, int last_run,const char* outDir,const char* save_name)
{
    //Loop over the runs
    vector< string > tmpFiles;
    for( int run = first_run; run < last_run + 1; ++run ){
        //Seperate string into run directories
        string runStr(TString::Format("%08d", run));
        string runStrParts[4];
        for(int i=0; i<4; ++i) runStrParts[i]=runStr.substr(i*2, 2);

        //Set the subrun path
        string subpath;
        if(!m_fullpath) subpath = m_realdata? "grid/minerva/ana/numibeam/":"grid/central_value/minerva/ana/";
        else subpath = "";

        TString inGlob(TString::Format("%s/%s%s/%s/%s/%s/%s/%s_*%s_*_%s*.root",
                                       inDirBase,
                                       subpath.c_str(),
                                       minerva_release.c_str(),
                                       runStrParts[0].c_str(),
                                       runStrParts[1].c_str(),
                                       runStrParts[2].c_str(),
                                       runStrParts[3].c_str(),
                                       m_realdata ? "MV" : "SIM",
                                       runStr.c_str(),
                                       tag));
        cout << "glob :" << inGlob.Data() << endl;
        glob_t g;
        glob(inGlob.Data(), 0, 0, &g);

        for(int i=0; i<(int)g.gl_pathc; ++i){
            string s(g.gl_pathv[i]);
            s.replace(0, n_pnfs, xrootd);
            tmpFiles.push_back(s);
        }
    }

    string runRangeStr(TString::Format("%08d-%08d", first_run, last_run));
    AddListOfFiles(runRangeStr.c_str(),tmpFiles);
    VerifyFilesandGetMeta();        
    Merge( outDir, save_name );
}

//Verify and grab meta tree info
//I know this is very... sequential, but this minimizes the number of times we open the file.  The overhead should be small
void MergeTool::VerifyFilesandGetMeta()
{
    TStopwatch ts;
   
    //Loop over every set of files we're supposed to merge
    map< string, vector< string > >::iterator mergeList;
    vector< string >::iterator filename;
    for( mergeList = m_listOfAllFiles.begin(); mergeList != m_listOfAllFiles.end(); mergeList++ )
    {
        //Set up the meta info storage
        m_metaInfo[mergeList->first]["POT_Total"] = 0;
        m_metaInfo[mergeList->first]["POT_Used"] = 0;
        m_metaInfo[mergeList->first]["nEntries_Reco"] = 0;
        if(!m_realdata) m_metaInfo[mergeList->first]["nEntries_Truth"] =0;
          
        //Loop over every file in this list
        for( filename = mergeList->second.begin(); filename != mergeList->second.end(); filename++ )
        {
            //Quality checks for the file itself
            TFile *f = TFile::Open(&(*filename->c_str()));//Open expects a const char*.  The extra & is needed due to ambiguity
            
            if(f->IsZombie()){
                f->Close();
                cout<<*filename<<" is a zombie. Skipping"<<endl;
                continue;
            }

            if(!m_realdata && !f->Get("Truth")){//Check this
                f->Close();
                cout<<*filename<<" has no Truth tree. Skipping"<<endl;
                continue;
            }

            //Quality checks for the meta tree 
            TTree* meta=(TTree*)f->Get("Meta");
            if(!meta){
                f->Close();
                cout<<*filename<<" has no Meta tree. Skipping"<<endl;
                continue;
            }

            if(meta->GetEntries() != 1){
                f->Close();
                cout<<*filename<<" Meta tree should have one entry. Skipping"<<endl;
                continue;
            }

            if(!meta->GetBranch("POT_Total") || !meta->GetBranch("POT_Used")){
                f->Close();
                cout<<*filename<<" is missing a POT_Total or POT_Used branch. Skipping"<<endl;
                continue;
            }

            if(!meta->GetBranch((string("nEntries_")+m_anatree).c_str())){
                f->Close();
                cout<<*filename<<" is missing "<<string("nEntries_")+m_anatree<<". Skipping"<<endl;
                continue;
            }

            if(!m_realdata && !meta->GetBranch("nEntries_Truth")){
                f->Close();
                cout<<*filename<<" is missing nEntries_Truth. Skipping"<<endl;
                continue;
            }

            //If the file passed all the checks, add it to the good list
            m_listOfGoodFiles[mergeList->first].push_back(*filename);

            cout<<"Getting POT for "<<*filename<<" and summing"<<endl;

            //With the meta, we can now grab info from it
            meta->GetEntry(0);
                                
            TLeaf* lPOT_Total = meta->GetLeaf("POT_Total");
            if(lPOT_Total) m_metaInfo[mergeList->first]["POT_Total"] += lPOT_Total->GetValue();

            TLeaf* lPOT_Used = meta->GetLeaf("POT_Used");
            m_metaInfo[mergeList->first]["POT_Used"]  += lPOT_Used->GetValue();

            TLeaf* lEntries_Reco = meta->GetLeaf((string("nEntries_")+m_anatree).c_str());
            m_metaInfo[mergeList->first]["nEntries_Reco"] += lEntries_Reco->GetValue();

            if(!m_realdata){
                TLeaf* lEntries_Truth = meta->GetLeaf("nEntries_Truth");
                m_metaInfo[mergeList->first]["nEntries_Truth"] += lEntries_Truth->GetValue();
            }          

            delete meta;//Added 210117
            f->Close();
            delete f;
            //sleep(1);//Just in case
        }
        cout<<mergeList->first<<endl;
        cout << "Found " << m_listOfGoodFiles[mergeList->first].size() << " good files out of " << mergeList->second.size() << endl;
        cout << "POT totals: Total = " << m_metaInfo[mergeList->first]["POT_Total"] << " Used = " << m_metaInfo[mergeList->first]["POT_Used"] << endl;

    }
    ts.Stop();
    cout << "----Finished verifying and getting POT.  Time: "<<endl;
    ts.Print();
    cout << endl; 
}

void MergeTool::Merge(const char* outDir, const char* save_name)
{
    TStopwatch tsall;

    map< string, vector< string > >::iterator mergeList;
    vector< string >::iterator filename;
    for( mergeList = m_listOfGoodFiles.begin(); mergeList != m_listOfGoodFiles.end(); mergeList++ )
    {
        if(mergeList->second.size()==0){
            cout << "No good files in "<<mergeList->first<<", nothing to do..." << endl;
            continue;
        }
        
        TStopwatch ts;

        //Open the output file
        TString output=TString::Format("%s/%s_%s.root",outDir,save_name,mergeList->first.c_str());
        TFile * fout = new TFile(output, "RECREATE");
        cout << "Merging ana tree" << endl;
        fout->cd(); // Just in case the surrounding lines get separated
      
        TChain * inChain = new TChain(m_anatree.c_str());
        TChain * inChainTruth = new TChain("Truth");

        //Add the files to the chain
        for( filename = mergeList->second.begin(); filename != mergeList->second.end(); filename++ )
        {
            cout << "Attempting to add " << *filename << endl;
            inChain->Add(&(*filename->c_str()));
            if(!m_realdata) inChainTruth->Add(&(*filename->c_str()));
        }
 
        fout->cd(); // Just in case the surrounding lines get separated

        if (inChain->Merge(fout, 32000, "keep SortBasketsByBranch") <= 0) {
          remove(output);
          cout << "Error: reco tree is not filled due to failure of input file access : No output produced !!!" << endl;
          exit(1);
        }
    
        //Merge the truth tree
        if(!m_realdata){
            cout << "Merging truth tree" << endl;
            fout->cd();
            const char* selection = "";
            if (m_nontracker_truthcut) {
                selection = tracker_cut.c_str();
            }
            TTree* outTreeTruth=inChainTruth->CopyTree(selection);
            outTreeTruth->Write();
            fout = outTreeTruth->GetCurrentFile();
            //update Total_Truth_Entries if tracker seletion applied.
            m_metaInfo[mergeList->first]["nEntries_Truth"]=outTreeTruth->GetEntries();
        }
      
        //Put in the meta tree
        cout<<"putting in the new meta tree"<<endl;    
        fout->cd();

        TTree* newMetaTree=new TTree("Meta", "");
        newMetaTree->Branch("POT_Used", &m_metaInfo[mergeList->first]["POT_Used"]);
        newMetaTree->Branch("POT_Total", &m_metaInfo[mergeList->first]["POT_Total"]);
        newMetaTree->Branch("Total_Reco_Entries",&m_metaInfo[mergeList->first]["nEntries_Reco"]);
        if(!m_realdata) newMetaTree->Branch("Total_Truth_Entries", &m_metaInfo[mergeList->first]["nEntries_Truth"]);

        newMetaTree->Fill();
        newMetaTree->Write();
        fout = newMetaTree->GetCurrentFile();
        fout->Close();

        ts.Stop();
        cout << mergeList->first <<" merging time:" << endl;
        ts.Print();

        delete fout;//Added 210117
    }
    tsall.Stop();
    cout << "Full merging time:" << endl;
    tsall.Print();
}

int main(int argc, char *argv[])
{
    /////////////////////////////////////////////////////////////////////
    //Environment variables
    ////////////////////////////////////////////////////////////////////i/
    
    char const * user_name = getenv("USER");
    if(!user_name){
        std::cout << "[WARNING]: Environment variable \"USER\" not found." << std::endl;
    }
    string username(user_name);
    
    char const * ana_tree = getenv("ANATREENAME");
    if(!ana_tree){
        std::cout << "[WARNING]: Environment variable \"ANATREE\" not set. To set see requirements file." << std::endl;
    }
    string anatree(ana_tree);// t -- set tree
    
    char const * ana_tool = getenv("ANATOOLNAME");
    if(!ana_tool){
        std::cout << "[WARNING]: Environment variable \"ANATOOLNAME\" not set. To set see requirements file." << std::endl;
    }
    string anatool(ana_tool);
    
    char const * minerva_release = getenv("MINERVA_RELEASE");
    if(!minerva_release){
        std::cerr << "[ERROR]: environment variable \"MINERVA_RELEASE\" not set. "
        "Cannot determine source tree location." << std::endl;
        return 1;
    }
    
    string ana_save_name = "minerva";// s -- savename

    /////////////////////////////////////////////////////////////////////
    //Options
    ////////////////////////////////////////////////////////////////////i/
    
    bool meta_data_check = true;
    bool real_data       = false;

    int merge = -1;

    bool full_path  = false;
    bool re_opt_i   = false; string infile;  
    bool re_opt_n   = false; string run_s;   
    bool re_opt_e   = false; string extPlaylist;
    bool re_opt_o   = false; string outfile; 
    bool is_per_dir = false;
    bool cut_tracker = false;
    
    char cc;
    while ((cc = getopt(argc, argv, "i:o:t:s:n:h:pa:m:e:cdfv:r")) != -1) {
        switch (cc){
            case 'v': minerva_release = optarg; break;
            case 'i': infile = optarg;  re_opt_i = true; break;
            case 'o': outfile = optarg; re_opt_o = true; break;
            case 't': anatree = optarg; break;
            case 'n': run_s = optarg; re_opt_n = true; break;
            case 's': ana_save_name = optarg; break;
            case 'm': merge = atoi(optarg); break;
            case 'e': extPlaylist = optarg; re_opt_e = true; break;
            case 'p': is_per_dir = true; break;
            case 'a': anatool = optarg; break;
            case 'c': meta_data_check = false; break;
            case 'd': real_data = true; break;
            case 'f': full_path = true; break;
            case 'r': cut_tracker =true; break;
            case 'h':
            //cout << argv[0] << endl
                cout << "|************************************************** Example *******************************************************|" << endl
                << "| To merge from runs 13200 to 13260 of an analysis output into a single root file do:                              |" << endl
                << "|                                                                                                                  |" << endl
                << "| MergeTool.exe -i /pnfs/minerva/persistant/persistent/users/dcoplowe/CC1P1Pi_PL13C_111216 -n 13200-13260          |" << endl
                << "|                                                                                                                  |" << endl
                << "| This will save the output into the -i directory                                                                  |" << endl
                << "|                                                                                                                  |" << endl
                << "| To merge individual runs do:                                                                                     |" << endl
                << "|                                                                                                                  |" << endl
                << "| MergeTool.exe -i /pnfs/minerva/persistant/persistent/users/dcoplowe/CC1P1Pi_PL13C_111216 -n 13200-13260 -m=1     |" << endl
                << "|                                                                                                                  |" << endl
                << "| To combine the output of -m=1 do:                                                                                |" << endl
                << "|                                                                                                                  |" << endl
                << "| MergeTool.exe -i /pnfs/minerva/persistant/persistent/users/dcoplowe/CC1P1Pi_PL13C_111216 -n 13200-13260 -m=2     |" << endl
                << "|                                                                                                                  |" << endl
                << "|******************************************************************************************************************|" << endl
                << " " << endl
                << " " << endl
                << "|*********************** Run Options ****************************|" << endl
                << "| Default is to merge files in the root directory of your analy- |" << endl
                << "| sis. The number of runs also needs to be defined and can be    |" << endl
                << "| either a single run or run range (see below).                  |" << endl
                << "|                                                                |" << endl
                << "| Options:                                                       |" << endl
                << "|   -i   :  Set input file dir if \"-per\" is defined only the   |" << endl
                << "|        :  root directory of your analysis in your persistent   |" << endl
                << "|        :  directory is required.                               |" << endl
                << "|        :                                                       |" << endl
                << "|   -n   :  Run or run range: start-end e.g 13200-13250          |" << endl
                << "|        :  will run over 50 runs from 13200 to 13250.           |" << endl
                << "|        :                                                       |" << endl
                << "|   -o   :  Set output file directory.                           |" << endl
                << "|        :                                                       |" << endl
                << "|  -per  :  Assume in/out files are in the users persistent      |" << endl
                << "|        :  directory. In such cases \"-i\" becomes the root     |" << endl
                << "|        :  directory name. (e.g. -per -i CC1P1Pi_PL13C_290916   |" << endl
                << "|        :  to merge files in /pnfs/minerva/persistent/users/    |" << endl
                << "|        :  dcoplowe/CC1P1Pi_PL13C_290916/                       |" << endl
                << "|        :                                                       |" << endl
                << "|   -e   :  Name of external playlist with paths of files to     |" << endl
                << "|        :  mergeSet name of analysis tree.                      |" << endl
                << "|        :                                                       |" << endl
                << "|   -a   :  Analysis tool name (Your analysis tool used to make  |" << endl
                << "|        :  the root files you want to merge)                    |" << endl
                << "         :  Default is currently " << anatool << endl
                << "|        :                                                       |" << endl
                << "|   -t   :  Set name of analysis tree.                           |" << endl
                << "         :  Default is currently " << anatree << endl
                << "|        :  If not set, this can be set in the PlotUtils requir- |" << endl
                << "|        :  ements files.                                        |" << endl
                << "|        :                                                       |" << endl
                << "|   -s   :  Set save name.                                       |" << endl
                << "|        :                                                       |" << endl
                << "|   -v   :  Set version of input anatuple files (vXrYpZ).        |" << endl
                << "|        :  If not specified, use current release.               |" << endl
                << "|        :                                                       |" << endl
                << "|   -m   :  Option 1 (-m 1): Combine each run into single root   |" << endl
                << "|        :  file.                                                |" << endl
                << "|        :  Option 2 (-m 2): Combine output of (-m 1) into a     |" << endl
                << "|        :  single root file.                                    |" << endl
                << "|        :                                                       |" << endl
                << "| -check :  Merge without checking POT in Meta tree is good      |" << endl
                << "|        :  (exists).                                            |" << endl
                << "|        :                                                       |" << endl
                << "| -data  :  Run on a data file (merge without Truth tree). Def-  |" << endl
                << "|        :  ault is to assume you are running on MC.             |" << endl
                << "|        :                                                       |" << endl
                << "| -full  :  Add full directory path. This enables you to merge   |" << endl
                << "|        :  special runs where the sub directories are not       |" << endl
                << "|        :  grid/central_value/minerva/ana. When \"-full\" is      |" << endl
                << "|        :  defined simply set the full path up to and including |" << endl
                << "|        :  ana using \"-i\".                                      |" << endl
                << "|        :                                                       |" << endl
                << "| -help  :  Print this.                                          |" << endl
                << "|        :                                                       |" << endl
                << "|****************************************************************|" << endl;
                return 1; break;
            default: return 1;
        }
    }
    
    //--need some input files
    if(!(re_opt_i || re_opt_n || re_opt_e)){
        cout << "|============ Minimum Requirements to run ============|" << endl;
        cout << "|                                                     |" << endl;
        cout << "|    -i     Set input file dir name                   |" << endl;
        cout << "|    -n     Number of runs to merge                   |" << endl;
        cout << "|    -help  For more options.                         |" << endl;
        cout << "|_____________________________________________________|" << endl;
        return 0;
    }

    //--parses the range of runs    
    int first_run = -999;
    int last_run =  -999;
    TString run_ts = run_s;
    if(run_ts.Contains("-",TString::kExact)){
        TString tmp_first( run_ts(0,run_ts.First("-")) );
        first_run = tmp_first.Atoi();
        
        TString tmp_last( run_ts(run_ts.First("-") + 1, run_ts.Length()) );
        last_run = tmp_last.Atoi();
    }
    else{
        first_run = run_ts.Atoi();
        last_run = first_run;
    }

    //--set the base directory is in persistent
    if(is_per_dir){
        string per_dir = "/pnfs/minerva/persistent/users/" + username + "/";
        
        string tmp_in(infile);
        infile = per_dir + tmp_in;
        
        string tmp_out(outfile);
        outfile = per_dir + tmp_out;
    }
    
    if(!re_opt_o) outfile = infile;
    
    cout << "Merge = " << merge << endl;
    
    cout << "|---------------------------------- Inputs ----------------------------------" << endl;
    if(meta_data_check) cout << "|                   Checking POT info is good. (Switch off using -check)" << endl;
    if(is_per_dir || full_path || real_data) cout << "| Option(s) called: " << endl;
    if(is_per_dir) cout << "|                   (-per)   In/out files are in persistent." << endl;
    if(real_data) cout << "|                   (-real) Merging real data files." << endl;
    if(full_path) cout << "|                   (-full) User defining full path." << endl;
    if(re_opt_e ) cout << "|                   (-e) Name of playlist " << extPlaylist <<endl;
    cout << "| Input  (-i): " << infile << endl;
    cout << "| N Runs (-n): " << run_s << endl;
    cout << "| Output (-o): " << outfile << endl;
    cout << "| Analysis Tree Name (-t): " << anatree << endl;
    cout << "| Analysis Tool Name (-a): " << anatool << endl;
    cout << "| Optional Save Name (-s): " << ana_save_name << endl;
    cout << "|--------------------------------- Running ----------------------------------" << endl;

    MergeTool * merger = new MergeTool();
    merger->SetCheckMetaData(meta_data_check);
    merger->SetMinervaRelease(minerva_release);
    merger->SetIsRealData(real_data);
    merger->SetFullPath(full_path);
    merger->SetAnaTree(anatree);
    merger->SetAnaTool(anatool);
    merger->SetTrackerCut(cut_tracker);
    if(merge == 2){
        cout<< "Only merging merged runs" << endl;
        merger->CombineMergedFiles(infile.c_str(), anatool.c_str(), ana_save_name.c_str());
    }
    else if(merge == 1){
        cout << "Merging sub-runs for each run" << endl;
        merger->MergeEachRun(infile.c_str(), anatool.c_str(), first_run, last_run, outfile.c_str(), ana_save_name.c_str());
    }
    else if(re_opt_e == true){
        cout << "Merging subruns according to external playlist" << endl;
        merger->SingleMergeFromExternalPlaylist(extPlaylist.c_str(), outfile.c_str(), ana_save_name.c_str());
    }
    else{
        cout << "Merging sub-runs for each run into a single root file" << endl;
        merger->SingleMergeRuns(infile.c_str(), anatool.c_str(), first_run, last_run, outfile.c_str(), ana_save_name.c_str());
    }
    delete merger;
    cout << "|-------------------------- Finished merging files --------------------------" << endl;
    return 0;
}
