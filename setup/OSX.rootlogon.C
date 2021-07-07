{
  if( gSystem->Getenv("PLOTUTILSROOT") )
  {
    gInterpreter->AddIncludePath( gSystem->ExpandPathName("$PLOTUTILSROOT") );
    string newpath = string(gROOT->GetMacroPath()) + ":" + string(gSystem->ExpandPathName("$PLOTUTILSROOT")) + "/PlotUtils";
    gROOT->SetMacroPath( newpath.c_str() );
    
    //    gSystem->Load( "libCintex.so" );  // needed to process the dictionaries for the objects
    //Cintex::Enable();
    gSystem->Load( gSystem->ExpandPathName("$PLOTUTILSROOT/libplotutils.so") );
    std::cout << "tried to load libplotutils.so" << endl;
    gInterpreter->ExecuteMacro("PlotStyle.C");
  }

}

