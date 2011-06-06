/** 
* @file RootOutputDesc.cxx
* @brief definition of the class RootOutputDesc
*
*  $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/RootOutputDesc.cxx,v 1.5 2009/09/12 16:00:40 heather Exp $
*  Original author: Heather Kelly heather@lheapop.gsfc.nasa.gov
*/
#ifndef RootOutputDesc_cxx
#define RootOutputDesc_cxx

#include "RootOutputDesc.h"
#include "rootUtil/RuUtil.h"
#include "facilities/Util.h"
#include "TROOT.h"
#include <vector>
#include <iostream>


RootOutputDesc::RootOutputDesc
( const std::string& outputFile, 
   const std::string & tree, 
   int compressionLevel, const std::string& treeTitle, bool verbose )
 : m_fileName(outputFile), m_treeName(tree), m_compressionLevel(compressionLevel), m_tree(0), 
   m_treeTitle(treeTitle), m_verbose(verbose), m_eventCounter(0), m_updated(false)
 {
	 openFile();
 }
    
RootOutputDesc::~RootOutputDesc() 
 {
  if (m_tree)
	   this->getCurrentFile()->Close();
   
  m_tree = 0;
 }

bool RootOutputDesc::openFile() {

	bool stat = true;
    // Save the current directory for the ntuple writer service 
    TDirectory * saveDir = gDirectory ;	
 if (m_tree)
   {
    getCurrentFile()->Close();
    delete m_tree;
    m_tree = 0;
   }

    facilities::Util::expandEnvVar(&m_fileName);
     // Create the new ROOT file
    TFile *f = TFile::Open(m_fileName.c_str(), "RECREATE");
    if (!f->IsOpen()) {
		stat = false;
		std::cout << "RootOutputDesc Error opening ROOT file " << m_fileName << std::endl;
        return stat;
    }
    f->cd();
    f->SetCompressionLevel(m_compressionLevel);
    
    // tree of events
    m_tree = new TTree(m_treeName.c_str(), m_treeTitle.c_str());

    saveDir->cd();

    return stat;
}

bool RootOutputDesc::setupBranch(const std::string& branchName, const std::string &classname, void* branchAddr, 
                                 int bufSize, int splitLevel) {
    
    TBranch *br = m_tree->Branch(branchName.c_str(),classname.c_str(), branchAddr, bufSize, splitLevel);
    if (0 == br) {
        return false;
    }
    return true;
}


TFile* RootOutputDesc::getCurrentFile() {
	if (m_tree)
		return m_tree->GetCurrentFile();
	return 0;
}
    
bool RootOutputDesc::fillTree(int autoSaveInterval)
 {
  bool stat = true;
  try
   {
    rootUtil::AutoCd cd(getCurrentFile()) ;
    if (getCurrentFile()->TestBits(TFile::kWriteError))
     { std::cerr <<"TestBits Failed" << std::endl;  throw ; }
    Int_t numBytes = m_tree->Fill() ;
    if (numBytes < 0) {
        std::cerr << "Filling TTree Error for " << m_treeName << std::endl;
        throw;
    }
    ++m_eventCounter;
    if ( m_eventCounter % autoSaveInterval == 0 )
     {
      if ( m_tree->AutoSave("flushbaskets") == 0 )
       { std::cerr << "AutoSave Failed" << std::endl;
         throw ; 
       }
     }
   }
  catch(...)
   {   
    std::cerr << "Failed to write the event to " << m_treeName << " Tree" << std::endl ;   
    std::cerr << "Exiting..." << std::endl ;   
    std::cerr.flush() ;   
    exit(1) ;   
   } 
  return stat;
 }

bool RootOutputDesc::closeFile()
 {
    // Purpose and Method:  Writes the ROOT file at the end of the run.
    //    The TObject::kWriteDelete parameter is used in the Write method
    //    Used rather than TObject::kOverwrite - supposed to be safer but slower
    //    since ROOT will periodically write to the ROOT file when the bufSize
    //    is filled.  Writing would create 2 copies of the same tree to be
    //    stored in the ROOT file, if we did not specify kOverwrite or kWriteDelete
    try {
        if (m_tree) {
            TDirectory *saveDir = gDirectory;
            TFile *f = m_tree->GetCurrentFile();
            f->cd();
            m_tree->BuildIndex("m_runId", "m_eventId");
            f->Write(0, TObject::kWriteDelete);
            f->Close();
            saveDir->cd();
            return true;
        }
    }
    catch(...) { 
        std::cerr << "Failed to final write to " << m_tree->GetName() << std::endl; 
        std::cerr << "Exiting..." << std::endl; 
        std::cerr.flush(); 
        exit(1); 
    } 

    return false;
}

bool RootOutputDesc::fileExists( const std::string & filename )
 {
  bool fileExists = false ;
  TFile * file = TFile::Open(filename.c_str()) ;
  if (file)
   {
    if (!file->IsZombie()) 
     {
      file->Close() ;
      fileExists = true ;
     }
    delete file ;
   }
  return fileExists ;
 }  

#endif

