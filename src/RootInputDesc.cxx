/** 
* @file RootInputDesc.cxx
* @brief definition of the class RootInputDesc
*
*  $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/RootInputDesc.cxx,v 1.7 2007/08/10 02:46:43 heather Exp $
*  Original author: Heather Kelly heather@lheapop.gsfc.nasa.gov
*/

#include "RootInputDesc.h"
#include "facilities/Util.h"
#include "TROOT.h"
#include <vector>
#include <iostream>


// ctor for vanilla reading from a list of ROOT files
// Setup the TChain* in this case, as we may have multiple files
RootInputDesc::RootInputDesc
 ( const StringArrayProperty& fileList, 
   const std::string & tree, 
   const std::string & branch, bool verbose )
 : m_tree(tree), m_branch(branch), m_chain(0), m_dataObject(0), m_verbose(verbose), m_runEvtIndex(0)
 {
  // Set up the input file list
  m_numEvents = setFileList(fileList, m_verbose) ;
 }

 // ctor for Event Collections
 // Set up the TChain* as we will retrieve a TChain* from the PointerSkim via the EventCollectionMgr
 RootInputDesc::RootInputDesc
     ( TChain *t,
     const std::string & treeName,
     const std::string & branchName, bool verbose)
     : m_tree(treeName), m_branch(branchName), m_chain(t), m_dataObject(0), m_verbose(verbose), 
       m_runEvtIndex(0) {
       m_numEvents = setEventCollection();
     }

    
RootInputDesc::~RootInputDesc() 
 {
  if (m_chain)
   {
    m_chain->GetFile()->Close() ;
    delete m_chain ;
   }
  if (m_dataObject) delete m_dataObject ;

 }

 Long64_t RootInputDesc::setEventCollection( ) {
     Long64_t numEvents = 0;
     // Save the current directory for the ntuple writer service
     TDirectory * saveDir = gDirectory ;	

     if (!m_chain) return numEvents;

     if (m_dataObject) 
         delete m_dataObject ;

     m_dataObject = new TObject * ;
     *m_dataObject = 0 ;

     // Check that we have the right TTree
     // To do this we will first need to "rediscover" the directory name
     std::string treeName = std::string(m_chain->GetName());
     if (treeName != m_tree)
     {
         return numEvents ;
     }

     // To set the data pointer we need to find the name of the branch contained
     // where we note that in Glast each file contains one branch
     TObjArray * branches = m_chain->GetListOfBranches() ;
     std::string branchName = "" ;
     int idx ;
     for( idx = 0 ; idx < branches->GetEntries() ; idx++ )
     {
         TBranch * branch = dynamic_cast<TBranch*>(branches->At(idx)) ;
         if (branch)
         {
             branchName = std::string(branch->GetName()) ;
             if (branchName == m_branch) break ;
         }
     }
     if (branchName!=m_branch) 
     { return numEvents ; }

     // Set the data pointer to receive the data and we are ready to go
     m_chain->SetBranchAddress(branchName.data(),m_dataObject) ;

     // If necessary, set up the TChainIndex
     TVirtualIndex *resetIndex = m_chain->GetTreeIndex();
     m_chain->BuildIndex("m_runId", "m_eventId") ;
     m_runEvtIndex = m_chain->GetTreeIndex();
     m_chain->SetTreeIndex(resetIndex);

     numEvents = m_chain->GetEntries() ;

     // Finally restore the context we entered with
     saveDir->cd() ;

     return numEvents ;
 }

Long64_t RootInputDesc::setFileList( const StringArrayProperty & fileList, bool verbose )
 {
  Long64_t numEvents = -1;

  m_fileList = fileList ;
 
  // Save the current directory for the ntuple writer service
  TDirectory * saveDir = gDirectory ;	

  // Check to see if we are changing files (which means a TChain exists)
  if (m_chain)
   {
    m_chain->GetFile()->Close() ;
    delete m_chain ;
    m_chain = 0;
   }
   if (m_dataObject) 
       delete m_dataObject ;

  // Get a new instance of a TChain for processing this file
  m_chain = new TChain(m_tree.c_str()) ;
  m_dataObject = new TObject * ;
  *m_dataObject = 0 ;

  // Add the files to the TChain
  typedef std::vector<std::string>::const_iterator StringVecIter ;
  StringVecIter fileListItr ;
  int fileCount = 0 ;
  for ( fileListItr = m_fileList.value().begin() ; 
        fileListItr != m_fileList.value().end() ;
        fileListItr++ )
   {
    std::string fileName = *fileListItr ;
    facilities::Util::expandEnvVar(&fileName);
    if (fileExists(fileName))
     {
      int nf = m_chain->Add(fileName.c_str()) ;
      if (verbose) std::cout << "RootInputDesc::setFileList opening: " << fileName << std::endl;
      if (nf != ++fileCount) {
          if (verbose) 
              std::cout << "RootInputDesc::setFileList number of files opened " 
                        << nf << " != filecount " << fileCount << std::endl;
         return numEvents ;
      }
     }
    else return numEvents ;
   }

  // Make sure the file is open
  if (m_chain->GetFile()->IsOpen() != kTRUE)
   {
    if (verbose) 
        std::cout << "RootInputDesc::setFileList failed to open" << std::endl;
    return numEvents ;
   }

  // TODO : those branch names could depend
  // on the data type
  // If necessary, set up the TChainIndex
  if (!m_chain->GetTreeIndex()) 
   {
    m_chain->BuildIndex("m_runId", "m_eventId") ;
    m_runEvtIndex = m_chain->GetTreeIndex();
   }

  // Get a pointer to the TTree within this data file
  // To do this we will first need to "rediscover" the directory name
  std::string treeName = std::string(m_chain->GetTree()->GetName());
  if (treeName != m_tree)
   {
    if (verbose) 
        std::cout << "RootInputDesc::setFileList failed to find tree " 
                  << m_tree << std::endl;
    return numEvents ;
   }
        
  // To set the data pointer we need to find the name of the branch contained
  // where we note that in Glast each file contains one branch
  TObjArray * branches = m_chain->GetListOfBranches() ;
  std::string branchName = "" ;
  int idx ;
  for( idx = 0 ; idx < branches->GetEntries() ; idx++ )
   {
    TBranch * branch = dynamic_cast<TBranch*>(branches->At(idx)) ;
    if (branch)
     {
      branchName = std::string(branch->GetName()) ;
      if (branchName == m_branch) break ;
     }
   }
  if (branchName!=m_branch) { 
    if (verbose) 
        std::cout << "RootInputDesc::setFileList failed to find branch " 
                  << m_branch << std::endl;
     return numEvents ; 
   }

  // Set the data pointer to receive the data and we are ready to go
  m_chain->SetBranchAddress(branchName.data(),m_dataObject) ;

  numEvents = m_chain->GetEntries() ;

  // Finally restore the context we entered with
  saveDir->cd() ;

  return numEvents ;
 }
    
bool RootInputDesc::fileExists( const std::string & filename )
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

TObject * RootInputDesc::getEvent( int index )
 {
  TObject * dataPtr = 0 ;
 
  // Save the current directory for the ntuple writer service
  TDirectory *saveDir = gDirectory;	

  if (m_chain) 
   {
    int ret = m_chain->GetEntry(index) ;
    dataPtr = *m_dataObject ;
   } 

  // Restore context we entered with
  saveDir->cd();

  return dataPtr ;
 }



TObject * RootInputDesc::getEvent( int runNum, int evtNum )
 {
  TObject * dataPtr = 0 ;
 
  // Save the current directory for the ntuple writer service
  TDirectory * saveDir = gDirectory ;	

  if (m_chain) 
   {
    if (m_runEvtIndex) m_chain->SetTreeIndex(m_runEvtIndex);
    int ret = m_chain->GetEntryWithIndex(runNum,evtNum) ;
    dataPtr = *m_dataObject ;
   } 

  // Restore context we entered with
  saveDir->cd() ;

  return dataPtr ;
 }

 bool RootInputDesc::checkEventAvailability( Long64_t index ) {
  TDirectory * saveDir = gDirectory ;	
  if (!m_chain) return false;
  if (index < 0) return false;
  if (index > m_chain->GetEntries()) return false;
  return true;
 }

 bool RootInputDesc::checkEventAvailability( int runNum, int evtNum ) {
  TDirectory * saveDir = gDirectory ;	
  if (!m_chain) return false;

  if (m_runEvtIndex) m_chain->SetTreeIndex(m_runEvtIndex);
  Long64_t readInd = m_chain->GetEntryNumberWithIndex(runNum,evtNum);
  saveDir->cd();
  if ((readInd<0)||(readInd>=m_chain->GetEntries())) return false;
  return true;
 }

void RootInputDesc::clearEvent()
 {
  TObject * dataPtr = *m_dataObject ;
  if (dataPtr)
   { dataPtr->Clear() ; }
 }