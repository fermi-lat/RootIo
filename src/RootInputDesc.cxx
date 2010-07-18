/** 
* @file RootInputDesc.cxx
* @brief definition of the class RootInputDesc
*
*  $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/RootInputDesc.cxx,v 1.23 2010/04/07 14:09:06 heather Exp $
*  Original author: Heather Kelly heather@lheapop.gsfc.nasa.gov
*/

#include "RootInputDesc.h"
#include "facilities/Util.h"
#include "TROOT.h"
#include <vector>
#include <iostream>
#include "TTreeIndex.h"
#include "TChainIndex.h"


// ctor for vanilla reading from a list of ROOT files
// Setup the TChain* in this case, as we may have multiple files
RootInputDesc::RootInputDesc( const StringArrayProperty& fileList, 
                              const std::string&         tree, 
                              const std::string&         branch, 
                              TObject**                  branchPtr,
                              bool                       rebuildIndex, 
                              bool                       verbose )
   : m_tree(tree), 
     m_branch(branch), 
     m_chain(0), 
     m_dataObject(branchPtr), 
     m_myDataObject(true), 
     m_verbose(verbose),
     m_rebuildIndex(rebuildIndex), 
     m_runEvtIndex(0),
     m_chainIndex(0)
 {
    // Check ownership of data object
    if (branchPtr) m_myDataObject = false;

    // Set up the input file list
    m_numEvents = setFileList(fileList, m_verbose) ;

    return;
 }

 // ctor for Event Collections
 // Set up the TChain* as we will retrieve a TChain* from the CompositeEventList via the CelManager
 RootInputDesc::RootInputDesc(TChain*            t,
                              const std::string& treeName,
                              const std::string& branchName, 
                              TObject**          branchPtr,
                              bool               rebuildIndex, 
                              bool               verbose)
     : m_tree(treeName), 
       m_branch(branchName), 
       m_chain(t), 
       m_dataObject(branchPtr), 
       m_myDataObject(true), 
       m_verbose(verbose), 
       m_rebuildIndex(rebuildIndex), 
       m_runEvtIndex(0),
       m_chainIndex(0)
 {
    // Check ownership of data object
    if (branchPtr) m_myDataObject = false;

    m_numEvents = setEventCollection();

    return;
 }

    
RootInputDesc::~RootInputDesc() 
{
    if (m_chainIndex) {
        delete m_chainIndex;
        m_chainIndex = 0;
    }

    if (m_chain)
    {
/*
        if (TVirtualIndex* indexPtr = m_chain->GetTreeIndex())
        {
            if (TChainIndex* chainPtr = dynamic_cast<TChainIndex*>(indexPtr))
            {
                for (int idx = 0; idx < chainPtr->GetN(); idx++)
                {
//                    TVirtualIndex* iptr = chainPtr->fEntries[idx]->fTreeIndex;
                    int j = 0;
                }
            }
        }
*/
        m_chain->GetFile()->Close() ;
        delete m_chain ;
    }

    // If we own it we need to clean up
    if (m_myDataObject && m_dataObject) 
    {
        if (*m_dataObject) delete *m_dataObject;
        delete m_dataObject ;
    }
}

Long64_t RootInputDesc::setEventCollection( )
{
    Long64_t numEvents = -1;
    // Save the current directory for the ntuple writer service
    TDirectory * saveDir = gDirectory ;	

    if (!m_chain) return numEvents;

    // If we own the data object then clean up
    if (m_myDataObject)
    {
        if (m_dataObject) 
        {
            if (*m_dataObject) delete *m_dataObject;
            delete m_dataObject ;
            m_dataObject = 0;
        }

        m_dataObject = new TObject * ;
        *m_dataObject = 0 ;
    }

    // Check that we have the right TTree
    // To do this we will first need to "rediscover" the directory name
    std::string treeName = std::string(m_chain->GetName());
    if (treeName != m_tree)
    {
        std::cout << "RootInputDesc:setEventCollection failed to find "
                  << "tree " << m_tree << std::endl;
        return numEvents ;
    }

    // To set the data pointer we need to find the name of the branch contained
    // where we note that in Glast each file contains one branch
    TObjArray*  branches   = m_chain->GetListOfBranches() ;
    std::string branchName = "" ;
    for(int idx = 0 ; idx < branches->GetEntries() ; idx++ )
    {
        TBranch* branch = dynamic_cast<TBranch*>(branches->At(idx)) ;
        if (branch)
        {
            branchName = std::string(branch->GetName()) ;
            if (branchName == m_branch) break ;
        }
    }
    if (branchName!=m_branch) 
    { 
        std::cout << "RootInputDesc:setEventCollection failed to find "
                  << "branch " << m_branch << std::endl;
        return numEvents ; 
    }

    // Set the data pointer to receive the data and we are ready to go
    m_chain->SetBranchAddress(branchName.data(), m_dataObject) ;

    // If necessary, set up the TChainIndex
    if (m_rebuildIndex)
    {
        for (int iTree = 0; iTree < m_chain->GetNtrees(); iTree++) 
        {
            m_chain->LoadTree((m_chain->GetTreeOffset())[iTree]);
            m_chain->GetTree()->SetTreeIndex(0);
        }
        m_chain->LoadTree(0);
        m_chainIndex = new TChainIndex(m_chain,"m_runId","m_eventId");
        //m_chain->BuildIndex("m_runId", "m_eventId") ;
    }

    // This will force reading in the file headers to determine the number
    // of events
    numEvents = m_chain->GetEntries() ;

    // Finally restore the context we entered with
    saveDir->cd() ;

    return numEvents ;
}

 

Long64_t RootInputDesc::setFileList( const StringArrayProperty & fileList, bool verbose )
{
    Long64_t numEvents = -1;
 
    // bool stat = checkForEnvVar(fileList);
    try {
        std::vector<std::string> expandedList;
        facilities::Util::expandEnvVarList(fileList.value(),expandedList);
        m_fileList.setValue(expandedList);
    } catch(...) {
        std::cout << "RootInputDesc::setFileList call to expandEnvVarList failed, taking the JO"
                  << " parameter as is and continuing on" << std::endl;
        m_fileList = fileList;
    }
    // If the checkForEnvVar method fails, just set the local m_fileList to the fileList input
    //if (!stat) m_fileList = fileList;
 
    // Save the current directory for the ntuple writer service
    TDirectory * saveDir = gDirectory ;	

    // Check to see if we are changing files (which means a TChain exists)
    if (m_chain)
    {
        if (m_chainIndex) delete m_chainIndex;
        m_chainIndex = 0; 

        m_chain->GetFile()->Close() ;
        delete m_chain ;
        m_chain = 0;
    }

    // If we own the data object must set up
    if (m_myDataObject)
    {
        // Do we already have one?
        if (m_dataObject) 
        {
            if (*m_dataObject) delete *m_dataObject;
            delete m_dataObject ;
            m_dataObject = 0;
        }
        m_dataObject = new TObject * ;
        *m_dataObject = 0 ;
    }

    // Get a new instance of a TChain for processing this file
    m_chain = new TChain(m_tree.c_str()) ;

    // Add the files to the TChain
    typedef std::vector<std::string>::const_iterator StringVecIter ;
    for (StringVecIter fileListItr = m_fileList.value().begin() ; 
                       fileListItr != m_fileList.value().end() ;
                       fileListItr++ )
    {
        std::string fileName = *fileListItr ;
        facilities::Util::expandEnvVar(&fileName);
        if (fileExists(fileName, m_tree.c_str()))
        {
            int nf = m_chain->Add(fileName.c_str()) ;
            if (verbose) 
                std::cout << "RootInputDesc::setFileList opening: " << fileName 
                          << std::endl;
            if (nf <= 0) {
                std::cout << "RootInputDesc::setFileList failed to TChain::Add " 
                          << fileName << " returned " << nf 
                          << " Continuing on" << std::endl;
            //return numEvents ;
            }
        }
        else 
        {
            std::cout << "RootInputDesc::setFileList Failed to open and add " 
                      << fileName << " the file may have no TTree entries"
                      << " Continuing on" << std::endl;
            //return numEvents ;
        }
    }

    // Make sure the file is open
    if (!m_chain->GetFile()) 
    {
        std::cout << "RootInputDesc::setFileList, no TFile available, no events"
                  << std::endl;
        return numEvents;
    }
    if (m_chain->GetFile()->IsOpen() != kTRUE)
    {
        std::cout << "RootInputDesc::setFileList failed to open" << std::endl;
        return numEvents ;
    }

    // Get a pointer to the TTree within this data file
    // To do this we will first need to "rediscover" the directory name
    std::string treeName = std::string(m_chain->GetTree()->GetName());
    if (treeName != m_tree)
    {
        std::cout << "RootInputDesc::setFileList failed to find tree " 
                  << m_tree << std::endl;
        return numEvents ;
    }
        
    // To set the data pointer we need to find the name of the branch contained
    // where we note that in Glast each file contains one branch
    TObjArray*  branches   = m_chain->GetListOfBranches() ;
    std::string branchName = "" ;
    for(int idx = 0 ; idx < branches->GetEntries() ; idx++ )
    {
        TBranch * branch = dynamic_cast<TBranch*>(branches->At(idx)) ;
        if (branch)
        {
            branchName = std::string(branch->GetName()) ;
            if (branchName == m_branch) break ;
        }
    }
    if (branchName!=m_branch) 
    {
        std::cout << "RootInputDesc::setFileList failed to find branch " 
                  << m_branch << std::endl;
        return numEvents ; 
    }

    // Set the data pointer to receive the data and we are ready to go
    m_chain->SetBranchAddress(branchName.data(),m_dataObject) ;

    // TODO : those branch names could depend
    // on the data type
    // If necessary, set up the TChainIndex
    if (m_rebuildIndex) 
    {
        int iTree;
        for (iTree = 0; iTree < m_chain->GetNtrees(); iTree++) 
        {
            m_chain->LoadTree((m_chain->GetTreeOffset())[iTree]);
            if (m_chain->GetTree()->GetTreeIndex()) 
            {
                TVirtualIndex* temp = m_chain->GetTree()->GetTreeIndex();
                delete temp;
                m_chain->GetTree()->SetTreeIndex(0);
            }
        }
        m_chain->LoadTree(0);
        m_chainIndex = new TChainIndex(m_chain,"m_runId","m_eventId");
    }

    // This call forces the file headers to be read in
    numEvents = m_chain->GetEntries() ;

    // Finally restore the context we entered with
    saveDir->cd() ;

    return numEvents ;
 }
    
bool RootInputDesc::fileExists( const std::string & filename, const char* treeName )
{
    bool fileExists = false ;
    TFile * file = TFile::Open(filename.c_str()) ;
    if (file)
    {
        if (!file->IsZombie()) 
        {
            if (treeName != 0) 
            {
                TTree *t = (TTree*)file->Get(treeName);
                if (!t) return false;
                if (t->GetEntries() <= 0) return false;
            }
            file->Close() ;
            fileExists = true ;
        }
        delete file ;
    }
    return fileExists ;
}  

TObject * RootInputDesc::getEvent( Long64_t index )
{
    TObject * dataPtr = 0 ;
 
    // Save the current directory for the ntuple writer service
    TDirectory *saveDir = gDirectory;	

    if (m_chain) 
    {
        int ret = m_chain->GetEntry(index) ;
        if (ret <= 0) 
        {
            std::cout << "RootInputDesc::getEvent bad read" << std::endl;
            return 0;
        }
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
        if (!m_chainIndex) {
            m_chainIndex = new TChainIndex(m_chain,"m_runId","m_eventId");
            if (!m_chainIndex) {
                std::cout << "RootInputDes::checkEventAvailability "
                          << "Failed to create TChainIndex" << std::endl;
                return false;
            }
        }
        Long64_t ind  = m_chainIndex->GetEntryNumberWithIndex(runNum,evtNum) ;
        dataPtr = getEvent(ind);
    }
    saveDir->cd();
    return dataPtr;
}

 bool RootInputDesc::checkEventAvailability( Long64_t index ) 
 {
    //TDirectory * saveDir = gDirectory ;	
    if (!m_chain) return false;
    if (index < 0) return false;
    if (index > m_chain->GetEntries()) return false;
    return true;
}

 bool RootInputDesc::checkEventAvailability( int runNum, int evtNum ) 
 {
    TDirectory * saveDir = gDirectory ;	
    if (!m_chain) return false;

    if (!m_chainIndex) {
        m_chainIndex = new TChainIndex(m_chain,"m_runId","m_eventId");
        if (!m_chainIndex) {
            std::cout << "RootInputDes::checkEventAvailability "
                      << "Failed to create TChainIndex" << std::endl;
            return false;
        }
    }
    
    Long64_t readInd = m_chainIndex->GetEntryNumberWithIndex(runNum,evtNum);
    saveDir->cd();
    if ((readInd<0)||(readInd>=m_chain->GetEntries())) return false;
    return true;
}

unsigned int RootInputDesc::getIndexByEventID( int runNum, int evtNum ) 
 {
    // Save the current directory for the ntuple writer service
    TDirectory * saveDir = gDirectory ;	

    unsigned int noIndex = (unsigned)-1;

    if (m_chain) 
    {
        if (!m_chainIndex) {
            m_chainIndex = new TChainIndex(m_chain,"m_runId","m_eventId");
            if (!m_chainIndex) {
                std::cout << "RootInputDes::checkEventAvailability "
                          << "Failed to create TChainIndex" << std::endl;
                saveDir->cd();
                return noIndex;
            }
        }
        Long64_t ind  = m_chainIndex->GetEntryNumberWithIndex(runNum,evtNum) ;
        saveDir->cd();

        return ind;
    }
}


void RootInputDesc::clearEvent()
{
    if (TObject* dataPtr = *m_dataObject)
    { 
        dataPtr->Clear() ; 
    }
}

bool RootInputDesc::setBranchStatus(const std::string& branch, int status) {
    if (!m_chain) return false;
    TDirectory * saveDir = gDirectory ;	
    unsigned int foundFlag;
    m_chain->SetBranchStatus(branch.c_str(),status,&foundFlag);
    saveDir->cd();
    return ((foundFlag==0) ? false : true);
}
