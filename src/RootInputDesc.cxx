/** 
* @file RootInputDesc.cxx
* @brief definition of the class RootInputDesc
*
*  $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/RootInputDesc.cxx,v 1.1 2007/05/09 21:32:39 usher Exp $
*  Original author: Heather Kelly heather@lheapop.gsfc.nasa.gov
*/

#include "RootInputDesc.h"
#include "TROOT.h"

RootInputDesc::RootInputDesc(const StringArrayProperty& fileList, 
                       const std::string&         tree, 
                       const std::string&         branch) : 
    m_tree(tree), m_branch(branch), m_chain(0), m_dataObject(0)
{
    // Set up the input file list
    m_numEvents = setFileList(fileList);
}
    
RootInputDesc::~RootInputDesc() 
{
    if (m_chain)      delete m_chain;
    if (m_dataObject) delete m_dataObject;
}

int RootInputDesc::setFileList(const StringArrayProperty& fileList)
{
    int numEvents = 0;

    m_fileList = fileList;
 
    // Save the current directory for the ntuple writer service
    TDirectory *saveDir = gDirectory;	

    // Check to see if we are changing files (which means a TChain exists)
    if (m_chain)
    {
        // Close the current file
        m_chain->GetFile()->Close();

        delete m_chain;
        delete m_dataObject;
    }

    // Get a new instance of a TChain for processing this file
    m_chain      = new TChain(m_tree.c_str());
    m_dataObject = new TObject*;

    *m_dataObject = 0;

    // Add the files to the TChain
    typedef std::vector<std::string>::const_iterator StringVecIter;

    int fileCount = 0;
    for(StringVecIter fileListItr = m_fileList.value().begin(); 
        fileListItr != m_fileList.value().end(); fileListItr++)
    {
        std::string fileName = *fileListItr;
        if (fileExists(fileName))
        {
            int nf = m_chain->Add(fileName.c_str());

            if (nf != ++fileCount) return numEvents;
        }
        else return numEvents;
    }

    // Make sure the file is open
    if (m_chain->GetFile()->IsOpen() != kTRUE)
    {
        return numEvents;
    }

    // If necessary, set up the TChainIndex
    if (!m_chain->GetTreeIndex()) 
    {
        m_chain->BuildIndex("m_runId", "m_eventId");
    }

    // Get a pointer to the TTree within this data file
    // To do this we will first need to "discover" the directory name
    std::string treeName = std::string(m_chain->GetTree()->GetName());

    if (treeName != m_tree)
    {
        return numEvents;
    }
        
    // To set the data pointer we need to find the name of the branch contained
    // where we note that in Glast each file contains one branch
    TObjArray*  branches   = m_chain->GetListOfBranches();
    std::string branchName = "";

    for(int idx = 0; idx < branches->GetEntries(); idx++)
    {
        TBranch* branch = dynamic_cast<TBranch*>(branches->At(idx));

        if (branch)
        {
            branchName = std::string(branch->GetName());

            if (branchName == m_branch) break;
        }
    }

    if (branchName != m_branch) 
    {
        return numEvents;
    }

    // Set the data pointer to receive the data and we are ready to go
    m_chain->SetBranchAddress(branchName.data(), m_dataObject);

    // Change the file status to unchanged

    numEvents = m_chain->GetEntries();

    // Finally restore the context we entered with
    saveDir->cd();

    return numEvents;
}
    
bool RootInputDesc::fileExists(const std::string& filename)
{
    bool   fileExists = false;
    TFile* file       = TFile::Open(filename.c_str());

    if (file)
    {
        if (!file->IsZombie()) 
        {
            file->Close();
            fileExists = true;
        }

        delete file;
    }
    
    return fileExists;
}  

TObject* RootInputDesc::getEvent(int index)
{
    TObject* dataPtr = 0;
 
    // Save the current directory for the ntuple writer service
    TDirectory *saveDir = gDirectory;	

    if (m_chain) 
    {
        int ret = m_chain->GetEntry(index);
        dataPtr = *m_dataObject;
    }

    // Restore context we entered with
    saveDir->cd();

    return dataPtr;
}

TObject* RootInputDesc::getEvent(int runNum, int evtNum)
{
    TObject* dataPtr = 0;
 
    // Save the current directory for the ntuple writer service
    TDirectory *saveDir = gDirectory;	

    if (m_chain) 
    {
        int ret = m_chain->GetEntryWithIndex(runNum, evtNum);
        dataPtr = *m_dataObject;
    }

    // Restore context we entered with
    saveDir->cd();

    return dataPtr;
}

void RootInputDesc::clearEvent()
{
    if (TObject* dataPtr = *m_dataObject)
    {
        dataPtr->Clear();
    }
}
