/** 
* @file EventCollectionMgr.cxx
* @brief definition of the class EventCollectionMgr
*
*  $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/EventCollectionMgr.cxx,v 1.2 2007/08/09 17:17:08 heather Exp $
*  Original author: Heather Kelly heather@lheapop.gsfc.nasa.gov
*/
#ifndef EventCollectionMgr_cxx
#define EventCollectionMgr_cxx

#include "EventCollectionMgr.h"
#include "facilities/Util.h"
#include "metaRootData/Component.h"
#include "TROOT.h"
#include <iostream>


EventCollectionMgr::~EventCollectionMgr() 
{
    m_treeCol.clear();
    if (m_compChainCol) m_compChainCol->Delete();
    if (m_masterChain) delete m_masterChain;
}

// Writing Methods 

bool EventCollectionMgr::initWrite(const std::string &fileName, const std::string &options, bool verbose) {
    bool stat = true;
    m_fileNameWrite = fileName;
    facilities::Util::expandEnvVar(&m_fileNameWrite);
    m_outputOptions = options;
    m_verbose = verbose;
    return stat;
}

bool EventCollectionMgr::initOutputFile() {
    bool stat = true;
    m_fileWrite = m_pointerSkimWrite.makeFile(m_fileNameWrite.c_str(), m_outputOptions.c_str());
    if (!m_fileWrite->IsOpen()) {
        stat = false;
        std::cout << "EventCollectionMgr Error opening ROOT file " << m_fileNameWrite << std::endl;
        return stat;
    }
    return stat;
}


UInt_t EventCollectionMgr::addComponent(const std::string &compName, TTree *t) {
    unsigned int ret = m_pointerSkimWrite.addComponent(compName);
    m_treeCol.push_back(t);
    return ret;
}

bool EventCollectionMgr::fillEvent() {
    bool stat = true;
    try{
    TDirectory *saveDir = gDirectory;
    // Need to call PointerSkim::makeFile after the AddComponent calls
    if (!m_fileWrite)
        if (!initOutputFile()) throw;
    if (m_fileWrite->TestBits(TFile::kWriteError)) {
        throw;
    }
    m_fileWrite->cd();
    std::vector<TTree*>::iterator treeIt;
    Long64_t numBytes;
    for (treeIt=m_treeCol.begin(); treeIt != m_treeCol.end(); treeIt++) 
        numBytes = (*treeIt)->LoadTree(m_eventCounter);
    m_pointerSkimWrite.fillEvent(m_treeCol);
    saveDir->cd();
    ++m_eventCounter;
    } catch(...) {
        std::cerr << "Error filling Meta ROOT file" << std::endl;
        std::cerr.flush();
        throw;
    }
    return stat;
}

    bool EventCollectionMgr::fillMeta() {
        bool stat = true;
        try {
            TDirectory *saveDir = gDirectory;
            m_fileWrite->cd();
            m_pointerSkimWrite.fillMeta();
            m_fileWrite->Write(0,TObject::kWriteDelete);
            m_fileWrite->Close();
            saveDir->cd();
        } catch(...) {
            std::cerr << "Failed final write to meta ROOT file" << std::endl; 
            std::cerr.flush(); 
            throw;
        }
        return stat;
    }



/// Reading Methods

bool EventCollectionMgr::initRead(const std::string &fileName, bool verbose) {
    bool stat = true;
    m_fileNameRead = fileName;
    facilities::Util::expandEnvVar(&m_fileNameRead);
    m_verbose = verbose;
    m_fileRead = m_pointerSkimRead.openFile(m_fileNameRead.c_str());
    if (!m_fileRead->IsOpen()) {
        stat = false;
        std::cout << "EventCollectionMgr Error opening ROOT file " << m_fileNameRead << std::endl;
        return stat;
    }

    // Set up our master TChain and the component TChains
    if (!m_compChainCol) m_compChainCol = new TObjArray();
    m_masterChain = m_pointerSkimRead.buildLinks(m_compChainCol);
    TIter itr(m_compChainCol);
    TChain *curChain;
    while ( (curChain = (TChain*)itr.Next()) ) {
        curChain->SetBranchStatus("*", 1);
        m_chainIndexCol[curChain->GetName()] = curChain->GetTreeIndex();
    }

   
    return stat;
}

Long64_t EventCollectionMgr::getEventIndex(const std::string &treeName, Long64_t index) {
    // make sure the pointer skim index is the active TVirtualIndex
    setIndex();
    Long64_t retVal = m_masterChain->LoadTree(index);  // Causes all component TChains to be loaded correctly via PointerIndex
    if (retVal < 0) return retVal;
    TChain *compChain = getChainByType(treeName);
    if (!compChain) return -1;
    return (compChain->GetReadEntry());
}

int EventCollectionMgr::setIndex() {
    int numSet=0;
    std::map<std::string, TVirtualIndex*>::iterator mapIt;
    for (mapIt=m_chainIndexCol.begin();mapIt != m_chainIndexCol.end(); mapIt++) {
        TChain *curChain = getChainByType(mapIt->first);
        if (curChain) {
            curChain->SetTreeIndex(mapIt->second);
            numSet++;
        }
    }
    return numSet;
}



TChain* EventCollectionMgr::getChainByType(const std::string& treeName) {
    TIter itr(m_compChainCol);
    TChain *curChain;
    while ((curChain = (TChain*)itr.Next())) {
        if (curChain->GetName() == treeName) return curChain;
    }
    return 0;
}

#endif