/** 
* @file EventCollectionMgr.cxx
* @brief definition of the class EventCollectionMgr
*
*  $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/EventCollectionMgr.cxx,v 1.1 2007/07/26 16:40:57 heather Exp $
*  Original author: Heather Kelly heather@lheapop.gsfc.nasa.gov
*/
#ifndef EventCollectionMgr_cxx
#define EventCollectionMgr_cxx

#include "EventCollectionMgr.h"
#include "facilities/Util.h"
#include "metaRootData/Component.h"
#include "TROOT.h"
#include <iostream>



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



EventCollectionMgr::~EventCollectionMgr() 
{
    m_treeCol.clear();
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



/// Reading
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
    // prime the pump so the TTrees are loaded
    m_pointerSkimRead.read(0);
    return stat;
}

Long64_t EventCollectionMgr::getEventIndex(const std::string &treeName, Long64_t index) {
    // Load the event
    m_pointerSkimRead.read(index);
    UInt_t componentInd = m_pointerSkimRead.componentIndex(treeName);
    if (componentInd == 0xFFFFFFFF) return -1;  // No such component in this event collection
    const Component* comp = m_pointerSkimRead.getComponent(componentInd);
    if (!comp) return -1;
    return (comp->eventPointer().eventIndex());
}

TTree* EventCollectionMgr::getTree(const std::string &treeName) {
    UInt_t componentInd = m_pointerSkimRead.componentIndex(treeName);
    if (componentInd == 0xFFFFFFFF) return 0;  // No such component in this event collection
    const Component* comp = m_pointerSkimRead.getComponent(componentInd);
    if (!comp) return 0;
    return (comp->getTree());
}

#endif