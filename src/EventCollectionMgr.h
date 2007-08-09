/** 
* @file EventCollectionMgr.h
* @brief definition of the class EventCollectionMgr
*        This class is used to set up and handle the event collection (meta) root files
*
*  $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/EventCollectionMgr.h,v 1.1 2007/08/08 14:14:45 heather Exp $
*  Original author: Heather Kelly heather@lheapop.gsfc.nasa.gov
*/

#ifndef EventCollectionMgr_h
#define EventCollectionMgr_h

#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include <string>
#include <vector>
#include <map>
#include "metaRootData/PointerSkim.h"


class EventCollectionMgr
 {
  public :
  
      EventCollectionMgr() : m_fileNameWrite(""),m_fileNameRead(""),
          m_fileWrite(0),m_fileRead(0),m_verbose(false), m_eventCounter(0),
          m_compChainCol(0), m_masterChain(0)  
      { };

      ~EventCollectionMgr() ; 

      /// Writing 
      bool initWrite(const std::string &fileName="meta.root",
          const std::string &options="RECREATE", bool verbose=false);

      bool initOutputFile();

    UInt_t addComponent(const std::string &compName, TTree *t);

	bool fillEvent(); 
    
    bool fillMeta();

    /// Reading
    bool initRead(const std::string &fileName="meta.root", bool verbose=false);
    Long64_t getNumEntries() { return (m_pointerSkimRead.entries()); };
    Long64_t getEventIndex(const std::string &treeName, Long64_t index);
    TChain* getChainByType(const std::string &treeName);
    int setIndex();


  private :

	std::string m_fileNameWrite ; /// The name of the original output file opened for writing
    std::string m_fileNameRead;
    std::string m_outputOptions;
    TFile *m_fileWrite, *m_fileRead;
    PointerSkim m_pointerSkimWrite, m_pointerSkimRead;
    std::vector<TTree*> m_treeCol;
	bool m_verbose;  /// set the chattiness for debug statements
    Long64_t m_eventCounter;  // Count number of events filled to the TTree so far
    TObjArray *m_compChainCol; // List of component TChains for reading
    TChain *m_masterChain;  // Master TChain for reading
    std::map<std::string, TVirtualIndex*> m_chainIndexCol;
    
 } ;



#endif



