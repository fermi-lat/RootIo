/** 
* @file RootInputDesc.h
* @brief definition of the class RootInputDesc
*        This class is used to set up and handle the actual root IO
*
*  $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/RootIoSvc.cxx,v 1.27 2007/03/19 01:14:35 heather Exp $
*  Original author: Heather Kelly heather@lheapop.gsfc.nasa.gov
*/

#ifndef RootInputDesc_h
#define RootInputDesc_h

#include <vector>
#include <string>

#include "GaudiKernel/Property.h"

#include "TChain.h"
#include "TFile.h"
#include "TObject.h"

class RootInputDesc
{
public:
    RootInputDesc(const StringArrayProperty& fileList, 
               const std::string&         tree, 
               const std::string&         branch);
    ~RootInputDesc(); 

    // Methods to return information about the TTree/TChain being accessed
    const StringArrayProperty& getFileList()   const {return m_fileList;}
    const std::string&         getTreeName()   const {return m_tree;}
    const std::string&         getBranchName() const {return m_branch;}
    const int                  getNumEvents()  const {return m_numEvents;}

    TChain*                    getTChain()           {return m_chain;}
    TObject*                   getTObject()          {return *m_dataObject;}

    // Methods to handle reading and clearing events
    TObject*                   getEvent(int index);
    TObject*                   getEvent(int runNum, int evtNum);
    void                       clearEvent();

    // Method to change the list of files in this TChain
    int  setFileList(const StringArrayProperty& fileList);

private:
    bool fileExists(const std::string& filename);

    StringArrayProperty m_fileList;        // A list of files (fully qualified) associated to this TChain
    std::string         m_tree;            // The name of the tree being accessed
    std::string         m_branch;          // Branch name for this tree
    TChain*             m_chain;           // Pointer to the TChain
    TObject**           m_dataObject;      // A pointer to the pointer to the data
    int                 m_numEvents;       // Number of events in current TChain
};

#endif
