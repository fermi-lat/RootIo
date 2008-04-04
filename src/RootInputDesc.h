/** 
* @file RootInputDesc.h
* @brief definition of the class RootInputDesc
*        This class is used to set up and handle the actual root IO
*
*  $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/RootInputDesc.h,v 1.8 2008/03/24 15:25:43 heather Exp $
*  Original author: Heather Kelly heather@lheapop.gsfc.nasa.gov
*/

#ifndef RootInputDesc_h
#define RootInputDesc_h

#include "GaudiKernel/Property.h"
#include "TChain.h"
#include "TFile.h"
#include "TObject.h"
#include <string>

class RootInputDesc
 {
  public :
  
    RootInputDesc
     ( const StringArrayProperty & fileList, 
       const std::string & treeName, 
       const std::string & branchName, bool verbose=false ) ;

    RootInputDesc
     (  TChain *t,
       const std::string & treename,
       const std::string & branchName, bool verbose=false);

    ~RootInputDesc() ; 

    /// Methods to return information about the TTree/TChain being accessed
    const StringArrayProperty & getFileList() const { return m_fileList ; }
    const std::string & getTreeName() const { return m_tree ; }
    const std::string & getBranchName() const { return m_branch ; }
    const int getNumEvents() const { return m_numEvents ; }

    TChain * getTChain() { return m_chain ; }
    TObject * getTObject() { return *m_dataObject ; }

    /// Methods to handle reading and clearing events
    TObject * getEvent( int index ) ;
    TObject * getEvent( int runNum, int evtNum ) ;
    bool checkEventAvailability( Long64_t index );
    bool checkEventAvailability( int runNum, int evtNum );
    void clearEvent() ;

    /// Method to change the list of files in this TChain
    Long64_t  setFileList( const StringArrayProperty & fileList, bool verbose = false ) ;

    /// Setup to read from an event collection
    Long64_t setEventCollection( );

  private :

    bool fileExists( const std::string & filename ) ;
    //bool checkForEnvVar( const StringArrayProperty & fileList);

    StringArrayProperty m_fileList ; // A list of files (fully qualified) associated to this TChain
    std::string m_tree ;             // The name of the tree being accessed
    std::string m_branch ;           // Branch name for this tree
    TChain * m_chain ;               // Pointer to the TChain
    TTree * m_treePtr ;              // For use with EventCollections
    TObject * * m_dataObject ;       // A pointer to the pointer to the data
    Long64_t m_numEvents ;           // Number of events in current TChain
    bool m_verbose;
    TVirtualIndex *m_runEvtIndex;    /// Save the RunId/EventId Index in case other indices are in use such as CompositeEventList
    
 } ;



#endif



