/** 
* @file RootOutputDesc.h
* @brief definition of the class RootOutputDesc
*        This class is used to set up and handle the actual root output files
*
*  $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/RootOutputDesc.h,v 1.1 2007/07/26 16:40:57 heather Exp $
*  Original author: Heather Kelly heather@lheapop.gsfc.nasa.gov
*/

#ifndef RootOutputDesc_h
#define RootOutputDesc_h

//#include "GaudiKernel/Property.h"
#include "TTree.h"
#include "TFile.h"
#include <string>

class RootOutputDesc
 {
  public :
  
    RootOutputDesc  ( const std::string& outputFile, 
       const std::string & treeName, 
	   int compressionLevel, const std::string& treeTitle, bool verbose=false ) ;
    ~RootOutputDesc() ; 

    bool openFile();

    bool closeFile();

    bool setupBranch(const std::string& name, const std::string &classname, void* branchAddr, 
        int bufSize=64000, int splitLevel=1);


    // Methods to return information about the TTree being accessed
    const std::string & getTreeName() const { return m_treeName; }

    TTree * getTree() { return m_tree ; }

	TFile* getCurrentFile();

    Long64_t getEventCounter() const { return m_eventCounter; };

    bool getUpdated() const { return m_updated; };
    void setUpdated(bool u) { m_updated = u; };

	bool LoadTree(Long64_t ievent){ m_tree->LoadTree(ievent); };

	bool fillTree(int autoSaveInterval=1000); 

  private :

    bool fileExists( const std::string & filename ) ;

	std::string m_fileName ; /// The name of the original output file opened for writing
	std::string m_treeName;  /// Name of the TTree
	int m_compressionLevel;  /// Compression level for writing
    TTree * m_tree ;         /// Pointer to the TTree
	std::string m_treeTitle; /// Title to pass to TTree ctor
	bool m_verbose;  /// set the chattiness for debug statements
    Long64_t m_eventCounter;  // Count number of events filled to the TTree so far
    bool m_updated;
    
 } ;



#endif



