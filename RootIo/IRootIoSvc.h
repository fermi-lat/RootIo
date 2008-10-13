/** 
* @file IRootIoSvc.h
* @brief definition of the interface for IRootIoSvc
*
*  $Header: /nfs/slac/g/glast/ground/cvs/RootIo/RootIo/IRootIoSvc.h,v 1.16 2008/03/06 04:29:01 heather Exp $
*/

#ifndef _H_IRootIoSvc
#define _H_IRootIoSvc

//
// The IRootIoSvc was originally created to provide a "Runable" interface for Gaudi,
// so to control the event loop when we are reading in data from root files.
// There are three "clients" we know about : HepRepSvc, the classes xxReaderRootAlg
// (which we will call the "readers"), the classes xxWriterRootAlg (which we will
// call the "writers").
//
// HepRepSvc is our interface to the event display, so that users can do the following:
// * set the run/event or index to specify what event to read in next
// * reset at any time the input ROOT files (mc, digi, recon, gcr)
//
// The readers:
// * register the number of events in their respective TTree with RootIoSvc
// * register the TTree * with RootIoSvc
//   ( used to check validity of requested run/event pairs
//     or indices from HepRepSvc )
// * retrieve a run/event pair that may have be set via HepRepSvc
// * retrieve an index that may be requested via HepRepSvc
// * check to see if a new input file has been requested,
//   such as from HepRepSvc if a new input file is to be opened,
//   the reader gets the filename from RootIoSvc via the getXXFile() call
//
// The writers retrieve the AutoSaveInterval that can be set via JO for
// RootIoSvc, to specify that we save every N events to file.
//

// includes
#include "GaudiKernel/IInterface.h"
#include <TChain.h>
#include <string>
#include "rootUtil/CompositeEventList.h"

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IRootIoSvc("RootIoSvc",4,1) ; 

/** 
* \class IRootIoSvc
* \brief The RootIoSvc gaudi service interface
*
* \author Heather Kelly heather@lheapop.gsfc.nasa.gov
* 
* $Header: /nfs/slac/g/glast/ground/cvs/RootIo/RootIo/IRootIoSvc.h,v 1.16 2008/03/06 04:29:01 heather Exp $
*/

class  IRootIoSvc : virtual public IInterface
 {
  public:

    //====================
    // For HepRepSvc
    //====================
    
    virtual bool setRootFile
     ( const char * mc, const char * digi, 
       const char * rec, const char * gcr="" ) = 0 ;
       
    virtual bool setIndex( Long64_t ) = 0 ;
    virtual Long64_t index() = 0 ;
    
    virtual bool setRunEventPair( std::pair<int,int> ) = 0 ;
    virtual std::pair<int,int> runEventPair() = 0 ;
  
  
    //====================
    // For readers
    //====================
    
    // file list
    virtual bool setFileList( const std::string & type, const StringArrayProperty & fileList ) = 0 ;
    virtual StringArrayProperty getFileList( const std::string & type) const = 0 ;
    virtual bool appendFileList( StringArrayProperty & fileList, const std::string & fileName ) = 0 ;

    virtual StatusCode prepareRootInput
     ( const std::string & type, 
       const std::string & tree,
       const std::string & branch,
       const StringArrayProperty & fileList) = 0 ;

    virtual StatusCode closeInput(const std::string& type) = 0;
       
    virtual TObject * getNextEvent( const std::string & type ) = 0 ;
    
    // the number of events we want to read
    // based on ApplicationMgr.EvtMax and number available in files
    virtual Long64_t getEvtMax() = 0 ;
    virtual void setEvtMax( Long64_t max ) = 0 ;


    //====================
    // For writers
    //====================
    
    virtual TTree * prepareRootOutput
     ( const std::string & type,
       const std::string & fileName,
       const std::string & treeName,
       int compressionLevel,
       const std::string & treeTitle ) = 0 ;

    //[David] hidden to check if it is used or not
    //virtual TTree * getTree( const std::string & type ) = 0 ;

    virtual StatusCode setupBranch( const std::string & type, const std::string & branchName, 
        const std::string & classname, void * branchAddr, int bufSize=64000, int splitLevel =1 ) = 0 ;

    virtual StatusCode fillTree( const std::string & type ) = 0;

    virtual int getAutoSaveInterval() = 0 ;

    virtual StatusCode closeFile( const std::string & type ) = 0 ;

    virtual CompositeEventList* getCel() = 0;


    //====================
    // Gaudi machinery
    //====================

    static const InterfaceID & interfaceID()
     { return IID_IRootIoSvc ; }
  
 } ;


#endif  // _H_IRootIoSvc


