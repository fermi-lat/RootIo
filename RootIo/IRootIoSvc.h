/** 
* @file IRootIoSvc.h
* @brief definition of the interface for IRootIoSvc
*
*  $Header: /nfs/slac/g/glast/ground/cvs/RootIo/RootIo/IRootIoSvc.h,v 1.4 2005/01/25 19:12:16 heather Exp $
*/
#ifndef _H_IRootIoSvc
#define _H_IRootIoSvc

// includes
#include "GaudiKernel/IInterface.h"
#include <TChain.h>

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IRootIoSvc("RootIoSvc", 2 , 0); 

/** 
* \class IRootIoSvc
* \brief The RootIoSvc gaudi service interface
*
* \author Heather Kelly heather@lheapop.gsfc.nasa.gov
* 
* $Header: /nfs/slac/g/glast/ground/cvs/RootIo/RootIo/IRootIoSvc.h,v 1.4 2005/01/25 19:12:16 heather Exp $
*/
class  IRootIoSvc : virtual public IInterface {
public:
    
    virtual Long64_t getEvtMax() = 0;
    virtual void setRootEvtMax(Long64_t max) = 0;
    virtual void setRootTimeMax(unsigned int max) = 0;

    virtual Long64_t index() = 0;
    virtual bool setIndex(Long64_t i) = 0;
    virtual void setActualIndex(Long64_t i) = 0;

    virtual void registerRootTree(TChain *ch) = 0;
    virtual bool setRunEventPair(std::pair<int,int> ids) = 0;
    virtual std::pair<int,int> runEventPair() = 0;

    virtual bool useIndex() = 0;
    virtual bool useRunEventPair() = 0;

    virtual int getAutoSaveInterval() = 0;

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IRootIoSvc; }
  

};

#endif  // _H_IRootIoSvc
