/** 
* @file IRootIoSvc.h
* @brief definition of the interface for IRootIoSvc
*
*  $Header$
*/
#ifndef _H_IRootIoSvc
#define _H_IRootIoSvc

// includes
#include "GaudiKernel/IInterface.h"


// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IRootIoSvc("RootIoSvc", 1 , 0); 

/** 
* \class IRootIoSvc
* \brief The RootIoSvc gaudi service interface
*
* \author Heather Kelly heather@lheapop.gsfc.nasa.gov
* 
* $Header$
*/
class  IRootIoSvc : virtual public IInterface {
public:
    
    virtual int getEvtMax() = 0;
    virtual void setRootEvtMax(unsigned int max) = 0;
    virtual void setRootTimeMax(unsigned int max) = 0;

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_IRootIoSvc; }
  

};

#endif  // _H_IRootIoSvc
