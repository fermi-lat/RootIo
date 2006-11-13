/** 
* @file IRootIoSvc.h
* @brief definition of the interface for IRootIoSvc
*
*  $Header: /nfs/slac/g/glast/ground/cvs/RootIo/RootIo/IRootIoSvc.h,v 1.6 2006/03/30 20:50:11 heather Exp $
*/
#ifndef _H_IRootIoSvc
#define _H_IRootIoSvc

// includes
#include "GaudiKernel/IInterface.h"
#include <TChain.h>
#include <string>

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_IRootIoSvc("RootIoSvc", 3 , 0); 

/** 
* \class IRootIoSvc
* \brief The RootIoSvc gaudi service interface
*
* \author Heather Kelly heather@lheapop.gsfc.nasa.gov
* 
* $Header: /nfs/slac/g/glast/ground/cvs/RootIo/RootIo/IRootIoSvc.h,v 1.6 2006/03/30 20:50:11 heather Exp $
*/
class  IRootIoSvc : virtual public IInterface {
public:

    //virtual bool setRootFile(const char *mc, const char *digi, const char *rec) = 0;

    virtual bool setRootFile(const char *mc, const char *digi, 
                             const char *rec, const char *gcr) = 0;
			     
    virtual std::string getMcFile() const = 0;
    virtual std::string getDigiFile() const = 0;
    virtual std::string getReconFile() const = 0;
    virtual std::string getGcrFile() const = 0;
    virtual bool fileChange() const = 0;
    
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
