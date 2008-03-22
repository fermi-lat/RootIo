/** @file INtupleReaderSvc.h
    @brief declare abstract INtupleReaderSvc

    $Header: /nfs/slac/g/glast/ground/cvs/ntupleWriterSvc/ntupleWriterSvc/INtupleReaderSvc.h,v 1.17 2007/12/14 20:14:14 heather Exp $
*/
#ifndef _H_INtupleReaderSvc_
#define _H_INtupleReaderSvc_

#include "GaudiKernel/IInterface.h"
#include <string>

// Declaration of the interface ID ( interface id, major version, minor version) 
static const InterfaceID IID_INtupleReaderSvc("INtupleReaderSvc",  1 ,1); 

/*! @class INtupleReaderSvc
 @brief Proper Gaudi abstract interface class for the ntupleReaderSvc 
*/
class INtupleReaderSvc : virtual public IInterface
{  

public:

    /// setup before event processing - required for Gaudi services
    virtual StatusCode initialize ()=0;
    
    /// cleanup after event processing - required for all Gaudi services
    virtual StatusCode finalize ()=0;
 
    //! return false when no more events
    virtual bool nextEvent(bool checkTds=true) = 0;

    //! return false if event idx does not exist
    virtual bool getEvent(long long idx, bool checkTds=true) = 0;

    virtual bool getEvent(unsigned int run, unsigned int evt, bool checkTds=true) = 0;

    virtual long long numEvents() const = 0;

    virtual bool addFile(const std::string & fileName, const std::string & treeName) = 0;

    virtual bool clearChain() = 0;

    virtual std::string getItemType(const std::string& itemName) const = 0; 

    /** @brief Set the pointer to the value of an existing item.
    @param tupleName  Name of the tuple
    @param itemName   Name of the item (or perhaps blank, see below)
    @param pointer    point pointer that will be set  

    Expect to throw exception if not found.
    @return a string describing the type (client then must cast)

    Note that the RootTupleSvc subclass implements a variation that it the name is empty,
    it will set the pointer to the TTree.
    
    */
    virtual std::string getItem(const std::string& itemName, void*& pointer)const =0;

    /// Retrieve interface ID
    static const InterfaceID& interfaceID() { return IID_INtupleReaderSvc; }

};


#endif // _H_INtupleReaderSvc
