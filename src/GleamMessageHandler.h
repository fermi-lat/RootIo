/** 
* @file GleamMessageHandler.h
* @brief definition of the class GleamMessageHandler
*        This class handles ROOT message to check for errors
*
*  $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/GleamMessageHandler.h,v 1.1 2008/10/13 20:29:04 heather Exp $
*  Original author: Heather Kelly heather625@gmail.com
*/

#ifndef GleamMessageHandler_h
#define GleamMessageHandler_h

#include "TClass.h"
#include "TMessageHandler.h"

class GleamMessageHandler : public TMessageHandler
 {
  public :

    typedef enum {
                 ROOT_WARNING = 1001,
                 ROOT_ERROR = 1002,
                 ROOT_SYSERROR = 1003,
                 ROOT_FATAL = 1004
    } MessageIds;
  
/*
    GleamMessageHandler(const GleamMessageHandler& other) //: TMessageHandler(other) {  
    {
     m_warningCount = other.m_warningCount;
     m_errorCount = other.m_errorCount;
     m_syserrorCount = other.m_syserrorCount;
     m_fatalCount = other.m_fatalCount;
     m_warning = other.m_warning;
     m_error = other.m_error;
     m_syserror = other.m_syserror;
     m_fatal = other.m_fatal;
   };
*/

    GleamMessageHandler(const TClass* cl, Bool_t derived = kTRUE) :
       TMessageHandler(cl, derived),
       m_warningCount(0),m_errorCount(0),m_syserrorCount(0),m_fatalCount(0),
       m_warning(false),m_error(false),m_syserror(false),m_fatal(false) { };

    GleamMessageHandler(const char* cl, Bool_t derived = kTRUE) :
       TMessageHandler(cl, derived), 
       m_warningCount(0),m_errorCount(0),m_syserrorCount(0),m_fatalCount(0),
       m_warning(false),m_error(false),m_syserror(false),m_fatal(false) { };

    virtual ~GleamMessageHandler() {}; 

    virtual Bool_t Notify();

    virtual void Print(Option_t* option="") const {
        TMessageHandler::Print(option); 
    }

    Int_t getWarningCount() const { return m_warningCount; }
    Int_t getErrorCount() const { return m_errorCount; }
    Int_t getFatalCount() const { return m_fatalCount; }
    Int_t getSysErrorCount() const { return m_syserrorCount; }

    Bool_t getWarning()const { return m_warning; }
    Bool_t getError() const { return m_error; }
    Bool_t getSysError() const { return m_syserror; }
    Bool_t getFatal() const { return m_fatal; }

    void resetFlags() {
        m_warning = false;
        m_error = false;
        m_syserror = false;
        m_fatal = false;
    }

    protected :
        Int_t m_warningCount;
        Int_t m_errorCount;
        Int_t m_syserrorCount;
        Int_t m_fatalCount;
        Bool_t m_warning; 
        Bool_t m_error;
        Bool_t m_syserror;
        Bool_t m_fatal;


    //ClassDef(GleamMessageHandler,0) // Gleam Message Handler
 } ;

#endif



