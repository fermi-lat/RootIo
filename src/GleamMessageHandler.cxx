/** 
* @file GleamMessageHandler.cxx
* @brief definition of the class GleamMessageHandler
*
*  $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/GleamMessageHandler.cxx,v 1.1 2008/10/13 20:29:04 heather Exp $
*  Original author: Heather Kelly heather@lheapop.gsfc.nasa.gov
*/

#include "GleamMessageHandler.h"
#include <iostream>

Bool_t GleamMessageHandler::Notify() {

  switch (fMessId) {
    case ROOT_WARNING :
      std::cout << "GleamMessageHandler ROOT WARNING for Class: "
                << fClass->GetName() << std::endl;
      ++m_warningCount;
      m_warning = true;
      break;
    case ROOT_ERROR :
      std::cout << "GleamMessageHandler ROOT ERROR for Class: "
                << fClass->GetName() << std::endl;
      ++m_errorCount;
      m_error = true;
      break;
    case ROOT_SYSERROR :
      std::cout << "GleamMessageHandler ROOT SYSERROR for Class: "
                << fClass->GetName() << std::endl;
      ++m_syserrorCount;
      m_syserror = true;
      break;
    case ROOT_FATAL :
      std::cout << "GleamMessageHandler ROOT FATAL for Class: "
                << fClass->GetName() << std::endl;
      ++m_fatalCount;
      m_fatal = true;
      break;
    default :
      break;
  }
  return kTRUE;
}

