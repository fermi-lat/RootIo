/** 
* @file RootIo_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header$
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(RootIo) {
    DECLARE_ALGORITHM( mcRootWriterAlg );
    DECLARE_ALGORITHM( mcRootReaderAlg );
}
  