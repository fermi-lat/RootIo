/** 
* @file RootIo_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/Dll/RootIo_load.cxx,v 1.11 2004/11/24 14:16:32 chamont Exp $
*/

#include "GaudiKernel/DeclareFactoryEntries.h"

DECLARE_FACTORY_ENTRIES(RootIo) {
    DECLARE_ALGORITHM( FhSetAlg );
    DECLARE_ALGORITHM( FhDemoGetAlg );
    DECLARE_ALGORITHM( FhDemoCaloSetAlg );
    DECLARE_ALGORITHM( mcRootWriterAlg );
    DECLARE_ALGORITHM( mcRootReaderAlg );
    DECLARE_ALGORITHM( digiRootWriterAlg );
    DECLARE_ALGORITHM( digiRootReaderAlg );
    DECLARE_ALGORITHM( reconRootWriterAlg );
    DECLARE_ALGORITHM( reconRootReaderAlg );
    DECLARE_ALGORITHM( relationRootWriterAlg );
    DECLARE_ALGORITHM( relationRootReaderAlg );
    DECLARE_TOOL( FhTool );
    DECLARE_TOOL( RootIoRandom );
    DECLARE_SERVICE( RootIoSvc );


}
  
