/** 
* @file RootIo_load.cpp
* @brief This is needed for forcing the linker to load all components
* of the library.
*
*  $Header: /nfs/slac/g/glast/ground/cvs/RootIo/src/Dll/RootIo_load.cxx,v 1.14 2006/11/13 10:38:56 claval Exp $
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
    DECLARE_ALGORITHM( gcrSelectRootWriterAlg );
    DECLARE_TOOL( FhTool );
#if 0 //THB: not needed if no use of random numbers
    DECLARE_TOOL( RootIoRandom );
#endif
    DECLARE_SERVICE( RootIoSvc );
    DECLARE_SERVICE( ntupleReaderSvc );


}
  
