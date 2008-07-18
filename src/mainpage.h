/** @mainpage package RootIo
* @author Heather Kelly
*
* @section intro Introduction
* This package defines the Gaudi algorithms that read and write ROOT files.
* The following Gaudi algorithms are defined in this package:
* - FhSetAlg: prepare the common data for the file headers to be written.
* - mcRootWriterAlg: writes a Monte Carlo ROOT file using the MC data on the TDS.
* - mcRootReaderAlg: read a Monte Carlo ROOT file, and puts MC data on the TDS.
* - digiRootWriterAlg: writes a Digi ROOT file using digi data from the TDS.
* - digiRootReaderAlg: reads in a Digi ROOT file and puts Digi data on the TDS.
* - reconRootWriterAlg: writes a Recon ROOT file using recon data from the TDS.
* - reconRootReaderAlg: reads in a Recon ROOT file and puts Recon data on the TDS.
* - relationRootWriterAlg: writes a relation ROOT file using the relations stored on the TDS.  This algorithm must be run
   after all other ROOT writing algorithms.
* - relationRootReaderAlg: reads in the Relation ROOT file and puts relations on the TDS.  This algorithm must run after
*  all other ROOT reading algorithms.
*
* Services:  RootIoSvc
* Allows for random event access via the GUI (when using RootDisplay) and can control
* the event loop depending upon the min number of events in all read root files
*
* One can also find a tool, FhTool, which give access to the current file
* headers, and a few helper classes (FhSystemEnv, FhJobOptions, FhCmtConfig
* and FhChronoStatTable) which ease the interpretation and use of the
* headers data, as demonstrated by FhDemoCaloSetAlg. 
*
* <hr>
* @section jobOptions jobOptions
* @param RootIoSvc.MaxTreeSize
*  Defaults to 500 GB
*  Size in MB for Trees, when a Tree hits this size, the current ROOT file
*  is closed and a new one is opened for writing.  This parameter is set for
*  ALL Trees - there is no way to assign it on a per tree basis
* @param RootIoSvc.AutoSaveInterval
*  Defaults to 1000
*  Controls how many events to process before actually writing data to a ROOT
*  file - until then, data is held in a buffer that may be lost if the job
*  crashes.
* @param RootIoSvc.RebuildIndex
*  Defaults to true
*  Determines if the TTreeIndex is automatically rebuilt from scratch whether
*  or not one pre-exists in the input ROOT files.
* @param mcRootWriterAlg.mcRootFile
*  Name of the output MC ROOT file, default mc.root
* @param mcRootWriterAlg.splitMode 
*  Specifies the split level when writing the TTree, default 1
* @param mcRootWriterAlg.bufferSize
* @param mcRootWriterAlg.compressionLevel
* @param mcRootWriterAlg.clearOption
*  Default is "", set to "ALL" to clear full McEvent with each iteration - this will impact performance
* @param mcRootReaderAlg.mcRootFileList
*  List of input MC ROOT file(s)
* @param mcRootReaderAlg.clearOption
*  Default is "", set to "ALL" to clear full McEvent with each iteration - this will impact performance
* @param digiRootWriterAlg.digiRootFile
*  Name of the output Digi ROOT file, default digi.root
* @param digiRootWriter.splitMode
* @param digiRootWriter.bufferSize
* @param digiRootWriter.compressionLevel
* @param digiRootReader.digiRootFileList
*  List of input Digi ROOT file(s)
* @param reconRootWriterAlg.reconRootFile
*  Name of the output Recon ROOT file
* @param reconRootWriterAlg.splitMode
* @param reconRootWriterAlg.bufferSize
* @param reconRootWriterAlg.compressionLevel
* @param reconRootReaderAlg.reconRootFileList
*  List of input Recon ROOT file(s), default recon.root
* @param relationRootReaderAlg.relationRootFileList
*  List of input Relation ROOT file(s), defaults to relations.root
* @param relationRootWriterAlg.relationRootFile
*  Name of output relation ROOT file.
*
* @section Tests Tests and Demonstrations
*
* This package contains one main test application, test_RootIo.
* It can be used to test both reading and writing of Monte Carlo, digi, and 
* reconstruction ROOT files - depending upon the settings in the 
* src/test/jobOptions.txt file. 
* To set for writing, leave the writeOptions.txt include file uncommented and
* comment out the readOptions.txt include file.  Do the reverse to run the
* reading test.
* When the test routine is set for writing, a special algorithm is used to create
* fake test data, rather than forcing the test routine to use the full simulation
* to create the data.
*
* On top of that, there are two options files, called writeHeadersOptions.txt
* and readHeadersOptions.txt, which specifically focus on the file headers.
* They make use of FhDemoCaloSetAlg, which demonstrates how one can add
* private specific data to the file headers.
*
* <hr>
* @section notes release notes
* release.notes
* <hr>
* @section requirements requirements
* @verbinclude requirements
*
*/

