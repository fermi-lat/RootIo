/** @mainpage package RootIo
* @author Heather Kelly
*
* @section intro Introduction
* This package defines the Gaudi algorithms that read and write ROOT files.
* The following Gaudi algorithms are defined in this package:
* - mcRootWriterAlg: writes a Monte Carlo ROOT file using the MC data on the TDS  It accepts the following parameters:
* - mcRootReaderAlg: read a Monte Carlo ROOT file, and puts MC data on the TDS
* - digiRootWriterAlg:  writes a Digi ROOT file using digi data from the TDS.
* - digiRootReaderAlg:  reads in a Digi ROOT file and puts Digi data on the TDS
* - reconRootWriterAlg: writes a Recon ROOT file using recon data from the TDS.
* - reconRootReaderAlg: reads in a Recon ROOT file and puts Recon data on the TDS.
*
* <hr>
* @section jobOptions jobOptions
* @param mcRootWriterAlg.mcRootFile
*  Name of the output MC ROOT file, default mc.root
* @param mcRootWriterAlg.splitMode 
*  Specifies the split level when writing the TTree, default 1
* @param mcRootWriterAlg.bufferSize
* @param mcRootWriterAlg.compressionLevel
* @param mcRootReaderAlg.mcRootFile
*  Name of the input MC ROOT file
* @param digiRootWriterAlg.digiRootFile
*  Name of the output Digi ROOT file, default digi.root
* @param digiRootWriter.splitMode
* @param digiRootWriter.bufferSize
* @param digiRootWriter.compressionLevel
* @param digiRootReader.digiRootFile
*  Name of the input Digi ROOT file
* @param reconRootWriterAlg.reconRootFile
*  Name of the output Recon ROOT file
* @param reconRootWriterAlg.splitMode
* @param reconRootWriterAlg.bufferSize
* @param reconRootWriterAlg.compressionLevel
* @param reconRootReaderAlg.reconRootFile
*  Name of the input Recon ROOT file, default recon.root
*
* @section Tests Tests
* This package contains one test application, test_RootIo.
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
* <hr>
* @section notes release notes
* release.notes
* <hr>
* @section requirements requirements
* @verbinclude requirements
*
*/

