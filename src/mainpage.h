/** @mainpage package RootIo
  @author Heather Kelly

  @section intro Introduction
  This package defines the Gaudi algorithms that read and write ROOT files.
  The following Gaudi algorithms are defined in this package:
  - mcRootWriterAlg: writes a Monte Carlo ROOT file using the MC data on the TDS  It accepts the following parameters:
    -# mcRootFile - name of the output ROOT file
    -# splitMode - default (1)
    -# bufferSize - default (64000)
    -# compressionLevel - default (1) a value 0-9 where 9 is maximum compression
  - mcRootReaderAlg: read a Monte Carlo ROOT file, and puts MC data on the TDS
    -# mcRootFile - name of the input ROOT file
  - digiRootWriterAlg:  writes a Digi ROOT file using digi data from the TDS.
    -# digiRootFile - name of the output Digi ROOT file
    -# splitMode - default (1)
    -# bufferSize - default (64000)
    -# compressionLevel - default (1) a value [0,9] where 9 is maximum compression
  - digiRootReaderAlg:  reads in a Digi ROOT file and puts Digi data on the TDS
    -# digiRootFile - name of in the input Digi ROOT file

  <hr>
  @section notes release notes
  release.notes
  @section requirements requirements
  @include requirements

*/

