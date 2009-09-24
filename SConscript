# -*- python -*-
# $Header: /nfs/slac/g/glast/ground/cvs/RootIo/SConscript,v 1.5 2008/10/13 21:30:09 glastrm Exp $
# Authors: Heather Kelly <heather@milkyway.gsfc.nasa.gov>, David Chamont <chamont@poly.in2p3.fr>
# Version: RootIo-21-12-00-gr2
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('RootIoLib', depsOnly = 1)

RootIo = libEnv.SharedLibrary('RootIo', listFiles(['src/*.cxx','src/Dll/*.cxx']))

# Have yet to convert "macro_append RootIo_cppflags " $(RootPolicy_cppflags) $(root_packages_include)" from requirements file.  Not sure what the values of the variables are

progEnv.Tool('RootIoLib')

if baseEnv['PLATFORM'] == 'win32':
	progEnv.AppendUnique(CPPDEFINES = 'GLEAM')
	progEnv.AppendUnique(CPPDEFINES = '__i386')
	progEnv.AppendUnique(CPPDEFINES = 'EFC_FILTER')
	progEnv.AppendUnique(CPPDEFINES = '_WIN32')

test_RootIo = progEnv.GaudiProgram('test_RootIo', listFiles(['src/test/*.cxx']), test = 1)

progEnv.Tool('registerObjects', package = 'RootIo', libraries = [RootIo], testApps = [test_RootIo])
