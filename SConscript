# -*- python -*-
# $Header: /nfs/slac/g/glast/ground/cvs/RootIo/SConscript,v 1.25 2009/12/02 19:18:38 heather Exp $
# Authors: Heather Kelly <heather@milkyway.gsfc.nasa.gov>, David Chamont <chamont@poly.in2p3.fr>
# Version: RootIo-24-06-01
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('RootIoLib', depsOnly = 1)

RootIo = libEnv.SharedLibrary('RootIo', listFiles(['src/*.cxx','src/Dll/*.cxx']))

progEnv.Tool('RootIoLib')

if baseEnv['PLATFORM'] == 'win32':
	progEnv.AppendUnique(CPPDEFINES = ['GLEAM'])
	progEnv.AppendUnique(CPPDEFINES = ['__i386'])
	progEnv.AppendUnique(CPPDEFINES = ['EFC_FILTER'])
	progEnv.AppendUnique(CPPDEFINES = ['_WIN32'])

test_RootIo = progEnv.GaudiProgram('test_RootIo',
                                   listFiles(['src/test/*.cxx']), test = 1)

progEnv.Tool('registerTargets', package = 'RootIo',
             libraryCxts=[[RootIo,libEnv]],testAppCxts=[[test_RootIo,progEnv]],
             includes = listFiles(['RootIo/*.h']))
