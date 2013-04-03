# -*- python -*-
# $Header: /nfs/slac/g/glast/ground/cvs/RootIo/SConscript,v 1.64 2013/02/20 19:22:55 usher Exp $
# Authors: Heather Kelly <heather@milkyway.gsfc.nasa.gov>, David Chamont <chamont@poly.in2p3.fr>
# Version: RootIo-25-07-03

Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('addLinkDeps', package='RootIo', toBuild='component')

RootIo =libEnv.ComponentLibrary('RootIo',
				listFiles(['src/*.cxx']))

progEnv.Tool('RootIoLib')

if baseEnv['PLATFORM'] == 'win32':
	progEnv.AppendUnique(CPPDEFINES = ['GLEAM'])
	progEnv.AppendUnique(CPPDEFINES = ['__i386'])
	progEnv.AppendUnique(CPPDEFINES = ['EFC_FILTER'])
	progEnv.AppendUnique(CPPDEFINES = ['_WIN32'])

test_RootIo = progEnv.GaudiProgram('test_RootIo',
                                   listFiles(['src/test/*.cxx']), test = 1,
				   package='RootIo')

progEnv.Tool('registerTargets', package = 'RootIo',
             libraryCxts=[[RootIo,libEnv]],testAppCxts=[[test_RootIo,progEnv]],
             includes = listFiles(['RootIo/*.h']),
	     jo = listFiles(['src/*.txt', 'src/test/*.txt']))
