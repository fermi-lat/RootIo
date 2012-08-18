#$Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/RootIo/RootIoLib.py,v 1.2 2008/10/27 19:05:56 ecephas Exp $
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['RootIo'])	
        if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
	    env.Tool('findPkgPath', package = 'RootIo') 
    env.Tool('addLibrary', library = env['rootLibs'])
    env.Tool('addLibrary', library = env['rootGuiLibs'])
    env.Tool('commonRootDataLib')
    env.Tool('rootUtilLib')
    env.Tool('digiRootDataLib')
    env.Tool('RootConvertLib')
    env.Tool('TriggerLib')
    env.Tool('OnboardFilterTdsLib')
    env.Tool('LdfEventLib')
    env.Tool('AncillaryDataEventLib')
    env.Tool('ntupleWriterSvcLib')
    env.Tool('gcrSelectRootDataLib')
def exists(env):
    return 1;
