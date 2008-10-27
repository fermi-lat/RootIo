#$Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/RootIo/RootIoLib.py,v 1.1 2008/08/15 21:42:36 ecephas Exp $
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['RootIo'])	
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
