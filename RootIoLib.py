#$Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['RootIo'])	
    env.Tool('addLibrary', library = env['rootLibs'])
    env.Tool('addLibrary', library = env['rootGuiLibs'])

def exists(env):
    return 1;
