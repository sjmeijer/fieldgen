#
# Karin Lagergren
#
# To build debug executable: scons debug=1
# To make gprof executable: scons prof=1
# To clean away all files built by scons: scons -c
# "scons ." will build everything

cflags = ['-Wall']
if ARGUMENTS.get('debug', 0):
    cflags.append('-g')
if ARGUMENTS.get('prof', 0):
    cflags.append('-pg')
if ARGUMENTS.get('vf', 0):
    cflags.append('-DCALC_VFIELD')
if ARGUMENTS.get('chatty',0):
   cflags.append('-DCHATTY');
libpath = []
includepath = []
libs = '-lm'

Default("f3d_gretina")

env = Environment(
        CCFLAGS=[cflags, '-DGRETINA'], 
        CPPPATH = includepath,
        LIBS=libs,
        LIBPATH=libpath,
        RPATH=libpath,
        TARFLAGS = '-c -j -h')#create, bzip2, dereference symlinks
env_coax = Environment(
        CCFLAGS=cflags, 
        CPPPATH = includepath,
        LIBS=libs,
        LIBPATH=libpath,
        RPATH=libpath,
        TARFLAGS = '-c -j -h')#create, bzip2, dereference symlinks

fields3d = env.StaticObject(target = 'fields3d.o',
                            source = 'fields3d.c')
fields3d_coax = env_coax.StaticObject(target = 'fields3d_coax.o',
                                      source = 'fields3d.c')
signal_calc_util = env.StaticObject(target = 'signal_calc_util.o',
                                    source = 'signal_calc_util.c')
signal_calc_util_coax = env_coax.StaticObject(target = 'signal_calc_util_coax.o',
                                              source = 'signal_calc_util.c')
point = env.StaticObject(target = 'point.o',
                         source = 'point.c')
point_coax = env_coax.StaticObject(target = 'point_coax.o',
                                   source = 'point.c');

env.Program(target = 'f3d_gretina', 
            source = ['field_init_gretina.c', fields3d,
                      signal_calc_util,point])
env_coax.Program(target = 'f3d_coax', 
            source = ['field_init_coax.c', fields3d_coax,
                      signal_calc_util_coax,point_coax])
env_coax.Program(target = 'f3d_wcoax', 
            source = ['field_weird_coax.c', fields3d_coax,
                      signal_calc_util_coax,point_coax])
env.Alias('tar','fields3d.tar.bz2')



env.Tar(target = 'fields3d.tar.bz2',
        source = ['fields3d.c',
                  'signal_calc_util.c',
                  'point.c',                  
                  'field_init_gretina.c',
                  'field_init_coax.c',
                  'field_weird_coax.c',
                  'signal_calc_util.h',
                  'point.h',
                  'assym_detector.h',
                  'field_init.h',
                  'coax_geometry.dat',
                  'field_calc_setup.dat',
                  'field_calc_setup_A1.dat',
                  'field_calc_setup_A2.dat',
                  'field_calc_setup_B1.dat',
                  'field_calc_setup_B2.dat',
                  'field_calc_setup_coax.dat',
                  'field_calc_setup_wc.dat',
                  'geometry_setup.dat',
                  'weird_coax_geometry.dat',
                  'weird_coax_geometry2.dat',
                  'SConstruct'])
