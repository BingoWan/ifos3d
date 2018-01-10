import os

inc_path = ['-isystem', 'libcseife']
lib_path = ['./../libcseife']
env = Environment(ENV = os.environ)

# compiler flags
#env['F90'] = ['mpiifort']
#env['CC'] = ['mpiicc']
#env['LINK'] = ['mpiicc', '-wd11021']
#env.Append(CCFLAGS = ['-Wall', '-O3'] + inc_path, F90FLAGS = ['-O2'])
#env.Append(LIBPATH = lib_path)
#env.Append(LIBS = ['m', 'cseife'])

env['F90'] = ['mpif90']
env['CC'] = ['mpicc']
env.Append(CCFLAGS = ['-O2', '-Wall', '-OPT:IEEE_arith=2'] + inc_path, F90FLAGS = ['-O2', '-OPT:IEEE_arith=2'] )
env.Append(LIBPATH = lib_path)
env.Append(LIBS = ['m', 'cseife', 'm_slave'])

# set the sub directories (key, value), where value is the name of directory#
# please make sure the source code are in src subdirectory
dirlist = [
   ('lib', 'lib'),
   ('bin', 'bin'),
   ('src', 'src'),
]
dirs = dict(dirlist)

model_scr = ARGUMENTS.get('model_scr')    
SConscript('src/SConscript', variant_dir = 'build', duplicate = 0, exports = 'env model_scr')
