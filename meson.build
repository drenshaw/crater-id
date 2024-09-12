project('c++ crater-id', 'cpp',
  default_options : [
    'c_std=c17', 
    'cpp_std=c++17', 
    'warning_level=1',
    # 'b_sanitize=address,undefined,memory',
    'buildtype=debug',
    'debug=true',
  ],
  version : '1.0.0',
  license : 'MIT',
)

add_global_arguments('-DSOME_TOKEN=value', 
                     '-DDEBUG',
                     '-DBOOST_LOG_DYN_LINK',
                     language : 'cpp',
)

incdir = include_directories('include', '/usr/local/include', '/usr/include')
cpp = meson.get_compiler('cpp')
# dir_base = meson.current_source_dir()
# dir_lib = join_paths(dir_base, 'lib')

subdir('include')
subdir('src')
subdir('test')


pkg_mod = import('pkgconfig')
pkg_mod.generate(
  # libraries : cidlib,
  version : '1.0',
  name : 'libcid',
  filebase : 'cid',
  # dependencies: lib_vtk,
  description : 'A C++ library for crater identification.'
)

