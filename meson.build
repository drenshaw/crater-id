project('c++ crater-id', 'cpp',
  default_options : ['c_std=c17', 'cpp_std=c++17'],
  version : '1.0.0',
  license : 'MIT')
add_global_arguments('-DSOME_TOKEN=value', language : 'cpp')
# glib_dep = dependency('glib-2.0')

incdir = include_directories('include')

subdir('include')
subdir('src')
subdir('test')

pkg_mod = import('pkgconfig')
pkg_mod.generate(
                #  libraries : cidlib,
                 version : '1.0',
                 name : 'libcid',
                 filebase : 'cid',
                 description : 'A C++ library for crater identification.'
)