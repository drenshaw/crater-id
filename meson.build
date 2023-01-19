project('c++ crater-id', 'cpp',
  default_options : ['c_std=c17', 'cpp_std=c++17'],
  version : '1.0.0',
  license : 'MIT')
add_global_arguments('-DSOME_TOKEN=value', language : 'cpp')
# glib_dep = dependency('glib-2.0')

incdir = include_directories('include')
cpp = meson.get_compiler('cpp')
# dir_base = meson.current_source_dir()
# dir_lib = join_paths(dir_base, 'lib')

# mylib_inc = include_directories('include', 'src')
# lib_vtk = cpp.find_library('vtk-9.2', dirs: ['/usr/local/lib','/usr/local','/usr/local/include'], required: false)
# lib_glfw = cpp.find_library('glfw3', dirs: dir_lib, required: true)

subdir('include')
subdir('src')
subdir('test')

pkg_mod = import('pkgconfig')
pkg_mod.generate(
                #  libraries : cidlib,
                 version : '1.0',
                 name : 'libcid',
                 filebase : 'cid',
                #  dependencies: lib_vtk,
                 description : 'A C++ library for crater identification.'
)