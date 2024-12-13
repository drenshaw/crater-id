cid_sources = [
  'conics.cpp', 
  'quadrics.cpp',
  'camera.cpp',
  'io.cpp', 
  'crater-id.cpp', 
  'combinatorics.cpp', 
  'math_utils.cpp',
  'visuals.cpp',
  'navigation.cpp',
]
cid_main = ['main.cpp']
boost_dep = dependency('boost', required: false)
eigen_dep = dependency('eigen3', required: true)
libspatial_dep = dependency('libspatialindex', required: false)
cv_dep = dependency('opencv4', required: true)
# vtk_dep = dependency('VTK', method: 'cmake',
#                    components: ['CommonCore', 'CommonColor', 'CommonExecutionModel', 'FiltersCore', 'RenderingCore', 'InteractionStyle', 'RenderingUI'],
#                    modules: ['VTK::CommonCore', 'VTK::CommonColor', 'VTK::CommonExecutionModel', 'VTK::FiltersCore', 'VTK::RenderingCore', 'VTK::InteractionStyle', 'VTK::RenderingUI'])
# -DVTK_DIR="<path to my vtk build>"
vtk_dep = dependency('VTK', required: false)
cpp = meson.get_compiler('cpp')
igl_dep = cpp.find_library('libigl', dirs: '/home/dqr0509/dev/libigl/lib', required: false)
# igl_dep = dependency('igl', required: false)
all_deps = [
  # boost_dep, 
  eigen_dep, 
  libspatial_dep, 
  cv_dep, 
  vtk_dep, 
  boost_dep,
  igl_dep,
]
# executable( 
#   'demo', 
#   cid_sources + cid_main, 
#   include_directories:incdir, 
#   dependencies : [
#     boost_dep, 
#     eigen_dep, 
#     libspatial_dep, 
#     vtk_dep, 
#     cv_dep, 
#     # igl_dep
#   ],
#   override_options : ['cpp_std=c++17'])
cidlib = shared_library(
  'crater-id',
  cid_sources,
  include_directories : incdir,
  dependencies : all_deps,
  install : true
)
