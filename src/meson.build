cid_sources = [
  'main.cpp', 
  'conics.cpp', 
  'quadrics.cpp',
  'camera.cpp',
  'io.cpp', 
  'crater-id.cpp', 
  'combinatorics.cpp', 
  'vector_math.cpp',
  'visuals.cpp',
]
boost_dep = dependency('boost')
eigen_dep = dependency('eigen3')
libspatial_dep = dependency('libspatialindex')
cv_dep = dependency('opencv4')
# -DVTK_DIR="<path to my vtk build>"
vtk_dep = dependency('VTK', required: true)
cpp = meson.get_compiler('cpp')
igl_dep = cpp.find_library('libigl', dirs: '/home/dqr0509/dev/libigl/lib', required: false)
# igl_dep = dependency('igl', required: true)
all_deps = [boost_dep, eigen_dep, libspatial_dep, cv_dep, vtk_dep]
executable( 
  'demo', 
  cid_sources, 
  include_directories:incdir, 
  dependencies : [boost_dep, eigen_dep, libspatial_dep, vtk_dep, cv_dep, igl_dep],
  override_options : ['cpp_std=c++17'])
cidlib = shared_library(
  'crater-id',
  cid_sources,
  include_directories : incdir,
  dependencies : all_deps,
  install : true
)
