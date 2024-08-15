vtk_dep = dependency('VTK', required: true)
gtest_proj = subproject('gtest')
# gtest_dep = dependency('gtest', required: true)
gtest_dep = gtest_proj.get_variable('gtest_main_dep')
# gtest_dep = dependency('gtest', fallback : ['gtest', 'gtest_dep'])
exe = executable('gtest_project', 'cidtest.cpp',
  dependencies : [gtest_dep],
  install : true)

test('basic', exe)

testexe = executable(
  'test_cid', 
  sources:'cidtest.cpp',
  include_directories : incdir,
  dependencies: [all_deps, vtk_dep, ],
  link_with: cidlib,
  # source: cid_sources,
  # link_with : cidlib,
  # link_with : libcid,
)
test( 
  'test_cid', 
  testexe,
  is_parallel:true
)
