# vtk_dep = dependency('VTK', required: true)
# gtest_proj = subproject('gtest')
# gtest_dep = gtest_proj.get_variable('gtest_main_dep')
# # gtest_dep = dependency('gtest', required: true)
# gtest_dep = gtest_proj.get_variable('gtest_main_dep')
# # gtest_dep = dependency('gtest', fallback : ['gtest', 'gtest_dep'])
# exe = executable('gtest_project', 'cidtest.cpp',
#   dependencies : [gtest_dep],
#   install : true)

# test('basic', exe)

gtest_dep = dependency('gtest', required : false)
if not gtest_dep.found()
  error('MESON_SKIP_TEST: gtest not installed.')
endif
# gtest_nomain = dependency('gtest', main : false, method : 'system')

# # e = executable('testprog', 'test.cc', dependencies : [gtest_dep])
# # test('gtest test', e, protocol : 'gtest')

# e = executable('testprog_nomain', 'cidtest.cpp', dependencies : [gtest_nomain])
# test('gtest nomain test', e, protocol : 'gtest')

# catch_dep = dependency('catch2', required : false)
# if not catch_dep.found()
#   catch_lib = subproject('catch2')
#   catch_dep = catch_lib.get_variable('catch2_dep')
# endif

# app_files = files('conics.cpp', 'math_utils.cpp')
main_file = files('test_main.cpp')
test_files = files([
  'test_math.cpp',
  'test_io.cpp',
  'test_conics.cpp', 
  'test_quadrics.cpp',
  'test_camera.cpp',
  'test_combinatorics.cpp',
  'test_invariants.cpp',
])

testexe = executable(
  'test_cid', 
  [main_file, test_files],
  include_directories : incdir,
  dependencies: [
    # catch_dep, 
    gtest_dep,
    cv_dep,
    # gtest_nomain,
  ],
  link_with : cidlib,
  # link_with : libcid,
)
test( 
  'test_cid', 
  testexe,
  is_parallel:true
)

test('math', testexe, suite: 'MathTest')
test('conic', testexe, suite: 'ConicTest')
