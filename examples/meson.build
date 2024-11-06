
# app_files = files('conics.cpp', 'math_utils.cpp')
main_file = files('main.cpp')

mainexe = executable(
  'example_cid', 
  [main_file],
  include_directories : incdir,
  dependencies: [
    cv_dep,
    eigen_dep,
  ],
  link_with : cidlib,
)
test( 
  'main_cid', 
  testexe,
  is_parallel:true
)

test('math', testexe, suite: 'MathTest')
test('conic', testexe, suite: 'ConicTest')
