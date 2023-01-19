cid_sources = ['main.cpp', 'conics.cpp', 'io.cpp', 'crater-id.cpp', 'combinatorics.cpp', 'vector_math.cpp']
boost_dep = dependency('boost')
eigen_dep = dependency('eigen3')
vtk_dep = dependency('vtk', required: false)
executable( 'demo', 
            cid_sources, 
            include_directories:incdir, 
            dependencies : [boost_dep, eigen_dep, vtk_dep],
            override_options : ['cpp_std=c++17'])
# cidlib = shared_library('crater-id',
#                         cid_sources,
#                         include_directories : inc,
#                         # dependencies : glib_dep,
#                         install : true)