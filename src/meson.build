cid_sources = ['main.cpp', 'conics.cpp']
boost_dep = dependency('boost')
executable( 'demo', 
            cid_sources, 
            include_directories:incdir, 
            dependencies : boost_dep,
            override_options : ['cpp_std=c++17'])
# cidlib = shared_library('crater-id',
#                         cid_sources,
#                         include_directories : inc,
#                         # dependencies : glib_dep,
#                         install : true)