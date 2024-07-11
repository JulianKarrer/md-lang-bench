lib_sources = [  # All source files (excluding headers)
    'lj_direct_summation.cpp',
    'atoms.cpp',
    'verlet.cpp'
]

lib = static_library(
    'lib',  # Unique name
    lib_sources,
    dependencies : [eigen]
)

lib_incdirs = include_directories('.')
