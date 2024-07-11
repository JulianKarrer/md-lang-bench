# https://mesonbuild.com/
project(
    'moleculardynamics', 
    ['c', 'cpp'], # Project type. We need a C++ compiler. The C compiler is required for MPI.
    default_options : ['cpp_std=c++17'], # Yes, we need C++17, at least for constexpr
    version : '0.1'
)

eigen = dependency('eigen3', version: '>=3.4.0')

subdir('src')

executable(
	'bench',
	'main.cpp',
	include_directories : [lib_incdirs],
	link_with : [lib],
	dependencies : [eigen]
)