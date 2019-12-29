from distutils.core import setup, Extension

setup(ext_modules=[Extension(
    'priors',
    sources=['priors.c', '_funcs.c'],
    include_dirs=['.']
)])
