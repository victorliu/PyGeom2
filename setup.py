from distutils.core import setup, Extension
from sys import platform
import os

if platform == 'darwin':
	os.environ['LDFLAGS'] = '-framework Cocoa -framework OpenGL -framework IOKit -framework CoreFoundation -framework CoreVideo';

PyGeom2 = Extension('PyGeom2',
	include_dirs = ['gl3w/include','/usr/local/include'],
	libraries = ['glfw3'],
	library_dirs = ['/usr/local/lib'],
	sources = ['PyGeom2.c', 'gl3w/src/gl3w.c', 'nanovg/src/nanovg.c']
)

setup (name = 'PyGeom2',
	version = '1.0',
	description = 'Interactive 2D geometry visualization tool',
ext_modules = [PyGeom2])