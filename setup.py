from distutils.core import setup, Extension
from sys import platform
import os

libraries = [];
if platform == 'darwin':
	libraries.append('glfw');
	os.environ['LDFLAGS'] = '-framework Cocoa -framework OpenGL -framework IOKit -framework CoreFoundation -framework CoreVideo';
elif platform == 'win32':
	libraries.append('glfw3');
	libraries.append('opengl32');
elif platform == 'linux':
	libraries.append('glfw3');
	libraries.append('GL');

PyGeom2 = Extension('PyGeom2',
	include_dirs = ['gl3w/include','/usr/local/include'],
	libraries = libraries,
	library_dirs = ['/usr/local/lib'],
	sources = ['PyGeom2.c', 'gl3w/src/gl3w.c', 'nanovg/src/nanovg.c', 'fonts/dejavu_sans_mono.c']
)

setup (name = 'PyGeom2',
	version = '1.0',
	description = 'Interactive 2D geometry visualization tool',
ext_modules = [PyGeom2])