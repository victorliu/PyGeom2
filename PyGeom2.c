#include <Python.h>
#include "structmember.h"

#include <GL/gl3w.h>
#include <GLFW/glfw3.h>

#include "nanovg/src/nanovg.h"
#define NANOVG_GL3_IMPLEMENTATION
#include "nanovg/src/nanovg_gl.h"

//#include "imgui_impl_glfw_gl3.cpp"
static int PyBoolean_converter(PyObject *arg, int *result){
	if(PyBool_Check(arg)){
		*result = (Py_True == arg);
		return 1;
	}
	PyErr_Format(PyExc_TypeError, "expected boolean");
	return 0;
}
static int color_converter(PyObject *arg, float *result){
	Py_ssize_t n = 0;
	int i;
	PyObject *fast = NULL;
	static const char *errmsg = "expected color as a length 3 or 4 sequence";
	if(
		(NULL != (fast = PySequence_Fast(arg, errmsg))) &&
		(((n = PySequence_Fast_GET_SIZE(fast)) == 3) || (4 == n))
	){
		float maxelem = 0;
		for(i = 0; i < n; ++i){
			PyObject *elem = PySequence_Fast_GET_ITEM(fast, i);
			float f = PyFloat_AsDouble(elem);
			result[i] = f;
			if(f > maxelem){
				maxelem = f;
			}
			if(-1 == result[i]){
				if(PyErr_Occurred()){
					PyErr_Print();
					return 0;
				}
			}
			if(result[i] < 0){
				PyErr_Format(PyExc_ValueError, "color components must be non-negative");
				return 0;
			}
		}
		maxelem = 1. / maxelem;
		if(4 == n){
			result[0] *= maxelem;
			result[1] *= maxelem;
			result[2] *= maxelem;
			result[3] *= maxelem;
		}else{
			result[0] *= maxelem;
			result[1] *= maxelem;
			result[2] *= maxelem;
			result[3] = 1;
		}
		return 1;
	}
	PyErr_SetString(PyExc_TypeError, errmsg);
	return 0;
}
static int halfwidth_converter(PyObject *arg, double *result){
	Py_ssize_t n = 0;
	int i;
	PyObject *fast = NULL;
	static const char *errmsg = "expected halfwidths as a length 2 sequence";
	if(
		(NULL != (fast = PySequence_Fast(arg, errmsg))) &&
		((n = PySequence_Fast_GET_SIZE(fast)) == 2)
	){
		for(i = 0; i < 2; ++i){
			PyObject *elem = PySequence_Fast_GET_ITEM(fast, i);
			result[i] = PyFloat_AsDouble(elem);
			if(result[i] < 0){
				PyErr_Format(PyExc_ValueError, "halfwidth must be non-negative");
				return 0;
			}
		}
		return 1;
	}
	PyErr_SetString(PyExc_TypeError, errmsg);
	return 0;
}
/*
Interface:

G = PyGeom2.new()

p = G.point(x, y, name = 'asdf', moveable = true)

G.circle(center = p, radius = 2)
G.line(p, q)

G.show(grid = 1)

*/

struct{
	double view_center[2];
	double view_scale;
	double view_angle;
	double view_grid; // if nonzero, then draw grid
	
	int window_size[3]; // 3rd element is minimum of first two
	
	int mouse_button_down[3];
	double mouse_down_pos[2];
	double mouse_down_center[2];
	int moused_point;
	
	double point_size; // positive = pixel space, negative = coordinate space
	double line_width; // same convention as above
	float default_color[4];
} viz_state;

void init_viz_state(){
	viz_state.view_center[0] = 0;
	viz_state.view_center[1] = 0;
	viz_state.view_scale = 1;
	viz_state.view_angle = 0;
	viz_state.view_grid = 0;
	
	viz_state.window_size[0] = 800;
	viz_state.window_size[1] = 600;
	viz_state.window_size[2] = 600;
	
	viz_state.mouse_button_down[0] = 0;
	viz_state.mouse_button_down[1] = 0;
	viz_state.mouse_button_down[2] = 0;
	viz_state.mouse_down_pos[0] = 0;
	viz_state.mouse_down_pos[1] = 0;
	viz_state.mouse_down_center[0] = 0;
	viz_state.mouse_down_center[1] = 0;
	viz_state.moused_point = -1;
	
	viz_state.point_size = 4;
	viz_state.line_width = 2;
	
	viz_state.default_color[0] = 1.00f;
	viz_state.default_color[1] = 0.75f;
	viz_state.default_color[2] = 0.00f;
	viz_state.default_color[3] = 1.00f;
}

void viz_screen_to_coord(const double *screen, double *coord){
	const double mindim = viz_state.window_size[2];
	coord[0] = ((screen[0] - 0.5*viz_state.window_size[0])/mindim - viz_state.view_center[0]) / viz_state.view_scale;
	coord[1] = ((screen[1] - 0.5*viz_state.window_size[1])/mindim - viz_state.view_center[1]) / viz_state.view_scale;
	if(0 != viz_state.view_angle){
		const double cs = cos(viz_state.view_angle/180*M_PI);
		const double sn = sin(viz_state.view_angle/180*M_PI);
		double t = cs*coord[0] + sn*coord[1];
		coord[1] = cs*coord[1] - sn*coord[0];
		coord[0] = t;
	}
}
void viz_coord_to_screen(const double *coord, double *screen){
	// translate center to origin
	// Apply scale and rotation
	double tc[2] = {
		coord[0]*viz_state.view_scale,
		coord[1]*viz_state.view_scale
	};
	if(0 != viz_state.view_angle){
		const double cs = cos(viz_state.view_angle/180*M_PI);
		const double sn = sin(viz_state.view_angle/180*M_PI);
		double t = cs*tc[0] - sn*tc[1];
		tc[1] = cs*tc[1] + sn*tc[0];
		tc[0] = t;
	}
	tc[0] += viz_state.view_center[0];
	tc[1] += viz_state.view_center[1];
	screen[0] = viz_state.window_size[0]*0.5*(tc[0]+1);
	screen[1] = viz_state.window_size[1]*0.5*(tc[1]+1);
}



////////////////////////////////////////////////////////

typedef struct{
	PyObject_HEAD
	double v[2];
} Vector;
static void Vector_dealloc(Vector *self){
	Py_TYPE(self)->tp_free((PyObject *) self);
}
static PyMemberDef Vector_members[] = {
	{"x", T_DOUBLE, offsetof(Vector, v[0]), 0, "x coordinate"},
	{"y", T_DOUBLE, offsetof(Vector, v[1]), 0, "y coordinate"},
	{NULL}  /* Sentinel */
};
static PyMethodDef Vector_methods[] = {
	//{"name", (PyCFunction)Noddy_name, METH_NOARGS, "Return the name, combining the first and last name" },
	{NULL}  /* Sentinel */
};

static PyObject* Vector_add(PyObject *o1, PyObject *o2); // v + v
static PyObject* Vector_sub(PyObject *o1, PyObject *o2); // v - v = p
static PyObject* Vector_neg(PyObject *o1); // negation
static PyObject* Vector_pos(PyObject *o1); // cast to point
static PyObject* Vector_mul(PyObject *o1); // scalar multiply or dot product
static PyObject* Vector_div(PyObject *o1); // scalar divide
static PyObject* Vector_xor(PyObject *o1); // cross product
static PyObject* Vector_invert(PyObject *o1); // rot90

static PyObject* Point_add(PyObject *o1, PyObject *o2); // p + v
static PyObject* Point_sub(PyObject *o1, PyObject *o2); // p - v


static PyNumberMethods Vector_number_methods = {
	.nb_add = Vector_add
};
static PyTypeObject VectorType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "PyGeom2.Vector",
	.tp_doc = "Vector object",
	.tp_basicsize = sizeof(Vector),
	.tp_itemsize = 0,
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
	//.tp_new = Point_new,
	//.tp_init = (initproc)Point_init,
	.tp_dealloc = (destructor)Vector_dealloc,
	.tp_members = Vector_members,
	.tp_methods = Vector_methods,
};

static PyObject* Vector_add(PyObject *o1, PyObject *o2){
	//if(PyT
	return NULL;
}

////////////////////////////////////////////////////////


typedef struct{
	PyObject_HEAD
	PyObject *g;
	int i;
	/*
	double p[2];
	PyObject *name;
	int moveable;*/
} Point;

static void Point_dealloc(Point *self){
	Py_XDECREF(self->g);
	Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyMethodDef Point_methods[] = {
	//{"name", (PyCFunction)Noddy_name, METH_NOARGS, "Return the name, combining the first and last name" },
	{NULL}  /* Sentinel */
};
static PyTypeObject PointType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "PyGeom2.Point",
	.tp_doc = "Point object",
	.tp_basicsize = sizeof(Point),
	.tp_itemsize = 0,
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    //.tp_new = Point_new,
    //.tp_init = (initproc)Point_init,
    .tp_dealloc = (destructor)Point_dealloc,
    .tp_methods = Point_methods,
};

////////////////////////////////////////////////////////

#define GPOINT_MOVEABLE 0x1
typedef struct{
	double p[2];
	PyObject *name;
	unsigned flags;
} gpoint;

typedef struct{
	Point **p;
	Py_ssize_t n;
} vertexlist;
void vertexlist_init(vertexlist *vlist){
	vlist->p = NULL;
	vlist->n = 0;
}
void vertexlist_free(vertexlist *vlist){
	if(NULL != vlist->p){ free(vlist->p); }
}
static int vertexlist_converter(PyObject *arg, vertexlist *result){
	Py_ssize_t n = 0;
	int i;
	PyObject *fast = NULL;
	static const char *errmsg = "expected vertex list of length at least 3";
	if(
		(NULL != (fast = PySequence_Fast(arg, errmsg))) &&
		((n = PySequence_Fast_GET_SIZE(fast)) >= 3)
	){
		result->n = n;
		result->p = (Point**)malloc(sizeof(Point*)*n);
		for(i = 0; i < n; ++i){
			PyObject *elem = PySequence_Fast_GET_ITEM(fast, i);
			if(!PyObject_TypeCheck(elem, &PointType)){
				PyErr_Format(PyExc_TypeError, "element %d is not a Point object", i);
				return 0;
			}
			result->p[i] = (Point*)elem;
		}
		return 1;
	}
	PyErr_SetString(PyExc_TypeError, errmsg);
	return 0;
}

typedef struct{
	PyObject_HEAD
	NVGcontext* vg;
	double t;
	
	int np;
	int np_alloc;
	gpoint *p;
	
	// We record all moveable points in the order in which they are created.
	int nmp;
	int nmp_alloc;
	int *mp;
	
	int nmp_cur;
} GraphicsContext;

static PyObject *GraphicsContext_save(GraphicsContext *ctx, PyObject *args);
static PyObject *GraphicsContext_restore(GraphicsContext *ctx, PyObject *args);
static PyObject *GraphicsContext_translate(GraphicsContext *ctx, PyObject *args);
static PyObject *GraphicsContext_rotate(GraphicsContext *ctx, PyObject *args);
static PyObject *GraphicsContext_scale(GraphicsContext *ctx, PyObject *args);

static PyObject *GraphicsContext_point(GraphicsContext *ctx, PyObject *args, PyObject *kwds);
static PyObject *GraphicsContext_line(GraphicsContext *ctx, PyObject *args, PyObject *kwds);
static PyObject *GraphicsContext_arc(GraphicsContext *ctx, PyObject *args, PyObject *kwds);

static PyObject *GraphicsContext_text(GraphicsContext *ctx, PyObject *args, PyObject *kwds);

static PyObject *GraphicsContext_circle(GraphicsContext *ctx, PyObject *args, PyObject *kwds);
static PyObject *GraphicsContext_ellipse(GraphicsContext *ctx, PyObject *args, PyObject *kwds);
static PyObject *GraphicsContext_rect(GraphicsContext *ctx, PyObject *args, PyObject *kwds);
static PyObject *GraphicsContext_polygon(GraphicsContext *ctx, PyObject *args, PyObject *kwds);

static PyMemberDef GraphicsContext_members[] = {
	{"time", T_DOUBLE, offsetof(GraphicsContext, t), 0, "time"},
	{NULL}  /* Sentinel */
};
static PyMethodDef GraphicsContext_methods[] = {
	{"save"     , (PyCFunction)GraphicsContext_save     , METH_NOARGS, "Save state" },
	{"restore"  , (PyCFunction)GraphicsContext_restore  , METH_NOARGS, "Restore state" },
	{"translate", (PyCFunction)GraphicsContext_translate, METH_VARARGS, "Apply a translation" },
	{"rotate"   , (PyCFunction)GraphicsContext_rotate   , METH_VARARGS, "Apply a rotation" },
	{"scale"    , (PyCFunction)GraphicsContext_scale    , METH_VARARGS, "Apply a scaling" },
	{"point", (PyCFunction)GraphicsContext_point, METH_VARARGS | METH_KEYWORDS, "Create a new point" },
	{"line" , (PyCFunction)GraphicsContext_line , METH_VARARGS | METH_KEYWORDS, "Create a line segment between two points" },
	{"arc"  , (PyCFunction)GraphicsContext_arc  , METH_VARARGS | METH_KEYWORDS, "Create a circular arc" },
	{"text", (PyCFunction)GraphicsContext_text, METH_VARARGS | METH_KEYWORDS, "Create a text string" },
	{"circle" , (PyCFunction)GraphicsContext_circle , METH_VARARGS | METH_KEYWORDS, "Create a circle" },
	{"ellipse", (PyCFunction)GraphicsContext_ellipse, METH_VARARGS | METH_KEYWORDS, "Create an ellipse" },
	{"rect"   , (PyCFunction)GraphicsContext_rect   , METH_VARARGS | METH_KEYWORDS, "Create a rectangle" },
	{"polygon", (PyCFunction)GraphicsContext_polygon, METH_VARARGS | METH_KEYWORDS, "Create a polygon" },
	{NULL}  /* Sentinel */
};
static PyTypeObject GraphicsContextType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "PyGeom2.GraphicsContext",
	.tp_doc = "Graphics context object",
	.tp_basicsize = sizeof(GraphicsContext),
	.tp_itemsize = 0,
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
	//.tp_new = Point_new,
	//.tp_init = (initproc)Point_init,
	//.tp_dealloc = (destructor)Point_dealloc,
	.tp_members = GraphicsContext_members,
	.tp_methods = GraphicsContext_methods,
};

static PyObject *GraphicsContext_save(GraphicsContext *ctx, PyObject *args){
	nvgSave(ctx->vg);
	Py_RETURN_NONE;
}
static PyObject *GraphicsContext_restore(GraphicsContext *ctx, PyObject *args){
	nvgRestore(ctx->vg);
	Py_RETURN_NONE;
}
static PyObject *GraphicsContext_translate(GraphicsContext *ctx, PyObject *args){
	double x, y;
	if(!PyArg_ParseTuple(args, "dd:translate", &x, &y)){ return NULL; }
	nvgTranslate(ctx->vg, x, y);
	Py_RETURN_NONE;
}
static PyObject *GraphicsContext_rotate(GraphicsContext *ctx, PyObject *args){
	double angle;
	if(!PyArg_ParseTuple(args, "d:rotate", &angle)){ return NULL; }
	nvgRotate(ctx->vg, angle*M_PI/180);
	Py_RETURN_NONE;
}
static PyObject *GraphicsContext_scale(GraphicsContext *ctx, PyObject *args){
	double x, y = 0;
	if(!PyArg_ParseTuple(args, "d|d:scale", &x, &y)){ return NULL; }
	if(0 == y){ y = x; }
	nvgScale(ctx->vg, x, y);
	Py_RETURN_NONE;
}
static PyObject *GraphicsContext_point(GraphicsContext *ctx, PyObject *args, PyObject *kwds){
	Point *p = (Point*)PointType.tp_alloc(&PointType, 0);
	if(NULL == p){
		return NULL;
	}
	p->g = (PyObject*)ctx;
	
	PyObject *name = NULL;
	int moveable = 0;
	double x = 0, y = 0;
	NVGcolor color;
	color.rgba[0] = viz_state.default_color[0];
	color.rgba[1] = viz_state.default_color[1];
	color.rgba[2] = viz_state.default_color[2];
	color.rgba[3] = viz_state.default_color[3];
	static char *kwlist[] = {"x", "y", "name", "moveable", "color", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "dd|O!O&O&:point", kwlist,
		&x, &y,
		&PyUnicode_Type, &name,
		&PyBoolean_converter, &moveable,
		&color_converter, &color.rgba[0]
	)){ return NULL; }
	
	if(ctx->np >= ctx->np_alloc){
		ctx->np_alloc *= 2;
		ctx->p = (gpoint*)realloc(ctx->p, sizeof(gpoint)*ctx->np_alloc);
	}
	p->i = ctx->np;
	ctx->p[p->i].name = name;
	ctx->p[p->i].flags = (moveable ? GPOINT_MOVEABLE : 0);
	ctx->np++;
	if(NULL != name){
		Py_INCREF(name);
	}
	if(moveable){
		if(ctx->nmp_cur >= ctx->nmp && ctx->nmp_cur < ctx->np){
			if(ctx->nmp >= ctx->nmp_alloc){
				ctx->nmp_alloc *= 2;
				ctx->mp = (int*)realloc(ctx->p, sizeof(int)*ctx->nmp_alloc);
			}
			ctx->mp[ctx->nmp] = p->i;
			ctx->nmp++;
			
			ctx->p[p->i].p[0] = x;
			ctx->p[p->i].p[1] = y;
		}else{
			x = ctx->p[ctx->mp[ctx->nmp_cur]].p[0];
			y = ctx->p[ctx->mp[ctx->nmp_cur]].p[1];
		}
		ctx->nmp_cur++;
	}else{
		ctx->p[p->i].p[0] = x;
		ctx->p[p->i].p[1] = y;
	}
	//printf("point(%f, %f)\n", p->p[0], p->p[1]);
	nvgBeginPath(ctx->vg);
	double pr = 0;
	if(viz_state.point_size > 0){
		const double mindim = viz_state.window_size[2];
		pr = viz_state.point_size/(viz_state.view_scale*mindim);
	}else if(viz_state.point_size < 0){
		pr = -viz_state.point_size;
	}
	nvgCircle(ctx->vg, x, y, pr);
	nvgFillColor(ctx->vg, color);
	nvgFill(ctx->vg);
	return (PyObject *)p;
}
static PyObject *GraphicsContext_line(GraphicsContext *ctx, PyObject *args, PyObject *kwds){
	Point *pfrom = NULL, *pto = NULL;
	NVGcolor color;
	color.rgba[0] = viz_state.default_color[0];
	color.rgba[1] = viz_state.default_color[1];
	color.rgba[2] = viz_state.default_color[2];
	color.rgba[3] = viz_state.default_color[3];
	static char *kwlist[] = {"from", "to", "stroke", /*"arrow_from", "arrow_to", "dashing",*/ NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|O&:line", kwlist,
		&PointType, (PyObject**)&pfrom,
		&PointType, (PyObject**)&pto,
		&color_converter, &color.rgba[0]
	)){ return NULL; }
	
	nvgBeginPath(ctx->vg);
	nvgMoveTo(ctx->vg, ctx->p[pfrom->i].p[0], ctx->p[pfrom->i].p[1]);
	nvgLineTo(ctx->vg, ctx->p[  pto->i].p[0], ctx->p[  pto->i].p[1]);
	double lw = 0;
	if(viz_state.line_width > 0){
		const double mindim = viz_state.window_size[2];
		lw = viz_state.line_width/(viz_state.view_scale*mindim);
	}else if(viz_state.line_width < 0){
		lw = -viz_state.line_width;
	}
	nvgStrokeWidth(ctx->vg, lw);
	nvgStrokeColor(ctx->vg, color);
	nvgStroke(ctx->vg);
	
	Py_RETURN_NONE;
}
static PyObject *GraphicsContext_arc(GraphicsContext *ctx, PyObject *args, PyObject *kwds){
	Point *pcenter = NULL;
	double radius = 0;
	double start = 0;
	double stop = 360;
	int clockwise = 0;
	NVGcolor color;
	color.rgba[0] = viz_state.default_color[0];
	color.rgba[1] = viz_state.default_color[1];
	color.rgba[2] = viz_state.default_color[2];
	color.rgba[3] = viz_state.default_color[3];
	static char *kwlist[] = {"center", "radius", "start", "stop", "clockwise", "stroke", /*"arrow_from", "arrow_to", "dashing",*/ NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!d|ddO&:arc", kwlist,
		&PointType, (PyObject**)&pcenter,
		&radius, &start, &stop,
		&PyBoolean_converter, &clockwise,
		&color_converter, &color.rgba[0]
	)){ return NULL; }
	
	start *= M_PI/180;
	stop  *= M_PI/180;
	
	nvgBeginPath(ctx->vg);
	nvgArc(ctx->vg, ctx->p[pcenter->i].p[0], ctx->p[pcenter->i].p[1], radius, start, stop, clockwise ? NVG_CW : NVG_CCW);
	double lw = 0;
	if(viz_state.line_width > 0){
		const double mindim = viz_state.window_size[2];
		lw = viz_state.line_width/(viz_state.view_scale*mindim);
	}else if(viz_state.line_width < 0){
		lw = -viz_state.line_width;
	}
	nvgStrokeWidth(ctx->vg, lw);
	nvgStrokeColor(ctx->vg, color);
	nvgStroke(ctx->vg);
	
	Py_RETURN_NONE;
}
static PyObject *GraphicsContext_text(GraphicsContext *ctx, PyObject *args, PyObject *kwds){
	Point *pbase = NULL;
	PyObject *str = NULL;
	double size = 1;
	NVGcolor color;
	color.rgba[0] = viz_state.default_color[0];
	color.rgba[1] = viz_state.default_color[1];
	color.rgba[2] = viz_state.default_color[2];
	color.rgba[3] = viz_state.default_color[3];
	static char *kwlist[] = {"base", "string", "size", "color", /*"arrow_from", "arrow_to", "dashing",*/ NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|dO&:text", kwlist,
		&PointType, (PyObject**)&pbase,
		&PyBytes_Type, (PyObject**)&str,
		&size,
		&color_converter, &color.rgba[0]
	)){ return NULL; }
	
	nvgFontSize(ctx->vg, size);
	nvgFillColor(ctx->vg, color);
	const char *cstr = PyBytes_AsString(str);
	nvgText(ctx->vg, ctx->p[pbase->i].p[0], ctx->p[pbase->i].p[1], cstr, NULL);
	
	Py_RETURN_NONE;
}
static PyObject *GraphicsContext_circle(GraphicsContext *ctx, PyObject *args, PyObject *kwds){
	Point *pcenter = NULL;
	double radius = 0;
	NVGcolor stroke_color;
	stroke_color.rgba[0] = viz_state.default_color[0];
	stroke_color.rgba[1] = viz_state.default_color[1];
	stroke_color.rgba[2] = viz_state.default_color[2];
	stroke_color.rgba[3] = viz_state.default_color[3];
	NVGcolor fill_color = { .rgba={0} };
	static char *kwlist[] = {"center", "radius", "stroke", "fill", /*"dashing",*/ NULL};;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!d|O&O&:circle", kwlist,
		&PointType, (PyObject**)&pcenter,
		&radius,
		&color_converter, &stroke_color.rgba[0],
		&color_converter, &fill_color.rgba[0]
	)){ return NULL; }
	
	nvgBeginPath(ctx->vg);
	nvgCircle(ctx->vg, ctx->p[pcenter->i].p[0], ctx->p[pcenter->i].p[1], radius);
	double lw = 0;
	if(viz_state.line_width > 0){
		const double mindim = viz_state.window_size[2];
		lw = viz_state.line_width/(viz_state.view_scale*mindim);
	}else if(viz_state.line_width < 0){
		lw = -viz_state.line_width;
	}
	nvgStrokeWidth(ctx->vg, lw);
	nvgFillColor(ctx->vg, fill_color);
	nvgFill(ctx->vg);
	nvgStrokeColor(ctx->vg, stroke_color);
	nvgStroke(ctx->vg);
	
	Py_RETURN_NONE;
}
static PyObject *GraphicsContext_ellipse(GraphicsContext *ctx, PyObject *args, PyObject *kwds){
	Point *pcenter = NULL;
	double halfwidth[2] = {0,0};
	NVGcolor stroke_color;
	stroke_color.rgba[0] = viz_state.default_color[0];
	stroke_color.rgba[1] = viz_state.default_color[1];
	stroke_color.rgba[2] = viz_state.default_color[2];
	stroke_color.rgba[3] = viz_state.default_color[3];
	NVGcolor fill_color = { .rgba={0} };
	static char *kwlist[] = {"center", "halfwidths", "stroke", "fill", /*"dashing",*/ NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!O&|O&O&:ellipse", kwlist,
		&PointType, (PyObject**)&pcenter,
		&halfwidth_converter, &halfwidth[0],
		&color_converter, &stroke_color.rgba[0],
		&color_converter, &fill_color.rgba[0]
	)){ return NULL; }
	
	nvgBeginPath(ctx->vg);
	nvgEllipse(ctx->vg, ctx->p[pcenter->i].p[0], ctx->p[pcenter->i].p[1], halfwidth[0], halfwidth[1]);
	double lw = 0;
	if(viz_state.line_width > 0){
		const double mindim = viz_state.window_size[2];
		lw = viz_state.line_width/(viz_state.view_scale*mindim);
	}else if(viz_state.line_width < 0){
		lw = -viz_state.line_width;
	}
	nvgStrokeWidth(ctx->vg, lw);
	nvgFillColor(ctx->vg, fill_color);
	nvgFill(ctx->vg);
	nvgStrokeColor(ctx->vg, stroke_color);
	nvgStroke(ctx->vg);
	
	Py_RETURN_NONE;
}
static PyObject *GraphicsContext_rect(GraphicsContext *ctx, PyObject *args, PyObject *kwds){
	Point *pcenter = NULL;
	double halfwidth[2] = {0,0};
	NVGcolor stroke_color;
	stroke_color.rgba[0] = viz_state.default_color[0];
	stroke_color.rgba[1] = viz_state.default_color[1];
	stroke_color.rgba[2] = viz_state.default_color[2];
	stroke_color.rgba[3] = viz_state.default_color[3];
	NVGcolor fill_color = { .rgba={0} };
	static char *kwlist[] = {"center", "halfwidths", "stroke", "fill", /*"dashing",*/ NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!O&|O&O&:rect", kwlist,
		&PointType, (PyObject**)&pcenter,
		&halfwidth_converter, &halfwidth[0],
		&color_converter, &stroke_color.rgba[0],
		&color_converter, &fill_color.rgba[0]
	)){ return NULL; }
	
	const double *c = &(ctx->p[pcenter->i].p[0]);
	nvgBeginPath(ctx->vg);
	nvgMoveTo(ctx->vg, c[0] - halfwidth[0], c[1] - halfwidth[1]);
	nvgLineTo(ctx->vg, c[0] + halfwidth[0], c[1] - halfwidth[1]);
	nvgLineTo(ctx->vg, c[0] + halfwidth[0], c[1] + halfwidth[1]);
	nvgLineTo(ctx->vg, c[0] - halfwidth[0], c[1] + halfwidth[1]);
	nvgClosePath(ctx->vg);
	double lw = 0;
	if(viz_state.line_width > 0){
		const double mindim = viz_state.window_size[2];
		lw = viz_state.line_width/(viz_state.view_scale*mindim);
	}else if(viz_state.line_width < 0){
		lw = -viz_state.line_width;
	}
	nvgStrokeWidth(ctx->vg, lw);
	nvgFillColor(ctx->vg, fill_color);
	nvgFill(ctx->vg);
	nvgStrokeColor(ctx->vg, stroke_color);
	nvgStroke(ctx->vg);
	
	Py_RETURN_NONE;
}
static PyObject *GraphicsContext_polygon(GraphicsContext *ctx, PyObject *args, PyObject *kwds){
	vertexlist vlist;
	vertexlist_init(&vlist);
	NVGcolor stroke_color;
	stroke_color.rgba[0] = viz_state.default_color[0];
	stroke_color.rgba[1] = viz_state.default_color[1];
	stroke_color.rgba[2] = viz_state.default_color[2];
	stroke_color.rgba[3] = viz_state.default_color[3];
	NVGcolor fill_color = { .rgba={0} };
	static char *kwlist[] = {"vertices", "stroke", "fill", /*"dashing",*/ NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O&|O&O&:rect", kwlist,
		&vertexlist_converter, &vlist,
		&color_converter, &stroke_color.rgba[0],
		&color_converter, &fill_color.rgba[0]
	)){ return NULL; }
	
	nvgBeginPath(ctx->vg);
	Py_ssize_t i;
	for(i = 0; i < vlist.n; ++i){
		const int j = vlist.p[i]->i;
		const double *c = &(ctx->p[j].p[0]);
		if(0 == i){
			nvgMoveTo(ctx->vg, c[0], c[1]);
		}else{
			nvgLineTo(ctx->vg, c[0], c[1]);
		}
	}
	nvgClosePath(ctx->vg);
	double lw = 0;
	if(viz_state.line_width > 0){
		const double mindim = viz_state.window_size[2];
		lw = viz_state.line_width/(viz_state.view_scale*mindim);
	}else if(viz_state.line_width < 0){
		lw = -viz_state.line_width;
	}
	nvgStrokeWidth(ctx->vg, lw);
	nvgFillColor(ctx->vg, fill_color);
	nvgFill(ctx->vg);
	nvgStrokeColor(ctx->vg, stroke_color);
	nvgStroke(ctx->vg);
	vertexlist_free(&vlist);
	Py_RETURN_NONE;
}

////////////////////////////////////////////////////////

/*
static void char_callback(GLFWwindow*, unsigned int c){
	ImGuiIO& io = ImGui::GetIO();
	if(c > 0 && c < 0x10000){
		io.AddInputCharacter((unsigned short)c);
	}
}
*/
static void glfw_framebuffer_size_callback(GLFWwindow* window, int w, int h){
	viz_state.window_size[0] = w;
	viz_state.window_size[1] = h;
	viz_state.window_size[2] = (w < h ? w : h);
	if(viz_state.window_size[2] < 1){ viz_state.window_size[2] = 1; }
	glViewport(0, 0, w, h);
}
static void glfw_mouse_button_callback(GLFWwindow* window, int button, int action, int mods){
	//const ImGuiIO& io = ImGui::GetIO();
	//if(io.WantCaptureMouse){ return; }
	if(!(0 <= button && button < 3)){ return; }
	if(action == GLFW_PRESS){
		viz_state.mouse_button_down[button] = 1;
		glfwGetCursorPos(window, &viz_state.mouse_down_pos[0], &viz_state.mouse_down_pos[1]);
		viz_state.mouse_down_pos[1] = viz_state.window_size[1] - viz_state.mouse_down_pos[1];
		viz_state.mouse_down_center[0] = viz_state.view_center[0];
		viz_state.mouse_down_center[1] = viz_state.view_center[1];
//printf("Clicked at glfw (%f, %f)\n", viz_state.mouse_down_pos[0], viz_state.mouse_down_pos[1]);
double pos[2];
viz_screen_to_coord(viz_state.mouse_down_pos, pos);
printf("  canvas point (%f, %f)\n", pos[0], pos[1]);
	}else if(action == GLFW_RELEASE){
		viz_state.mouse_button_down[button] = 0;
	}
}
static void glfw_cursor_callback(GLFWwindow* window, double x, double y){
	y = viz_state.window_size[1] - y;
	const double xy[2] = {x, y};
	double coord[2];
	viz_screen_to_coord(xy, coord);
	GraphicsContext *ctx = (GraphicsContext*)glfwGetWindowUserPointer(window);
	//const ImGuiIO& io = ImGui::GetIO();
	//if(io.WantCaptureMouse){ return; }
	if(viz_state.mouse_button_down[0]){
		if(viz_state.moused_point < 0){ // if no selected point, drag view
			const double mindim = viz_state.window_size[2];
			const double curscale = viz_state.view_scale / mindim;
			viz_state.view_center[0] = curscale * (x - viz_state.mouse_down_pos[0]) + viz_state.mouse_down_center[0];
			viz_state.view_center[1] = curscale * (y - viz_state.mouse_down_pos[1]) + viz_state.mouse_down_center[1];
		}else if(0 <= viz_state.moused_point && viz_state.moused_point < ctx->np && (ctx->p[viz_state.moused_point].flags & GPOINT_MOVEABLE)){
			ctx->p[viz_state.moused_point].p[0] = coord[0];
			ctx->p[viz_state.moused_point].p[1] = coord[1];
		}
	}else if(viz_state.mouse_button_down[1]){
	}else if(viz_state.mouse_button_down[2]){
	}else{
		int i;
		double dd_threshold;
		if(viz_state.point_size > 0){
			const double mindim = viz_state.window_size[2];
			const double pr = viz_state.point_size/(viz_state.view_scale*mindim);
			dd_threshold = 1.5*pr;
		}else{
			const double pr = viz_state.point_size;
			dd_threshold = 1.5*pr;
		}
		dd_threshold *= dd_threshold;
		
		viz_state.moused_point = -1;
		for(i = 0; i < ctx->nmp; ++i){
			const int j = ctx->mp[i];
			const double d[2] = {
				coord[0] - ctx->p[j].p[0],
				coord[1] - ctx->p[j].p[1]
			};
			const double dd = d[0]*d[0]+d[1]*d[1];
			if(dd < dd_threshold){
				viz_state.moused_point = j;
				break;
			}
		}
	}
}
static void glfw_mouse_scroll_callback(GLFWwindow* window, double x, double y){
	//g_MouseWheel += (float)y; // Use fractional mouse wheel, 1.0 unit 5 lines.
	//const ImGuiIO& io = ImGui::GetIO();
	//if(io.WantCaptureMouse){ return; }
	//arcball.mouseScroll(x, y);
	if(0 != y){
		double view_scale_factor = 1;
		if(y > 0){
			view_scale_factor = 1./0.95;
		}else{
			view_scale_factor = 0.95;
		}
		viz_state.view_scale *= view_scale_factor;
	}
}
static void glfw_error_callback(int error, const char* description){
	PyErr_Format(PyExc_RuntimeError, "GLFW internal error %d: %s", error, description);
}
static PyObject *PyGeom2_show(PyTypeObject *type, PyObject *args, PyObject *kwds){
	GLFWwindow* window;
	NVGcontext* vg = NULL;
	GraphicsContext *ctx = (GraphicsContext*)GraphicsContextType.tp_alloc(&GraphicsContextType, 0);
	if(NULL == ctx){ return NULL; }
	
	double grid = 0;
	PyObject *func;
	{
		static char *kwlist[] = {"function", "grid", NULL};
		if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|d", kwlist,
			&func,
			&grid
		)){ return NULL; }
		if(!PyCallable_Check(func)){
			PyErr_SetString(PyExc_TypeError, "parameter must be callable");
			return NULL;
		}
		Py_INCREF(func);
	}
	
	glfwSetErrorCallback(glfw_error_callback);
	if(!glfwInit()){
		return PyErr_Format(PyExc_RuntimeError, "Failed to initalize GLFW");
	}
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if __APPLE__
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
	init_viz_state();
	viz_state.view_grid = grid;
	
	window = glfwCreateWindow(viz_state.window_size[0], viz_state.window_size[1], "PyGeom2", NULL, NULL);
	if(!window){
		glfwTerminate();
		return PyErr_Format(PyExc_RuntimeError, "Failed to create GLFW window");
	}
	glfwSetWindowUserPointer(window, (void*)ctx);
	
	//glfwSetCharCallback(window, char_callback);
	//glfwSetKeyCallback(window, key_callback);
	glfwSetFramebufferSizeCallback(window, glfw_framebuffer_size_callback);
	glfwSetCursorPosCallback(window, glfw_cursor_callback);
	glfwSetMouseButtonCallback(window, glfw_mouse_button_callback);
	glfwSetScrollCallback(window, glfw_mouse_scroll_callback);
	
	glfwMakeContextCurrent(window);
	
	gl3wInit();
	
	vg = nvgCreateGL3(NVG_ANTIALIAS | NVG_STENCIL_STROKES | NVG_DEBUG);
	
	nvgFontFace(vg, "sans");
	
	//ImGui_ImplGlfwGL3_Init(window, false);
	glfwSwapInterval(1);
	
	Py_INCREF(ctx);
	ctx->vg = vg;
	ctx->t = 0;
	ctx->np = 0;
	ctx->np_alloc = 128;
	ctx->p = (gpoint*)malloc(sizeof(gpoint)*ctx->np_alloc);
	ctx->nmp = 0;
	ctx->nmp_alloc = 128;
	ctx->mp = (int*)malloc(sizeof(int)*ctx->nmp_alloc);
	
	while(!glfwWindowShouldClose(window)){
		ctx->np = 0;
		ctx->nmp_cur = 0;
		//ImGui_ImplGlfwGL3_NewFrame();
		
		glClearColor(0.2f, 0.2f, 0.2f, 1.0f );
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		
		int winWidth, winHeight;
		int fbWidth, fbHeight;
		
		glfwGetWindowSize(window, &winWidth, &winHeight);
		glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
		const float pxRatio = (float)fbWidth / (float)winWidth;
		glViewport(0, 0, fbWidth, fbHeight);
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
		nvgBeginFrame(vg, winWidth, winHeight, pxRatio);
//printf("window size: %d, %d\n", winWidth, winHeight);
		// Set up base transform
		float mindim = (winHeight < winWidth ? winHeight : winWidth);
		
		nvgReset(vg);
		nvgTranslate(vg, 0.5*winWidth, 0.5*winHeight);
		nvgScale(vg, mindim, -mindim);
		nvgTranslate(vg, viz_state.view_center[0], viz_state.view_center[1]);
		nvgScale(vg, viz_state.view_scale, viz_state.view_scale);
		nvgRotate(vg, viz_state.view_angle * M_PI/180);
		if(0){
			float xform[6];
			nvgCurrentTransform(vg, xform);
			printf("cur xform: %f, %f, %f\n           %f, %f, %f\n", xform[0], xform[2], xform[4], xform[1], xform[3], xform[5]);
		}
		PyObject *arglist = Py_BuildValue("(O)", ctx);
		if(NULL == arglist){
			// Handle error
		}else{
			PyObject *result = PyObject_CallObject(func, arglist);
			Py_DECREF(arglist);
			if(NULL == result){
				PyErr_Print();
				break;
			}else{
				Py_DECREF(result);
			}
		}
		// Draw selection for point
		if(viz_state.moused_point >= 0){
			const double x = ctx->p[viz_state.moused_point].p[0];
			const double y = ctx->p[viz_state.moused_point].p[1];
			nvgBeginPath(vg);
			double pr = 0;
			if(viz_state.point_size > 0){
				const double mindim = viz_state.window_size[2];
				pr = viz_state.point_size/(viz_state.view_scale*mindim);
			}else if(viz_state.point_size < 0){
				pr = -viz_state.point_size;
			}
			nvgCircle(vg, x, y, 1.5*pr);
			double lw = 0;
			if(viz_state.line_width > 0){
				const double mindim = viz_state.window_size[2];
				lw = viz_state.line_width/(viz_state.view_scale*mindim);
			}else if(viz_state.line_width < 0){
				lw = -viz_state.line_width;
			}
			nvgStrokeWidth(ctx->vg, lw);
			nvgStrokeColor(vg, nvgRGBA(255,192,0,255));
			nvgStroke(vg);
		}
		// Current matrix is: projection*modelview, but projection is simply identity
		//   modelview:
		//     [ s*cs  -s*sn  0 0 ] [ 1 0 0 -cx ]
		//     [ s*sn   s*cs  0 0 ] [ 0 1 0 -cy ]
		//     [  0      0    1 0 ] [ 0 0 1  0  ]
		//     [  0      0    0 1 ] [ 0 0 0  1  ]
		// ...
		//RenderGui();

		//ImGui::Render();
		nvgEndFrame(vg);
		
		glfwSwapBuffers(window);
		glfwPollEvents();
	}
	
	free(ctx->mp);
	
	Py_DECREF(func);
	Py_DECREF(ctx);
	//ImGui_ImplGlfwGL3_Shutdown();
	glfwDestroyWindow(window);
	
	nvgDeleteGL3(vg);
	glfwTerminate();
	
	Py_RETURN_NONE;
}

static PyMethodDef PyGeom2_methods[] = {
	{"show", (PyCFunction)PyGeom2_show, METH_VARARGS | METH_KEYWORDS, "Show the visualization window" },
	{NULL}  /* Sentinel */
};

#if PY_MAJOR_VERSION >= 3
static PyModuleDef PyGeom2module = {
	PyModuleDef_HEAD_INIT,
	.m_name = "PyGeom2",
	.m_doc = "Interactive 2D geometry visualization tool.",
	.m_size = -1,
	.m_methods = PyGeom2_methods,
};
#endif

#if PY_MAJOR_VERSION >= 3
#define INIT_RETVAL m
#else
#define INIT_RETVAL
#endif

PyMODINIT_FUNC
#if PY_MAJOR_VERSION >= 3
PyInit_PyGeom2(void)
#else
initPyGeom2(void)
#endif
{
	PyObject *m = NULL;
	if(PyType_Ready(&PointType) < 0){ return INIT_RETVAL; }
	if(PyType_Ready(&VectorType) < 0){ return INIT_RETVAL; }
	if(PyType_Ready(&GraphicsContextType) < 0){ return INIT_RETVAL; }
#if PY_MAJOR_VERSION >= 3
	m = PyModule_Create(&PyGeom2module);
#else
	m = Py_InitModule3("PyGeom2", PyGeom2_methods, "Interactive 2D geometry visualization tool.");
#endif

	if(m == NULL){ return INIT_RETVAL; }
	Py_INCREF(&PointType);
	PyModule_AddObject(m, "Point", (PyObject*)&PointType);
	PyModule_AddObject(m, "Vector", (PyObject*)&VectorType);
	PyModule_AddObject(m, "GraphicsContext", (PyObject*)&GraphicsContextType);
	return INIT_RETVAL;
}
