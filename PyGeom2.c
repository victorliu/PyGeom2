#include <Python.h>
#include "structmember.h"

#if PY_MAJOR_VERSION >= 3
#define Py_TPFLAGS_CHECKTYPES 0
#endif

#include <GL/gl3w.h>
#include <GLFW/glfw3.h>
#include <float.h>

#include "nanovg/src/nanovg.h"
#define NANOVG_GL3_IMPLEMENTATION
#include "nanovg/src/nanovg_gl.h"

#include "fonts/dejavu.h"

static double tween(double x0){
	static const double initial_slope = 0.1;
	static const double a = 2*(1-initial_slope);
	double x = x0;
	if(x0 > 0.5){ x = 1-x; }
	double y = a*(x+initial_slope)*x;
	if(x0 > 0.5){ y = 1-y; }
	return y;
}

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

typedef struct{
	double center[2];
	double scale;
	double angle; // CCW, in degrees
} view_state;

typedef struct{
	view_state view;
	view_state view_targ;
	view_state view_prev;
	
	int tween_num_timesteps;
	int tween_timesteps_remaining;
	
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

void init_viz_state(viz_state *v){
	v->view.center[0] = 0;
	v->view.center[1] = 0;
	v->view.scale = 1;
	v->view.angle = 0;
	
	v->view_targ.center[0] = 0;
	v->view_targ.center[1] = 0;
	v->view_targ.scale = 1;
	v->view_targ.angle = 0;
	
	v->view_prev.center[0] = 0;
	v->view_prev.center[1] = 0;
	v->view_prev.scale = 1;
	v->view_prev.angle = 0;
	
	v->tween_num_timesteps = 20;
	v->tween_timesteps_remaining = 0;
	
	v->view_grid = 0;
	
	v->window_size[0] = 800;
	v->window_size[1] = 600;
	v->window_size[2] = 600;
	
	v->mouse_button_down[0] = 0;
	v->mouse_button_down[1] = 0;
	v->mouse_button_down[2] = 0;
	v->mouse_down_pos[0] = 0;
	v->mouse_down_pos[1] = 0;
	v->mouse_down_center[0] = 0;
	v->mouse_down_center[1] = 0;
	v->moused_point = -1;
	
	v->point_size = -4;
	v->line_width = -2;
	
	v->default_color[0] = 1.00f;
	v->default_color[1] = 0.75f;
	v->default_color[2] = 0.00f;
	v->default_color[3] = 1.00f;
}

double viz_state_get_coord_size(const viz_state *v, double s){
	if(s < 0){
		const double r = v->view.scale/v->window_size[2];
		return -s*r;
	}
	return s;
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

static PyObject* Vector_norm(Vector *v, PyObject *args){
	return PyFloat_FromDouble(hypot(v->v[0], v->v[1]));
}
static PyObject* Vector_angle(Vector *v, PyObject *args){
	return PyFloat_FromDouble(180.*(atan2(v->v[1], v->v[0])/M_PI));
}
static PyObject* Vector_normalize(Vector *v, PyObject *args);
static PyMethodDef Vector_methods[] = {
	{"angle", (PyCFunction)Vector_angle, METH_NOARGS, "Return the angle of the vector from the x-axis" },
	{"norm", (PyCFunction)Vector_norm, METH_NOARGS, "Return the Euclidean length of the vector" },
	{"normalize", (PyCFunction)Vector_normalize, METH_NOARGS, "Normalizes and returns the vector" },
	{NULL}  /* Sentinel */
};

static PyObject* Vector_add(PyObject *o1, PyObject *o2); // v + v
static PyObject* Vector_sub(PyObject *o1, PyObject *o2); // v - v = p
static PyObject* Vector_neg(PyObject *o1); // negation
static PyObject* Vector_mul(PyObject *o1, PyObject *o2); // scalar multiply or dot product
static PyObject* Vector_div(PyObject *o1, PyObject *o2); // scalar divide
static PyObject* Vector_xor(PyObject *o1, PyObject *o2); // cross product
static PyObject* Vector_invert(PyObject *o1); // rot90

static PyObject* Point_add(PyObject *o1, PyObject *o2); // p + v
static PyObject* Point_sub(PyObject *o1, PyObject *o2); // p - v, p - p


static PyNumberMethods Vector_number_methods = {
	.nb_add = Vector_add,
	.nb_subtract = Vector_sub,
	.nb_multiply = Vector_mul,
	.nb_divide = Vector_div,
	.nb_true_divide = Vector_div,
	.nb_negative = Vector_neg,
	.nb_xor = Vector_xor,
	.nb_invert = Vector_invert,
};
PyObject* Vector_repr(Vector *v){
	char str[64];
	snprintf(str, 64, "[%f, %f]", v->v[0], v->v[1]);
	return PyUnicode_FromString(str);
}
static PyTypeObject VectorType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "PyGeom2.Vector",
	.tp_doc = "Vector object",
	.tp_basicsize = sizeof(Vector),
	.tp_itemsize = 0,
	.tp_as_number = &Vector_number_methods,
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES,
	//.tp_new = Point_new,
	//.tp_init = (initproc)Point_init,
	.tp_repr = (reprfunc)Vector_repr,
	.tp_dealloc = (destructor)Vector_dealloc,
	.tp_members = Vector_members,
	.tp_methods = Vector_methods,
};
static int Vector_check(PyObject *obj){
	return PyObject_TypeCheck(obj, &VectorType);
}

static PyObject* Vector_add(PyObject *o1, PyObject *o2){ // v + v
	if(!Vector_check(o1)){
		return PyErr_Format(PyExc_TypeError, "expected Point on LHS of plus operator");
	}
	if(!Vector_check(o2)){
		return PyErr_Format(PyExc_TypeError, "expected Point on RHS of plus operator");
	}
	Vector *lhs = (Vector*)o1;
	Vector *rhs = (Vector*)o2;
	Vector *v = (Vector*)VectorType.tp_alloc(&VectorType, 0);
	if(NULL == v){
		PyErr_Format(PyExc_MemoryError, "could not create a new Vector");
	}
	v->v[0] = lhs->v[0] + rhs->v[0];
	v->v[1] = lhs->v[1] + rhs->v[1];
	return (PyObject*)v;
}

static PyObject* Vector_sub(PyObject *o1, PyObject *o2){ // v - v
	if(!Vector_check(o1)){
		return PyErr_Format(PyExc_TypeError, "expected Vector on LHS of minus operator");
	}
	if(!Vector_check(o2)){
		return PyErr_Format(PyExc_TypeError, "expected Vector on RHS of minus operator");
	}
	Vector *lhs = (Vector*)o1;
	Vector *rhs = (Vector*)o2;
	Vector *v = (Vector*)VectorType.tp_alloc(&VectorType, 0);
	if(NULL == v){
		PyErr_Format(PyExc_MemoryError, "could not create a new Vector");
	}
	v->v[0] = lhs->v[0] - rhs->v[0];
	v->v[1] = lhs->v[1] - rhs->v[1];
	return (PyObject*)v;
}
static PyObject* Vector_neg(PyObject *o1){
	if(!Vector_check(o1)){
		return PyErr_Format(PyExc_TypeError, "expected Vector");
	}
	Vector *u = (Vector*)o1;
	Vector *v = (Vector*)VectorType.tp_alloc(&VectorType, 0);
	if(NULL == v){
		PyErr_Format(PyExc_MemoryError, "could not create a new Vector");
	}
	v->v[0] = -u->v[0];
	v->v[1] = -u->v[1];
	return (PyObject*)v;
}
static PyObject* Vector_mul(PyObject *o1, PyObject *o2){
	if(!Vector_check(o2)){
		return PyErr_Format(PyExc_TypeError, "expected Vector on RHS of multiplication operator");
	}
	Vector *rhs = (Vector*)o2;
	if(Vector_check(o1)){ // scalar dot product
		Vector *lhs = (Vector*)o1;
		return PyFloat_FromDouble(lhs->v[0]*rhs->v[0] + lhs->v[1]*rhs->v[1]);
	}else if(PyFloat_Check(o1)){
		double s = PyFloat_AsDouble(o1);
		Vector *v = (Vector*)VectorType.tp_alloc(&VectorType, 0);
		if(NULL == v){
			PyErr_Format(PyExc_MemoryError, "could not create a new Vector");
		}
		v->v[0] = s*rhs->v[0];
		v->v[1] = s*rhs->v[1];
		return (PyObject*)v;
	}else{
		return PyErr_Format(PyExc_TypeError, "expected number or Vector on LHS of multiplication operator");
	}
}
static PyObject* Vector_div(PyObject *o1, PyObject *o2){
	if(!Vector_check(o1)){
		return PyErr_Format(PyExc_TypeError, "expected Vector on LHS of division operator");
	}
	Vector *lhs = (Vector*)o1;
	if(PyFloat_Check(o2)){
		double s = PyFloat_AsDouble(o2);
		Vector *v = (Vector*)VectorType.tp_alloc(&VectorType, 0);
		if(NULL == v){
			PyErr_Format(PyExc_MemoryError, "could not create a new Vector");
		}
		v->v[0] = lhs->v[0]/s;
		v->v[1] = lhs->v[1]/s;
		return (PyObject*)v;
	}else{
		return PyErr_Format(PyExc_TypeError, "expected number on RHS of division operator");
	}
}
static PyObject* Vector_xor(PyObject *o1, PyObject *o2){
	if(!Vector_check(o1)){
		return PyErr_Format(PyExc_TypeError, "expected Vector on LHS of ^ operator");
	}
	Vector *lhs = (Vector*)o1;
	if(!Vector_check(o2)){
		return PyErr_Format(PyExc_TypeError, "expected Vector on RHS of ^ operator");
	}
	Vector *rhs = (Vector*)o2;
	return PyFloat_FromDouble(lhs->v[0]*rhs->v[1] - lhs->v[1]*rhs->v[0]);
}
static PyObject* Vector_invert(PyObject *o1){
	if(!Vector_check(o1)){
		return PyErr_Format(PyExc_TypeError, "expected Vector");
	}
	Vector *u = (Vector*)o1;
	Vector *v = (Vector*)VectorType.tp_alloc(&VectorType, 0);
	if(NULL == v){
		PyErr_Format(PyExc_MemoryError, "could not create a new Vector");
	}
	v->v[0] = -u->v[1];
	v->v[1] =  u->v[0];
	return (PyObject*)v;
}
static PyObject* Vector_normalize(Vector *v, PyObject *args){
	double ilen = 1./hypot(v->v[0], v->v[1]);
	v->v[0] *= ilen;
	v->v[1] *= ilen;
	Vector *r = (Vector*)VectorType.tp_alloc(&VectorType, 0);
	if(NULL == r){
		PyErr_Format(PyExc_MemoryError, "could not create a new Vector");
	}
	r->v[0] = v->v[0];
	r->v[1] = v->v[1];
	return (PyObject*)r;
}
////////////////////////////////////////////////////////


typedef struct{
	PyObject_HEAD
	PyObject *g;
	int i;
} Point;

static void Point_dealloc(Point *self){
	Py_XDECREF(self->g);
	Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyMethodDef Point_methods[] = {
	//{"name", (PyCFunction)Noddy_name, METH_NOARGS, "Return the name, combining the first and last name" },
	{NULL}  /* Sentinel */
};

static PyNumberMethods Point_number_methods = {
	.nb_add = Point_add,
	.nb_subtract = Point_sub,
};
static PyObject* Point_repr(Point *p);

static PyTypeObject PointType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	.tp_name = "PyGeom2.Point",
	.tp_doc = "Point object",
	.tp_basicsize = sizeof(Point),
	.tp_itemsize = 0,
	.tp_repr = (reprfunc)Point_repr,
	.tp_as_number = &Point_number_methods,
	.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES,
	//.tp_new = Point_new,
	//.tp_init = (initproc)Point_init,
	.tp_dealloc = (destructor)Point_dealloc,
	.tp_methods = Point_methods,
};
static int Point_check(PyObject *obj){
	return PyObject_TypeCheck(obj, &PointType);
}


////////////////////////////////////////////////////////

#define GPOINT_MOVEABLE 0x1
#define GPOINT_VISIBLE  0x2
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
			if(!Point_check(elem)){
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
	viz_state viz;
	double t;
	
	int np;
	int np_alloc;
	gpoint *p;
	
	// We record all moveable points in the order in which they are created.
	int nmp;
	int nmp_alloc;
	int *mp;
	
	int nmp_cur;
	
	double bound[4];
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

static Point* GraphicsContext_point_add(GraphicsContext *ctx, double x, double y){
	Point *p = (Point*)PointType.tp_alloc(&PointType, 0);
	if(NULL == p){
		return NULL;
	}
	p->g = (PyObject*)ctx;
	if(ctx->np >= ctx->np_alloc){
		ctx->np_alloc *= 2;
		ctx->p = (gpoint*)realloc(ctx->p, sizeof(gpoint)*ctx->np_alloc);
	}
	p->i = ctx->np;
	ctx->p[p->i].p[0] = x;
	ctx->p[p->i].p[1] = y;
	ctx->p[p->i].name = NULL;
	ctx->p[p->i].flags = 0;
	ctx->np++;
	return p;
}

static PyObject* Point_repr(Point *p){
	char str[64];
	double x, y;
	const int i = p->i;
	if(NULL == p->g){ return NULL; }
	GraphicsContext *ctx = (GraphicsContext*)p->g;
	x = ctx->p[i].p[0];
	y = ctx->p[i].p[1];
	snprintf(str, 64, "(%f, %f)", x, y);
	return PyUnicode_FromString(str);
}
static PyObject* Point_add(PyObject *o1, PyObject *o2){ // p + v
	if(!Point_check(o1)){
		return PyErr_Format(PyExc_TypeError, "expected Point on LHS of plus operator");
	}
	Point *lhs = (Point*)o1;
	if(Vector_check(o2)){
		Point *p;
		Vector *rhs = (Vector*)o2;
		GraphicsContext *lg = (GraphicsContext*)lhs->g;
		p = GraphicsContext_point_add(lg,
			lg->p[lhs->i].p[0] + rhs->v[0],
			lg->p[lhs->i].p[1] + rhs->v[1]
		);
		if(NULL == p){
			return PyErr_Format(PyExc_MemoryError, "could not create a new Point");
		}
		return (PyObject*)p;
	}else{
		return PyErr_Format(PyExc_TypeError, "expected a Vector on RHS of plus operator");
	}
}
static PyObject* Point_sub(PyObject *o1, PyObject *o2){ // p - v, p - p
	if(!Point_check(o1)){
		return PyErr_Format(PyExc_TypeError, "expected Point on LHS of minus operator");
	}
	Point *lhs = (Point*)o1;
	if(Point_check(o2)){
		Point *rhs = (Point*)o2;
		Vector *v = (Vector*)VectorType.tp_alloc(&VectorType, 0);
		if(NULL == v){
			PyErr_Format(PyExc_MemoryError, "could not create a new Vector");
		}
		GraphicsContext *lg = (GraphicsContext*)lhs->g;
		GraphicsContext *rg = (GraphicsContext*)rhs->g;
		v->v[0] = lg->p[lhs->i].p[0] - rg->p[rhs->i].p[0];
		v->v[1] = lg->p[lhs->i].p[1] - rg->p[rhs->i].p[1];
		return (PyObject*)v;
	}else if(Vector_check(o2)){
		Point *p;
		Vector *rhs = (Vector*)o2;
		GraphicsContext *lg = (GraphicsContext*)lhs->g;
		p = GraphicsContext_point_add(lg,
			lg->p[lhs->i].p[0] - rhs->v[0],
			lg->p[lhs->i].p[1] - rhs->v[1]
		);
		if(NULL == p){
			PyErr_Format(PyExc_MemoryError, "could not create a new Point");
		}
		return (PyObject*)p;
	}else{
		return PyErr_Format(PyExc_TypeError, "expected Point or Vector on RHS of minus operator");
	}
}

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
	PyObject *name = NULL;
	Point *pp = NULL;
	int moveable = 0;
	int visible = 0;
	double x = 0, y = 0;
	NVGcolor color;
	color.rgba[0] = ctx->viz.default_color[0];
	color.rgba[1] = ctx->viz.default_color[1];
	color.rgba[2] = ctx->viz.default_color[2];
	color.rgba[3] = ctx->viz.default_color[3];
	static char *kwlist[] = {"x", "y", "name", "moveable", "visible", "color", NULL};
	static char *kwlist2[] = {"point", "name", "moveable", "visible", "color", NULL};
	if(PyArg_ParseTupleAndKeywords(args, kwds, "O!|O!O&O&O&:point", kwlist2,
		&PointType, &pp,
		&PyUnicode_Type, &name,
		&PyBoolean_converter, &moveable,
		&PyBoolean_converter, &visible,
		&color_converter, &color.rgba[0]
	)){
		x = ctx->p[pp->i].p[0];
		y = ctx->p[pp->i].p[1];
	}else{
		PyErr_Clear();
		if(!PyArg_ParseTupleAndKeywords(args, kwds, "dd|O!O&O&O&:point", kwlist,
			&x, &y,
			&PyUnicode_Type, &name,
			&PyBoolean_converter, &moveable,
			&PyBoolean_converter, &visible,
			&color_converter, &color.rgba[0]
		)){ return NULL; }
	}
	Point *p;
	if(moveable){
		visible = 1;
		if(ctx->nmp_cur >= ctx->nmp){
			if(ctx->nmp >= ctx->nmp_alloc){
				ctx->nmp_alloc *= 2;
				ctx->mp = (int*)realloc(ctx->mp, sizeof(int)*ctx->nmp_alloc);
			}
			p = GraphicsContext_point_add(ctx, x, y);
			ctx->mp[ctx->nmp] = p->i;
			ctx->nmp++;
			
			ctx->p[p->i].p[0] = x;
			ctx->p[p->i].p[1] = y;
		}else{
			x = ctx->p[ctx->mp[ctx->nmp_cur]].p[0];
			y = ctx->p[ctx->mp[ctx->nmp_cur]].p[1];
			p = GraphicsContext_point_add(ctx, x, y);
		}
		ctx->p[p->i].flags |= GPOINT_MOVEABLE;
		ctx->nmp_cur++;
	}else{
		p = GraphicsContext_point_add(ctx, x, y);
	}
	
	if(visible){
		ctx->p[p->i].flags |= GPOINT_VISIBLE;
	}
	
	ctx->p[p->i].name = name;
	if(NULL != name){
		Py_INCREF(name);
	}

	if(x < ctx->bound[0]){ ctx->bound[0] = x; }
	if(y < ctx->bound[1]){ ctx->bound[1] = y; }
	if(x > ctx->bound[2]){ ctx->bound[2] = x; }
	if(y > ctx->bound[3]){ ctx->bound[3] = y; }
	
	//printf("point(%f, %f)\n", p->p[0], p->p[1]);
	if(visible){
		nvgBeginPath(ctx->vg);
		const double pr = viz_state_get_coord_size(&ctx->viz, ctx->viz.point_size);
		nvgCircle(ctx->vg, x, y, pr);
		nvgFillColor(ctx->vg, color);
		nvgFill(ctx->vg);
	}
	return (PyObject *)p;
}
static PyObject *GraphicsContext_line(GraphicsContext *ctx, PyObject *args, PyObject *kwds){
	Point *pfrom = NULL, *pto = NULL;
	NVGcolor color;
	color.rgba[0] = ctx->viz.default_color[0];
	color.rgba[1] = ctx->viz.default_color[1];
	color.rgba[2] = ctx->viz.default_color[2];
	color.rgba[3] = ctx->viz.default_color[3];
	static char *kwlist[] = {"from", "to", "stroke", /*"arrow_from", "arrow_to", "dashing",*/ NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|O&:line", kwlist,
		&PointType, (PyObject**)&pfrom,
		&PointType, (PyObject**)&pto,
		&color_converter, &color.rgba[0]
	)){ return NULL; }
	
	nvgBeginPath(ctx->vg);
	nvgMoveTo(ctx->vg, ctx->p[pfrom->i].p[0], ctx->p[pfrom->i].p[1]);
	nvgLineTo(ctx->vg, ctx->p[  pto->i].p[0], ctx->p[  pto->i].p[1]);
	const double lw = viz_state_get_coord_size(&ctx->viz, ctx->viz.line_width);
	nvgStrokeWidth(ctx->vg, lw);
	nvgStrokeColor(ctx->vg, color);
	nvgStroke(ctx->vg);
	
	Py_RETURN_NONE;
}
static PyObject *GraphicsContext_arc(GraphicsContext *ctx, PyObject *args, PyObject *kwds){
	Point *pc = NULL;
	double radius = 0;
	double start = 0;
	double stop = 360;
	Point *pd = NULL;
	int clockwise = 0;
	NVGcolor color;
	color.rgba[0] = ctx->viz.default_color[0];
	color.rgba[1] = ctx->viz.default_color[1];
	color.rgba[2] = ctx->viz.default_color[2];
	color.rgba[3] = ctx->viz.default_color[3];
	static char *kwlist[] = {"center", "radius", "start", "stop", "clockwise", "stroke", /*"arrow_from", "arrow_to", "dashing",*/ NULL};
	static char *kwlist2[] = {"from", "to", "bulge", "stroke", /*"arrow_from", "arrow_to", "dashing",*/ NULL};
	if(PyArg_ParseTupleAndKeywords(args, kwds, "O!d|ddO&:arc", kwlist,
		&PointType, (PyObject**)&pc,
		&radius, &start, &stop,
		&PyBoolean_converter, &clockwise,
		&color_converter, &color.rgba[0]
	)){
		start *= M_PI/180;
		stop  *= M_PI/180;
		
		nvgBeginPath(ctx->vg);
		nvgArc(ctx->vg, ctx->p[pc->i].p[0], ctx->p[pc->i].p[1], radius, start, stop, clockwise ? NVG_CW : NVG_CCW);
		const double lw = viz_state_get_coord_size(&ctx->viz, ctx->viz.line_width);
		nvgStrokeWidth(ctx->vg, lw);
		nvgStrokeColor(ctx->vg, color);
		nvgStroke(ctx->vg);
		Py_RETURN_NONE;
	}else{
		PyErr_Clear();
		if(PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|dO&:arc", kwlist2,
			&PointType, (PyObject**)&pc,
			&PointType, (PyObject**)&pd,
			&radius,
			&color_converter, &color.rgba[0]
		)){
			nvgBeginPath(ctx->vg);
			nvgMoveTo(ctx->vg, ctx->p[pc->i].p[0], ctx->p[pc->i].p[1]);
			nvgArcBulgeTo(ctx->vg, ctx->p[pd->i].p[0], ctx->p[pd->i].p[1], radius);
			const double lw = viz_state_get_coord_size(&ctx->viz, ctx->viz.line_width);
			nvgStrokeWidth(ctx->vg, lw);
			nvgStrokeColor(ctx->vg, color);
			nvgStroke(ctx->vg);
			Py_RETURN_NONE;
		}
		return NULL;
	}
}
static PyObject *GraphicsContext_text(GraphicsContext *ctx, PyObject *args, PyObject *kwds){
	Point *pbase = NULL;
	PyObject *str = NULL;
	double size = 1;
	NVGcolor color;
	color.rgba[0] = ctx->viz.default_color[0];
	color.rgba[1] = ctx->viz.default_color[1];
	color.rgba[2] = ctx->viz.default_color[2];
	color.rgba[3] = ctx->viz.default_color[3];
	static char *kwlist[] = {"base", "string", "size", "color", /*"arrow_from", "arrow_to", "dashing",*/ NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!|dO&:text", kwlist,
		&PointType, (PyObject**)&pbase,
		&PyBytes_Type, (PyObject**)&str,
		&size,
		&color_converter, &color.rgba[0]
	)){ return NULL; }
	size = viz_state_get_coord_size(&ctx->viz, size);
	
	nvgFontSize(ctx->vg, size);
	nvgFontFace(ctx->vg, "mono");
	nvgFontBlur(ctx->vg, 0);
	nvgTextAlign(ctx->vg, NVG_ALIGN_LEFT|NVG_ALIGN_BASELINE);
	nvgFillColor(ctx->vg, color);
	const char *cstr = PyBytes_AsString(str);
	float bounds[4];
	nvgTextBounds(ctx->vg, ctx->p[pbase->i].p[0], ctx->p[pbase->i].p[1], cstr, NULL, bounds);
	//printf("%f, %f (%s) %f %f %f %f\n", ctx->p[pbase->i].p[0], ctx->p[pbase->i].p[1], cstr, bounds[0], bounds[1], bounds[2], bounds[3]);
	nvgText(ctx->vg, ctx->p[pbase->i].p[0], ctx->p[pbase->i].p[1], cstr, NULL);
	
	Py_RETURN_NONE;
}
static PyObject *GraphicsContext_circle(GraphicsContext *ctx, PyObject *args, PyObject *kwds){
	Point *pcenter = NULL;
	double r = 0;
	NVGcolor stroke_color;
	stroke_color.rgba[0] = ctx->viz.default_color[0];
	stroke_color.rgba[1] = ctx->viz.default_color[1];
	stroke_color.rgba[2] = ctx->viz.default_color[2];
	stroke_color.rgba[3] = ctx->viz.default_color[3];
	NVGcolor fill_color = { .rgba={0} };
	static char *kwlist[] = {"center", "radius", "stroke", "fill", /*"dashing",*/ NULL};;
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O!d|O&O&:circle", kwlist,
		&PointType, (PyObject**)&pcenter,
		&r,
		&color_converter, &stroke_color.rgba[0],
		&color_converter, &fill_color.rgba[0]
	)){ return NULL; }
	
	double x = ctx->p[pcenter->i].p[0];
	double y = ctx->p[pcenter->i].p[1];
	if(x-r < ctx->bound[0]){ ctx->bound[0] = x-r; }
	if(y-r < ctx->bound[1]){ ctx->bound[1] = y-r; }
	if(x+r > ctx->bound[2]){ ctx->bound[2] = x+r; }
	if(y+r > ctx->bound[3]){ ctx->bound[3] = y+r; }
	
	nvgBeginPath(ctx->vg);
	nvgCircle(ctx->vg, x, y, r);
	const double lw = viz_state_get_coord_size(&ctx->viz, ctx->viz.line_width);
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
	stroke_color.rgba[0] = ctx->viz.default_color[0];
	stroke_color.rgba[1] = ctx->viz.default_color[1];
	stroke_color.rgba[2] = ctx->viz.default_color[2];
	stroke_color.rgba[3] = ctx->viz.default_color[3];
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
	const double lw = viz_state_get_coord_size(&ctx->viz, ctx->viz.line_width);
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
	stroke_color.rgba[0] = ctx->viz.default_color[0];
	stroke_color.rgba[1] = ctx->viz.default_color[1];
	stroke_color.rgba[2] = ctx->viz.default_color[2];
	stroke_color.rgba[3] = ctx->viz.default_color[3];
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
	const double lw = viz_state_get_coord_size(&ctx->viz, ctx->viz.line_width);
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
	stroke_color.rgba[0] = ctx->viz.default_color[0];
	stroke_color.rgba[1] = ctx->viz.default_color[1];
	stroke_color.rgba[2] = ctx->viz.default_color[2];
	stroke_color.rgba[3] = ctx->viz.default_color[3];
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
	const double lw = viz_state_get_coord_size(&ctx->viz, ctx->viz.line_width);
	nvgStrokeWidth(ctx->vg, lw);
	nvgFillColor(ctx->vg, fill_color);
	nvgFill(ctx->vg);
	nvgStrokeColor(ctx->vg, stroke_color);
	nvgStroke(ctx->vg);
	vertexlist_free(&vlist);
	Py_RETURN_NONE;
}

////////////////////////////////////////////////////////

static void glfw_key_callback(GLFWwindow* window, int key, int scancode, int action, int mods){
	if(GLFW_PRESS == action){
		switch(key){
		case GLFW_KEY_ESCAPE:
			glfwSetWindowShouldClose(window, GL_TRUE); break;
		default: break;
		}
	}
}
/*
static void char_callback(GLFWwindow*, unsigned int c){
	ImGuiIO& io = ImGui::GetIO();
	if(c > 0 && c < 0x10000){
		io.AddInputCharacter((unsigned short)c);
	}
}
*/
static void glfw_framebuffer_size_callback(GLFWwindow* window, int w, int h){
	GraphicsContext *ctx = (GraphicsContext*)glfwGetWindowUserPointer(window);
	ctx->viz.window_size[0] = w;
	ctx->viz.window_size[1] = h;
	ctx->viz.window_size[2] = (w < h ? w : h);
	if(ctx->viz.window_size[2] < 1){ ctx->viz.window_size[2] = 1; }
	glViewport(0, 0, w, h);
}
static void glfw_mouse_button_callback(GLFWwindow* window, int button, int action, int mods){
	GraphicsContext *ctx = (GraphicsContext*)glfwGetWindowUserPointer(window);
	//const ImGuiIO& io = ImGui::GetIO();
	//if(io.WantCaptureMouse){ return; }
	if(!(0 <= button && button < 3)){ return; }
	if(action == GLFW_PRESS){
		ctx->viz.mouse_button_down[button] = 1;
		glfwGetCursorPos(window, &ctx->viz.mouse_down_pos[0], &ctx->viz.mouse_down_pos[1]);
		ctx->viz.mouse_down_pos[1] = ctx->viz.window_size[1] - ctx->viz.mouse_down_pos[1];
		ctx->viz.mouse_down_center[0] = ctx->viz.view.center[0];
		ctx->viz.mouse_down_center[1] = ctx->viz.view.center[1];
		float pos[2];
		float xform[6], xformi[6], r[2];
		const float src[2] = { ctx->viz.mouse_down_pos[0], ctx->viz.mouse_down_pos[1] };
		nvgCurrentTransform(ctx->vg, xform);
		nvgTransformInverse(xformi, xform);
		nvgTransformPoint(&pos[0], &pos[1], xformi, src[0], src[1]);
		printf("Click pos: (%f, %f)\n", pos[0], pos[1]);
		fflush(stdout);
	}else if(action == GLFW_RELEASE){
		ctx->viz.mouse_button_down[button] = 0;
	}
}
static void glfw_cursor_callback(GLFWwindow* window, double x, double y){
	GraphicsContext *ctx = (GraphicsContext*)glfwGetWindowUserPointer(window);
	y = ctx->viz.window_size[1] - y;
	float coord[2];
	float xform[6], xformi[6], r[2];
	nvgCurrentTransform(ctx->vg, xform);
	nvgTransformInverse(xformi, xform);
	nvgTransformPoint(&coord[0], &coord[1], xformi, x, y);
	
	//const ImGuiIO& io = ImGui::GetIO();
	//if(io.WantCaptureMouse){ return; }
	if(ctx->viz.mouse_button_down[0]){
		if(ctx->viz.moused_point < 0){ // if no selected point, drag view
			const double mindim = ctx->viz.window_size[2];
			const double curscale = ctx->viz.view.scale / mindim;
			ctx->viz.view.center[0] = curscale * (ctx->viz.mouse_down_pos[0] - x) + ctx->viz.mouse_down_center[0];
			ctx->viz.view.center[1] = curscale * (ctx->viz.mouse_down_pos[1] - y) + ctx->viz.mouse_down_center[1];
		}else{
			if(0 <= ctx->viz.moused_point && ctx->viz.moused_point < ctx->np && (ctx->p[ctx->viz.moused_point].flags & GPOINT_MOVEABLE)){
				ctx->p[ctx->viz.moused_point].p[0] = coord[0];
				ctx->p[ctx->viz.moused_point].p[1] = coord[1];
			}
		}
	}else if(ctx->viz.mouse_button_down[1]){
	}else if(ctx->viz.mouse_button_down[2]){
	}else{
		int i;
		double dd_threshold = viz_state_get_coord_size(&ctx->viz, ctx->viz.point_size);
		dd_threshold *= 1.5;
		dd_threshold *= dd_threshold;
		
		ctx->viz.moused_point = -1;
		for(i = 0; i < ctx->nmp; ++i){
			const int j = ctx->mp[i];
			const double d[2] = {
				coord[0] - ctx->p[j].p[0],
				coord[1] - ctx->p[j].p[1]
			};
			const double dd = d[0]*d[0]+d[1]*d[1];
			if(dd < dd_threshold){
				ctx->viz.moused_point = j;
				break;
			}
		}
	}
}
static void glfw_mouse_scroll_callback(GLFWwindow* window, double x, double y){
	GraphicsContext *ctx = (GraphicsContext*)glfwGetWindowUserPointer(window);
	//g_MouseWheel += (float)y; // Use fractional mouse wheel, 1.0 unit 5 lines.
	//const ImGuiIO& io = ImGui::GetIO();
	//if(io.WantCaptureMouse){ return; }
	//arcball.mouseScroll(x, y);
	if(0 != y){
		double view_scale_factor = 1;
		if(y > 0){
			view_scale_factor = 0.95;
		}else{
			view_scale_factor = 1./0.95;
		}
		ctx->viz.view.scale *= view_scale_factor;
		{ // apply translation to keep point under the cursor fixed
			glfwGetCursorPos(window, &x, &y);
			y = ctx->viz.window_size[1] - y;
			
			float xform[6], xformi[6], r[2];
			const float src[2] = { x, y };
			nvgCurrentTransform(ctx->vg, xform);
			nvgTransformInverse(xformi, xform);
			nvgTransformPoint(&r[0], &r[1], xformi, src[0], src[1]);
			ctx->viz.view.center[0] = r[0] - (r[0] - ctx->viz.view.center[0]) * view_scale_factor;
			ctx->viz.view.center[1] = r[1] - (r[1] - ctx->viz.view.center[1]) * view_scale_factor;
		}
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
#ifndef _WIN32 // don't require this on win32, and works with more cards
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#endif
	glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, 1);
	
	init_viz_state(&ctx->viz);
	ctx->viz.view_grid = grid;
	
	window = glfwCreateWindow(ctx->viz.window_size[0], ctx->viz.window_size[1], "PyGeom2", NULL, NULL);
	if(!window){
		glfwTerminate();
		return PyErr_Format(PyExc_RuntimeError, "Failed to create GLFW window");
	}
	glfwSetWindowUserPointer(window, (void*)ctx);
	
	//glfwSetCharCallback(window, char_callback);
	glfwSetKeyCallback(window, glfw_key_callback);
	glfwSetFramebufferSizeCallback(window, glfw_framebuffer_size_callback);
	glfwSetCursorPosCallback(window, glfw_cursor_callback);
	glfwSetMouseButtonCallback(window, glfw_mouse_button_callback);
	glfwSetScrollCallback(window, glfw_mouse_scroll_callback);
	
	glfwMakeContextCurrent(window);
	
	gl3wInit();
	
	//vg = nvgCreateGL3(NVG_ANTIALIAS | NVG_STENCIL_STROKES | NVG_DEBUG);
	vg = nvgCreateGL3(NVG_ANTIALIAS);
	
	nvgCreateFontMem(vg, "mono", &DejaVuSansMono[0], sizeof(DejaVuSansMono)-1, 0);
	
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
	
	int frame0 = 1;
	
	glfwSwapInterval(1);
	while(!glfwWindowShouldClose(window)){
		ctx->np = 0;
		ctx->nmp_cur = 0;
		ctx->bound[0] = DBL_MAX;
		ctx->bound[1] = DBL_MAX;
		ctx->bound[2] = -DBL_MAX;
		ctx->bound[3] = -DBL_MAX;
		//ImGui_ImplGlfwGL3_NewFrame();
		
		
		int winWidth, winHeight;
		int fbWidth, fbHeight;
		
		glfwGetWindowSize(window, &winWidth, &winHeight);
		glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
		const float pxRatio = (float)fbWidth / (float)winWidth;
		glViewport(0, 0, fbWidth, fbHeight);
	
		glEnable(GL_STENCIL_TEST);
		glClearColor(0.2f, 0.2f, 0.2f, 1.0f );
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
		nvgBeginFrame(vg, winWidth, winHeight, pxRatio);
//printf("window size: %d, %d\n", winWidth, winHeight);
		// Set up base transform
		
		if(ctx->viz.tween_timesteps_remaining >= 0){
			const double t = 1 - (double)ctx->viz.tween_timesteps_remaining / (double)ctx->viz.tween_num_timesteps;
			const double f = tween(t);
			ctx->viz.view.center[0] = (1-t)*ctx->viz.view_prev.center[0] + t*ctx->viz.view_targ.center[0];
			ctx->viz.view.center[1] = (1-t)*ctx->viz.view_prev.center[1] + t*ctx->viz.view_targ.center[1];
			ctx->viz.view.scale = (1-t)*ctx->viz.view_prev.scale + t*ctx->viz.view_targ.scale;
			ctx->viz.view.angle = (1-t)*ctx->viz.view_prev.angle + t*ctx->viz.view_targ.angle;
			ctx->viz.tween_timesteps_remaining--;
		}else{
			ctx->viz.view_prev = ctx->viz.view;
		}
		
		const double mindim = ctx->viz.window_size[2];
		
		nvgReset(vg);
		// Apply current transform
		/*
		nvgTranslate(vg, 0.5*winWidth, 0.5*winHeight);
		nvgScale(vg, mindim/ctx->viz.view.scale, mindim/ctx->viz.view.scale);
		nvgRotate(vg, -ctx->viz.view.angle * M_PI/180);
		nvgTranslate(vg, -ctx->viz.view.center[0], -ctx->viz.view.center[1]);
		*/
		if(0 == ctx->viz.view.angle){
			const float s = (float)ctx->viz.window_size[2]/(float)ctx->viz.view.scale;
			nvgTransform(vg,
				s, 0, 0, s,
				0.5*(float)winWidth  - s*(float)ctx->viz.view.center[0],
				0.5*(float)winHeight - s*(float)ctx->viz.view.center[1]
			);
		}else{
			const float s = (float)ctx->viz.window_size[2]/(float)ctx->viz.view.scale;
			const float angle = -nvgDegToRad(ctx->viz.view.angle);
			const float cs = cosf(angle);
			const float sn = sinf(angle);
			nvgTransform(vg,
				s*cs, s*sn, -s*sn, s*cs,
				0.5*(float)winWidth  - s*(cs*(float)ctx->viz.view.center[0]-sn*(float)ctx->viz.view.center[1]),
				0.5*(float)winHeight - s*(cs*(float)ctx->viz.view.center[1]+sn*(float)ctx->viz.view.center[0])
			);
		}
		/*
		nvgStrokeColor(vg, nvgRGBA(255,192,0,255));
		nvgBeginPath(vg);
		nvgCircle(vg, 10, 10, 2);
		nvgStroke(vg);
		nvgBeginPath(vg);
		nvgCircle(vg, 100, 10, 2);
		nvgStroke(vg);
		nvgBeginPath(vg);
		nvgCircle(vg, 10, 100, 2);
		nvgStroke(vg);
		*/
		if(0){
			float xform[6];
			nvgCurrentTransform(vg, xform);
			printf("cur xform: %f, %f, %f\n           %f, %f, %f\n", xform[0], xform[2], xform[4], xform[1], xform[3], xform[5]);
		}
		if(1){
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
		}
		if(frame0 && ctx->bound[0] < ctx->bound[2] && ctx->bound[1] < ctx->bound[3]){ // set initial view
			double dx = ctx->bound[2] - ctx->bound[0];
			double dy = ctx->bound[3] - ctx->bound[1];
			ctx->viz.view_targ.scale = 1.2*(dx > dy ? dx : dy);
			ctx->viz.view_targ.center[0] = 0.5*(ctx->bound[0] + ctx->bound[2]);
			ctx->viz.view_targ.center[1] = 0.5*(ctx->bound[1] + ctx->bound[3]);
			ctx->viz.view_targ.angle = 0;
			ctx->viz.tween_timesteps_remaining = ctx->viz.tween_num_timesteps;
		}
		// Draw selection for point
		if(ctx->viz.moused_point >= 0){
			const double x = ctx->p[ctx->viz.moused_point].p[0];
			const double y = ctx->p[ctx->viz.moused_point].p[1];
			nvgBeginPath(vg);
			const double pr = viz_state_get_coord_size(&ctx->viz, ctx->viz.point_size);
			nvgCircle(vg, x, y, 1.5*pr);
			const double lw = viz_state_get_coord_size(&ctx->viz, ctx->viz.line_width);
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
		frame0 = 0;
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
