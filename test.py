import PyGeom2

def viz(g):
	a=g.point(0,0,color=(1,0,0))
	b=g.point(0.1,0,moveable=True)
	c=g.point(0,0.1,moveable=True)
	g.text(base=c, string = "asdf", size=-18, color = (0,1,0))
	g.line(a,b)
	g.rect(c, (0.05,0.06))
	g.polygon(vertices = (a, b, c), fill = (1,2,3))
	#g.point(50,50)
	#g.point(100,50)
	#g.point(50,100)

PyGeom2.show(viz)
