/*
 * CHANGELOG: New functions for animations:
 *
 * clearplot()
 * plotdisplay()
 * plotbackground()
 * ------eg.------
 *
 * plotinit();
 * plotbackground();
 * FRAME 1 { plot, plot, plot, ... }
 * plotdisplay();
 *
 * clearplot();
 * plotbackground();
 * FRAME 1 { plot, plot, plot, ... }
 * plotdisplay();
 * 
 * ...
 *
 */
/* 
 * plot - plot on root window
 */

/* 
 * SYNOPSIS
 *
 * #include "plot.h"
 * 
 * void 
 * plotinit (double xmin, double xmax, double ymin, double ymax, 
 * 		const char *fgcolorstr, const char *bgcolorstr, 
 * 		const char *xlabel, const char *ylabel, 
 * 		const char *title,...);
 *
 * void 
 * plot (double x, double y, double xmin, double xmax, double ymin, 
 * 		double ymax, int state);
 *
 * void
 * plotexit(void);
 *
 */

#include <X11/Xlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <math.h>

Display *display;	
Pixmap pix;
GC gc;
int screen;
Window win;		
unsigned int display_width, display_height;
Colormap screen_colormap;

int xX (double x, double xmin, double xmax, int width) {
	double w = (double)width - 44.0;
	if (isinf(x)) return width/2;
	if (isnan(x)) return 0;
	if (x<xmin) return 22;
	if (x>xmax) return (width)/2;
	return (int)(22.0 + (w*x/(xmax-xmin) - w*xmin/(xmax-xmin))/2);
}
int yY (double y, double ymin, double ymax, int height) {
	double h = (double)height - 44.0;
	if (isinf(y)) return 0;
	if (isnan(y)) return height;
	if (y<ymin) return height-22;
	if (y>ymax) return 22;
	return (int)(22.0 + (y-ymax)*h/(ymin-ymax));
}

int p_xX (double x, double xmin, double xmax, int width) {
	double w = (double)width - 44.0;
	if (isinf(x)) return width-22;
	if (isnan(x)) return 0;
	if (x<xmin) return width/2;
	if (x>xmax) return width-22;
	return (int)(width + (w*x/(xmax-xmin) - w*xmin/(xmax-xmin)))/2;
}
int p_yY (double y, double ymin, double ymax, int height) {
	double h = ((double)height - 44.0)/2;
	if (isinf(y)) return 0;
	if (isnan(y)) return height-22;
	if (y<ymin) return height-22;
	if (y>ymax) return 22;
	return (int)(22 + (y-ymax)*h/(ymin-ymax));
}

int pu_yY (double y, double ymin, double ymax, int height) {
	if (isinf(y)) return height/2;
	if (isnan(y)) return height-22;
	if (y<ymin) return height-22;
	if (y>ymax) return height/2;
	return (int)(height/2 + (y-ymax)*(height/2-22)/(ymin-ymax));
}
/* mode 0: continue 1: newline 
 * style 0: points, 1: line: 2: dashed line, 3: long-dashed */
void plot (double x, double y, double xmin, double xmax, double ymin, 
		double ymax, int mode, int linestyle, int linewidth, 
		const char *colorstr) {
	static int init=1;
	int line_style = LineSolid;	
	int cap_style = CapButt;
	int join_style = JoinBevel; 
	int X,Y;
	static int oldx, oldy;
	XColor color;
	Status rc;
	// LineOnOffDash, LineDoubleDash
	
	line_style = LineOnOffDash;
	switch (linestyle) {
	case 1:
		line_style = LineSolid;
		break;
	case 2:
		line_style = LineOnOffDash;
		break;
	case 3:
		line_style = LineDoubleDash;
		break;
	default:
		line_style = LineSolid; /* 0: Points */
		break;
	}

	rc = XAllocNamedColor (display, screen_colormap, colorstr, &color, &color);
	if (rc == 0) {
		perror("colour fault.");
		exit (1);
	}
	XSetForeground(display, gc, color.pixel);
	
	X = xX(x,xmin,xmax,display_width);
	Y = yY(y,ymin,ymax,display_height);
	if (mode==1) init = 1;
	if (linestyle==0) {
		if (linewidth > 0) {
			//XSetFillStyle(display, gc, FillSolid);
			XFillArc(display, pix, gc, X-linewidth/2, Y-linewidth/2, linewidth, linewidth, 0, 360*64);
		} else {
			XDrawPoint(display, pix, gc, X, Y);
		}
		return;
	}
	if (init) {
		oldx = X;
		oldy = Y;
		init = 0;
	} else {
		XSetLineAttributes(display, gc, linewidth, line_style, 
			cap_style, join_style); 
		XDrawLine(display, pix, gc, oldx, oldy, X, Y);
		oldx = X; 
		oldy = Y;
	}	
	
	/* draw two intersecting lines, one horizontal and one vertical,
	 * which intersect at point "50,100".                           
	XDrawLine(display, pix, gc, 50, 0, 50, 200);
	XDrawLine(display, pix, gc, 0, 100, 200, 100);
	XSetWindowBackgroundPixmap (display, win, pix);
	XClearWindow(display,win);
	XFlush(display);
	*/
}

void plot_time(double x, double xmin, double xmax) {
	int line_style = LineSolid;
	int cap_style = CapButt;
	int join_style = JoinBevel; 
	int linewidth = 2;
	XColor color;
	Status rc;
	int X;

	rc = XAllocNamedColor (display, screen_colormap, "red", &color, &color);
	if (rc == 0) {
		perror("colour fault.");
		exit (1);
	}
	XSetForeground(display, gc, color.pixel);
	XSetLineAttributes(display, gc, linewidth, line_style, cap_style, join_style); 

	X = p_xX(x,xmin,xmax,display_width);
	XDrawLine(display, pix, gc, X, 22, X, display_height-22);
}
void plotinflexion(double x, double xmin, double xmax) {
	int line_style = LineSolid;	
	int cap_style = CapButt;
	int join_style = JoinBevel; 
	int linewidth = 4;
	int X;
	XColor color;
	Status rc;

	rc = XAllocNamedColor (display, screen_colormap, "red", &color, &color);
	if (rc == 0) {
		perror("colour fault.");
		exit (1);
	}
	XSetForeground(display, gc, color.pixel);
	XSetLineAttributes(display, gc, linewidth, line_style, cap_style, join_style); 
	
	X = xX(x,xmin,xmax,display_width);
	XDrawLine(display, pix, gc, X, 22, X, display_height-22);

	//X = p_xX(x,xmin,xmax,display_width);
	//XDrawLine(display, pix, gc, X, -10+display_height/2, X, display_height/2);
	
}

void plot_inflexionline(int *vec, int N) {
	int i;
	int line_style = LineSolid;	
	int cap_style = CapButt;
	int join_style = JoinBevel; 
	int linewidth = 4;
	int X;
	int line=0;
	int linestart;
	int Xold=-1;
	XColor color;
	Status rc;

	rc = XAllocNamedColor (display, screen_colormap, "red", &color, &color);
	if (rc == 0) {
		perror("colour fault.");
		exit (1);
	}
	XSetForeground(display, gc, color.pixel);
	XSetLineAttributes(display, gc, linewidth, line_style, cap_style, join_style); 

	for (i=0; i<N; i++) {
		if ((!(line)) && (vec[i])) {
			linestart = i;
			Xold = p_xX((double)i,0,(double)N,display_width);
			line = 1;
			linewidth=1;
			XSetLineAttributes(display, gc, linewidth, line_style, cap_style, join_style); 
			if (Xold != display_width/2)
				XDrawLine(display, pix, gc, Xold, 22, Xold, display_height-22);
		} 
		if ((line) && (!(vec[i]))) {
			line = 0;
			X = p_xX((double)i,0,(double)N,display_width);
			linewidth=4;
			XSetLineAttributes(display, gc, linewidth, line_style, cap_style, join_style); 
			XDrawLine(display, pix, gc, Xold, display_height/2, X, display_height/2);
			linewidth=1;
			XSetLineAttributes(display, gc, linewidth, line_style, cap_style, join_style); 
			//if (Xold != display_width/2)
				XDrawLine(display, pix, gc, X, 22, X, display_height-22);
		} 
		if ((i==N-1)&&(line)) {
			X = p_xX((double)i,0,(double)N,display_width);
			linewidth=4;
			XSetLineAttributes(display, gc, linewidth, line_style, cap_style, join_style); 
			XDrawLine(display, pix, gc, Xold, display_height/2, X, display_height/2);
		}
	}

}
void plotpressure(double x, double y, double xmin, double xmax, double ymin, double ymax, int first) {
	static int init = 1;
	int line_style = LineSolid;	
	int cap_style = CapButt;
	int join_style = JoinBevel; 
	static int oldx, oldy;
	int linewidth = 2;
	char pressure_range[100]; 
	int X,Y;

	XColor color;
	Status rc;

	snprintf(pressure_range, 99, "[%f]", ymax);
	
	rc = XAllocNamedColor (display, screen_colormap, "green", &color, &color);
	if (rc == 0) {
		perror("colour fault.");
		exit (1);
	}
	XSetForeground(display, gc, color.pixel);
	XSetLineAttributes(display, gc, linewidth, line_style, cap_style, join_style); 
	
	X = p_xX(x,xmin,xmax,display_width);
	Y = p_yY(y,ymin,ymax,display_height);
	if (first) init = 1;
	if (init) {
		XDrawString(display, pix, gc, display_width/2+180, 18, pressure_range, strlen(pressure_range));
		oldx = X;
		oldy = Y;
		init = 0;
	} else {
		XDrawLine(display, pix, gc, oldx, oldy, X, Y);
		oldx = X; 
		oldy = Y;
	}
}

void plotflow(double x, double y, double xmin, double xmax, double ymin, double ymax, int first) {
	static int init = 1;
	int line_style = LineSolid;	
	int cap_style = CapButt;
	int join_style = JoinBevel; 
	static int oldx, oldy;
	int linewidth = 2;
	char pressure_range[100]; 
	int X,Y;

	XColor color;
	Status rc;

	snprintf(pressure_range, 99, "[%f]", ymax);
	

	rc = XAllocNamedColor (display, screen_colormap, "green", &color, &color);
	if (rc == 0) {
		perror("colour fault.");
		exit (1);
	}
	XSetForeground(display, gc, color.pixel);
	XSetLineAttributes(display, gc, linewidth, line_style, cap_style, join_style); 
	
	X = p_xX(x,xmin,xmax,display_width);
	Y = pu_yY(y,ymin,ymax,display_height);
	if (first) init = 1;
	if (init) {
		XDrawString(display, pix, gc, display_width/2+50, 18+display_height/2, pressure_range, strlen(pressure_range));
		oldx = X;
		oldy = Y;
		init = 0;
	} else {
		XDrawLine(display, pix, gc, oldx, oldy, X, Y);
		oldx = X; 
		oldy = Y;
	}
}

void plotfflow(double x, double y, double xmin, double xmax, double ymin, double ymax, int first) {
	static int init = 1;
	int line_style = LineSolid;	
	int cap_style = CapButt;
	int join_style = JoinBevel; 
	static int oldx, oldy;
	int linewidth = 2;
	char pressure_range[100]; 
	int X,Y;

	XColor color;
	Status rc;

	snprintf(pressure_range, 99, "[%f]", ymax);
	

	rc = XAllocNamedColor (display, screen_colormap, "yellow", &color, &color);
	if (rc == 0) {
		perror("colour fault.");
		exit (1);
	}
	XSetForeground(display, gc, color.pixel);
	XSetLineAttributes(display, gc, linewidth, line_style, cap_style, join_style); 
	
	X = p_xX(x,xmin,xmax,display_width);
	Y = pu_yY(y,ymin,ymax,display_height);
	if (first) init = 1;
	if (init) {
		XDrawString(display, pix, gc, display_width/2+50, 18+display_height/2, pressure_range, strlen(pressure_range));
		oldx = X;
		oldy = Y;
		init = 0;
	} else {
		XDrawLine(display, pix, gc, oldx, oldy, X, Y);
		oldx = X; 
		oldy = Y;
	}
}

void clearplot(void) {
	XSetWindowBackgroundPixmap (display, win, pix);
	XClearWindow(display,win);
//	XFlush(display);
	
	XSetForeground(display, gc, BlackPixel(display, screen));
 	XFillRectangle(display, pix, gc, 0, 0, display_width, display_height);
	XSetWindowBackgroundPixmap (display, win, pix);
//	XClearWindow(display,win);
//	XFlush(display);
}
void plotdisplay(void) {
	XSetWindowBackgroundPixmap (display, win, pix);
	XClearWindow(display,win);
	XFlush(display);

}
void plotbackground (double xmin, double xmax, double ymin, double ymax, const char *fgcolorstr, 
		const char *bgcolorstr, const char *xlabel, const char *ylabel, 
		const char *title, int width, int height) {

	char *font_name ="10x20";
	const char *veloclabel = "velocity";
	const char *pressurelabel = "pressure gradient";
	const char *flowlabel = "flow";
	XFontStruct *font_info;
	unsigned int line_width = 1;
	int line_style = LineSolid;	
	int cap_style = CapButt;
	int join_style = JoinBevel; 
	char yrange[100]; 
	XColor fgcolor, bgcolor;
	Status rc;
	snprintf(yrange, 99, "[%f]", ymax);

	XSetLineAttributes(display, gc, line_width, line_style, cap_style, join_style); 
	/* Color */
	screen_colormap = DefaultColormap(display, screen);
	rc = XAllocNamedColor (display, screen_colormap, fgcolorstr, &fgcolor, &fgcolor);
	if (rc == 0) {
		perror("colour fault.");
		exit (1);
	}	
	rc = XAllocNamedColor (display, screen_colormap, bgcolorstr, &bgcolor, &bgcolor);
	if (rc == 0) {
		perror("colour fault.");
		exit (1);
	}
	XSetForeground(display, gc, fgcolor.pixel);
	XSetBackground(display, gc, bgcolor.pixel);

	/* Clear Background */
	XSetForeground(display, gc, bgcolor.pixel);
	XFillRectangle(display, pix, gc, 0, 0, display_width, display_height);
	XSetForeground(display, gc, fgcolor.pixel);

	/* Font */
	font_info = XLoadQueryFont(display, font_name);
	XSetFont(display, gc, font_info->fid);

	XDrawString(display, pix, gc, display_width-10*strlen(xlabel)-22, display_height, xlabel, strlen(xlabel));
//	XDrawString(display, pix, gc, 22, 18, ylabel, strlen(ylabel));
	XDrawString(display, pix, gc, 24, 18, veloclabel, strlen(veloclabel));
	XDrawString(display, pix, gc, 34+10*strlen(veloclabel), 18, yrange, strlen(yrange));
	XDrawString(display, pix, gc, display_width/2+2, 18, pressurelabel, strlen(pressurelabel));
	XDrawString(display, pix, gc, display_width/2+2, display_height/2+18, flowlabel, strlen(flowlabel));

	
	XDrawRectangle(display, pix, gc, 22, 22, display_width-44, display_height-44);
	line_style = LineSolid;
	XSetLineAttributes(display, gc, 1, line_style, cap_style, join_style); 
	// vert Mittellinie
	XDrawLine(display, pix, gc, display_width/2, 
			22, display_width/2, display_height-22);
	// rechts horizontal Trennlinie
	XDrawLine(display, pix, gc, display_width/2, 
			display_height/2, display_width-22, display_height/2);
	line_style = LineOnOffDash;
	XSetLineAttributes(display, gc, 1, line_style, cap_style, join_style); 
	// rechts oben Nulllinie
	XDrawLine(display, pix, gc, display_width/2, display_height/4+11, display_width-22, display_height/4+11);
	// rechts unten Nulllinie
	XDrawLine(display, pix, gc, display_width/2, display_height*3/4-11, display_width-22, display_height*3/4-11);

	if ((xmin<0.0) && (xmax>0.0)) {
		line_style = LineOnOffDash;
		XSetLineAttributes(display, gc, line_width, line_style, cap_style, join_style); 
		XDrawLine(display, pix, gc, xX(0.0,xmin,xmax,display_width), 
				22, xX(0.0,xmin,xmax,display_width), display_height-22);
	}
	if ((ymin<0.0) && (ymax>0.0)) {
		line_style = LineOnOffDash;
		XSetLineAttributes(display, gc, line_width, line_style, cap_style, join_style); 
		XDrawLine(display, pix, gc, 22, yY(0.0,ymin,ymax,display_height), 
				display_width-22, yY(0.0,ymin,ymax,display_height));
	}
}

void plotinit(double xmin, double xmax, double ymin, double ymax, const char *fgcolorstr, 
		const char *bgcolorstr, const char *xlabel, const char *ylabel, 
		const char *title, int width, int height) {


	char *display_name = getenv("DISPLAY");
//	char *font_name ="10x20";
//	XFontStruct *font_info;
	unsigned long valuemask = 0;
	XGCValues values;
//	unsigned int line_width = 1;
//	int line_style = LineSolid;	
//	int cap_style = CapButt;
//	int join_style = JoinBevel; 
//	char xrange[100]; 
//	char yrange[100]; 
//	snprintf(xrange, 99, "[%f, %f]", xmin, xmax);
//	snprintf(yrange, 99, "[%f, %f]", ymin, ymax);

	display = XOpenDisplay(display_name);
	if (display == NULL) {
		fprintf(stderr, "cannot connect to X server '%s'\n",
			display_name);
		exit(1);
	}
	/* get the geometry of the default screen for our display. */
	screen = DefaultScreen(display);
	display_width = DisplayWidth(display, screen);
	display_height = DisplayHeight(display, screen);

	if ((width==0)||(height==0)) win = RootWindow(display, screen);
	else { 
		win = XCreateSimpleWindow(display,RootWindow(display,screen),
			0,0,width,height,0,BlackPixel(display,screen),WhitePixel(display,screen));
		XMapWindow(display,win);
	}
	gc = XCreateGC(display, win, valuemask, &values);
	if (gc < 0) {
		fprintf(stderr, "XCreateGC: \n");
	}
	//XSetForeground(display, gc, WhitePixel(display, screen));
	//XSetBackground(display, gc, BlackPixel(display, screen)); 
	//XSetFillStyle(display, gc, FillSolid);

	pix = XCreatePixmap(display, win, display_width, display_height, (unsigned int)DefaultDepth(display, screen));

	XSync(display, False);


}
void plotexit(void) {
	XSetWindowBackgroundPixmap (display, win, pix);
	XClearWindow(display,win);
	XFlush(display);
	XCloseDisplay(display);
}

