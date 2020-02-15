/* $Id: pp.c,v 1.10 2004/05/17 19:32:50 elmar Exp $ */

/*
 * pp - post processor
 *
 * USAGE
 * 	pp WIDTH HEIGHT DATAFILE [ DATAFILE2 ... ]
 *
 * DATAFILE FORMAT
 * 	binary array of doubles... (machine dependent)
 */


/*
o binary format self describtable
	header:
	width height origscale endian
	comments
	binary data
o multipics mousewheel
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <math.h>
#include <gtk/gtk.h>

#include "colourmap.h"
#include "readfile.h"

#define EVENT_METHOD(i, x) GTK_WIDGET_GET_CLASS(i)->x
int XSIZE=600, YSIZE=600;
int WXSIZE=600, WYSIZE=600;
int MAXSCALE=10;
int SCALE=1;
double ZMIN=0.0, ZMAX=1.0;
// double XMIN=0.0, XMAX=1.0, YMIN=0.0, YMAX=1.0;
int TFX,TFY,TTX,TTY; 	/* translate from/to x,y */
int OFFX=0, OFFY=0;	/* x,y offset */	
int FILES=1;		/* number of datasets */
int SET=0;
struct datafile **DATA;

GtkWidget 	*zminspinner, *zmaxspinner,*ezminspinner, *ezmaxspinner, *scalespinner;
GtkWidget	*colourmapcombo;
GtkWidget      *xframe, *yframe, *zframe; 
GtkWidget	*xlabel, *ylabel, *zlabel, *flabel, *zminlabel, *zmaxlabel, *setlabel;
GtkWidget	*area, *cbararea;
GtkWidget	*hrule, *vrule, *crule;
guchar *cbuf;	/* scaled char buffer */
guchar *vbuf;	/* visual buffer */
guchar cbarbuf[50*256];	/* colourmap legend buffer */
//guchar *blackbuf;	/* black background */
GdkRgbCmap *cmap=NULL;

static gboolean area_expose(GtkWidget *widget, GdkEvent *event, gpointer data) ;
double scale (double a, double b, double A, double B, double x) {
	return (x*(A-B)/(a-b)+(a*B-b*A)/(a-b));
}

void rescale(void) {
	int i,j,o;
        for (i=0; i<SCALE*YSIZE; i++)   {
                for (j=0; j<SCALE*XSIZE; j++)
                        vbuf[i*SCALE*XSIZE+j] = cbuf[i/SCALE*XSIZE+j/SCALE];
        }
}

void setruler(void) {
	double x1,x2;
	x1 = scale(OFFX,OFFX+XSIZE*SCALE,DATA[SET]->xmin,DATA[SET]->xmax,0);
	x2 = scale(OFFX,OFFX+XSIZE*SCALE,DATA[SET]->xmin,DATA[SET]->xmax,WXSIZE);
        gtk_ruler_set_range(GTK_RULER(hrule), x1, x2, 0, x1);
	x1 = scale(OFFY,OFFY+YSIZE*SCALE,DATA[SET]->ymax,DATA[SET]->ymin,0);
	x2 = scale(OFFY,OFFY+YSIZE*SCALE,DATA[SET]->ymax,DATA[SET]->ymin,WYSIZE);
        gtk_ruler_set_range(GTK_RULER(vrule), x1, x2, 0, x1);
}

void setcrule(void) {
	gtk_ruler_set_range(GTK_RULER(crule), ZMIN, ZMAX, 0, ZMIN);
}

void redraw(void) {
	static int s=0;
	GdkColor bg;
	if (s!=SCALE) rescale();
	s=SCALE;
//	gdk_draw_gray_image(area->window, area->style->fg_gc[GTK_STATE_NORMAL],0,0,WXSIZE,WYSIZE,GDK_RGB_DITHER_NONE,blackbuf, WXSIZE);
//	gdk_draw_gray_image(area->window, area->style->fg_gc[GTK_STATE_NORMAL],OFFX,OFFY,SCALE*XSIZE,SCALE*YSIZE,GDK_RGB_DITHER_NONE,vbuf, SCALE*XSIZE);
	gdk_draw_indexed_image(area->window, area->style->fg_gc[GTK_STATE_NORMAL],OFFX,OFFY,SCALE*XSIZE,SCALE*YSIZE,GDK_RGB_DITHER_NONE,vbuf, SCALE*XSIZE, cmap);
//      gdk_rgb_gc_set_background(gc, 0); // gc ?
	gtk_widget_modify_base(GTK_WIDGET(area), GTK_STATE_NORMAL, &bg);
	setruler();
}

void redrawcbar(void) {
	gdk_draw_indexed_image(cbararea->window, cbararea->style->fg_gc[GTK_STATE_NORMAL],0,0,50,256,GDK_RGB_DITHER_NONE,cbarbuf, 50, cmap);
}

void recalc() {
	int i;
	double x;
	
	for (i=0; i<XSIZE*YSIZE; i++) {
		x = scale(ZMIN,ZMAX,0,255,DATA[SET]->data[i]);
		if ((x>0.0) && (x<255.0)) {
			cbuf[i] = (guchar)x;
		} else if (x<=0.0) {
			cbuf[i] = 0;
		} else {
			cbuf[i] = 255;
		}
	}
}

/*
void readfile(char *filename, int N) {
	int i;
	FILE *fp;
	fp = fopen(filename,"r");
	if (fp == NULL) {
		fprintf(stderr,"error: cannot open file %s\n", filename);
		exit(1);
	}
	for (i=0; i<N; i++) {
		if (fread(&dbuf[SET*XSIZE*YSIZE+i],sizeof(double), 1, fp) != 1) {
			fprintf(stderr,"error: cannot read from file %s at position %d\n", filename, i);
			exit(1);
		}
	}
	recalc();
}
*/

double height(double x, double y) {
	return 1.0;
}
static gboolean area_expose(GtkWidget *widget, GdkEvent *event, gpointer data) {
	int i;
	redraw();
	redrawcbar();
	return TRUE;
}
static gboolean delete_event(GtkWidget * widget, GdkEvent * event, gpointer data) {
	gtk_main_quit();
	return FALSE;
}
static gboolean keypress_event(GtkWidget * widget, GdkEventKey * event, gpointer data) {
	if (event->length > 0)
		if (event->string[0] == 27)
			gtk_main_quit();
	return FALSE;
}
static gboolean scrollevent(GtkWidget *widget, GdkEventScroll *event, gpointer data) {
	char buf[40];
	if (event->direction == GDK_SCROLL_UP )
		SET++;
	else if (event->direction == GDK_SCROLL_DOWN )
		SET--;
	if (SET == FILES)	SET=0;
	if (SET == -1 ) 	SET += FILES;

	snprintf(buf,sizeof(buf),"%d",SET);
	gtk_label_set(GTK_LABEL(setlabel), buf);
	snprintf(buf,sizeof(buf),"%s",DATA[SET]->filename);
	gtk_label_set(GTK_LABEL(flabel), buf);
	gtk_frame_set_label(GTK_FRAME(xframe), DATA[SET]->xlabel);
	gtk_frame_set_label(GTK_FRAME(yframe), DATA[SET]->ylabel);
	gtk_frame_set_label(GTK_FRAME(zframe), DATA[SET]->zlabel);
	recalc();
	rescale();
	redraw();
	return TRUE;
}

static gboolean buttonpressevent(GtkWidget *widget, GdkEventButton *event, gpointer data) {
	int i,j;
	char xstring[20], ystring[20], zstring[20];
	snprintf(xstring,sizeof(xstring),"%g",scale(OFFX,OFFX+XSIZE*SCALE,DATA[SET]->xmin,DATA[SET]->xmax,event->x));
	snprintf(ystring,sizeof(ystring),"%g",scale(OFFY,OFFY+YSIZE*SCALE,DATA[SET]->ymax,DATA[SET]->ymin,event->y));
	snprintf(zstring,sizeof(zstring),"%g",height(event->x,event->y));
	gtk_label_set(GTK_LABEL(xlabel), xstring);
	gtk_label_set(GTK_LABEL(ylabel), ystring);
	TFX = event->x;
	TFY = event->y;

	i = (event->x - OFFX)/SCALE;
	j = (event->y - OFFY)/SCALE;
	if ((i<0) || (i>XSIZE) || (j<0) || (j>YSIZE))
		snprintf(zstring,sizeof(zstring),"out");
	else 
		snprintf(zstring,sizeof(zstring),"%g", DATA[SET]->data[j*XSIZE+i]);
	gtk_label_set(GTK_LABEL(zlabel), zstring);
	
	return TRUE;

}

static gboolean buttonreleaseevent(GtkWidget *widget, GdkEventButton *event, gpointer data) {
	static int oldox=0, oldoy=0;
	TTX = event->x;
	TTY = event->y;
	OFFX += TTX-TFX;
	OFFY += TTY-TFY;
	if ((oldox!=OFFX)||(oldoy!=OFFY)) {
		rescale();
		redraw();
	}
	oldox = OFFX;
	oldoy = OFFY;
	return TRUE;
}

static gboolean colourmapcallback(GtkWidget *widget, GdkEventButton *event, gpointer data) {
	const char *c;
	c = gtk_entry_get_text(GTK_ENTRY(GTK_COMBO(colourmapcombo)->entry));
	if (strlen(c)>2)
		gdk_rgb_cmap_free(cmap);
	if (!(strcmp(c,"jet"))) cmap = gdk_rgb_cmap_new(cmap_jet, 256);
	else if (!(strcmp(c,"autumn"))) cmap = gdk_rgb_cmap_new(cmap_autumn, 256);
//	else if (!(strcmp(c,"blackwhite"))) cmap = gdk_rgb_cmap_new(cmap_blackwhite, 256);
	else if (!(strcmp(c,"contourmap123"))) cmap = gdk_rgb_cmap_new(cmap_contourmap123, 256);
	else if (!(strcmp(c,"contourmap15"))) cmap = gdk_rgb_cmap_new(cmap_contourmap15, 256);
	else if (!(strcmp(c,"contourmap3"))) cmap = gdk_rgb_cmap_new(cmap_contourmap3, 256);
	else if (!(strcmp(c,"contourmap31"))) cmap = gdk_rgb_cmap_new(cmap_contourmap31, 256);
	else if (!(strcmp(c,"contourmap61"))) cmap = gdk_rgb_cmap_new(cmap_contourmap61, 256);
	else if (!(strcmp(c,"contourmap7"))) cmap = gdk_rgb_cmap_new(cmap_contourmap7, 256);
	else if (!(strcmp(c,"cool"))) cmap = gdk_rgb_cmap_new(cmap_cool, 256);
	else if (!(strcmp(c,"coolmap"))) cmap = gdk_rgb_cmap_new(cmap_coolmap, 256);
	else if (!(strcmp(c,"coolmap2"))) cmap = gdk_rgb_cmap_new(cmap_coolmap2, 256);
	else if (!(strcmp(c,"copper"))) cmap = gdk_rgb_cmap_new(cmap_copper, 256);
	else if (!(strcmp(c,"graymap"))) cmap = gdk_rgb_cmap_new(cmap_graymap, 256);
	else if (!(strcmp(c,"graymap2"))) cmap = gdk_rgb_cmap_new(cmap_graymap2, 256);
	else if (!(strcmp(c,"hot"))) cmap = gdk_rgb_cmap_new(cmap_hot, 256);
	else if (!(strcmp(c,"hotmap"))) cmap = gdk_rgb_cmap_new(cmap_hotmap, 256);
	else if (!(strcmp(c,"hotmap2"))) cmap = gdk_rgb_cmap_new(cmap_hotmap2, 256);
	else if (!(strcmp(c,"hotmap3"))) cmap = gdk_rgb_cmap_new(cmap_hotmap3, 256);
	else if (!(strcmp(c,"hsvmap"))) cmap = gdk_rgb_cmap_new(cmap_hsvmap, 256);
//	else if (!(strcmp(c,"ocean"))) cmap = gdk_rgb_cmap_new(cmap_ocean, 256);
	else if (!(strcmp(c,"rainbow"))) cmap = gdk_rgb_cmap_new(cmap_rainbow, 256);
	else if (!(strcmp(c,"spring"))) cmap = gdk_rgb_cmap_new(cmap_spring, 256);
	else if (!(strcmp(c,"summer"))) cmap = gdk_rgb_cmap_new(cmap_summer, 256);
	else if (!(strcmp(c,"wackymap"))) cmap = gdk_rgb_cmap_new(cmap_wackymap, 256);
	else if (!(strcmp(c,"wackymap2"))) cmap = gdk_rgb_cmap_new(cmap_wackymap2, 256);
//	else if (!(strcmp(c,"whiteblack"))) cmap = gdk_rgb_cmap_new(cmap_whiteblack, 256);
	else if (!(strcmp(c,"winter"))) cmap = gdk_rgb_cmap_new(cmap_winter, 256);
	else {
		cmap = gdk_rgb_cmap_new(cmap_graymap, 256);
	}
	redrawcbar();
	redraw();
}

static void zminmaxcallback(GtkWidget *widget) {
	char zminstring[20], zmaxstring[20];
	gdouble zmin0, zmax0;
	gdouble ezmin, ezmax;
	zmin0 = gtk_spin_button_get_value(GTK_SPIN_BUTTON(zminspinner));
	zmax0 = gtk_spin_button_get_value(GTK_SPIN_BUTTON(zmaxspinner));
	ezmin = gtk_spin_button_get_value(GTK_SPIN_BUTTON(ezminspinner));
	ezmax = gtk_spin_button_get_value(GTK_SPIN_BUTTON(ezmaxspinner));
	ZMIN = zmin0*pow(10.0,ezmin);
	ZMAX = zmax0*pow(10.0,ezmax);

	snprintf(zminstring,sizeof(zminstring),"%g",ZMIN);
	snprintf(zmaxstring,sizeof(zmaxstring),"%g",ZMAX);
	gtk_label_set(GTK_LABEL(zminlabel), zminstring);
	gtk_label_set(GTK_LABEL(zmaxlabel), zmaxstring);
	recalc();
	rescale();
	redraw();
	redrawcbar();
}
static void scalecallback(GtkWidget *widget) {
	int s;
	s = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(scalespinner));
	SCALE = s;
	redraw();
}

void usage(char **argv) {
	fprintf(stderr,"usage: %s datafile [datafile2 ...]\n",argv[0]);
	exit(1);
}

int main(int args, char **argv) {
	GtkWidget      *window;
	GtkWidget      *button;
	GtkWidget	*redrawbox, *redrawbutton;
	GtkWidget	*zminbox,*zmaxbox,*ezminbox,*ezmaxbox,*scalebox;
	GtkWidget      *box1, *box2, *box3, *middlebox, *cbox;
	GtkWidget      *separator, *vseparator, *bseparator, *mseparator;
	GtkWidget      *label; 
	GtkWidget	*fframe, *setframe;
	GtkWidget      *quitbox;
	GtkWidget	*picbox, *picarea, *table;

	GtkWidget	*image;
	GdkPixbuf	*origpixbuf;
	GdkPixbuf	*scaledpixbuf;

	GList 	*colourmaplist;
	GtkAdjustment	*adj;
	GdkCursor 	*cursor = NULL;
	char *var;
        char zminstring[20], zmaxstring[20];
	int i,j;

	if (var=getenv("SCALE")) SCALE=atoi(var);
	if (SCALE>MAXSCALE)	SCALE=MAXSCALE;
	/*
	if (var=getenv("zmin")) ZMIN=atof(var);
	if (var=getenv("zmax")) ZMAX=atof(var);
	*/

	if (args < 2)	usage(argv);
	FILES = args-1;

	// memory for the datafile structure
	DATA = (struct datafile **)malloc(FILES*(ssize_t)sizeof(struct datafile*));
	if (DATA == NULL) { perror("DATA"); exit(1); }
	for (i=0; i<FILES; i++) {
		printf("malloc DATA[%d]: %d",i,sizeof(struct datafile));
		DATA[i] = (struct datafile *)malloc(sizeof(struct datafile));
		if (DATA[i] == NULL) { perror("DATA structure"); exit(1); }
		printf(".\n");
		readfile(DATA[i], argv[i+1]);
	}
	SET = 0;
	XSIZE = DATA[SET]->cols;
	YSIZE = DATA[SET]->rows;

	// memory for picture buffers
	vbuf = (guchar *)malloc(MAXSCALE*MAXSCALE*XSIZE*YSIZE*(size_t)sizeof(guchar));
	if (vbuf == NULL) { perror("vbuf"); exit(1); }
	cbuf = (guchar *)malloc(XSIZE*YSIZE*(size_t)sizeof(guchar));
	if (cbuf == NULL) { perror("cbuf"); exit(1); }

	// dimension test, otherwise overflow of [cv]buf
	for (i=0; i<FILES; i++) {
		if ( (DATA[i]->cols != XSIZE) || (DATA[i]->rows != YSIZE) ) {
			fprintf(stderr,"error: datasets must have same dimensions. Set [%d] %dx%d != %dx%d\n",i,DATA[i]->cols,DATA[i]->rows,XSIZE,YSIZE);
			exit(1);
		}
		printf("File %s:\n%s\n",DATA[i]->filename,DATA[i]->header);
	}

	recalc();
	WXSIZE=SCALE*XSIZE;
	WYSIZE=SCALE*YSIZE;

	for (i=0; i<50*256; i++) cbarbuf[i] = 255-i/50;

	gtk_init(&args, &argv);
	cmap = gdk_rgb_cmap_new(cmap_jet, 256);
	window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	g_signal_connect(G_OBJECT(window), "delete_event", G_CALLBACK(delete_event), NULL);
	gtk_container_set_border_width(GTK_CONTAINER(window), 0);

//	g_signal_connect(G_OBJECT(window), "configure_event", G_CALLBACK(configure_event), NULL);
	g_signal_connect(G_OBJECT(window), "key_press_event", G_CALLBACK(keypress_event), NULL);

	box1 = gtk_vbox_new(FALSE, 0);
	box2 = gtk_hbox_new(FALSE,0);

	/* zmin zmax spinner */
	ezmaxbox = gtk_vbox_new(FALSE,0);
	label = gtk_label_new ("E");
	gtk_misc_set_alignment (GTK_MISC (label), 0, 0.5);
	gtk_box_pack_start (GTK_BOX (ezmaxbox), label, FALSE, TRUE, 0);
	//                                         (init, lower, upper, stepinc, pageinc, unused)
	adj = (GtkAdjustment *) gtk_adjustment_new (0.0, -100, 100, 1, 10, 0.0);
	ezmaxspinner = gtk_spin_button_new (adj, 1.0, 0);
	gtk_box_pack_start(GTK_BOX (ezmaxbox), ezmaxspinner, FALSE, TRUE, 0);
	gtk_box_pack_end(GTK_BOX(box2), ezmaxbox, FALSE, FALSE, 0);
	gtk_widget_show(label);
	gtk_widget_show(ezmaxspinner);
	gtk_widget_show(ezmaxbox);
	g_signal_connect(G_OBJECT(adj),"value_changed",G_CALLBACK(zminmaxcallback),NULL);

	zmaxbox = gtk_vbox_new(FALSE,0);
	label = gtk_label_new ("zmax");
	gtk_misc_set_alignment (GTK_MISC (label), 0, 0.5);
	gtk_box_pack_start (GTK_BOX (zmaxbox), label, FALSE, TRUE, 0);
	//                                         (init, lower, upper, stepinc, pageinc, unused)
	adj = (GtkAdjustment *) gtk_adjustment_new (ZMAX, -10.0, 10.0, 0.01, 0.1, 0.0);
	zmaxspinner = gtk_spin_button_new (adj, 1.0, 2);
	gtk_spin_button_set_wrap (GTK_SPIN_BUTTON (zmaxspinner), TRUE);
	gtk_spin_button_set_snap_to_ticks (GTK_SPIN_BUTTON (zmaxspinner), FALSE);
	gtk_box_pack_start(GTK_BOX (zmaxbox), zmaxspinner, FALSE, TRUE, 0);
	gtk_box_pack_end(GTK_BOX(box2), zmaxbox, FALSE, FALSE, 0);
	gtk_widget_show(label);
	gtk_widget_show(zmaxspinner);
	gtk_widget_show(zmaxbox);
	g_signal_connect(G_OBJECT(adj),"value_changed",G_CALLBACK(zminmaxcallback),NULL);

	ezminbox = gtk_vbox_new(FALSE,0);
	label = gtk_label_new ("E");
	gtk_misc_set_alignment (GTK_MISC (label), 0, 0.5);
	gtk_box_pack_start (GTK_BOX (ezminbox), label, FALSE, TRUE, 0);
	adj = (GtkAdjustment *) gtk_adjustment_new (0.0, -100, 100, 1, 10, 0.0);
	ezminspinner = gtk_spin_button_new (adj, 1.0, 0);
	gtk_box_pack_start(GTK_BOX (ezminbox), ezminspinner, FALSE, TRUE, 0);
	gtk_box_pack_end(GTK_BOX(box2), ezminbox, FALSE, FALSE, 0);
	gtk_widget_show(label);
	gtk_widget_show(ezminspinner);
	gtk_widget_show(ezminbox);
	g_signal_connect(G_OBJECT(adj),"value_changed",G_CALLBACK(zminmaxcallback),NULL);

	zminbox = gtk_vbox_new(FALSE,0);
	label = gtk_label_new ("zmin");
	gtk_misc_set_alignment (GTK_MISC (label), 0, 0.5);
	gtk_box_pack_start (GTK_BOX (zminbox), label, FALSE, TRUE, 0);
	adj = (GtkAdjustment *) gtk_adjustment_new (ZMIN, -10.0, 10.0, 0.01, 0.1, 0.0);
	zminspinner = gtk_spin_button_new (adj, 1.0, 2);
	gtk_spin_button_set_wrap (GTK_SPIN_BUTTON (zminspinner), TRUE);
	gtk_spin_button_set_snap_to_ticks (GTK_SPIN_BUTTON (zmaxspinner), TRUE);
	gtk_box_pack_start(GTK_BOX (zminbox), zminspinner, FALSE, TRUE, 0);
	gtk_box_pack_end(GTK_BOX(box2), zminbox, FALSE, FALSE, 0);
	gtk_widget_show(label);
	gtk_widget_show(zminspinner);
	gtk_widget_show(zminbox);
	g_signal_connect(G_OBJECT(adj),"value_changed",G_CALLBACK(zminmaxcallback),NULL);

	/* vseparator */
	vseparator = gtk_vseparator_new();
	gtk_box_pack_end(GTK_BOX(box2), vseparator, FALSE, FALSE, 0);
	gtk_widget_show(vseparator);

	/* scale spinner */	
	scalebox = gtk_vbox_new(FALSE,0);
	label = gtk_label_new ("scale");
	gtk_misc_set_alignment (GTK_MISC (label), 0, 0.5);
	gtk_box_pack_start (GTK_BOX (scalebox), label, FALSE, TRUE, 0);
	adj = (GtkAdjustment *) gtk_adjustment_new (SCALE, 1, MAXSCALE, 1, 1, 0.0);
	scalespinner = gtk_spin_button_new (adj, 1.0, 2);
	gtk_spin_button_set_wrap (GTK_SPIN_BUTTON (scalespinner), FALSE);
	gtk_box_pack_start(GTK_BOX (scalebox), scalespinner, FALSE, TRUE, 0);
	gtk_box_pack_end(GTK_BOX(box2), scalebox, FALSE, FALSE, 0);
	gtk_widget_show(label);
	gtk_widget_show(scalespinner);
	gtk_widget_show(scalebox);
	g_signal_connect(G_OBJECT(adj),"value_changed",G_CALLBACK(scalecallback),NULL);


	/* quitbox */
	quitbox = gtk_vbox_new(FALSE, 0);
	button = gtk_button_new_with_label("Quit");
	g_signal_connect_swapped(G_OBJECT(button), "clicked", G_CALLBACK(gtk_main_quit), G_OBJECT(window));
	gtk_box_pack_start(GTK_BOX(quitbox), button, TRUE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(box2), quitbox, FALSE, FALSE, 0);
	gtk_widget_show(button);
	gtk_widget_show(quitbox);

	/* colourmap combo */
	colourmapcombo = gtk_combo_new();
	//gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(colourmapcombo)->entry),"greyscale");
	gtk_box_pack_start(GTK_BOX(box2), colourmapcombo, FALSE, FALSE, 0);
	colourmaplist = NULL;
	colourmaplist = g_list_append(colourmaplist, "jet");
	colourmaplist = g_list_append(colourmaplist, "autumn");
	colourmaplist = g_list_append(colourmaplist, "rainbow");
//	colourmaplist = g_list_append(colourmaplist, "blackwhite");
	colourmaplist = g_list_append(colourmaplist, "contourmap123");
	colourmaplist = g_list_append(colourmaplist, "contourmap15");
	colourmaplist = g_list_append(colourmaplist, "contourmap3");
	colourmaplist = g_list_append(colourmaplist, "contourmap31");
	colourmaplist = g_list_append(colourmaplist, "contourmap7");
	colourmaplist = g_list_append(colourmaplist, "cool");
	colourmaplist = g_list_append(colourmaplist, "coolmap");
	colourmaplist = g_list_append(colourmaplist, "coolmap2");
	colourmaplist = g_list_append(colourmaplist, "copper");
	colourmaplist = g_list_append(colourmaplist, "graymap");
	colourmaplist = g_list_append(colourmaplist, "graymap2");
	colourmaplist = g_list_append(colourmaplist, "hot");
	colourmaplist = g_list_append(colourmaplist, "hotmap");
	colourmaplist = g_list_append(colourmaplist, "hotmap2");
	colourmaplist = g_list_append(colourmaplist, "hotmap3");
	colourmaplist = g_list_append(colourmaplist, "hsvmap");
//	colourmaplist = g_list_append(colourmaplist, "ocean");
	colourmaplist = g_list_append(colourmaplist, "rainbow");
	colourmaplist = g_list_append(colourmaplist, "spring");
	colourmaplist = g_list_append(colourmaplist, "summer");
	colourmaplist = g_list_append(colourmaplist, "wackymap");
	colourmaplist = g_list_append(colourmaplist, "wackymap2");
//	colourmaplist = g_list_append(colourmaplist, "whiteblack");
	colourmaplist = g_list_append(colourmaplist, "winter");
	gtk_combo_set_popdown_strings(GTK_COMBO(colourmapcombo),colourmaplist);
	gtk_widget_show(colourmapcombo);
//	g_signal_connect(G_COBJECT(GTK_COMBO(colourmapcombo->entry), "activate",G_CALLBACK(colourmapcallback),(gpointer)my))
	g_signal_connect(G_OBJECT(GTK_COMBO(colourmapcombo)->entry),"changed",G_CALLBACK(colourmapcallback),NULL);

	/* set label */
	setframe = gtk_frame_new("set");
	setlabel = gtk_label_new("0");
	gtk_container_add(GTK_CONTAINER(setframe),setlabel);
	gtk_box_pack_start(GTK_BOX(box2), setframe, FALSE, FALSE, 0);
	gtk_widget_show(setframe);
	gtk_widget_show(setlabel);

	/* Pack box2 into box1 (the vbox remember ? :) */
	gtk_box_pack_start(GTK_BOX(box1), box2, FALSE, FALSE, 0);
	gtk_widget_show(box2);

	/* A separator for the bottom. */
	separator = gtk_hseparator_new();
	gtk_widget_set_size_request(separator, WXSIZE, 1);
	gtk_box_pack_start(GTK_BOX(box1), separator, FALSE, TRUE, 0);
	gtk_widget_show(separator);
	
	/* middle box */
	middlebox = gtk_hbox_new(FALSE, 0);
	gtk_box_pack_start(GTK_BOX(box1), middlebox, FALSE, TRUE, 0);

	cbox = gtk_vbox_new(FALSE,0);
	gtk_box_pack_start(GTK_BOX(middlebox), cbox, FALSE, TRUE, 0);
	gtk_widget_show(cbox);

	/* zmin/zmax colourmap label */
	snprintf(zminstring,sizeof(zminstring),"%g",ZMIN);
	snprintf(zmaxstring,sizeof(zmaxstring),"%g",ZMAX);

	zmaxlabel = gtk_label_new(zminstring);
	gtk_box_pack_start(GTK_BOX(cbox), zmaxlabel, FALSE, FALSE, 0);
	gtk_widget_show(zmaxlabel);

	cbararea = gtk_drawing_area_new();
	gtk_widget_set_size_request(GTK_WIDGET(cbararea), 50, 256);
	gtk_box_pack_start(GTK_BOX(cbox), cbararea, FALSE, TRUE, 0);
	gtk_widget_show(cbararea);

	zminlabel = gtk_label_new(zmaxstring);
	gtk_box_pack_start(GTK_BOX(cbox), zminlabel, FALSE, FALSE, 0);
	gtk_widget_show(zminlabel);

	gtk_label_set(GTK_LABEL(zminlabel), zminstring);
	gtk_label_set(GTK_LABEL(zmaxlabel), zmaxstring);


	mseparator = gtk_vseparator_new();
//	gtk_widget_set_size_request(separator, WXSIZE, 1);
	gtk_box_pack_start(GTK_BOX(middlebox), mseparator, FALSE, TRUE, 0);
	gtk_widget_show(mseparator);

	/* picbox */
	picbox = gtk_hbox_new(FALSE, 0);

	picarea = gtk_event_box_new();
	gtk_box_pack_start(GTK_BOX(picbox),picarea,TRUE,FALSE,0);

	table = gtk_table_new(3,2,FALSE);
//	gtk_box_pack_start(GTK_BOX(picbox),table,TRUE,FALSE,0);
	gtk_container_add(GTK_CONTAINER(picarea),table);
//	g_signal_connect(G_OBJECT(picarea), "configure_event", G_CALLBACK(picarea_configure), NULL);

	/* drawing area */
//	area = gtk_event_box_new();
	area = gtk_drawing_area_new();
	gtk_widget_set_size_request(GTK_WIDGET(area), WXSIZE, WYSIZE);
	gtk_table_attach(GTK_TABLE(table), area, 1, 2, 1, 2, GTK_EXPAND | GTK_FILL, GTK_FILL, 0, 0);
	gtk_widget_set_events(area, GDK_SCROLL | GDK_BUTTON_PRESS_MASK | GDK_BUTTON_RELEASE_MASK | GDK_POINTER_MOTION_MASK | GDK_POINTER_MOTION_HINT_MASK);
//	gtk_widget_set_events(area, GDK_BUTTON_PRESS_MASK);
	g_signal_connect(G_OBJECT(area),"scroll_event",G_CALLBACK(scrollevent),NULL);
	g_signal_connect(G_OBJECT(area),"button_press_event",G_CALLBACK(buttonpressevent),NULL);
	g_signal_connect(G_OBJECT(area),"button_release_event",G_CALLBACK(buttonreleaseevent),NULL);
	g_signal_connect(G_OBJECT(area),"expose_event",G_CALLBACK(area_expose),NULL);
//	g_signal_connect(G_OBJECT(area), "configure_event", G_CALLBACK(area_configure), NULL);

	/* cursor ??
	cursor = gdk_cursor_new(GDK_CROSS);
	gdk_window_set_cursor(GDK_WINDOW(image), cursor);
//	gdk_cursor_destroy(cursor);
	*/

        hrule = gtk_hruler_new();
        gtk_ruler_set_metric(GTK_RULER(hrule), GTK_PIXELS);
        //gtk_ruler_set_range(GTK_RULER(hrule), YMIN, YMAX, 0, SCALE*YSIZE);
        g_signal_connect_swapped(G_OBJECT(area), "motion_notify_event", G_CALLBACK(EVENT_METHOD(hrule, motion_notify_event)), G_OBJECT(hrule));
        gtk_table_attach(GTK_TABLE(table), hrule, 1, 2, 0, 1, GTK_EXPAND | GTK_SHRINK | GTK_FILL, GTK_FILL, 0, 0);

        vrule = gtk_vruler_new();
        gtk_ruler_set_metric(GTK_RULER(vrule), GTK_PIXELS);
        //gtk_ruler_set_range(GTK_RULER(vrule), XMAX, XMIN, 0, SCALE*YSIZE);
        g_signal_connect_swapped(G_OBJECT(area), "motion_notify_event", G_CALLBACK(EVENT_METHOD(vrule, motion_notify_event)), G_OBJECT(vrule));
        //gtk_table_attach(GTK_TABLE(table), vrule, 0, 1, 1, 2, GTK_FILL, GTK_EXPAND | GTK_SHRINK | GTK_FILL, 0, 0);
        gtk_table_attach(GTK_TABLE(table), vrule, 0, 1, 1, 2, GTK_EXPAND, GTK_SHRINK | GTK_FILL | GTK_FILL, 0, 0);

	setruler();

	gtk_box_pack_start(GTK_BOX(middlebox), picbox, TRUE, FALSE,0);

	/* bottom separator */
	bseparator = gtk_hseparator_new();
	gtk_widget_set_size_request(bseparator, WXSIZE, 1);
	gtk_box_pack_start(GTK_BOX(box1), bseparator, FALSE, TRUE, 0);
	gtk_widget_show(bseparator);

	/* status line box */
	box3 = gtk_hbox_new(FALSE,0);
	gtk_box_pack_start(GTK_BOX(box1), box3, FALSE, TRUE, 0);
	gtk_widget_show(box3);

	/* coordinates */
	zframe = gtk_frame_new(DATA[SET]->zlabel);
	zlabel = gtk_label_new("000.000000");
	gtk_container_add(GTK_CONTAINER(zframe),zlabel);
	gtk_box_pack_end(GTK_BOX(box3), zframe, FALSE, FALSE, 0);
	gtk_widget_show(zframe);
	gtk_widget_show(zlabel);

	yframe = gtk_frame_new(DATA[SET]->ylabel);
	ylabel = gtk_label_new("000.000000");
	gtk_container_add(GTK_CONTAINER(yframe),ylabel);
	gtk_box_pack_end(GTK_BOX(box3), yframe, FALSE, FALSE, 0);
	gtk_widget_show(yframe);
	gtk_widget_show(ylabel);

	xframe = gtk_frame_new(DATA[SET]->xlabel);
	xlabel = gtk_label_new("000.000000");
	gtk_container_add(GTK_CONTAINER(xframe),xlabel);
	gtk_box_pack_end(GTK_BOX(box3), xframe, FALSE, FALSE, 0);
	gtk_widget_show(xframe);
	gtk_widget_show(xlabel);

	/* filename */
	fframe = gtk_frame_new("file");
	flabel = gtk_label_new(DATA[SET]->filename);
	gtk_container_add(GTK_CONTAINER(fframe),flabel);
	gtk_box_pack_end(GTK_BOX(box3), fframe, FALSE, FALSE, 0);
	gtk_widget_show(fframe);
	gtk_widget_show(flabel);

	gtk_container_add(GTK_CONTAINER(window), box1);

	gtk_widget_show(area);	
	gtk_widget_show(hrule);	
	gtk_widget_show(vrule);	
	gtk_widget_show(table);
	gtk_widget_show(picarea);
	gtk_widget_show(picbox);
	gtk_widget_show(middlebox);

	gtk_widget_show(box1);
	gtk_widget_show(window);

	gtk_main();
	return 0;
}
