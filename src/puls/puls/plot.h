void plot_time(double x, double xmin, double xmax); 
void plotpressure(double x, double y, double xmin, double xmax, double ymin, double ymax, int first);
void plotflow(double x, double y, double xmin, double xmax, double ymin, double ymax, int first);
void plotfflow(double x, double y, double xmin, double xmax, double ymin, double ymax, int first);
void plotinflexion(double x, double xmin, double xmax);
void plot_inflexionline(int *, int);
void plot (double x, double y, double xmin, double xmax, 
		double ymin, double ymax, int state, int linestyle, int linewidth, const char *color);
void plotinit (double xmin, double xmax, double ymin, double ymax, 
		const char *fgcolorstr, const char *bgcolorstr, 
		const char *xlabel, const char *ylabel, const char *title, int width, int height);
void plotbackground (double xmin, double xmax, double ymin, double ymax, 
		const char *fgcolorstr, const char *bgcolorstr, 
		const char *xlabel, const char *ylabel, const char *title, int width, int height);
void plotdisplay(void);
void clearplot(void);
void plotexit(void);
