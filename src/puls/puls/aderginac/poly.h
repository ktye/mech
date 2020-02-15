ex polyintegrate(ex p, symbol x) {
	ex poly;
	for (int i=p.ldegree(x); i<=p.degree(x); ++i) {
		poly = poly + p.coeff(x,i)*pow(x,i+1)/(i+1); 
	}
	return poly;
}

ex innerex(ex p, ex q, ex w, symbol x) {
	ex a = 0;
	ex b = 1;
//	ex w;
//	w=x;
	ex i, res;
	i = polyintegrate(expand(p*q*w),x);
	return i.subs(x==b)-i.subs(x==a);
	
}
numeric inner(ex p, ex q, ex w, symbol x) {
	ex res = innerex(p, q, w, x);
	if (is_a<numeric>(res))
		return ex_to<numeric>(res);
	else {
		cerr << "error in inner: res = " << " is not numeric" << endl;
		exit(1);
	}
}

/*
 * ansatzpolynoms
 * p0: prod(x-zero_i) * weight0;
 */
lst polyansatz(lst zeros, ex weight, ex weight0, symbol x, int N){
	lst l;
	ex f = 1;
	for (lst::const_iterator i = zeros.begin(); i != zeros.end(); ++i)
		f = f*(x - *i);
	l.append(f*weight0);

	for (int i=1; i<=N; i++) {
		lst eqns, vars;
		ex F;
		ex a,p;
		ex unscaled;
		lst xlst;
		for (int k=1; k<=i; k++) {
			xlst.append(pow(x,k));
		}
		matrix X(1,i,xlst);
		a = symbolic_matrix(i,1,"a");
		p = 1 + (X*a).evalm()[0];
		p = expand(p*f);
		eqns.remove_all();
		vars.remove_all();
		for (int j=1; j<=i; j++) {
			F = l[j-1];
			ex gl = innerex(F,p,weight,x)==0;
			eqns.append(gl);
		}
		for (int j=1; j<=i; j++) {
			vars.append( a.op(j-1) );
		}
		l.append(p.subs(lsolve(eqns,vars)));
	}

	for (int i=0; i<=N; i++)
		l[i] = l[i]/sqrt(innerex(l[i],l[i],weight,x));
	
	return l;
}
