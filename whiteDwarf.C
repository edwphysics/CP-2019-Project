//IMPLEMENTATION OF THE RUNGE-KUTTA 4-TH ORDER
//WHITE DWARF - Relativistic Case
//y = pressure, z = mass

/*
	COMMENTS: Results for three initial values of pressure
	ADDITIONS: A bucle for evey new p0, prints 3 graphs for every p0
	ISSUES:
*/

//GLOBAL CONSTANTS AND VARIABLES
int n = 8000, l = 1;
double alpha = 1.473, beta = 52.46; //DE Constants
double DeltaT = 1e-30;				//Sliding constant to avoid r=0 divergence
  
double f(double y[], double z[], double t, int j){
	//GENERALIZED VELOCITY
	double v[l];
	v[0] = -alpha*pow(sqrt(sqrt(y[0])),3)*z[0]/pow(t,2);
	return v[j];
}

double g(double y[], double z[], double t, int j){
	//GENERALIZED VELOCITY
	double v[l];
	v[0] = beta*pow(t,2)*pow(sqrt(sqrt(y[0])),3);
	return v[j];
}

double RungeKuttaG(double y[], double z[], double t, double dt, int j){
	double d1[l], d2[l], d3[l], d4[l];
	double c1[l], c2[l], c3[l], c4[l];
	double D1[l], D2[l], D3[l], D4[l];
	double C1[l], C2[l], C3[l], C4[l];

	d1[j] = f(y, z, t, j);
	c1[j] = g(y, z, t, j);

	D2[j] = y[j] + dt*d1[j]/2;
	C2[j] = z[j] + dt*c1[j]/2;
	d2[j] = f(D2, C2, t+dt/2, j);
	c2[j] = g(D2, C2, t+dt/2, j);

	D3[j] = y[j] + dt*d2[j]/2;
	C3[j] = z[j] + dt*c2[j]/2;
	d3[j] = f(D3, C3, t+dt/2, j);
	c3[j] = g(D3, C3, t+dt/2, j);

	D4[j] = y[j] + dt*d3[j];
	C4[j] = z[j] + dt*c3[j];
	d4[j] = f(D4, C4, t + dt, j);
	c4[j] = g(D4, C4, t + dt, j);

	d1[j] = y[j] + dt*(d1[j]+2*(d2[j]+d3[j])+d4[j])/6;
	return d1[j];
}

double RungeKuttaF(double y[], double z[], double t, double dt, int j){
	double d1[l], d2[l], d3[l], d4[l];
	double c1[l], c2[l], c3[l], c4[l];
	double D2[l], D3[l], D4[l];
	double C2[l], C3[l], C4[l];

	d1[j] = f(y, z, t, j);
	c1[j] = g(y, z, t, j);

	D2[j] = y[j] + dt*d1[j]/2;
	C2[j] = z[j] + dt*c1[j]/2;
	d2[j] = f(D2, C2, t+dt/2, j);
	c2[j] = g(D2, C2, t+dt/2, j);

	D3[j] = y[j] + dt*d2[j]/2;
	C3[j] = z[j] + dt*c2[j]/2;
	d3[j] = f(D3, C3, t+dt/2, j);
	c3[j] = g(D3, C3, t+dt/2, j);

	D4[j] = y[j] + dt*d3[j];
	C4[j] = z[j] + dt*c3[j];
	d4[j] = f(D4, C4, t + dt, j);
	c4[j] = g(D4, C4, t + dt, j);

	c1[j] = z[j] + dt*(c1[j]+2*(c2[j]+c3[j])+c4[j])/6;
	return c1[j];
}

void principal(double p0){
	//FLAG & CRITERION TO DISPLAY THE RESULT
	int signal;
	bool signalb = false;
	double criterion = 1e-22;

	//VARIABLES
	double y1[n+1], y[l];
	double z1[n+1], z[l];

	//SET UP TIME STEP AND INITIAL VALUES
	double dt = 15100./n;
	y1[0] = p0;
	y[0] = p0;

	z1[0] = 0.;
	z[0] = 0.;

	double T[n+1];

	//PERFORM THE 4TH-ORDER RUNGE-KUTTA INTEGRATION
	for (int i = 0; i < n; ++i)
	{
		double t = dt*i + DeltaT;
		T[i] = t;

		for (int j = 0; j < l; ++j)
		{
			y[j] = RungeKuttaG(y, z, t, dt, j);
			z[j] = RungeKuttaF(y, z ,t, dt, j);
		}

		y1[i+1] = y[0];
		z1[i+1] = z[0];		

		//FIND THE FINAL RADIUS AND MASS VALUE
		if ((y[i+1] <= criterion) && (signalb == false) && (i>10))
		{
			signal = i+1;
			signalb = true;
		}
	}

	//PRINT RESULT
	printf("%.1e\t%lf\t%lf\n", p0, T[signal], z[signal]);

	//GRAPHS
	TCanvas *c1 = new TCanvas();
	c1->Divide(2,1);
	c1->cd(1);
	TGraphErrors *g1 = new TGraphErrors(signal, T, y1);
	g1->SetTitle("Pression; r(km); #bar{p}");
	g1->SetLineColor(kRed);
	TAxis* a = g1->GetXaxis();
    a->SetMaxDigits(5);
	g1->Draw("al");
	c1->cd(2);
	TGraphErrors *g2 = new TGraphErrors(signal, T, z1);
	g2->SetTitle("Mass; r(km); #bar{#Mu}");
	g2->SetLineColor(kBlue);
	TAxis* b = g2->GetXaxis();
    b->SetMaxDigits(5);
	g2->Draw("al");
}

void whiteDwarfV5(){
	double p0 = 1e-16; //Initial Pressure

	printf("p0\tR\t\tMass\n");
	for (int i = 0; i < 3; ++i)
	{
		principal(p0);
		p0 *= 1e1;
	}
}