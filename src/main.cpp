#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>

#ifdef _WIN32
	#define popen _popen
    #define pclose _pclose
#endif

double U = 7.272e-5;
const double d_time = 0.00500;

const double rad2Deg = 180 / 3.14159265358979323846;
const double deg2Rad = 3.14159265358979323846 / 180;

double fi = 56.0000037975 * deg2Rad;
double lam = 36.9999992699 * deg2Rad;
double h = 101.455;

double a = 6378137.0; 
double e2 = 0.00669437999014;
double Re = a; 
double g = 9.81; 


double f_b[3];
double w_b[3];

double Cb_ll[9] = { 1, 0, 0,
					0, 1, 0,
					0, 0, 1 };


struct Data {
	double v_ll[3] = { 0.0, 0.0, 0 }; // Ve, Vn, Vup
	double rpy[3] = { 0, 0, 0 }; // roll, pitch, yaw 
	double coords[3] = { fi, lam, h };

	double f_ll[3] = { 0, 0, 0 };
	double w_ll[3] = { 0, 0, 0 };
};

struct GraphData {
	GraphData() { 
		time.reserve(maxSize);
		v_e.reserve(maxSize);
		v_n.reserve(maxSize);
		pitch.reserve(maxSize);
		roll.reserve(maxSize);
		yaw.reserve(maxSize);
		fi.reserve(maxSize);
		lam.reserve(maxSize);
		h.reserve(maxSize);
	}
	const int maxSize = 727000;

	std::vector<double> time;
	std::vector<double> v_e;
	std::vector<double> v_n;
	std::vector<double> pitch;
	std::vector<double> roll;
	std::vector<double> yaw;
	std::vector<double> fi;
	std::vector<double> lam;
	std::vector<double> h;
};

void initGnuplot(FILE*& gp)
{
	gp = popen("gnuplot -persist", "w");
	if (!gp) {
		std::cerr << "Cannot open gnuplot. Make sure it's installed and in PATH" << std::endl;
		exit(-1);
	}

	fprintf(gp, "set terminal qt size 1400,1000\n");
	fprintf(gp, "set multiplot layout 4,2\n");
	fprintf(gp, "set style data lines\n");
	fprintf(gp, "set grid\n");  
	fprintf(gp, "set autoscale\n");  
}

void closeGnuplot(FILE*& gp)
{
	fprintf(gp, "unset multiplot\n");
	pclose(gp);
}

void drawParam(FILE* gp, std::string title, const std::vector<double>& etParam, const std::vector<double>& calcParam, const std::vector<double>& time)
{
	static int counter = 0;
	std::string etBlock = "$etalon" + std::to_string(counter);
	std::string calcBlock = "$calc" + std::to_string(counter);
	counter++;

	fprintf(gp, "%s << EOD\n", etBlock.c_str());
	for (size_t i = 0; i < etParam.size(); ++i) {
		fprintf(gp, "%lf %lf\n", time[i], etParam[i]);
	}
	fprintf(gp, "EOD\n");

	fprintf(gp, "%s << EOD\n", calcBlock.c_str());
	for (size_t i = 0; i < calcParam.size(); ++i) {
		fprintf(gp, "%lf %lf\n", time[i], calcParam[i]);
	}
	fprintf(gp, "EOD\n");

	fprintf(gp, "set title '%s'\n", title.c_str());
	fprintf(gp, "set xlabel 'time'\n");
	fprintf(gp, "set ylabel 'deg'\n");
	fprintf(gp, "plot %s using 1:2 with lines title 'etalon', "
		"%s using 1:2 with lines title 'calc'\n",
		etBlock.c_str(), calcBlock.c_str());
}

void drawGraphs(const GraphData& etalon, const GraphData& calc)
{
	FILE* gp;

	initGnuplot(gp);

	drawParam(gp, "Roll", etalon.roll, calc.roll, etalon.time); 
	drawParam(gp, "Pitch", etalon.pitch, calc.pitch, etalon.time); 
	drawParam(gp, "Yaw", etalon.yaw, calc.yaw, etalon.time); 
	drawParam(gp, "Latitude", etalon.fi, calc.fi, etalon.time); 
	drawParam(gp, "Longitude", etalon.lam, calc.lam, etalon.time); 
	drawParam(gp, "Altitude", etalon.h, calc.h, etalon.time); 
	drawParam(gp, "V_e", etalon.v_e, calc.v_e, etalon.time); 
	drawParam(gp, "V_n", etalon.v_n, calc.v_n, etalon.time); 

	closeGnuplot(gp);
}

void toLL(const double* f_b, double* f_ll, const double* Cb_ll)
{
	f_ll[0] = Cb_ll[0] * f_b[0] + Cb_ll[1] * f_b[1] + Cb_ll[2] * f_b[2];
	f_ll[1] = Cb_ll[3] * f_b[0] + Cb_ll[4] * f_b[1] + Cb_ll[5] * f_b[2];
	f_ll[2] = Cb_ll[6] * f_b[0] + Cb_ll[7] * f_b[1] + Cb_ll[8] * f_b[2];
}

void mxMpy(const double* A, const double* B, double* C)
{
	C[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];		C[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];		C[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];
	C[3] = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];		C[4] = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];		C[5] = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];
	C[6] = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];		C[7] = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];		C[8] = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];
}

void updateRadii(double& R_fi, double& R_lam, const double& sin_fi)
{
	double denom = sqrt(1 - e2 * sin_fi * sin_fi);
	R_lam = a / denom; 
	R_fi = a * (1 - e2) / (denom * denom * denom); 
}

void calcPoissonEquation(const double* w_ig)
{
	double w_b_mat[9] = { 0, -w_b[2], w_b[1],
						w_b[2], 0, -w_b[0],
						-w_b[1], w_b[0], 0 };
	double w_ig_mat[9] = { 0, -w_ig[2], w_ig[1],
						w_ig[2], 0, -w_ig[0],
						-w_ig[1], w_ig[0], 0 };

	double firstSummand[9], secondSummand[9];
	mxMpy(Cb_ll, w_b_mat, firstSummand); 
	mxMpy(w_ig_mat, Cb_ll, secondSummand); 

	for (int i = 0; i < 9; ++i)
		Cb_ll[i] += (firstSummand[i] - secondSummand[i]) * d_time;
}

void calcRpy(const double* Cb_ll, double* rpy)
{
	double val = Cb_ll[6];
	val = fmax(-1.0, fmin(1.0, val));

	rpy[1] = atan2(Cb_ll[7], Cb_ll[8]) * rad2Deg; 	// pitch
	rpy[0] = -asin(val) * rad2Deg;					// roll
	rpy[2] = atan2(Cb_ll[3], Cb_ll[0]) * rad2Deg;	// yaw

	rpy[2] = std::fmod(rpy[2] + 360.0, 360) * -1 + 360;
}

void navigationAlgorithm(double* v_ll, double* rpy, double* coords)
{
	double f_ll[3] = { 0, 0, 0 };
	double w_ll[3] = { 0, 0, 0 };

	double sin_fi = sin(coords[0]);
	double cos_fi = cos(coords[0]);

	double R_fi, R_lam;
	updateRadii(R_fi, R_lam, sin_fi);

	double w_ie[3] = { 0, U * cos_fi, U * sin_fi };

	double w_eg[3];
	w_eg[0] = -v_ll[1] / (R_fi + coords[2]);
	w_eg[1] = v_ll[0] / (R_lam + coords[2]);
	w_eg[2] = v_ll[0] / (R_lam + coords[2]) * (sin_fi / cos_fi);

	double w_ig[3];
	for (int i = 0; i < 3; ++i)
		w_ig[i] = w_ie[i] + w_eg[i];

	calcPoissonEquation(w_ig);
	calcRpy(Cb_ll, rpy);

	toLL(f_b, f_ll, Cb_ll);

	double cor[3];
	cor[0] = (2 * w_ie[1] + w_eg[1]) * v_ll[2] - (2 * w_ie[2] + w_eg[2]) * v_ll[1];
	cor[1] = (2 * w_ie[2] + w_eg[2]) * v_ll[0] - (2 * w_ie[0] + w_eg[0]) * v_ll[2];
	cor[2] = (2 * w_ie[0] + w_eg[0]) * v_ll[1] - (2 * w_ie[1] + w_eg[1]) * v_ll[0];

	v_ll[0] += (f_ll[0] - cor[0]) * d_time;
	v_ll[1] += (f_ll[1] - cor[1]) * d_time;
	v_ll[2] += (f_ll[2] - 9.81 - cor[2]) * d_time;

	coords[0] += (v_ll[1] / (R_fi + coords[2])) * d_time;
	coords[1] += (v_ll[0] / ((R_lam + coords[2]) * cos_fi)) * d_time;

}

void saveData(GraphData& data, double time, double* rpy, double* v_ll, double* coords)
{
	data.time.push_back(time);
	data.roll.push_back(rpy[0]);
	data.pitch.push_back(rpy[1]);
	data.yaw.push_back(rpy[2]);
	data.v_e.push_back(v_ll[0]);
	data.v_n.push_back(v_ll[1]);
	data.fi.push_back(coords[0]);
	data.lam.push_back(coords[1]);
	data.h.push_back(coords[2]);
}

int main()
{
	//initGnuplot(gp);

	GraphData etalon, calc;
	
	std::ifstream in("Ref_Data_200Hz.txt");
	if (!in.is_open()) {
		std::cout << "Cannot open file Ref_Data_200Hz.txt" << std::endl;
		return -1;
	}
	
	int it = 0;
	while (!in.eof()) {
		double time;
		double v_ll[3] = { 0.0, 0.0, 0 }; // Ve, Vn, Vup
		double rpy[3] = { 0, 0, 0 }; // roll, pitch, yaw
		double coords[3] = { fi, lam, h };

		in >> time;
		in >> w_b[0] >> w_b[1] >> w_b[2] >> f_b[0] >> f_b[1] >> f_b[2];

		for (int i = 0; i < 3; ++i)
			w_b[i] *= deg2Rad;

		in >> rpy[0] >> rpy[1] >> rpy[2] >> v_ll[0] >> v_ll[1];
		in >> coords[0] >> coords[1] >> coords[2];

		saveData(etalon, time, rpy, v_ll, coords);

		navigationAlgorithm(v_ll, rpy, coords);

		saveData(calc, time, rpy, v_ll, coords);

		if (it++ > etalon.maxSize)
			break;
	}

	std::cout << "end";
	drawGraphs(etalon, calc);
}