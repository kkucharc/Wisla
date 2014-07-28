#include<stdio.h>
#include<ncurses.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include <xmmintrin.h> 

char ch, file_name[25];

float wisla10Dat[36][90];
float sola10Dat[36][90];
float przemsza10Dat[36][90];

int i, rok, j;
int aproach;
char line[1000];
float test;
int k, x, u, w;
float ro, epsilon;

int startStateDivision;
int howManyInflows;
int howManyTimeSteps;
int howManyControls;

int *elementsArray;

typedef struct control {
	float rG1;
	float rG2;
	float rT1;
	float rT2;
	float rT3;
	float r3;
} control;

typedef struct state {
	float SG;
	float ST;
} state;

typedef struct inflows {
	float iG;
	float iT;
	float i1;
} inflows;

state *statesCombinations;
control *controlCombinations;
inflows **disturbanceCombinations;
float **Jprev, **Jactual;
control **Policy;

float rand_FloatRange(float a, float b);
void wczytaj(char * str, float *tab);
void wczytajHist(char * str, int typ);
float maxy(float num1, float num2);
void generateStartStates(int howMany);
void generateControls(int rG1div, int rG2div, int rT1div, int rT2, int rT3div,
		int r3div);
void generateInflows(int howMany, int timeStep);
float g(state x, control u, inflows w);
void prepareJ(int howManyTimeSteps, int startStateDivision);
void addJ(int nrOfstateToAdd, int timeStep, float bestQ);
void preparePolicy(int howManyTimeSteps, int startStateDivision);
float norm(float **matrix, float **odejmowana, int x, int y);
float J(int timeStep, state x, control u, inflows w);
float koszt(float rG2, float ik1, float r1, float r2, float r3, float rG1,
		float rT2);
int convertNewFloatStateTo(state s);
void readMatrix(float **matrix, int x, int y);

void main(int argc, char **argv) {
	int rG1division, rG2division, rT1division, rT2division, rT3division,
			r3division;
	srand (time(NULL));wczytajHist
	("WISLA10.DAT", 1); //goczalkowice
	wczytajHist("SOLA10.DAT", 2); //sola
	wczytajHist("PRZEMSZA.DAT", 3); //przemsza

	aproach = 1;
	ro = 1;
	epsilon = 430;
	howManyInflows = 3; //wielokrotnosc 3
	howManyTimeSteps = 36;

	rG1division = 7;
	rG2division = 5;
	rT1division = 7;
	rT2division = 5;
	rT3division = 6;
	r3division = 7;
	startStateDivision = 7;
	howManyControls = rG1division * rG2division * rT1division * rT2division
			* rT3division * r3division;

	preparePolicy(howManyTimeSteps, startStateDivision);
	prepareJ(howManyTimeSteps, startStateDivision);
	generateStartStates(startStateDivision);
	generateControls(rG1division, rG2division, rT1division, rT2division,
			rT3division, r3division);
	generateInflows(howManyInflows, howManyTimeSteps);

	while (ro > 0.01) {
		for (k = 0; k < howManyTimeSteps; k++) {
			for (x = 0; x < startStateDivision * startStateDivision; x++) {
				state startState, startStateAtTheBeginig;
				float best_Q;
				control bestControl;
				startStateAtTheBeginig = startState = statesCombinations[x];
				best_Q = 999999999.9;
				for (u = 0; u < howManyControls; u++) {
					float Q;
					control startControl;
					startControl = controlCombinations[u];
					Q = 0;

					for (w = 0; w < howManyInflows; w++) {
						inflows startDisturbance;
						float g_J;
						startDisturbance = disturbanceCombinations[k][w];
						switch (aproach) {
						case 1:	//minmax

							g_J = maxy(
									g(startState, startControl,
											startDisturbance),
									J(k, startState, startControl,
											startDisturbance));
							Q = maxy(Q, g_J);

							break;

						case 2: // discount
							g_J = g(startState, startControl, startDisturbance)
									+ (float) 0.1
											* J(k, startState, startControl,
													startDisturbance);
							Q = Q + g_J;
							break;
						default:

							break;
						}

					}

					if (Q < best_Q) {
						best_Q = Q;
						bestControl = startControl;

					}

				}
				//dodanie polityki
				Policy[k][x] = bestControl;
				
				// dodanie do macierzy J
				addJ(convertNewFloatStateTo(statesCombinations[x]), k, best_Q); // x instead of startState

			}
		}
		// r = || J - Jprev ||
		ro = norm(Jactual, Jprev, startStateDivision * startStateDivision,
				howManyTimeSteps);

		// Jprev = J
		for (i = 0; i < howManyTimeSteps; i++)
			for (j = 0; j < startStateDivision * startStateDivision; j++)
				Jprev[i][j] = Jactual[i][j];
	}

	getch();
}

void prepareJ(int howManyTimeSteps, int startStateDivision) {
	int a, b;
	Jprev = (float**) malloc(howManyTimeSteps * sizeof(float*));
	for (a = 0; a < howManyTimeSteps; a++) {
		Jprev[a] = (float*) malloc(
				startStateDivision * startStateDivision * sizeof(float));
	}

	for (a = 0; a < howManyTimeSteps; a++)
		for (b = 0; b < startStateDivision * startStateDivision; b++)
			Jprev[a][b] = 0;

	Jactual = (float**) malloc(howManyTimeSteps * sizeof(float*));
	for (a = 0; a < howManyTimeSteps; a++) {
		Jactual[a] = (float*) malloc(
				startStateDivision * startStateDivision * sizeof(float));
	}

	for (a = 0; a < howManyTimeSteps; a++)
		for (b = 0; b < startStateDivision * startStateDivision; b++)
			Jactual[a][b] = 0;
}

void preparePolicy(int howManyTimeSteps, int startStateDivision) {
	int a;
	Policy = (control**) malloc(howManyTimeSteps * sizeof(control*));
	for (a = 0; a < howManyTimeSteps; a++) {
		Policy[a] = (control*) malloc(
				startStateDivision * startStateDivision * sizeof(control));
	}

}

void addJ(int nrOfstateToAdd, int timeStep, float bestQ) {
	Jactual[timeStep][nrOfstateToAdd] = bestQ;
}

float g(state x, control u, inflows w) {
	return koszt(u.rG2, w.i1, u.rG2 + w.i1, u.rT1 - u.r3, u.r3, u.rG1, u.rT2);
}

float J(int timeStep, state x, control u, inflows w) {
	state newState;

	newState.SG = x.SG + (w.iG - u.rG1 - u.rG2 + u.rT3) * (10 * 24 * 60 * 60);
	newState.ST = x.ST + (w.iT - u.rT1 - u.rT2 - u.rT3) * (10 * 24 * 60 * 60);
	if (newState.SG < 0 || newState.ST < 0)
		return 0;
	else if (convertNewFloatStateTo(newState) > 48)
		return 0.0; //przepelnienie

	return Jprev[timeStep][convertNewFloatStateTo(newState)];

}

int convertNewFloatStateTo(state s) { // pass division as parameter- to do
	int foundState, x, y;
	x = (int) ((s.SG / 1000000) - 20) / (140 / 7);
	y = (int) ((s.ST / 1000000) - 5) / (125 / 7);
	foundState = x * 7 + y;
	if (foundState > 48) {
		return 49;
	}

	return foundState;
}

void generateStartStates(int howMany) {
	int x1;
	int x2;
	int pom;
	pom = 0;
	statesCombinations = (state *) malloc(howMany * howMany * sizeof(state));
	for (x1 = 0; x1 < howMany; x1++) {
		for (x2 = 0; x2 < howMany; x2++) {
			statesCombinations[pom].SG = (float) ((20 + (140 / howMany) * x1)
					* 1000000); //10days interval
			statesCombinations[pom].ST = (float) ((5 + (125 / howMany) * x2)
					* 1000000);
			pom++;
		}
	}
}

void generateInflows(int howMany, int timeStep) {
	int t, x;
	int a, b, p, k;
	inflows inf;
	disturbanceCombinations = (inflows **) malloc(timeStep * sizeof(inflows*));
	for (a = 0; a < timeStep; a++)
		disturbanceCombinations[a] = (inflows *) malloc(
				howMany * sizeof(inflows));

	t = p = howMany / 3;
	if (t > 19)
		p = 20; //bo mamy tylko dla 20 lat 

	for (k = 0; k < timeStep; k++) {
		x = 0;
		for (a = 0; a < t; a++)
			for (b = 0; b < t; b++)
				for (p = 0; p < t; p++) {

					inf.iT = wisla10Dat[k][a];
					inf.iG = sola10Dat[k][b];
					inf.iT = przemsza10Dat[k][p];

					disturbanceCombinations[k][x] = inf;
					x++;
				}
	}
}

void generateControls(int rG1div, int rG2div, int rT1div, int rT2div,
		int rT3div, int r3div) {
	int rG1div0, rG2div0, rT1div0, rT2div0, rT3div0, r3div0;
	int pom = 0;
	controlCombinations = (control *) malloc(
			rG1div * rG2div * rT1div * rT2div * rT3div * r3div
					* sizeof(control));

	for (rG1div0 = 0; rG1div0 < rG1div; rG1div0++)
		for (rG2div0 = 0; rG2div0 < rG2div; rG2div0++)
			for (rT1div0 = 0; rT1div0 < rT1div; rT1div0++)
				for (rT2div0 = 0; rT2div0 < rT2div; rT2div0++)
					for (rT3div0 = 0; rT3div0 < rT3div; rT3div0++)
						for (r3div0 = 0; r3div0 < r3div; r3div0++) {
							controlCombinations[pom].rG1 = (float) (0
									+ (5 / rG1div) * rG1div0);
							controlCombinations[pom].rG2 = (float) (0
									+ (22 / rG2div) * rG2div0);
							controlCombinations[pom].rT1 = (float) (0
									+ (335 / rT1div) * rT1div0);
							controlCombinations[pom].rT2 = (float) (0
									+ (20 / rT2div) * rT2div0);
							controlCombinations[pom].rT3 = (float) (0
									+ (3.9 / rT3div) * rT3div0);
							controlCombinations[pom].r3 = (float) (2.5
									+ (9.5 / r3div) * r3div0);

							pom++;
						}

}

void wczytajHist(char * str, int typ) {
	FILE *fp;

	fp = fopen(str, "r"); // read mode 
	if (fp == NULL) {
		perror("Error while opening the file GOCZAL10.DEM \n");
		exit (EXIT_FAILURE);
	}
	//skip first line
	while (fgets(line, sizeof line, fp) != NULL) /* read a line from a file */{
		fprintf(stdout, "%s", line); //print the file contents on stdout.
		break;
	}

	if (typ == 1) {
		for (rok = 0; rok < 90; rok++) {
			fscanf(fp, "%d", &test); //to remove year

			for (i = 0; i < 36; i++) {
				fscanf(fp, "%f", &wisla10Dat[i][rok]);
			}
		}

	}

	else if (typ == 2) {

		for (rok = 0; rok < 90; rok++) {
			fscanf(fp, "%d", &test); //to remove year

			for (i = 0; i < 36; i++) {
				fscanf(fp, "%f", &sola10Dat[i][rok]);

			}
		}

	}

	else {
		for (rok = 0; rok < 20; rok++) {
			for (i = 0; i < 36; i++) {
				fscanf(fp, "%f", &przemsza10Dat[i][rok]);
			}
		}
	}

	fclose(fp);
}

float koszt(float rG2, float ik1, float r1, float r2, float r3, float rG1,
		float rT2) {

	return maxy(0, 20 - rG2 - ik1) + maxy(0, 5 - r2) + maxy(0, 10 - r3)
			+ maxy(0, 30 - r1 - r2) + maxy(0, 5 - rG1) + maxy(0, 8 - rT2);

}

float maxy(float num1, float num2) {
	if (num1 > num2)
		return num1;
	return num2;

}

void readMatrix(float **matrix, int y, int x) {
	int i, j = 0;
	printf("\nCzytam macierz:\n");
	for (i = 0; i < x; i++) {
		for (j = 0; j < y; j++) {
			printf("%.2f ", matrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

float norm(float **matrix, float **odejmowana, int y, int x) {
	float sum;
	float **a;
	int i, j;
	sum = 0;
	a = (float**) malloc(x * sizeof(float*));

	for (i = 0; i < x; i++) {
		a[i] = (float*) malloc(y * sizeof(float));
	}

	//odejmowanie
	for (i = 0; i < x; i++)
		for (j = 0; j < y; j++) {
			a[i][j] = matrix[i][j] - odejmowana[i][j];
		}

	//potega
	for (i = 0; i < x; i++)
		for (j = 0; j < y; j++) {
			sum += (a[i][j] * a[i][j]);
		}

	sum = sqrt(sum);

	return sum;
}
