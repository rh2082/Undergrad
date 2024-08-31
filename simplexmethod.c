//import i/o and other necessary libraries
#include <stdio.h> // include access to the printf statement and file I/O
#include <math.h> //include maths functions eg sqrt

//declare a coordinates structure which will be used for all the 2D points
typedef struct coords {
	double x0;
	double x1;
} coords;

//declare function F, being the Rosenbrock function
double F(double x0, double x1)
{
	return 100*(x1 - (x0*x0))*(x1 - (x0*x0))+ (1-x0)*(1-x0);
}

//this function expands the triangle
double expand(double pstarx, double pbarx)
{
	printf("Transformation = expand \n");
	double x;
	x = (2*pstarx) - pbarx;
	return x;
}

//this function contracts the triangle
double contract(double phx, double pbarx)
{
	printf("Transformation = contract \n");
	double x;
	x = ((phx + pbarx)/2);
	return x;
}

int minimum(double y0, double y1, double y2, double ybar)
{
	//this checks whether or not the minimum has been reached
	//it applies the condition provided in the instructions
	//then returns 1 if this condition is met
	//and 0 if it is not
	double t0, t1, t2, currentVal;
	t0 = ((y0 - ybar)*(y0-ybar))/2;
	t1 = ((y1 - ybar)*(y1-ybar))/2;
	t2 = ((y2 - ybar)*(y2-ybar))/2;
	currentVal = sqrt(t0 + t1 + t2);
	printf("Minimum test %e\n", currentVal);
	if (currentVal < 10e-8){
		return 1;
		}
	else{
		return 0;
		}
}

int downhillSimplex()
{
//initially p0 = (0,0), p1=(2,0) and p2=(0,2)
	//declare coordinate points, and find their respective y values
	coords p0, p1, p2;
	double y0, y1, y2;
	p0.x0 = 0.9;
	p0.x1 = 0.9;
	p1.x0 = 2;
	p1.x1 = 0;
	p2.x0 = 0;
	p2.x1 = 2;
	y0 = F(p0.x0, p0.x1);
	y1 = F(p1.x0, p1.x1);
	y2 = F(p2.x0, p2.x1);
	
	//print to check all is well
	printf("y0 = %e\n y1 = %e\n y2 = %e\n", y0, y1, y2);
	
	int count;
	int min;
	count = 1;
	min = 0;
	
	//sort variables
	double yl, ym, yh;
	coords pl, ph, pm;
	int h, m, l;
	
	//initialise pbar, p*
	coords pbar, pstar;
	double ybar, ystar;
	
	//initialise p** value
	coords p2star;
	double y2star;
	
	//set up a loop to count the iterations
	//and a flag to indicate if the minimum has been reached
	//these are the two conditions that could end the loop
	//if the minimum is reached then the loop will exit early
	//if it is not the loop will continue for 1000 iterations
	//start loop
	//count starts at 0 so limit of 999 is needed for 1000 iterations
	while (count < 20){
		//just to keep track of what count the loop is on
		printf("________________%d__________________\n", count);
		printf("initial values:\n");
		printf("p0 (%e,%e),\n p1 (%e,%e),\n p2(%e,%e)\n", p0.x0, p0.x1, p1.x0, p1.x1, p2.x0, p2.x1);

		//now need to order the y values by size
		//first need to find y values again
		y0 = F(p0.x0, p0.x1);
		y1 = F(p1.x0, p1.x1);
		y2 = F(p2.x0, p2.x1);
		
		//do so by sorting them
		//begin by taking one as the smallest arbitrarily
		
		yl = y0;
		l = 0;
		pl.x0=p0.x0;
		pl.x1=p0.x1;
		if (y1 < yl)
		{
			yl = y1;
			l = 1;
			pl.x0 = p1.x0;
			pl.x1 = p1.x1;
		}
		if (y2 < yl)
		{
			yl = y2;
			l = 2;
			pl.x0=p2.x0;
			pl.x1=p2.x1;
		}
		//the value yl holds is now the minimum y value
		yh = y0;
		h = 0;
		ph.x0=p0.x0;
		ph.x1=p0.x1;
		if (y1 > yh)
		{
			yh = y1;
			h = 1;
			ph.x0=p1.x0;
			ph.x1=p1.x1;
		}
		if (y2 > yh)
		{
			yh = y2;
			h = 2;
			ph.x0=p2.x0;
			ph.x1=p2.x1;
		}
		//yh now holds the maximum y value
		//just need to assign ym the remaining value
		if (yl == y0)
		{
			if (yh == y1)
			{
				ym = y2;
				m = 2;
				pm.x0 = p2.x0;
				pm.x1 = p2.x1;
			}
			else{
				ym = y1;
				m = 1;
				pm.x0 = p1.x0;
				pm.x1 = p1.x1;
			}
		}

		
	
		//print to check all is well
		printf("yl = %e\n ym = %e\n yh = %e\n", yl, ym, yh);
		//p0, p1, p2 and their respective y values are now in magnitude order
	
		//now that the p and y values are obtained
		//need to find pbar, p* and their y values
		
		//pbar is the midpoint between ym and yl so obtain this
		pbar.x0 = ((pm.x0 + pl.x0)/2);
		pbar.x1 = ((pm.x1 + pl.x1)/2);
		ybar = F(pbar.x0, pbar.x1);
		//p* is found by the equation given
		pstar.x0 = 2*(pbar.x0) - ph.x0;
		pstar.x1 = 2*(pbar.x1) - ph.x1;
		ystar = F(pstar.x0, pstar.x1);
		printf("p* = (%e,%e) \n", pstar.x0, pstar.x1);
		//print to check all is well
	
		printf("ybar = %e\n y* = %e\n ", ybar, ystar);
	
		//finally, all the initial values have been found
		//we can now move into the main flowchart for the downhill simplex
		
		//if y* < yl, we need to expand our triangle
		//determine p** in the expand function
		//then find y**
		if (ystar < yl){
			//y* is less than yl
			//determine p** for expansion case
			p2star.x0=expand(pstar.x0, pbar.x0);
			p2star.x1=expand(pstar.x1, pbar.x1);
			//determine y**
			y2star = F(p2star.x0, p2star.x1);
			//now we want to check how y** relates to yl
			//if y** < yl, replace ph with p**
			if (y2star < yl){
				ph.x0 = p2star.x0;
				ph.x1 = p2star.x1;
				printf("ph = p**\n");
				//we have now reached the bottom of the flowchart
				//need to check if minimum has been reached
				min = minimum(y0, y1, y2, ybar);
				}
			else{
				//as y** is not less than yl
				//we want to replace ph with p*
				ph.x0 = pstar.x0;
				ph.x1 = pstar.x1;
				printf("ph = p*\n");
				//we have now reached the bottom of the flowchart
				//need to check if minimum has been reached
				min = minimum(y0, y1, y2, ybar);
				}
			}
		//ystar is not less than yl
		else{
			//now we check how y* relates to ym
			if (ystar > ym){
				//as y* is greater than yl
				//we now have another condition
				//we want to check how y* relates to yh
				//if y* is not greater than yh we want to replace ph with p*
				if (ystar < yh){
					ph.x0 = pstar.x0;
					ph.x1 = pstar.x1;
					}
				if (ystar == yh){
					ph.x0 = pstar.x0;
					ph.x1 = pstar.x1;
					}
				//now we can move on as this branch is over
				//we need to find p** and y** in a contraction case
				p2star.x0 = contract(ph.x0, pbar.x0);
				p2star.x1 = contract(ph.x1, pbar.x1);
				y2star = F(p2star.x0, p2star.x1);
				//this has now happened
				//finally we check how y** relates to yh
				//if y** is greater than yh we need to replace all pi's
				if (y2star > yh){
					p0.x0 = (p0.x0+pl.x0)/2;
					p0.x1 = (p0.x1+pl.x1)/2;
					p1.x0 = (p1.x0+pl.x0)/2;
					p1.x1 = (p1.x1+pl.x1)/2;
					p2.x0 = (p2.x0+pl.x0)/2;
					p2.x1 = (p2.x1+pl.x1)/2;
					printf("replace pis\n");
					//we have now reached the bottom of the flowchart
					//need to check if minimum has been reached
					min = minimum(y0, y1, y2, ybar);
					}
				//if y** is not greater than yh we need to replace ph by p**
				else{
					ph.x0 = p2star.x0;
					ph.x1 = p2star.x1;
					printf("ph = p**\n");
					//we have now reached the bottom of the flowchart
					//need to check if minimum has been reached
					min = minimum(y0, y1, y2, ybar);
					}
				}
			else{
				//as y* is not greater than ym
				//we want to replace ph with p*
				ph.x0 = pstar.x0;
				ph.x1 = pstar.x1;
				printf("ph = p*\n");
				//we have now reached the bottom of the flowchart
				//need to check if minimum has been reached
				min = minimum(y0, y1, y2, ybar);
				}
			}
		
		//update coordinate values
		if (h == 0){
			p0.x0 = ph.x0;
			p0.x1 = ph.x1;
			}
		if (h == 1){
			p1.x0 = ph.x0;
			p1.x1 = ph.x1;
			}
		if (h == 2){
			p2.x0 = ph.x0;
			p2.x1 = ph.x1;
			}
		
		count = count + 1;
		if (min == 1){
			break;
			}
		
		
		
		
	
		printf("final values\n");
		printf("p0 = (%e,%e)\n", p0.x0, p0.x1);
		printf("p1 = (%e,%e)\n", p1.x0, p1.x1);
		printf("p2 = (%e,%e)\n", p2.x0, p2.x1);
		printf("ph = (%e, %e), pbar = (%e,%e)\n", ph.x0, ph.x1, pbar.x0, pbar.x1);
		printf("p* (%e,%e), p** (%e,%e)\n", pstar.x0, pstar.x1, p2star.x0, p2star.x1);
	}
	
			
	if (min == 1){
		printf("Minimum was reached in %d iterations \n", count);
		printf("Coordinates reached were:\n");
		printf("p0 = (%e,%e)\n", p0.x0, p0.x1);
		printf("p1 = (%e,%e)\n", p1.x0, p1.x1);
		printf("p2 = (%e,%e)\n", p2.x0, p2.x1);
		}
	else{
		printf("Minimum was not reached within 1000 iterations\n");
		}
}

//declare main function, where program will default to starting
//this is where the if statements will mostly reside
int main()
{
	print("i am running")
	downhillSimplex();
}