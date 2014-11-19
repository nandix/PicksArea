/**
 * Program: Pick's Theorem
 *
 * Description:
 * 			This program uses computational geometry methods to implement
 * 			pick's theorem. It computes the area of a simple convex polygon 
 * 			by counting the number of integer (x,y) points ("Lattice Points")
 * 			lie inside a polygon and lie on the edge of a polygon using the
 * 			formula:
 *
 * 						Area(P) = I(P) + B(P)/2 − 1
 *
 *			where I(P) is the number of points inside and B(P) is the number
 *			of points lying on an edge.
 *
 * 			Some of the code used here was originally used in a programming 
 * 			team practice on 10/22/2014.
 */

#include <fstream>
#include <vector>
#include <iostream>
#include <iomanip>
#include <map>
#include <string>

using namespace std;

typedef long int LI;
typedef long double LD;

ofstream debug;
const bool DEBUGGING = false;

struct point
{
	LD x,y;
};

struct vect
{
	LD i,j;	
};

struct poly{
	vector<point> polygon;
};

void print_point(point p, ofstream &fout ){

	fout << "( " << p.x << ", " << p.y << " ) ";

	return;
}

bool colinear(point p1, point p2, point p3)
{
	double result;
	
	result = p1.x * (p2.y - p3.y);
	result += p2.x *( p3.y - p1.y);
	result += p3.x * (p1.y - p2.y);

	if (result == 0)
		return true;
	else
		return false; 
}


bool straddle (point p1, point p2, point p3, point p4)
{
	vect line, v, w;
	double cross1;
	double cross2;

	line.i = p2.x - p1.x;
	line.j = p2.y - p1.y;

	v.i = p3.x - p2.x;
	v.j = p3.y - p2.y;

	w.i = p4.x - p2.x;
	w.j = p4.y - p2.y;

	cross1 = line.i * v.j - line.j * v.i;
	cross2 = line.i * w.j - line.j * w.i;

	if(cross1 > 0 && cross2 > 0)
		return false;

	if(cross1 < 0 && cross2 < 0)
		return false;

	if( cross1 == 0 && cross2 == 0)
	{
		if(!colinear (p2, p3, p4))
			return false;
	}
	return true;
}

bool intersect(point p1, point p2, point p3, point p4)
{
	if(!straddle (p1, p2, p3, p4))
		return false;

	if(!straddle (p3, p4, p1, p2))
		return false;

	return true;
}

bool inside_poly(vector<point> points, point checkPoint)
{
	point op;
	op.x = 1000;
	op.y = 1000;
	LI ints = 0;

	for(int i = 0 ; i < points.size() ;i++)
	{
		int iwrap = (i+1)%points.size();

		// If the point is colinear with points[i] and points[iwrap]
		// then it is a boundary point, and don't count it.
		if( colinear( checkPoint, points[i], points[iwrap] ))
			return false;

		// If the line of intersection is colinear with a vertex,
		// add 1 to the outside point x-coord.
		if(colinear(op, checkPoint, points[i]))
		{	
			op.x += 1;
		}
		if(colinear(op, checkPoint, points[iwrap]))
		{	
			op.x += 1;
		}

		if(intersect(op, checkPoint, points[i], points[iwrap]))
			ints++;
	}


	if(ints % 2 == 0)
		return false;
	else
		return true;
}


// Returns true if point bet is colinear between p1 and p2
bool between_colinear(point p1, point p2, point bet){


	// Calculate the slope from p1 to bet and p2 to bet
	// If the point is between and colinear, the slopes will be identical
	// with opposite sign
	// DOES NOT TAKE CARE OF 0 SLOPE CASE
	
	// print_point(p1, debug);
	// print_point(p2, debug);
	// print_point(bet, debug);
	// debug << endl;
	bool ret_val;

	int num_slope1 = (bet.y - p1.y);
	int denom_slope1 = (bet.x - p1.x);
	int num_slope2 = (bet.y - p2.y);
	int denom_slope2 = (bet.x - p2.x);

	// So we don't double count verticies, see if our point is the first 
	// vertex passed.
	if( denom_slope1 == 0 && num_slope1 == 0 ){
		
		if( DEBUGGING ){
			debug << "Vertex: ";
			print_point(bet, debug);
			debug << endl;
		}	

		return true;
	}

	// So we don't accidentally double count verticies, make sure the 
	// second point (p2) is not the same point as bet.
	if( !(denom_slope2 == 0 && num_slope2 == 0)){
		// If the numerator is zero, see if the point lies between the two
		if( num_slope1 == 0 && num_slope2 == 0 ){
			ret_val = !((denom_slope1 > 0) == (denom_slope2 > 0 ));
			if( DEBUGGING && ret_val ){
				debug << "Zero Numerator: ";
				print_point(bet, debug);
				debug << endl;
			}
			return ret_val;
		}
		// If the denominator is zero, see if the point lies between the two
		if( denom_slope1 == 0 && denom_slope2 == 0 ){
			ret_val = !((num_slope1 > 0) == (num_slope2 > 0 ));

			if( DEBUGGING && ret_val ){
				debug << "Zero Numerator: ";
				print_point(bet, debug);
				debug << endl;
			}

			return ret_val;
		}
	}


	// An attempt to use slopes to find colinear. Taking another apporach
	// double slope1 = num_slope1/denom_slope1;
	// double slope2 = num_slope2/denom_slope2;

	// if( slope2 == -slope1 ){
	// 	return true;
	// }
	// 
	// To determine if a point (without a zero rise or run) is on the line,
	// make a bounding box and see if the point is inside the box. If the 
	// point lies inside the box, check if it is colinear.
	vector<point> bound;
	point temp;

	// Create the bounding box
	bound.push_back(p1);

	temp.x = p1.x;
	temp.y = p2.y;
	bound.push_back(temp);

	bound.push_back(p2);

	temp.x = p2.x;
	temp.y = p1.y;
	bound.push_back(temp);

	if( inside_poly(bound, bet) ){
		if( colinear(p1, p2, bet)){
			debug << "Inside & colin: ";
			print_point(bet, debug);
			debug << endl;
			return true;
		}
	}

	return false;

}

double area( vector<point> p ){
	double result = 0;

	for( int i = 0; i < p.size(); i++ ){
		int j = (i + 1) % p.size();
		result += p[i].x * p[j].y;
		result -= p[i].y * p[j].x;
	}

	return result/2;
}

// Function to count the number of points inside the shape
vector<point> I(LI min_x, LI min_y, LI max_x, LI max_y, poly p, LI &count){

	vector<point> inside_points;
	count = 0;
	point temp;

	// Check if each point inside the range is inside the shape
	for( LI x = min_x; x <= max_x; x++ )
		for( LI y = min_y; y <= max_y; y++ ){

			temp.x = x;
			temp.y = y;

			if( inside_poly(p.polygon, temp) ){
				count ++;
				inside_points.push_back(temp);
			}

		}

	return inside_points;
}

vector<point> B(LI min_x, LI min_y, LI max_x, LI max_y, poly p, LI &count){
	
	count = 0;
	point temp;

	vector<point> boundary_points;

	// Check if each point inside the range lies on an edge
	for( LI x = min_x; x <= max_x; x++ )
		for( LI y = min_y; y <= max_y; y++ ){

			temp.x = x;
			temp.y = y;

			for( int i = 0; i < p.polygon.size(); i++ )
				if( between_colinear(p.polygon[i], p.polygon[(i+1)%p.polygon.size()], temp) ){
					count ++;
					boundary_points.push_back(temp);
				}

		}

	return boundary_points;
}


int main(int argc, char** argv)
{

	debug.open("DEBUG.txt");

	ifstream fin;
	ofstream fout;
	poly shape;
	point temp;

	vector<point> boundary_points;
	vector<point> inside_points;

	LI min_x, min_y, max_x, max_y;
	LI n_inside, n_boundary;

	min_x = min_y = -100;
	max_x = max_y = 100;

	if( argc != 2 ){
		cout << "Usage: ./picks <file_prefix>" << endl;
		return 1;
	}

	string in = argv[1];
	in += ".in";

	string out = argv[1];
	out += ".out";

	fin.open(in.c_str());
	fout.open(out.c_str());

	if(!fin){
		cout << "Could not open input file " << in << endl;
		return 1;
	}
	else
		fout << "Reading input from file " << in << endl;


	// Read in all the points in the polygon
	while(fin >> temp.x >> temp.y){
		shape.polygon.push_back(temp);
	}


	
	inside_points = I(min_x, min_y, max_x, max_y, shape, n_inside);
	boundary_points = B(min_x, min_y, max_x, max_y, shape, n_boundary);

	fout << "Area = " << n_inside + n_boundary/2.0 - 1 << endl;
	fout << "I(P) = " << n_inside << endl;
	fout << "B(P) = " << n_boundary << endl;
	fout << "Interior points" << endl;


	for(int i = 0; i < inside_points.size(); i++ ){
		fout << "     " << inside_points[i].x << " " << inside_points[i].y << endl;
	}
	
	fout << "Boundary points" << endl;
	for(int i = 0; i < boundary_points.size(); i++ ){
		fout << "     " << boundary_points[i].x << " " << boundary_points[i].y << endl;
	}
	




	fin.close();
	fout.close();

	return 0;
}
