/*
Author: Aoun Hussain
Class: ECE 6122-A
Last Date Modified: Fri, Sep 24, 2021
Description: This is the cpp source file for the solution of problem 3, lab 1.
             It takes a single command line argument that represents the initial
             reflection’s x location along the AB segment and calculates the
             number of times the beam is reflected off an internal surface of the
             white cell before exiting.
*/

/*
 References:
 https://www.geeksforgeeks.org/program-for-point-of-intersection-of-two-lines/
 https://www.wikihow.com/Find-the-Equation-of-a-Perpendicular-Line
*/

#include <cctype>
#include <algorithm>
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
using namespace std;

bool is_numeric(const string &strIn, long double &xLoc)
{
    /*
    Checks if the argument is valid integer or float and returns true/false accordingly
    */

    bool bRC = false;     // flags for valid number, integer, and floats
    bool fFlt = false;
    bool fNeg = false;

    for (int i = 0; i < strIn.length(); i ++)
    {
        if ((i == 0) && (strIn.at(i) == '-'))
        {
            bRC = true;
            fNeg = true;
        }
        else if ((i == 1) && (strIn.at(i) == '.'))
        {
            bRC = true;
            fFlt = true;
        }
        else if ((i == 2) && (strIn.at(i) == '.') && (fNeg))
        {
            bRC = true;
            fFlt = true;
        }
        else if (::isdigit(strIn.at(i)))
        {
            bRC = true;
        }
        else
        {
            bRC = false;
            break;
        }
    }

    if ((bRC) and (fFlt))
    {
        xLoc = stof(strIn);
        return true;
    }
    else if (bRC)
    {
        xLoc = stoi(strIn);
        return true;
    }
    else
    {
        return false;
    }

}

class ECEPoint
{
    /*
     This class ECEPoint is a structure class
     to store x and y coordinates of the point
     of a line segment and calculate distance
     between two points.
    */

public:
    long double m_x = 0.0;
    long double m_y = 0.0;

    double distance(const ECEPoint &inPt)
    {
        /*
         function to calculate distance between two points
        */
        return sqrt((m_x - inPt.m_x) * (m_x - inPt.m_x) + (m_y - inPt.m_y) * (m_y - inPt.m_y));
    }
};

class ECELineSegment
{
    /*
    This class ECELineSegment is a structure class of a line segment
     which helps us to store three sides of the triangular cell
     and perform functions such as intersection and reflection of arrays,
     associated with the line segments of the cell.
    */

public:
    ECEPoint m_pt1, m_pt2;
    long double m_length = 0.0;
    long double m_normalX = 0.0;
    long double m_normalY = 0.0;
    int num;

    ECELineSegment() {}

    ECELineSegment(const ECEPoint  &pt1, const ECEPoint &pt2) : m_pt1(pt1), m_pt2(pt2)
    {
        /*
        function to assign points to the line segment - overloaded constructor
        */

        m_length = m_pt1.distance(m_pt2);
    }

    int getNumber()
    {
        /*
        function to return the side number
        */

        return num;
    }

    void setNumber(int n)
    {
        /*
        function to set the side number
        */

        num = n;
    }

    /*
    functions to return the vector from two points of  a line segment in x and y direction
    */

    long double getDeltaX() const { return (m_pt2.m_x - m_pt1.m_x); }
    long double getDeltaY() const { return (m_pt2.m_y - m_pt1.m_y); }

    void setNormal(long double nX, long double nY)
    {
        /*
        function to calculate the normal of the line segments directed inside the cell
        */

        long double length = sqrt(nX * nX + nY * nY);
        if (length <= 0.0)
        {
            length = 1.0;
        }
        m_normalX = nX / length;
        m_normalY = nY / length;
    }

    /*
    functions to set/change the points of the line segments
    */

    void setPoint1(const ECEPoint &pt) { m_pt1 = pt; }
    void setPoint2(const ECEPoint &pt) { m_pt2 = pt; }

    bool getReflectionPt(const ECELineSegment &inLine, ECEPoint &reflectPt)
    {
        /*
        function to return the intersection point between two line segments
        */
        // https://www.geeksforgeeks.org/program-for-point-of-intersection-of-two-lines/

        // Line AB represented as a1x + b1y = c1
        long double a1 = m_pt2.m_y - m_pt1.m_y;
        long double b1 = m_pt1.m_x - m_pt2.m_x;
        long double c1 = a1 * (m_pt1.m_x) + b1 * (m_pt1.m_y);

        // Line CD represented as a2x + b2y = c2
        long double a2 = inLine.m_pt2.m_y - inLine.m_pt1.m_y;
        long double b2 = inLine.m_pt1.m_x - inLine.m_pt2.m_x;
        long double c2 = a2 * (inLine.m_pt1.m_x) + b2 * (inLine.m_pt1.m_y);

        long double determinant = a1 * b2 - a2 * b1;

        if (determinant == 0)
        {
            // The lines are parallel.
            return false;
        }
        else
        {
            // on intersection
            reflectPt.m_x = (b2 * c1 - b1 * c2)/determinant;
            reflectPt.m_y = (a1 * c2 - a2 * c1)/determinant;

            return true;
        }
    }

    bool getReflectedLine(const ECELineSegment &inLine, ECELineSegment &reflectLine, ECEPoint &reflectPt)
    {
        /*
        function to return the reflected array
        */

        if (getReflectionPt(inLine, reflectPt))
        {
            // Check if point is inside of the segment
            if (reflectPt.distance(m_pt1) > m_length || reflectPt.distance(m_pt2) > m_length)
            {
                return false;
            }
            // The point is inside the segment so check if the y value is less than the cutout section
            if (reflectPt.m_y <= 0.01)
            {
                return false;
            }
            else // Determine the reflected ray
            {
                // r = d - 2(d * n)n;
                // d is the inLine
                // n is the normal vector of this side
                // r is the reflectLine

                long double vectorX = 100 * (inLine.getDeltaX() - 2 * ((inLine.getDeltaX() * m_normalX) + (inLine.getDeltaY() * m_normalY)) * m_normalX);
                long double vectorY = 100 * (inLine.getDeltaY() - 2 * ((inLine.getDeltaX() * m_normalX) + (inLine.getDeltaY() * m_normalY)) * m_normalY);

                reflectLine.m_pt1.m_x = reflectPt.m_x;
                reflectLine.m_pt1.m_y = reflectPt.m_y;

                reflectLine.m_pt2.m_x = reflectPt.m_x + vectorX;
                reflectLine.m_pt2.m_y = reflectPt.m_y + vectorY;

                return true;
            }
        }
        else
        {
            return false;
        }
    }

};


void refCalc(long double &xLoc, unsigned long &refNum)
{
    /*
     Calculates the number of reflections until the line segment exits the triangular cell
     by calling the functions of the above declared classes
     A(-10, root300), B(10, root300), C(0,0)
    */

    long double yLoc = sqrt(20.0 * 20.0 - 10.0 * 10.0);
    bool bStillInside = false;

    ECELineSegment sides[3]; // These are the three side of the triangle

    // Initialize the 3 sides

    ECEPoint A = ECEPoint();
    A.m_x = -10;
    A.m_y = sqrt(300);

    ECEPoint B = ECEPoint();
    B.m_x = 10;
    B.m_y = sqrt(300);

    ECEPoint C = ECEPoint();

    // side 0 is A-B
    sides[0] = {A, B};
    sides[0].setNormal(sides[0].getDeltaX(), sides[0].getDeltaY());
    sides[0].setNumber(0);

    // side 1 is B-C
    sides[1] = {B, C};
    sides[1].setNormal(sides[1].getDeltaX(), sides[1].getDeltaY());
    sides[1].setNumber(1);

    // side 2 is C-A
    sides[2] = {C, A};
    sides[2].setNormal(sides[2].getDeltaX(), sides[2].getDeltaY());
    sides[2].setNumber(2);

    ECELineSegment incidentRay = { { 0.0, 0.0 }, { xLoc, yLoc } };
    ECELineSegment reflectedRay;
    ECEPoint intersectPt = { xLoc, yLoc };

    unsigned long nNumReflects{ 0 };
    int nReflectedSide = -1; // Top segment is the first that is 0 but for now its an invalid value

    //cout << 0.0 << "," << 0.0 << std::endl;
    do
    {
        bStillInside = false;
        for (int ii = 0; ii < 3; ii++)
        {
            if (nReflectedSide < 0 || ii != nReflectedSide)
            {
                if (sides[ii].getReflectedLine(incidentRay, reflectedRay, intersectPt))
                {
                    // update the values and swap the incident and reflected rays
                    //cout << intersectPt.m_x << "," << intersectPt.m_y << std::endl;
                    nReflectedSide = sides[ii].getNumber();  //updates the last side reflected
                    incidentRay = {reflectedRay.m_pt1, reflectedRay.m_pt2};   //swaps the incident and reflected ray
                    nNumReflects++;
                    bStillInside = true;
                }
            }
        }

    } while (bStillInside);

    //cout << intersectPt.m_x << "," << intersectPt.m_y << std::endl;

    //std::cout << nNumReflects << std::endl;

    // assigning final value to output
    refNum = nNumReflects;
}


int main(int argc, char* argv[]){
    /*
    main function to call all functions and check for invalid inputs
     and write output to file
    */

    if (argc == 2)
    {
        long double xLoc{ 0 };
        unsigned long refNum{ 0 };
        string strInput(argv[1]);
        bool bIsValid = is_numeric(strInput, xLoc);
        if ((bIsValid) && (xLoc > -10) && (xLoc < 10))
        {
            // perform function with xLoc

            refCalc(xLoc, refNum);

            ofstream myfile ("output3.txt");
            if (myfile.is_open())
            {
                myfile << to_string(refNum) << "\n";
                myfile.close();
            }
            else cout << "Unable to open file";
        }
        else
        {
            cout << "Invalid Input" <<endl;
            ofstream myfile ("output3.txt");
            if (myfile.is_open())
            {
                myfile << "Invalid Input" << "\n";
                myfile.close();
            }
            else cout << "Unable to open file";
        }


    }
    else
    {
        cout << "Invalid Input" <<endl;
        ofstream myfile ("output3.txt");
        if (myfile.is_open())
        {
            myfile << "Invalid Input" << "\n";
            myfile.close();
        }
        else cout << "Unable to open file";
    }

    return 0;
}


