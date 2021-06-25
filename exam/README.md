Examination project: Akima sub-spline

Author: Jeppe Thøis Madsen, au648224

Student no: 201904744
44 mod 22 = 0 -> Akima sub-spline

PROBLEM
Implement the Akima sub-spline. See the section "Akima sub-spline interpolation" in the book for the inspiration.

DESCRIPTION OF SOLUTION
This examination project interpolates a number of 2-D-points with the Akima sub-spline and a simpler cubic spline and plots the two.
The Akima sub-spline is described in the course material. It prioritizes continuity of the second derivative at the expense of maximal differentiability. This is useful for reducing wiggling.

The project contains cspline which is the algorithm of a regular cubic spline, and aspline which is the Akima sub-spline -- the actual project.
They are plotted in plot.png which also shows the points restricting the splines and the function which generated these.
The points are chosen to be as similar as possible to Figure 1.2 in D.V.Fedorov's "Introduction to numerical methods".
