#include "object.hpp"

#include <cmath>
#include <cfloat>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <iostream>


bool Object::intersect(Ray ray, Intersection &hit) const 
{
    // Assert the correct values of the W coords of the origin and direction.
    // You can comment this out if you take great care to construct the rays
    // properly.
    ray.origin[3] = 1;
    ray.direction[3] = 0;

    Ray local_ray(i_transform * ray.origin, i_transform * ray.direction);
	//!!! USEFUL NOTES: to calculate depth in localIntersect(), if the intersection happens at
	//ray.origin + ray.direction * t, then t is just the depth
	//!!! USEFUL NOTES: Here direction might be scaled, so you mustn't renormalize it in
	//localIntersect(), or you will get a depth in local coordinate system,
	//which can't be compared with intersections with other objects
    if (localIntersect(local_ray, hit)) 
	{
        // Assert correct values of W.
        hit.position[3] = 1;
        hit.normal[3] = 0;
        
		// Transform intersection coordinates into global coordinates.
        hit.position = transform * hit.position;
        hit.normal = (n_transform * hit.normal).normalized();
        
		return true;
    }

    return false;
}


bool Sphere::localIntersect(Ray const &ray, Intersection &hit) const 
{
    //////////////////
    // YOUR CODE HERE 
	// For this part you are required to calculate if a ray has intersected your sphere.
	// Be sure to cover all the possible intersection scenarios(zero, one, and two points of intersection).
	// Test your result by comparing the output depth image of your algorithm with the provided example solution's results.

	// Here in local coordinate system, the sphere is centered on (0, 0, 0)
	//with radius 1.0
	//
	// NOTE: hit.depth is the current closest intersection depth, so don't
	// accept any intersection that happens further away than that.
	   // just like in the case of plane

   // solving the quadratic equation
    double radius = this->radius;
    Vector v = Vector(0,0,0) - ray.origin;
    double a = ray.direction.dot(ray.direction);
    double b = 2 * ray.direction.dot(ray.origin);
    double c = v.dot(v) - radius * radius;
    double discr = b * b - 4 * a * c;  
       // now check for the conditions
    if (discr < 0) { // if what is under the root is negative, return false
        return false;
     }
   else if (discr>0){  // 
        double d1 = (-b+sqrt(discr))/(2*a);
        double d2 =(-b-sqrt(discr))/(2*a);
        if(d1>0 && d2>0 && d2<hit.depth){ // grab the near
            hit.depth=d2;
            hit.normal=(ray.origin+d2*ray.direction).normalized();
            hit.position=ray.origin+d2*ray.direction;
            return true;
        }
        else if (d1>0 && d2<0 && d1<hit.depth){
            hit.depth=d1;
            hit.normal=(ray.origin+d1*ray.direction).normalized();
            hit.position=ray.origin+d1*ray.direction;
            return false;
        }
        else return false;
        
    }
      

     // Geometric solution, source: 
    Vector L = Vector(0,0,0) - ray.origin; // centre - origin
    double pro = L.dot(ray.direction);
      if(pro < 0)
      	return false;

    double t4 = L.dot(L) - pro * pro;
       if( t4 > radius * radius)
       	return false;

    double thc = sqrt(radius * radius - t4);
    double t0  = pro - thc;
    double t1  = pro + thc;
    
       
    // swap them if t0 is greater than t1
    if (t0 > t1) std::swap(t0, t1);

    if (t0 < 0) {
        t0 = t1;  // use t1 in this case 
    
    if (t0 < 0 )  // a case where both are negative
    	return false; 
      
    }
    if (t0 > 0 && t1 > 0 && t1 < hit.depth){
        hit.depth = t1;
        hit.normal = (ray.origin + t1 * ray.direction).normalized();
        hit.position = ray.origin + t1 * ray.direction;
        return true;

    } 

    if(t0 > 0 && t1 < 0 && t0 < hit.depth){
        hit.depth = t0;
        hit.normal=(ray.origin + t0 * ray.direction).normalized();
        hit.position=ray.origin + t0 * ray.direction;
        return false;
    }
       else return false;

      double t = t0;
      Vector p = ray.origin + ray.direction * t;
      Vector n = (p - Vector(0,0,0)).normalized();

    //double t = (-ray.origin[2]/ray.direction[2]);
    //Vector n(0,0,1);

    if (t >= hit.depth){
        return false;
    }
    else
    {
    	hit.normal = n;
        hit.depth = t;
        hit.position = ray.origin + ray.direction * t;
        return true;
        
    }
      
    
}


bool Plane::localIntersect(Ray const &ray, Intersection &hit) const
{
	// YOUR CODE HERE
	// Do not accept intersection while ray is inside the plane
	// Here in local coordinate system, the plane is at z = 0
	//
	// NOTE: hit.depth is the current closest intersection depth, so don't
	// accept any intersection that happens further away than that.
  
    double t = (-ray.origin[2])/(ray.direction[2]);
    Vector n(0,0,1);  // define vector n;

    if (t < 0 || ( t > hit.depth)){
        return false;
    }
    else
       {
        hit.normal = n;
        hit.depth =  t;
        hit.position = ray.origin + hit.depth * ray.direction;
        return true;
   }
}
  // helper function for checking a sign

// helper function for checking a sign on the three sides of the triangle.
int checkSign(Vector normal, Vector intersectPoint, Vector edgeStart, Vector edgeStop)
{
    Vector v  = edgeStop - edgeStart; 
    Vector vec = intersectPoint - edgeStart; 
    Vector result = v.normalized().cross(vec.normalized());
    double result1 = normal.dot(result);
    if (result1 >= 0)
    {
        return false;
    }
    else
    {
        return true;
    }
}
bool Mesh::intersectTriangle(Ray const &ray,Triangle const &tri, Intersection &hit) const
{
    // Extract vertex positions from the mesh data.
    Vector const &p0 = positions[tri[0].pi];
    Vector const &p1 = positions[tri[1].pi];
    Vector const &p2 = positions[tri[2].pi];

    //YOUR CODE HERE
    // Decide whether ray intersects the triangle (p0,p1,p2).
    // If so, fill in intersection information in hit and return true.
    // You may find it useful to use the routine implicitLineEquation()
    // to compute the result of the implicit line equation in 2D.
    //
    // NOTE: hit.depth is the current closest intersection depth, so don't
    // accept any intersection that happens further away than that.
    Vector e1 = p1 - p0;
    Vector e2 = p2 - p0;

    Vector normal = e1.cross(e2).normalized();  // normal vector
    double d = normal[0] * p0[0] + normal [1] * p0[1] + normal[2] * p0[2] ; // 
    double top = (ray.origin.dot(normal) - d) * -1;
    double bottom = ray.direction.dot(normal);
    double t = top/bottom;
    if (t > 0 && t < hit.depth)
    {
        Vector intersectionP = ray.origin + t * ray.direction;
          // introduce a helper function to check if the sign is negative or not
        if((!checkSign(normal ,intersectionP, p0, p1)) && (!checkSign(normal ,intersectionP, p1, p2)) && (!checkSign(normal ,intersectionP, p2, p0))){
            hit.position = ray.origin + t * ray.direction;
            hit.normal = normal;
            hit.depth = t;
            return true;
        }
    }
           return false;
  
}

bool Conic::localIntersect(Ray const &ray, Intersection &hit) const {
    // YOUR CODE HERE (creative license)
    // coming soon

    //const Point3D Cone::CAP_CENTRE = Point3D(0, 0, 0);
   // const Vector3D Cone::CAP_NORMAL = Vector3D(0, -1, 0);

    //int cflag = 0;

    double t_near = zMin, t_far = zMax;
    double t[2];

    Intersection near, far;
   // double t_near;
    double x_D = ray.direction[0];
    double x_O = ray.origin[0];
    double y_D = ray.direction[1];
    double y_O = ray.origin[1];
    double z_D = ray.direction[2];
    double z_O = ray.direction[2];
    
   double radius1 =  this->radius1;
   double radius2 =  this->radius2;
   double fac = (radius1 * radius1)/(radius2 * radius2);

   return false;
}
  

// Intersections!
bool Mesh::localIntersect(Ray const &ray, Intersection &hit) const
{
	// Bounding box check
	double tNear = -DBL_MAX, tFar = DBL_MAX;
	for (int i = 0; i < 3; i++) {
		if (ray.direction[i] == 0.0) {
			if (ray.origin[i] < bboxMin[i] || ray.origin[i] > bboxMax[i]) {
				// Ray parallel to bounding box plane and outside of box!
				return false;
			}
			// Ray parallel to bounding box plane and inside box: continue;
		}
		else {
			double t1 = (bboxMin[i] - ray.origin[i]) / ray.direction[i];
			double t2 = (bboxMax[i] - ray.origin[i]) / ray.direction[i];
			if (t1 > t2) std::swap(t1, t2); // Ensure t1 <= t2

			if (t1 > tNear) tNear = t1; // We want the furthest tNear
			if (t2 < tFar) tFar = t2; // We want the closest tFar

			if (tNear > tFar) return false; // Ray misses the bounding box.
			if (tFar < 0) return false; // Bounding box is behind the ray.
		}
	}
	// If we made it this far, the ray does intersect the bounding box.

	// The ray hits the bounding box, so check each triangle.
	bool isHit = false;
	for (size_t tri_i = 0; tri_i < triangles.size(); tri_i++) {
		Triangle const &tri = triangles[tri_i];

		if (intersectTriangle(ray, tri, hit)) {
			isHit = true;
		}
	}
	return isHit;
}

double Mesh::implicitLineEquation(double p_x, double p_y,
	double e1_x, double e1_y,
	double e2_x, double e2_y) const
{
	return (e2_y - e1_y)*(p_x - e1_x) - (e2_x - e1_x)*(p_y - e1_y);
}
