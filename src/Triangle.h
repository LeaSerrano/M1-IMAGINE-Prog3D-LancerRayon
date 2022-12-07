#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "Vec3.h"
#include "Ray.h"
#include "Plane.h"

struct RayTriangleIntersection{
    bool intersectionExists;
    float t;
    float w0,w1,w2;//coordonnées barycentrique, permettra de faire les calcul de shading où si le point est contenu dans le triangle w2 = Aire(T2)/Aire(T), p = w0*v0 + w1*v1 + w2*v2 -> coordonnées du point p, n = w0n0 + w1n1 + w2n2 -> normale du point p
    unsigned int tIndex;
    Vec3 intersection;
    Vec3 normal;
};

class Triangle {
private:
    Vec3 m_c[3] , m_normal;
    float area;
public:
    Triangle() {}
    Triangle( Vec3 const & c0 , Vec3 const & c1 , Vec3 const & c2 ) {
        m_c[0] = c0;
        m_c[1] = c1;
        m_c[2] = c2;
        updateAreaAndNormal();
    }
    void updateAreaAndNormal() {
        Vec3 nNotNormalized = Vec3::cross( m_c[1] - m_c[0] , m_c[2] - m_c[0] );
        float norm = nNotNormalized.length();
        m_normal = nNotNormalized / norm;
        area = norm / 2.f;
    }
    void setC0( Vec3 const & c0 ) { m_c[0] = c0; } // remember to update the area and normal afterwards!
    void setC1( Vec3 const & c1 ) { m_c[1] = c1; } // remember to update the area and normal afterwards!
    void setC2( Vec3 const & c2 ) { m_c[2] = c2; } // remember to update the area and normal afterwards!
    Vec3 const & normal() const { return m_normal; }
    Vec3 projectOnSupportPlane( Vec3 const & p ) const {
        Vec3 result;
        //TODO completer
        return result;
    }
    float squareDistanceToSupportPlane( Vec3 const & p ) const {
        float result;
        //TODO completer
        return result;
    }
    float distanceToSupportPlane( Vec3 const & p ) const { 
        return sqrt( squareDistanceToSupportPlane(p) ); 
    }
    bool isParallelTo( Line const & L ) const {
        bool result;
        //TODO completer
        /*Vec3 OD = L.direction() - L.origin();

        Vec3 P0P1 = m_c[1] - m_c[0];
        Vec3 P1P2 = m_c[2] - m_c[1];
        Vec3 P2P0 = m_c[0] - m_c[2];

        if ((OD[0]/P0P1[0]) == (OD[1]/P0P1[1])) {
            result = true;
        }
        else if ((OD[0]/P1P2[0]) == (OD[1]/P1P2[1])) {
            result = true;
        }
        else if ((OD[0]/P2P0[0]) == (OD[1]/P2P0[1])) {
            result = true;
        }

        result = false;*/

        return result;
    }
    Vec3 getIntersectionPointWithSupportPlane( Line const & L ) const {
        // you should check first that the line is not parallel to the plane!
        Vec3 result;
        //TODO completer
        return result;
    }
    void computeBarycentricCoordinates( Vec3 const & p , float & u0 , float & u1 , float & u2 ) const {
        //TODO Complete
    }

    RayTriangleIntersection getIntersection( Ray const & ray ) const {
        RayTriangleIntersection result;
        // 1) check that the ray is not parallel to the triangle:

        // 2) check that the triangle is "in front of" the ray:

        // 3) check that the intersection point is inside the triangle:
        // CONVENTION: compute u,v such that p = w0*c0 + w1*c1 + w2*c2, check that 0 <= w0,w1,w2 <= 1

        // 4) Finally, if all conditions were met, then there is an intersection! :

        //float area = m_normal.length()/2;

        if (Vec3::dot(m_normal, ray.direction()) != 0) {
            float D = Vec3::dot(m_normal, m_c[0]);
            float t = (D - Vec3::dot(m_normal, ray.origin()))/Vec3::dot(m_normal, ray.direction());

            if (t >= 0) {
                Vec3 P = ray.origin() + t*ray.direction();

                if (Vec3::dot(m_normal, Vec3::cross(m_c[1] - m_c[0], P - m_c[0])) > 0 && Vec3::dot(m_normal, Vec3::cross(m_c[2] - m_c[1], P - m_c[1])) > 0 && Vec3::dot(m_normal, Vec3::cross(m_c[0] - m_c[2], P - m_c[2])) > 0) {
                    result.intersection = P;
                    result.intersectionExists = true;
                    result.t = t;
                    result.normal = m_normal;
                    result.normal.normalize();

                    float areaP12 = Vec3::dot(Vec3::cross(m_c[2] - m_c[1], P - m_c[1]), m_normal);
                    float areaOP2 = Vec3::dot(Vec3::cross(m_c[0] - m_c[2], P - m_c[2]), m_normal);
                    float area01P = Vec3::dot(Vec3::cross(m_c[1] - m_c[0], P - m_c[0]), m_normal);
                    float area012 = Vec3::dot(Vec3::cross(m_c[1] - m_c[0], m_c[2] - m_c[0]), m_normal);

                    result.w0 = areaP12/area012;
                    result.w1 = areaOP2/area012;
                    result.w2 = area01P/area012;
                }

                else {
                    result.intersectionExists = false;
                }
            }
            else {
                result.intersectionExists = false;
            }
        }
        else {
            result.intersectionExists = false;
        }

        return result;
    }
};
#endif
