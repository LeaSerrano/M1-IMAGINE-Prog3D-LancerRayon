#ifndef Sphere_H
#define Sphere_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

struct RaySphereIntersection{
    bool intersectionExists;
    float t;
    float theta,phi;//coordonnées sphériques
    Vec3 intersection; 
    Vec3 secondintersection;//dans le cas où on fait de la réfraction, on a deux intersection
    Vec3 normal;//normale à une sphère (position du point - le centre normalisé)
};

static
Vec3 SphericalCoordinatesToEuclidean( Vec3 ThetaPhiR ) {
    return ThetaPhiR[2] * Vec3( cos(ThetaPhiR[0]) * cos(ThetaPhiR[1]) , sin(ThetaPhiR[0]) * cos(ThetaPhiR[1]) , sin(ThetaPhiR[1]) );
}
static
Vec3 SphericalCoordinatesToEuclidean( float theta , float phi ) {
    return Vec3( cos(theta) * cos(phi) , sin(theta) * cos(phi) , sin(phi) );
}

static
Vec3 EuclideanCoordinatesToSpherical( Vec3 xyz ) {
    float R = xyz.length();
    float phi = asin( xyz[2] / R );
    float theta = atan2( xyz[1] , xyz[0] );
    return Vec3( theta , phi , R );
}



class Sphere : public Mesh {
public:
    Vec3 m_center;
    float m_radius;

    Sphere() : Mesh() {}
    Sphere(Vec3 c , float r) : Mesh() , m_center(c) , m_radius(r) {}

    void build_arrays(){
        unsigned int nTheta = 20 , nPhi = 20;
        positions_array.resize(3 * nTheta * nPhi );
        normalsArray.resize(3 * nTheta * nPhi );
        uvs_array.resize(2 * nTheta * nPhi );
        for( unsigned int thetaIt = 0 ; thetaIt < nTheta ; ++thetaIt ) {
            float u = (float)(thetaIt) / (float)(nTheta-1);
            float theta = u * 2 * M_PI;
            for( unsigned int phiIt = 0 ; phiIt < nPhi ; ++phiIt ) {
                unsigned int vertexIndex = thetaIt + phiIt * nTheta;
                float v = (float)(phiIt) / (float)(nPhi-1);
                float phi = - M_PI/2.0 + v * M_PI;
                Vec3 xyz = SphericalCoordinatesToEuclidean( theta , phi );
                positions_array[ 3 * vertexIndex + 0 ] = m_center[0] + m_radius * xyz[0];
                positions_array[ 3 * vertexIndex + 1 ] = m_center[1] + m_radius * xyz[1];
                positions_array[ 3 * vertexIndex + 2 ] = m_center[2] + m_radius * xyz[2];
                normalsArray[ 3 * vertexIndex + 0 ] = xyz[0];
                normalsArray[ 3 * vertexIndex + 1 ] = xyz[1];
                normalsArray[ 3 * vertexIndex + 2 ] = xyz[2];
                uvs_array[ 2 * vertexIndex + 0 ] = u;
                uvs_array[ 2 * vertexIndex + 1 ] = v;
            }
        }
        triangles_array.clear();
        for( unsigned int thetaIt = 0 ; thetaIt < nTheta - 1 ; ++thetaIt ) {
            for( unsigned int phiIt = 0 ; phiIt < nPhi - 1 ; ++phiIt ) {
                unsigned int vertexuv = thetaIt + phiIt * nTheta;
                unsigned int vertexUv = thetaIt + 1 + phiIt * nTheta;
                unsigned int vertexuV = thetaIt + (phiIt+1) * nTheta;
                unsigned int vertexUV = thetaIt + 1 + (phiIt+1) * nTheta;
                triangles_array.push_back( vertexuv );
                triangles_array.push_back( vertexUv );
                triangles_array.push_back( vertexUV );
                triangles_array.push_back( vertexuv );
                triangles_array.push_back( vertexUV );
                triangles_array.push_back( vertexuV );
            }
        }
    }

    /*Vec3 reflection(Vec3 view, Vec3 normal) {
        return view - 2*Vec3::dot(view, normal) * normal;
    }

    float scattering(Ray ray, RaySphereIntersection intersection, Vec3 attenuation) {
        Vec3 reflected = reflection(ray.direction(), intersection.normal);
        Ray scattered = Ray(intersection.intersection, reflected);

        float result = Vec3::dot(scattered.direction(), intersection.normal);
        if (result > 0) {
            return result;
        }
    }*/


    RaySphereIntersection intersect(const Ray &ray) const {
        //TODO calcul l'intersection rayon sphere
        RaySphereIntersection intersection;

        float t, t1, t2;

        float a = Vec3::dot(ray.direction(), ray.direction());
        float b = Vec3::dot((2*ray.direction()), (ray.origin() - m_center));
        float c = Vec3::dot((ray.origin()-m_center), (ray.origin()-m_center)) - pow(m_radius, 2);

        float discriminant = pow(b, 2)-4*a*c;

        if (discriminant > 0 && 2*a != 0 && 4*a*c >= 0) {
            t1 = (-b + sqrt(discriminant))/(2*a);
            t2 = (-b - sqrt(discriminant))/(2*a);

            if (t1 > 0 && t2 > 0) {
                intersection.intersectionExists = true;
                intersection.t = std::min(t1, t2);
                intersection.intersection = ray.origin() + intersection.t*ray.direction();
                intersection.normal = intersection.intersection-m_center;
                intersection.normal.normalize();
            }
            else if (t1 > 0 && t2 < 0) {
                intersection.intersectionExists = true;
                intersection.t = t1;
                intersection.intersection = ray.origin() + intersection.t*ray.direction();
                intersection.normal = intersection.intersection-m_center;
                intersection.normal.normalize();
            }
            else if (t1 < 0 && t2 > 0) {
                intersection.intersectionExists = true;
                intersection.t = t2;
                intersection.intersection = ray.origin() + intersection.t*ray.direction();
                intersection.normal = intersection.intersection-m_center;
                intersection.normal.normalize();
            }
            else {
                intersection.intersectionExists = false;
            }
        }
        else {
            intersection.intersectionExists = false;
        }

        return intersection;
    }
};
#endif
