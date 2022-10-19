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


    RaySphereIntersection intersect(const Ray &ray) const {
        RaySphereIntersection intersection;
        //TODO calcul l'intersection rayon sphere
        int distance = sqrt(pow((ray.direction()[0] - ray.origin()[0]), 2) + pow((ray.direction()[1] - ray.origin()[1]), 2) + pow((ray.direction()[2] - ray.origin()[2]), 2));
        Vec3 P;
        float a, b, k, t1, t2;
        for (int i = 0; i < distance; i++) {
            P = ray.origin() + i*ray.direction();
            if (Vec3::dot((P-m_center), (P-m_center)) == pow(m_radius, 2)) {
                a = Vec3::dot(ray.direction(), ray.direction());
                b = Vec3::dot(2*ray.direction(), ray.origin()-m_center);
                k = Vec3::dot(ray.origin()-m_center, ray.origin()-m_center)-pow(m_radius, 2);

                t1 = (-b+sqrt(pow(b, 2)-4*a*k))/2*a;
                t2 = (-b-sqrt(pow(b, 2)-4*a*k))/2*a;

                if (t1 > 0 && t2 > 0) {
                    if (t1 > t2) {
                        intersection.t = t2;
                    }
                    else {
                        intersection.t = t1;
                    }
                }
                else if (t1 > 0 && t2 < 0) {
                    intersection.t = t1;
                }
                else if (t1 < 0 && t2 > 0) {
                    intersection.t = t2;
                }

                intersection.intersectionExists = true;
                intersection.intersection = Vec3(ray.origin()[0] + i*ray.direction()[0], ray.origin()[1] + i*ray.direction()[1], ray.origin()[2] + i*ray.direction()[2]);
                //std::cout << intersection << std::endl;
                /*Vec3 SphericalCoordinates = SphericalCoordinatesToEuclidean(intersection.intersection);
                intersection.theta = SphericalCoordinates[0];
                intersection.phi = SphericalCoordinates[1];
                intersection.normal = intersection.intersection - m_center;*/
                std::cout << intersection.intersection << std::endl;
                return intersection;
            }
        }
        intersection.intersectionExists = false;
        return intersection;
    }
};
#endif
