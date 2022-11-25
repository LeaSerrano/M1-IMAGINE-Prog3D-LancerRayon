#ifndef SQUARE_H
#define SQUARE_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

struct RaySquareIntersection{
    bool intersectionExists;
    float t;
    float u,v;
    Vec3 intersection;
    Vec3 normal;
    Vec3 bouncedir;
};


class Square : public Mesh {
public:
    Vec3 m_normal;
    Vec3 m_bottom_left;//point en bas à gauche
    Vec3 m_right_vector;//vecteur vers la droite
    Vec3 m_up_vector;//vecteur vers le haut
     Mat3 transformMat; //ajout ludovic
    bool hasTransform = false;


    Square() : Mesh() {}
    Square(Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
           float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) : Mesh() {
        setQuad(bottomLeft, rightVector, upVector, width, height, uMin, uMax, vMin, vMax);
    }

    void setQuad( Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
                  float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) {
        m_right_vector = rightVector;
        m_up_vector = upVector;
        m_normal = Vec3::cross(rightVector , upVector);
        m_bottom_left = bottomLeft;

        m_normal.normalize();
        m_right_vector.normalize();
        m_up_vector.normalize();

        m_right_vector = m_right_vector*width;
        m_up_vector = m_up_vector*height;

        vertices.clear();
        vertices.resize(4);
        vertices[0].position = bottomLeft;                                      vertices[0].u = uMin; vertices[0].v = vMin;
        vertices[1].position = bottomLeft + m_right_vector;                     vertices[1].u = uMax; vertices[1].v = vMin;
        vertices[2].position = bottomLeft + m_right_vector + m_up_vector;       vertices[2].u = uMax; vertices[2].v = vMax;
        vertices[3].position = bottomLeft + m_up_vector;                        vertices[3].u = uMin; vertices[3].v = vMax;
        vertices[0].normal = vertices[1].normal = vertices[2].normal = vertices[3].normal = m_normal;
        triangles.clear();
        triangles.resize(2);
        triangles[0][0] = 0;
        triangles[0][1] = 1;
        triangles[0][2] = 2;
        triangles[1][0] = 0;
        triangles[1][1] = 2;
        triangles[1][2] = 3;

    }

    Vec3 bottomLeft() const { 
        return this->vertices[0].position; 
    }

    Vec3 bottomRight() const { 
        return this->vertices[1].position; 
    }

    Vec3 upRight() const { 
        return this->vertices[2].position; 
    }

    Vec3 upLeft() const { 
        return this->vertices[3].position; 
    }
    
    Vec3 normal() const { 
        return Vec3::cross((bottomRight() - bottomLeft()) , (upLeft() - bottomLeft())); 
    }

    RaySquareIntersection intersect(const Ray &ray) const {
        //TODO calculer l'intersection rayon quad
        RaySquareIntersection intersection;

        Vec3 bottomLeft = this->bottomLeft();
        Vec3 bottomRight = this->bottomRight();
        Vec3 upLeft = this->upLeft();
        Vec3 normal = this->normal();

        float t = Vec3::dot(bottomLeft-ray.origin(), normal)/Vec3::dot(ray.direction(), normal);
        Vec3 P = ray.origin() + t*ray.direction();

        Vec3 RL = bottomLeft - bottomRight;
        Vec3 RP = P - bottomRight;
        float u = (Vec3::dot(RL, RP))/RL.length();

        Vec3 UB = bottomLeft - upLeft;
        Vec3 UP = P - upLeft;
        float v = (Vec3::dot(UB, UP))/UB.length();

        if (t >= 0 && u < RL.norm() && u > 0 && v < UB.norm() && v > 0) {
            intersection.intersectionExists = true;
            intersection.t = t;
            intersection.u = u;
            intersection.v = v;
            intersection.intersection = P;
            intersection.normal = normal;
            intersection.normal.normalize();
        }
        else {
            intersection.intersectionExists = false;
        }

        return intersection;

    }

};
#endif // SQUARE_H
