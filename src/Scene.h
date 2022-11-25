#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <string>
#include "Mesh.h"
#include "Sphere.h"
#include "Square.h"
#include <algorithm>

#include <GL/glut.h>


enum LightType {
    LightType_Spherical,
    LightType_Quad
};

struct Light {
    Vec3 material;
    bool isInCamSpace;
    LightType type;

    Vec3 pos;
    float radius;

    Mesh quad;

    float powerCorrection;

    Light() : powerCorrection(1.0) {}

};

struct RaySceneIntersection{
    bool intersectionExists;
    unsigned int typeOfIntersectedObject; //dans quelle liste aller cherche liste des mesh, square, sphere
    unsigned int objectIndex;// indice dans la liste des mesh, square ...
    float t; //permettre de comparer toutes les intersection
    RayTriangleIntersection rayMeshIntersection;
    RaySphereIntersection raySphereIntersection;
    RaySquareIntersection raySquareIntersection;
    RaySceneIntersection() : intersectionExists(false) , t(FLT_MAX) {}

    Vec3 normal() {
        if (typeOfIntersectedObject == 0) {
            return rayMeshIntersection.normal;
        }
        else if (typeOfIntersectedObject == 1) {
            return raySphereIntersection.normal;
        }
        else if (typeOfIntersectedObject == 2) {
            return raySquareIntersection.normal;
        }
    }

    Vec3 intersection() {
        if (typeOfIntersectedObject == 0) {
            return rayMeshIntersection.intersection;
        }
        else if (typeOfIntersectedObject == 1) {
            return raySphereIntersection.intersection;
        }
        else if (typeOfIntersectedObject == 2) {
            return raySquareIntersection.intersection;
        }
    }
};



class Scene {
    std::vector< Mesh > meshes;
    std::vector< Sphere > spheres;
    std::vector< Square > squares;
    std::vector< Light > lights;

public:


    Scene() {
    }

    void draw() {
        // iterer sur l'ensemble des objets, et faire leur rendu :
        for( unsigned int It = 0 ; It < meshes.size() ; ++It ) {
            Mesh const & mesh = meshes[It];
            mesh.draw();
        }
        for( unsigned int It = 0 ; It < spheres.size() ; ++It ) {
            Sphere const & sphere = spheres[It];
            sphere.draw();
        }
        for( unsigned int It = 0 ; It < squares.size() ; ++It ) {
            Square const & square = squares[It];
            square.draw();
        }

    }


    Material getMaterial (RaySceneIntersection intersect) {
        if (intersect.typeOfIntersectedObject == 0) {
            return meshes[intersect.objectIndex].material;
        }
        else if (intersect.typeOfIntersectedObject == 1) {
            return spheres[intersect.objectIndex].material;
        }
        else if (intersect.typeOfIntersectedObject == 2) {      
            return squares[intersect.objectIndex].material;
        }
    }


    Vec3 getColor(RaySceneIntersection intersect) {
        if (intersect.typeOfIntersectedObject == 0) {
            return meshes[intersect.objectIndex].material.diffuse_material;
        }
        else if (intersect.typeOfIntersectedObject == 1) {
            return spheres[intersect.objectIndex].material.diffuse_material;
        }
        else if (intersect.typeOfIntersectedObject == 2) {      
            return squares[intersect.objectIndex].material.diffuse_material;
        }
    }


    RaySceneIntersection computeIntersection(Ray const & ray, float znear) {
        //TODO calculer les intersections avec les objets de la scene et garder la plus proche

        //on apelle les intersect des sphere, des square

        //position = origin + t*dIJ -> pour avoir la sphere la plus proche
        //on veut tous les t >= 0 et on veut savoir lesquels sont avant lesquels
        //on va comparer les t de toutes les intersections pour savoir quel objet est le plus proche
        //result.intersectionExists = true;

        RaySceneIntersection result;
        RaySphereIntersection raySphere;
        RaySquareIntersection raySquare;

        for (int i = 0; i < spheres.size(); i++) {
            raySphere = spheres[i].intersect(ray);
            if(raySphere.t >= 0 && raySphere.t < result.t && raySphere.intersectionExists && raySphere.t > znear){
                result.raySphereIntersection = raySphere;
                result.intersectionExists = true;
                result.objectIndex = i;
                result.typeOfIntersectedObject = 1;
                result.t = raySphere.t;
            }
        }

        for (int i = 0; i < squares.size(); i++) {
            raySquare = squares[i].intersect(ray);
            if(raySquare.t >= 0 && raySquare.t < result.t && raySquare.intersectionExists && raySquare.t > znear){
                result.raySquareIntersection = raySquare;
                result.intersectionExists = true;
                result.objectIndex = i;
                result.typeOfIntersectedObject = 2;
                result.t = raySquare.t;
            }
        }

        return result;
    }

    Vec3 displayPhongIllumination(RaySceneIntersection intersect, Ray ray, float znear) {
        Material material;
        Vec3 color;

        Vec3 ptIntersect = intersect.intersection();
        Vec3 normal = intersect.normal();
        Vec3 view = -1 * ray.direction();
        view.normalize(); 
        Vec3 lightVec;

        Vec3 ambientVec, diffuseVec, specularVec;  
        
        material = getMaterial(intersect);

        for (int n = 0; n < lights.size(); n++) {
            ambientVec = material.ambient_material;

            lightVec = lights[n].pos - ptIntersect;
            lightVec.normalize();

            float diffuse = Vec3::dot(normal, lightVec);
            diffuseVec = material.diffuse_material * diffuse;
        

            Vec3 reflection = 2*(Vec3::dot(normal, lightVec)) * normal - lightVec;
            reflection.normalize();

            float specular = pow(std::max(0.f,Vec3::dot(reflection, view)), material.shininess);
            specularVec = material.specular_material * specular;

            color += Vec3::compProduct(lights[n].material, ambientVec + diffuseVec + specularVec);

        }

        return color;
    }

    
    Vec3 displayPhongIlluminationWithShadow(RaySceneIntersection intersect, Ray ray, float znear) {
        Material material;
        Vec3 color;
        bool isShadow;

        Vec3 ptIntersect = intersect.intersection();
        Vec3 normal = intersect.normal();
        Vec3 view = -1 * ray.direction();
        view.normalize(); 

        Vec3 ambientVec, diffuseVec, specularVec;  
        
        material = getMaterial(intersect);

        if (intersect.typeOfIntersectedObject == 1) {
        for (int n = 0; n < lights.size(); n++) {
            ambientVec = Vec3::compProduct(lights[n].material, material.ambient_material);


            Vec3 lightVec = lights[0].pos - intersect.intersection();
            lightVec.normalize();
        
            Ray rayon = Ray(intersect.intersection(), lightVec);
            RaySceneIntersection rayIntersection = computeIntersection(rayon, znear);

            if (rayIntersection.intersectionExists) {
                return Vec3(0., 0., 0.);
            }
            else {
                float diffuse = Vec3::dot(normal, lightVec);
                diffuseVec = Vec3::compProduct(lights[n].material, material.diffuse_material) * diffuse;


                Vec3 reflection = 2*(Vec3::dot(normal, lightVec)) * normal - lightVec;
                reflection.normalize();

                float specular = pow(std::max(0.f,Vec3::dot(reflection, view)), material.shininess);
                specularVec = Vec3::compProduct(lights[n].material, material.specular_material) * specular;

                color = ambientVec + diffuseVec + specularVec;
            }

        }
        }

        return color;
    }


    Vec3 rayTraceRecursive( Ray ray , int NRemainingBounces, float znear) {
        //TODO RaySceneIntersection raySceneIntersection = computeIntersection(ray);
        //une fois qu'on a le pt intersection, on va ensuite reappeler raytrace pour lancer un rayon avec la direction le rayon réfléchi et l'origine notre point d'intersection jusqu'à un nb (on peut dire 5)
        Vec3 color;
        RaySceneIntersection raySceneIntersection = computeIntersection(ray, znear);

        if (!raySceneIntersection.intersectionExists) {
            return Vec3(0, 0, 0);
        }

        if (NRemainingBounces == 0) {
            return color;
        }

        color = displayPhongIllumination(raySceneIntersection, ray, znear);
        //color = displayPhongIlluminationWithShadow(raySceneIntersection, ray, znear);

        /*Vec3 reflection = ray.direction() - 2*(Vec3::dot(raySceneIntersection.normal(), ray.direction())) * raySceneIntersection.normal();
        Ray reflectionRay = Ray(raySceneIntersection.intersection(), reflection);*/

        NRemainingBounces--;
        rayTraceRecursive(ray, NRemainingBounces, znear);
        
        return color;
    }



    Vec3 rayTrace( Ray const & rayStart , float znear) {
        //TODO appeler la fonction recursive
        //pour i, j de l'image du rendu
        //colorIJ = rayTraceR(r(i, j))
        Vec3 color;

       /*RaySceneIntersection intersect = computeIntersection(rayStart, znear);
       if (intersect.intersectionExists) {
            color = displayPhongIllumination(intersect, rayStart, znear);
        }

        else if (!intersect.intersectionExists) {
            return getColor(intersect);
        }*/

        color = rayTraceRecursive(rayStart, 5, znear);

        return color;
    }

    /*void setup_single_sphere() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        //Sphère rouge
        {
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(2. , 0. , -4);
            s.m_radius = 2.5f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 0.2,0.2,0.2 );
            s.material.shininess = 20;
        }
        //Sphère verte
        {
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-2. , 0. , -2.5);
            s.m_radius = 3.5f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 0.,1.,0. );
            s.material.specular_material = Vec3( 0.2,0.2,0.2 );
            s.material.shininess = 20;
        }
    }*/

    /*void setup_single_square() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        {
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.0 ,0.,0. );
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }
    }*/

    void setup_cornell_box(){
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0.0, 1.5, 0.0 );
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;

        }

        { //Back Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.,1.,0. );
            s.material.specular_material = Vec3( 1.,1.,1. );
            s.material.shininess = 16;
        }

        { //Left Wall

            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
        }

        { //Right Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,0.,1. );
            s.material.specular_material = Vec3( 0.0,1.0,0.0 );
            s.material.shininess = 16;
        }

        { //Floor
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        { //Ceiling
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.,0.,1. );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        { //Front Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(180);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }


        { //GLASS Sphere

            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.4;
        }


        { //MIRRORED Sphere
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3( 0.,0.,1. );
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }

        /*{ 
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(0., 2.2, 0.);
            s.m_radius = 0.2f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3( 0., 0., 0);
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }*/

    }

};



#endif
