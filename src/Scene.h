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

    Vec3 getNormal() {
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

    Vec3 getIntersection() {
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

class AABB {
    public :
    Vec3 min, max;

    AABB(std::vector< Vec3 > meshes) {
        float xmin = FLT_MAX, ymin = FLT_MAX, zmin = FLT_MAX, xmax = -FLT_MAX, ymax = -FLT_MAX, zmax = -FLT_MAX;

        for (int elt = 0; elt < meshes.size(); elt ++) {

            if (meshes[elt][0] < xmin) {
                xmin = meshes[elt][0];
            }
            if (meshes[elt][0] >= xmax) {
                xmax = meshes[elt][0];
            }

            if (meshes[elt][1] < ymin) {
                ymin = meshes[elt][1];
            }
            if (meshes[elt][1] >= ymax) {
                ymax = meshes[elt][1];
            }

            if (meshes[elt][2] < zmin) {
                zmin = meshes[elt][2];
            }
            if (meshes[elt][2] >= zmax) {
                zmax = meshes[elt][2];
            }
        }

        float epsilon = 0.1;

        min = Vec3(xmin-epsilon, ymin-epsilon, zmin-epsilon);
        max = Vec3(xmax+epsilon, ymax+epsilon, zmax+epsilon);
    }

    AABB(AABB* aabb1, AABB* aabb2) {
        if (aabb1 != nullptr && aabb2 != nullptr) {
            min[0] = std::min(aabb1->min[0], aabb2->min[0]);
            max[0] = std::max(aabb1->max[0], aabb2->max[0]);

            min[1] = std::min(aabb1->min[1], aabb2->min[1]);
            max[1] = std::max(aabb1->max[1], aabb2->max[1]);

            min[2] = std::min(aabb1->min[2], aabb2->min[2]);
            max[2] = std::max(aabb1->max[2], aabb2->max[2]);
        }
        else if (aabb1 == nullptr) {
            min[0] = aabb2->min[0];
            max[0] = aabb2->max[0];

            min[1] = aabb2->min[1];
            max[1] = aabb2->max[1];

            min[2] = aabb2->min[2];
            max[2] = aabb2->max[2];
        }
        else if (aabb2 == nullptr) {
            min[0] = aabb1->min[0];
            max[0] = aabb1->max[0];

            min[1] = aabb1->min[1];
            max[1] = aabb1->max[1];
            
            min[2] = aabb1->min[2];
            max[2] = aabb1->max[2];
        }
    }
};


class KdTreeNode {
    public :
    KdTreeNode *backNode = nullptr, *frontNode = nullptr;
    std::vector<Vec3> point;

    int dimSplit;
    float splitDistance;
    bool isLeaf = false;

    AABB* aabb = nullptr;

    AABB* getAABB() {
        return aabb;
    }

    KdTreeNode(std::vector<Vec3> objectsList, int depth = 0) {
        if (objectsList.empty()) {
            return;
        }
        else if(objectsList.size() == 1 || depth > 2) {
            point = objectsList;
            isLeaf = true;
            aabb = new AABB(objectsList);

            return;
        }

        this->splitDistance = std::round(objectsList.size() / 2.0);

        dimSplit = depth%3;

        std::sort(objectsList.begin(), objectsList.end(), [&](Vec3 a, Vec3 b) {
            return a[dimSplit] < b[dimSplit];
        });

        std::vector<Vec3> objectsListBack;
        std::vector<Vec3> objectsListFront;

        for (int i = 0; i < splitDistance; i++) {
            objectsListBack.push_back(objectsList[i]);
        }

        for (int j = splitDistance; j < objectsList.size(); j++) {
            objectsListFront.push_back(objectsList[j]);
        }

        backNode = new KdTreeNode(objectsListBack, depth+1);
        frontNode = new KdTreeNode(objectsListFront, depth+1);

        if (backNode != nullptr || frontNode != nullptr) {
            aabb = new AABB(backNode->getAABB(), frontNode->getAABB());
        }

    };

    float intersection(Ray &ray) {
        if (ray.direction()[0] == 0 && ray.origin()[0] < aabb->min[0] && ray.origin()[0] > aabb->max[0]) {
            return -1;
        }
       
        float tMinX = (aabb->min[0] - ray.origin()[0]) / ray.direction()[0];
        float tMaxX = (aabb->max[0] - ray.origin()[0]) / ray.direction()[0];

        float tMinY = (aabb->min[1] - ray.origin()[1]) / ray.direction()[1];
        float tMaxY = (aabb->max[1] - ray.origin()[1]) / ray.direction()[1];

        float tMinZ = (aabb->min[2] - ray.origin()[2]) / ray.direction()[2];
        float tMaxZ = (aabb->max[2] - ray.origin()[2]) / ray.direction()[2];


        if (tMinX > tMaxX) {
            std::swap(tMinX, tMaxX);
        }
        if (tMinY > tMaxY) {
            std::swap(tMinY, tMaxY);
        }
        if (tMinZ > tMaxZ) {
            std::swap(tMinZ, tMaxZ);
        }


        float tStartXY = std::max(tMinX, tMinY);
        float tEndXY = std::min(tMaxX, tMaxY);

        float tMinXY = std::min(ray.origin()[0], ray.origin()[1]);

        if (tStartXY > tEndXY) {
            return -1;
        }
        if (tEndXY < tMinXY) {
            return -1;
        }

        if(tStartXY > tMinXY) {
            return tStartXY;
        }
        else {
            return tEndXY;
        }


    }

    float traverseKdTree(Ray &ray, float t_start, float t_end) {

        if(isLeaf) {
            return intersection(ray);
        }

        float t = (splitDistance - ray.origin()[this->dimSplit]) / ray.direction()[this->dimSplit];

        if (t <= t_start) {
            return backNode->traverseKdTree(ray, t_start, t_end);
        }
        else if (t >= t_end) {
            return frontNode->traverseKdTree(ray, t_start, t_end);
        }
        else {
            float hit = frontNode->traverseKdTree(ray, t_start, t);

            if (hit <= t) {
                return hit;
            }

            return backNode->traverseKdTree(ray, t, t_end);
        }

    }
};


struct BoundingBox {
    Vec3 min;
    Vec3 max;

    std::vector<Mesh> triangleList;
};


class Scene {
    std::vector< Mesh > meshes;
    std::vector< Sphere > spheres;
    std::vector< Square > squares;
    std::vector< Light > lights;

    KdTreeNode *node = nullptr;

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

    
    RaySceneIntersection computeIntersection(Ray & ray, float znear) {
        //TODO calculer les intersections avec les objets de la scene et garder la plus proche

        //on apelle les intersect des sphere, des square

        //position = origin + t*dIJ -> pour avoir la sphere la plus proche
        //on veut tous les t >= 0 et on veut savoir lesquels sont avant lesquels
        //on va comparer les t de toutes les intersections pour savoir quel objet est le plus proche
        //result.intersectionExists = true;

        RaySceneIntersection result;
        RayTriangleIntersection rayMesh;
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

            //Sans Kd-Tree
         /*for (int i = 0; i < meshes.size(); i++) {
            rayMesh = meshes[i].intersect(ray);

            if(rayMesh.t >= 0 && rayMesh.t < result.t && rayMesh.intersectionExists && rayMesh.t > znear){
                result.rayMeshIntersection = rayMesh;
                result.intersectionExists = true;
                result.objectIndex = i;
                result.typeOfIntersectedObject = 0;
                result.t = rayMesh.t;
            }
            }*/

            //Kd-Tree
            if (node == nullptr) {
                std::vector<Vec3> listForKdTree;

                for (int i = 0; i < meshes.size(); i++) {
                    for (int j = 0; j < meshes[i].vertices.size(); j++) {
                        listForKdTree.push_back(meshes[i].vertices[meshes[i].triangles[j][0]].position);
                    }
                }

                node = new KdTreeNode(listForKdTree);
            }

            float tMinX = (node->aabb->min[0] - ray.origin()[0]) / ray.direction()[0];
            float tMaxX = (node->aabb->max[0] - ray.origin()[0]) / ray.direction()[0];

            float tMinY = (node->aabb->min[1] - ray.origin()[1]) / ray.direction()[1];
            float tMaxY = (node->aabb->max[1] - ray.origin()[1]) / ray.direction()[1];

            float tStart = std::max(tMinX, tMinY);
            float tEnd = std::min(tMaxX, tMaxY);

            float t = node->traverseKdTree(ray, tStart, tEnd);

            if (t < result.t) {
                for (int i = 0; i < meshes.size(); i++) {
                rayMesh = meshes[i].intersect(ray);

                    if(rayMesh.t >= 0 && rayMesh.t < result.t && rayMesh.intersectionExists && rayMesh.t > znear){
                        result.rayMeshIntersection = rayMesh;
                        result.intersectionExists = true;
                        result.objectIndex = i;
                        result.typeOfIntersectedObject = 0;
                        result.t = rayMesh.t;
                    }
                }
            }

        return result;
    }

    Vec3 displayPhongIllumination(RaySceneIntersection intersect, Ray ray, float znear) {
        Vec3 color;

        Vec3 ptIntersect = intersect.getIntersection();
        Vec3 normal = intersect.getNormal();
        Vec3 view = -1 * ray.direction();
        view.normalize(); 
        Vec3 lightVec;

        Vec3 diffuseVec, specularVec;  
        
        Material material = getMaterial(intersect);

        Vec3 ambientVec = material.ambient_material;

        //for (int n = 0; n < lights.size(); n++) {
            lightVec = lights[0].pos - ptIntersect;
            lightVec.normalize();

            float diffuse = Vec3::dot(normal, lightVec);
            diffuseVec = material.diffuse_material * diffuse;
        

            Vec3 reflection = 2*(Vec3::dot(normal, lightVec)) * normal - lightVec;
            reflection.normalize();

            float specular = pow(std::max(0.f,Vec3::dot(reflection, view)), material.shininess);
            specularVec = material.specular_material * specular;

            color += Vec3::compProduct(lights[0].material, ambientVec + diffuseVec + specularVec);

        //}

        return color;
    }

     Vec3 displayShadow(RaySceneIntersection raySceneIntersection, Ray ray, float znear) {
        Vec3 color;
        Vec3 intersectPt = raySceneIntersection.getIntersection();

        Vec3 L = lights[0].pos - intersectPt;
        L.normalize();

        Ray rayS = Ray(intersectPt, L);
        RaySceneIntersection shadow = computeIntersection(rayS, 0.01);

        if (shadow.intersectionExists && shadow.t <= L.length()) {
            color = Vec3(0, 0, 0);
        }
        else {
            color = displayPhongIllumination(raySceneIntersection, ray, znear);
        }

        return color;
    }

    Vec3 displaySoftShadow(RaySceneIntersection raySceneIntersection, Ray ray,float znear) {
        Vec3 color;
        Vec3 intersectPt = raySceneIntersection.getIntersection();

        int nbLumieres = 10;
        int cpt = 0;

        for (int i = 0; i < nbLumieres; i++) {

            //on prend x et z aléatoires entre -2 et 2
            double randomX = -0.5 + ((double)rand()/RAND_MAX) * (0.5 + 0.5);
            double randomZ = -0.5 + ((double)rand()/RAND_MAX) * (0.5 + 0.5);

            Vec3 randomLight = Vec3(lights[i].pos[0] + randomX, lights[i].pos[1], lights[i].pos[2] + randomZ);
        
            Vec3 L = randomLight - intersectPt;
            L.normalize();

            Ray rayS = Ray(intersectPt, L);
            RaySceneIntersection shadow = computeIntersection(rayS, 0.00001);

            if (shadow.intersectionExists && shadow.t <= L.length() && shadow.typeOfIntersectedObject == 1) {
                color = displayPhongIllumination(raySceneIntersection, ray, znear);
                cpt++;
            }
            else {
                color = displayPhongIllumination(raySceneIntersection, ray, znear);
            }
        }
        
        if (cpt > 0) {
            color /= cpt;
        }

        return color;
    }


    Vec3 rayTraceRecursive( Ray ray , int NRemainingBounces, float znear) {
        //TODO RaySceneIntersection raySceneIntersection = computeIntersection(ray);
        //une fois qu'on a le pt intersection, on va ensuite reappeler raytrace pour lancer un rayon avec la direction le rayon réfléchi et l'origine notre point d'intersection jusqu'à un nb (on peut dire 5)
        Vec3 color;

        RaySceneIntersection raySceneIntersection = computeIntersection(ray, znear);

        Material m = getMaterial(raySceneIntersection);

        if (!raySceneIntersection.intersectionExists) {
            return Vec3(0, 0, 0);
        }

        if (NRemainingBounces == 0) {
            return color;
        }
        

        if (m.type == Material_Mirror || m.type == Material_Glass) {
            color = m.diffuse_material;
        }
        else {
            //return getColor(raySceneIntersection);
            //return displayPhongIllumination(raySceneIntersection, ray, znear);
            //return displayShadow(raySceneIntersection, ray, znear);
            return displaySoftShadow(raySceneIntersection, ray, 0.0);
        }

        if (raySceneIntersection.typeOfIntersectedObject == 0) {
            color = Vec3(raySceneIntersection.rayMeshIntersection.w0, raySceneIntersection.rayMeshIntersection.w1, 1 - raySceneIntersection.rayMeshIntersection.w0 - raySceneIntersection.rayMeshIntersection.w1);
        }

        if(m.type == Material_Mirror) {

            Vec3 reflection = ray.direction() - 2*(Vec3::dot(ray.direction(), raySceneIntersection.getNormal())) * raySceneIntersection.getNormal();
            Ray reflectionRay = Ray(raySceneIntersection.getIntersection(), reflection);

            if (Vec3::dot(reflectionRay.direction(), raySceneIntersection.getNormal()) > 0) {
                return rayTraceRecursive(reflectionRay, NRemainingBounces-1, 0.005);
            }

        }
        else if (m.type == Material_Glass) {

            float cosI = std::clamp(-1.0f, 1.0f, Vec3::dot(ray.direction(), raySceneIntersection.getNormal()));
            float etaI = 1; float etaT = m.index_medium;
            Vec3 n = raySceneIntersection.getNormal();

            if (cosI < 0) {
                cosI = -cosI;
            }
            else {
                std::swap(etaI, etaT);
                n = -1*raySceneIntersection.getNormal();
            }

            float eta = etaI/etaT;
            float k = 0.5 - eta * eta * (0 - cosI * cosI);

            if ( k > 0) {
                Vec3 refraction = eta * ray.direction() + (eta * cosI - sqrtf(k)) * n;
                Ray refractionRay = Ray(raySceneIntersection.getIntersection(), refraction);
                return rayTraceRecursive(refractionRay, NRemainingBounces-1, 0.001);
            }
        }

    }


    Vec3 rayTrace( Ray const & rayStart , float znear) {
        //TODO appeler la fonction recursive
        //pour i, j de l'image du rendu
        //colorIJ = rayTraceR(r(i, j))
        Vec3 color;

        color = rayTraceRecursive(rayStart, 5, znear);

        return color;
    }

    void setup_single_sphere() {
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
    }

    void setup_single_square() {
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
    }

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
            s.material.diffuse_material = Vec3( 1.,1.,1. );
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
            s.material.diffuse_material = Vec3( 0.0,1.0,0.0 );
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
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
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
            s.material.type = Material_Glass;
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
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,1.,1. );
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }

        { 
            meshes.resize(meshes.size() + 1);
            Mesh & m = meshes[meshes.size() - 1];
            m.loadOFF("data/tetrahedron.off");
            m.centerAndScaleToUnit();
            m.translate(Vec3(0., 0.5, 0));
            m.build_arrays();
            m.material.diffuse_material = Vec3(1, 0, 1);
            m.material.specular_material = Vec3(1, 0, 1);
            m.material.shininess = 16;
        }

    }

    void setup_cornell_box_personnalised(){
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
            //light.type = LightType_Spherical;
            light.type = LightType_Quad;
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
            s.material.diffuse_material = Vec3(0.72, 0.69, 0.55);
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
            s.material.diffuse_material = Vec3(0.15, 0.22, 0.27);
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
            s.material.diffuse_material = Vec3(0.94, 0.83, 0.57);
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
            s.material.diffuse_material = Vec3(0.9, 0.9, 0.9);
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
            s.material.diffuse_material = Vec3(0.52, 0.62, 0.72);
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
            /*s.material.diffuse_material = Vec3(0.23, 0.32, 0.41);
            s.material.specular_material = Vec3(0.23, 0.32, 0.41);*/
            s.material.diffuse_material = Vec3(1, 1, 1);
            s.material.specular_material = Vec3(1, 1, 1);
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
            s.material.diffuse_material = Vec3(0.94, 0.58, 0.34);
            s.material.specular_material = Vec3(0.94, 0.58, 0.34);
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }
        
        /*{ 
            meshes.resize(meshes.size() + 1);
            Mesh & m = meshes[meshes.size() - 1];
            m.loadOFF("data/tetrahedron.off");
            m.centerAndScaleToUnit();
            m.build_arrays();
            m.material.diffuse_material = Vec3(1, 0, 1);
            m.material.specular_material = Vec3(1, 0, 1);
            m.material.shininess = 16;
        }*/

    }

};



#endif
