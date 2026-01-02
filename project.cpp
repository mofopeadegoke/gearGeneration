#include <iostream>
#include <fstream>
#include <string>
#include <math.h>   
#include <vector>

struct Point3D {
    float x, y, z;
};


struct Triangle{
    Point3D v1, v2, v3, normal;
};

enum GearType {
    SPUR,      // Straight teeth (default)
    HELICAL,   // Angled teeth for smoother meshing
    BEVEL,     // Conical gears for changing axis direction
    WORM       // Screw-type gear for high reduction ratios
};

struct GearParams {
    int numTeeth;           // How many teeth the gear has
    float innerRadius;      // Radius of the gear body/hub
    float pitchRadius;      // Middle of the tooth (where gears mesh)
    float outerRadius;      // Tip of the tooth
    float thickness;        // How thick the gear is (Z direction)
    float toothWidth;       // What fraction of space each tooth takes (0-1)
    int segments;           // How many points per tooth (smoothness)
    GearType type;          // Type of gear (SPUR, HELICAL, etc.)
    float helixAngle;       // Angle for helical gears (degrees)
    float pressureAngle;    // Pressure angle for involute teeth (usually 20°)
    bool useInvolute;       // Use involute tooth profile?
    bool roundedTeeth;      // Add rounding to tooth tips and valleys?
    float roundingRadius;   // Radius for tooth rounding
};

void writeTriangle(std::ofstream &file, 
                   float normalX, float normalY, float normalZ,
                   float v1x, float v1y, float v1z,
                   float v2x, float v2y, float v2z,
                   float v3x, float v3y, float v3z) {
    file << "    facet normal " << normalX << " " << normalY << " " << normalZ << "\n";
    file << "        outer loop\n";
    file << "        vertex " << v1x << " " << v1y << " " << v1z << "\n";
    file << "        vertex " << v2x << " " << v2y << " " << v2z << "\n";
    file << "        vertex " << v3x << " " << v3y << " " << v3z << "\n";
    file << "        endloop\n";
    file << "    endfacet\n";
}

void writeTriangle(std::ofstream &file, const Triangle &triangle) {
    writeTriangle(file,
                  triangle.normal.x, triangle.normal.y, triangle.normal.z,
                  triangle.v1.x, triangle.v1.y, triangle.v1.z,
                  triangle.v2.x, triangle.v2.y, triangle.v2.z,
                  triangle.v3.x, triangle.v3.y, triangle.v3.z);
}
Point3D calculateNormal(Point3D v1, Point3D v2, Point3D v3) {
    Point3D u = {v2.x - v1.x, v2.y - v1.y, v2.z - v1.z};
    Point3D v = {v3.x - v1.x, v3.y - v1.y, v3.z - v1.z};
    
    Point3D normal = {
        u.y * v.z - u.z * v.y,
        u.z * v.x - u.x * v.z,
        u.x * v.y - u.y * v.x
    };
    
    float length = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    if (length > 0) {
        normal.x /= length;
        normal.y /= length;
        normal.z /= length;
    }
    
    return normal;
}

float distance(Point3D p1, Point3D p2) {
    float dx = p2.x - p1.x;
    float dy = p2.y - p1.y;
    float dz = p2.z - p1.z;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}
Point3D normalize(Point3D vector) {
    float length = std::sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
    if (length == 0) return {0, 0, 0};
    return {vector.x / length, vector.y / length, vector.z / length};
}

Point3D crossProduct(Point3D a, Point3D b) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

std::vector<Point3D> generateCirclePoints(float radius, int segments, float z) {
    std::vector<Point3D> points;
    for (int i = 0; i < segments; ++i) {
        float angle = 2.0f * M_PI * i / segments;
        points.push_back({radius * cos(angle), radius * sin(angle), z});
    }
    for (auto i : points) {
        std::cout << "Point: (" << i.x << ", " << i.y << ", " << i.z << ")\n";
    }
    return points;
}



std::vector<Point3D> generateGearProfile(GearParams params) {
    std::vector<Point3D> profile;
    
    float anglePerTooth = 2.0f * M_PI / params.numTeeth;
    
    float toothArcAngle = anglePerTooth * params.toothWidth;
    float valleyArcAngle = anglePerTooth * (1.0f - params.toothWidth);
    
    for (int tooth = 0; tooth < params.numTeeth; ++tooth) {
        float toothStartAngle = tooth * anglePerTooth;
        

        for (int i = 0; i < params.segments; ++i) {
            float angle = toothStartAngle + (toothArcAngle * i) / params.segments;
            profile.push_back({
                params.outerRadius * cos(angle),
                params.outerRadius * sin(angle),
                0 
            });
        }
        
        float valleyStartAngle = toothStartAngle + toothArcAngle;
        for (int i = 0; i < params.segments; ++i) {
            float angle = valleyStartAngle + (valleyArcAngle * i) / params.segments;
            profile.push_back({
                params.innerRadius * cos(angle),
                params.innerRadius * sin(angle),
                0
            });
        }
    }
    
    return profile;
}

void generateGear(std::ofstream &file, GearParams params) {
    std::vector<Point3D> bottomProfile = generateGearProfile(params);
    
    std::vector<Point3D> topProfile;
    for (const auto& point : bottomProfile) {
        topProfile.push_back({point.x, point.y, params.thickness});
    }
    
    int numPoints = bottomProfile.size();
    
    for (int i = 0; i < numPoints; ++i) {
        int next = (i + 1) % numPoints;
        
        Point3D b1 = bottomProfile[i];
        Point3D b2 = bottomProfile[next];
        Point3D t1 = topProfile[i];
        Point3D t2 = topProfile[next];
        Point3D normal1 = calculateNormal(b1, b2, t1);
        
        writeTriangle(file,
                      normal1.x, normal1.y, normal1.z,
                      b1.x, b1.y, b1.z,
                      b2.x, b2.y, b2.z,
                      t1.x, t1.y, t1.z);
        
        Point3D normal2 = calculateNormal(t1, b2, t2);
        writeTriangle(file,
                      normal2.x, normal2.y, normal2.z,
                      t1.x, t1.y, t1.z,
                      b2.x, b2.y, b2.z,
                      t2.x, t2.y, t2.z);
    }
    Point3D bottomCenter = {0, 0, 0};
    for (int i = 0; i < numPoints; ++i) {
        int next = (i + 1) % numPoints;
        writeTriangle(file,
                      0, 0, -1,
                      bottomCenter.x, bottomCenter.y, bottomCenter.z,
                      bottomProfile[next].x, bottomProfile[next].y, bottomProfile[next].z,
                      bottomProfile[i].x, bottomProfile[i].y, bottomProfile[i].z);
    }
    
    Point3D topCenter = {0, 0, params.thickness};
    for (int i = 0; i < numPoints; ++i) {
        int next = (i + 1) % numPoints;
        
        writeTriangle(file,
                      0, 0, 1,
                      topCenter.x, topCenter.y, topCenter.z,
                      topProfile[i].x, topProfile[i].y, topProfile[i].z,
                      topProfile[next].x, topProfile[next].y, topProfile[next].z);
    }
}
void generateGearWithHole(std::ofstream &file, GearParams params, float holeRadius) {
    std::vector<Point3D> bottomProfile = generateGearProfile(params);
    std::vector<Point3D> topProfile;
    for (const auto& point : bottomProfile) {
        topProfile.push_back({point.x, point.y, params.thickness});
    }
    
    int numPoints = bottomProfile.size();
    
    for (int i = 0; i < numPoints; ++i) {
        int next = (i + 1) % numPoints;
        
        Point3D b1 = bottomProfile[i];
        Point3D b2 = bottomProfile[next];
        Point3D t1 = topProfile[i];
        Point3D t2 = topProfile[next];
        
        Point3D normal1 = calculateNormal(b1, b2, t1);
        writeTriangle(file,
                      normal1.x, normal1.y, normal1.z,
                      b1.x, b1.y, b1.z,
                      b2.x, b2.y, b2.z,
                      t1.x, t1.y, t1.z);
        
        Point3D normal2 = calculateNormal(t1, b2, t2);
        writeTriangle(file,
                      normal2.x, normal2.y, normal2.z,
                      t1.x, t1.y, t1.z,
                      b2.x, b2.y, b2.z,
                      t2.x, t2.y, t2.z);
    }
    int holeSegments = 32;
    std::vector<Point3D> bottomHole = generateCirclePoints(holeRadius, holeSegments, 0);
    std::vector<Point3D> topHole = generateCirclePoints(holeRadius, holeSegments, params.thickness);
    
    for (int i = 0; i < holeSegments; ++i) {
        int next = (i + 1) % holeSegments;
        
        Point3D b1 = bottomHole[i];
        Point3D b2 = bottomHole[next];
        Point3D t1 = topHole[i];
        Point3D t2 = topHole[next];
        Point3D normal1 = calculateNormal(b1, t1, b2);
        writeTriangle(file,
                      normal1.x, normal1.y, normal1.z,
                      b1.x, b1.y, b1.z,
                      t1.x, t1.y, t1.z,
                      b2.x, b2.y, b2.z);
        
        Point3D normal2 = calculateNormal(t1, t2, b2);
        writeTriangle(file,
                      normal2.x, normal2.y, normal2.z,
                      t1.x, t1.y, t1.z,
                      t2.x, t2.y, t2.z,
                      b2.x, b2.y, b2.z);
    }
    
    for (int i = 0; i < numPoints; ++i) {
        int next = (i + 1) % numPoints;  
        int holeIdx = (i * holeSegments) / numPoints;
        int holeNext = (holeIdx + 1) % holeSegments;
        
        Point3D outerP1 = bottomProfile[i];
        Point3D outerP2 = bottomProfile[next];
        Point3D innerP1 = bottomHole[holeIdx];
        Point3D innerP2 = bottomHole[holeNext];
        
        writeTriangle(file,
                      0, 0, -1,
                      outerP1.x, outerP1.y, outerP1.z,
                      innerP1.x, innerP1.y, innerP1.z,
                      outerP2.x, outerP2.y, outerP2.z);
        
        writeTriangle(file,
                      0, 0, -1,  
                      outerP2.x, outerP2.y, outerP2.z,
                      innerP1.x, innerP1.y, innerP1.z,
                      innerP2.x, innerP2.y, innerP2.z);
    }
    
    for (int i = 0; i < numPoints; ++i) {
        int next = (i + 1) % numPoints;
        
        int holeIdx = (i * holeSegments) / numPoints;
        int holeNext = (holeIdx + 1) % holeSegments;
        
        Point3D outerP1 = topProfile[i];
        Point3D outerP2 = topProfile[next];
        Point3D innerP1 = topHole[holeIdx];
        Point3D innerP2 = topHole[holeNext];
        
        writeTriangle(file,
                      0, 0, 1, 
                      outerP1.x, outerP1.y, outerP1.z,
                      outerP2.x, outerP2.y, outerP2.z,
                      innerP1.x, innerP1.y, innerP1.z);
        
        writeTriangle(file,
                      0, 0, 1,
                      outerP2.x, outerP2.y, outerP2.z,
                      innerP2.x, innerP2.y, innerP2.z,
                      innerP1.x, innerP1.y, innerP1.z);
    }
}
void generateDisk(std::ofstream &file, 
                  float innerRadius, 
                  float outerRadius, 
                  float z, 
                  int segments) {
    
    std::vector<Point3D> innerCircle = generateCirclePoints(innerRadius, segments, z);
    std::vector<Point3D> outerCircle = generateCirclePoints(outerRadius, segments, z);
    
    for (int i = 0; i < segments; ++i) {
        int next = (i + 1) % segments;
        
        Point3D inner1 = innerCircle[i];
        Point3D inner2 = innerCircle[next];
        Point3D outer1 = outerCircle[i];
        Point3D outer2 = outerCircle[next];
        
        writeTriangle(file,
                      0, 0, 1, 
                      inner1.x, inner1.y, inner1.z,
                      outer1.x, outer1.y, outer1.z,
                      outer2.x, outer2.y, outer2.z);
        
        writeTriangle(file,
                      0, 0, 1, 
                      inner1.x, inner1.y, inner1.z,
                      outer2.x, outer2.y, outer2.z,
                      inner2.x, inner2.y, inner2.z);
    }
}

void generateCylinder(std::ofstream &file, 
                     float radius, 
                     float height, 
                     int segments) {
    
    std::vector<Point3D> bottomCircle = generateCirclePoints(radius, segments, 0);
    std::vector<Point3D> topCircle = generateCirclePoints(radius, segments, height);
    
    for (int i = 0; i < segments; ++i) {
        int next = (i + 1) % segments; 
        
        Point3D b1 = bottomCircle[i];    
        Point3D b2 = bottomCircle[next]; 
        Point3D t1 = topCircle[i];         
        Point3D t2 = topCircle[next];   
        Point3D normal1 = calculateNormal(b1, b2, t1);
        writeTriangle(file,
                      normal1.x, normal1.y, normal1.z,
                      b1.x, b1.y, b1.z,
                      b2.x, b2.y, b2.z,
                      t1.x, t1.y, t1.z);
        
        Point3D normal2 = calculateNormal(t1, b2, t2);
        writeTriangle(file,
                      normal2.x, normal2.y, normal2.z,
                      t1.x, t1.y, t1.z,
                      b2.x, b2.y, b2.z,
                      t2.x, t2.y, t2.z);
    }
    Point3D bottomCenter = {0, 0, 0};
    for (int i = 0; i < segments; ++i) {
        int next = (i + 1) % segments;
        
        Point3D p1 = bottomCircle[i];
        Point3D p2 = bottomCircle[next];

        writeTriangle(file,
                      0, 0, -1,
                      bottomCenter.x, bottomCenter.y, bottomCenter.z,
                      p2.x, p2.y, p2.z, 
                      p1.x, p1.y, p1.z);
    }
    
    Point3D topCenter = {0, 0, height};
    for (int i = 0; i < segments; ++i) {
        int next = (i + 1) % segments;
        
        Point3D p1 = topCircle[i];
        Point3D p2 = topCircle[next];
        writeTriangle(file,
                      0, 0, 1,
                      topCenter.x, topCenter.y, topCenter.z,
                      p1.x, p1.y, p1.z,
                      p2.x, p2.y, p2.z);
    }
}

float calculateToothWidth(float radius,
                          float baseWidth,
                          float tipWidth,
                          float currentRadius)
{
    if (currentRadius < 0) currentRadius = 0;
    if (currentRadius > radius) currentRadius = radius;
    float t = currentRadius / radius;

    // Linear interpolation
    return baseWidth + t * (tipWidth - baseWidth);
}


int main() {
    std::ofstream stlFile("cube.stl");
    if (!stlFile.is_open()) {
        std::cerr << "Error: Could not create STL file!" << std::endl;
        return 1;
    }
     stlFile << "solid gear\n";
    GearParams params;
    params.numTeeth = 8;      
    params.innerRadius = 30.0f;
    params.outerRadius = 50.0f;
    params.thickness = 10.0f;  
    params.toothWidth = 0.4f;  
    params.segments = 7;       
    
    std::cout << "Generating gear...\n";
    std::cout << "Number of teeth: " << params.numTeeth << "\n";
    std::cout << "Inner radius: " << params.innerRadius << "\n";
    std::cout << "Outer radius: " << params.outerRadius << "\n";
    std::cout << "Thickness: " << params.thickness << "\n";
    
    generateGear(stlFile, params);
    
    stlFile << "endsolid gear\n";
    stlFile.close();
    std::cout << "Gear created successfully! (gear.stl)\n\n";

    std::ofstream gearWithHoleFile("gearhole.stl");
    if (!gearWithHoleFile.is_open()) {
        std::cerr << "Error: Could not create gear_with_hole.stl!" << std::endl;
        return 1;
    }
    
    gearWithHoleFile << "solid gear_with_hole\n";
    
    std::cout << "Generating gear with mounting hole...\n";
    std::cout << "Hole radius: 10\n";
    
    generateGearWithHole(gearWithHoleFile, params, 10.0f);
    
    gearWithHoleFile << "endsolid gear_with_hole\n";
    gearWithHoleFile.close();
    std::cout << "Gear with hole created successfully! (gear_with_hole.stl)\n";
    
    return 0;
}