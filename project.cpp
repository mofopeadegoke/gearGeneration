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
    
    // Normalize the normal vector
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

void generateDisk(std::ofstream &file, 
                  float innerRadius, 
                  float outerRadius, 
                  float z, 
                  int segments) {
    
    // Step 1: Generate two circles at the same Z height
    std::vector<Point3D> innerCircle = generateCirclePoints(innerRadius, segments, z);
    std::vector<Point3D> outerCircle = generateCirclePoints(outerRadius, segments, z);
    
    // Step 2: Create the flat disk/ring surface with triangles
    // Connect inner and outer circles with 2 triangles per segment
    for (int i = 0; i < segments; ++i) {
        int next = (i + 1) % segments;  // Wrap around to 0 at the end
        
        Point3D inner1 = innerCircle[i];
        Point3D inner2 = innerCircle[next];
        Point3D outer1 = outerCircle[i];
        Point3D outer2 = outerCircle[next];
        
        // Triangle 1: inner-current → outer-current → outer-next
        // Counter-clockwise from above (normal points up)
        writeTriangle(file,
                      0, 0, 1,  // Normal points up
                      inner1.x, inner1.y, inner1.z,
                      outer1.x, outer1.y, outer1.z,
                      outer2.x, outer2.y, outer2.z);
        
        // Triangle 2: inner-current → outer-next → inner-next
        // Counter-clockwise from above (normal points up)
        writeTriangle(file,
                      0, 0, 1,  // Normal points up
                      inner1.x, inner1.y, inner1.z,
                      outer2.x, outer2.y, outer2.z,
                      inner2.x, inner2.y, inner2.z);
    }
}

void generateCylinder(std::ofstream &file, 
                     float radius, 
                     float height, 
                     int segments) {
    
    // Step 1: Generate circle points at bottom (z=0) and top (z=height)
    std::vector<Point3D> bottomCircle = generateCirclePoints(radius, segments, 0);
    std::vector<Point3D> topCircle = generateCirclePoints(radius, segments, height);
    
    // Step 2: Create side faces (connecting bottom and top circles)
    // Each segment creates 2 triangles forming a rectangle
    for (int i = 0; i < segments; ++i) {
        int next = (i + 1) % segments;  // Wrap around to 0 at the end
        
        Point3D b1 = bottomCircle[i];      // Bottom current
        Point3D b2 = bottomCircle[next];   // Bottom next
        Point3D t1 = topCircle[i];         // Top current
        Point3D t2 = topCircle[next];      // Top next
        
        // Triangle 1: bottom-current → bottom-next → top-current
        Point3D normal1 = calculateNormal(b1, b2, t1);
        writeTriangle(file,
                      normal1.x, normal1.y, normal1.z,
                      b1.x, b1.y, b1.z,
                      b2.x, b2.y, b2.z,
                      t1.x, t1.y, t1.z);
        
        // Triangle 2: top-current → bottom-next → top-next
        Point3D normal2 = calculateNormal(t1, b2, t2);
        writeTriangle(file,
                      normal2.x, normal2.y, normal2.z,
                      t1.x, t1.y, t1.z,
                      b2.x, b2.y, b2.z,
                      t2.x, t2.y, t2.z);
    }
    
    // Step 3: Create bottom cap (triangles from center to edge)
    Point3D bottomCenter = {0, 0, 0};
    for (int i = 0; i < segments; ++i) {
        int next = (i + 1) % segments;
        
        Point3D p1 = bottomCircle[i];
        Point3D p2 = bottomCircle[next];
        
        // Normal points down (negative Z)
        writeTriangle(file,
                      0, 0, -1,
                      bottomCenter.x, bottomCenter.y, bottomCenter.z,
                      p2.x, p2.y, p2.z,  // Note: reversed order for correct winding
                      p1.x, p1.y, p1.z);
    }
    
    // Step 4: Create top cap (triangles from center to edge)
    Point3D topCenter = {0, 0, height};
    for (int i = 0; i < segments; ++i) {
        int next = (i + 1) % segments;
        
        Point3D p1 = topCircle[i];
        Point3D p2 = topCircle[next];
        
        // Normal points up (positive Z)
        writeTriangle(file,
                      0, 0, 1,
                      topCenter.x, topCenter.y, topCenter.z,
                      p1.x, p1.y, p1.z,
                      p2.x, p2.y, p2.z);
    }
}

int main() {
    // Create and open the STL file
    std::ofstream stlFile("cube.stl");
    
    // Check if file was opened successfully
    if (!stlFile.is_open()) {
        std::cerr << "Error: Could not create STL file!" << std::endl;
        return 1;
    }
    
    // Write the STL header
    stlFile << "solid cube\n";


    generateDisk(stlFile, 5.0f, 10.0f, 0.0f, 36);


    // Write the STL footer
    stlFile << "endsolid cube\n";
    
    // Close the file
    stlFile.close();
    
    std::cout << "STL file 'cube.stl' created successfully!" << std::endl;
    
    return 0;
}