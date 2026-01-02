#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

enum GearType {
    SPUR,
    HELICAL,
    BEVEL,
    WORM
};

struct GearParams {
    int numTeeth;
    float innerRadius;
    float pitchRadius;
    float outerRadius;
    float thickness;
    float toothWidth;
    int segments;
    GearType type;
    float helixAngle;
    float pressureAngle;
    bool useInvolute;
    bool roundedTeeth;
    float roundingRadius;
    float holeRadius;  // NEW: 0 = no hole, >0 = hole with this radius
};

struct Point3D {
    float x, y, z;
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

std::vector<Point3D> generateCirclePoints(float radius, int segments, float z) {
    std::vector<Point3D> points;
    points.reserve(segments);
    for (int i = 0; i < segments; ++i) {
        float angle = 2.0f * M_PI * i / segments;
        points.push_back({radius * std::cos(angle), radius * std::sin(angle), z});
    }
    return points;
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
    if (length > 1e-6f) {
        normal.x /= length;
        normal.y /= length;
        normal.z /= length;
    }
    return normal;
}

Point3D rotateZ(Point3D point, float angle) {
    float cosA = std::cos(angle);
    float sinA = std::sin(angle);
    return {
        point.x * cosA - point.y * sinA,
        point.x * sinA + point.y * cosA,
        point.z
    };
}

Point3D applyHelixAngle(Point3D point, float helixAngle, float thickness) {
    float rotationAngle = (point.z / thickness) * helixAngle * M_PI / 180.0f;
    return rotateZ(point, rotationAngle);
}

std::vector<Point3D> generateInvoluteToothProfile(GearParams params, float toothStartAngle) {
    std::vector<Point3D> profile;
    float anglePerTooth = 2.0f * M_PI / params.numTeeth;
    float toothArcAngle = anglePerTooth * params.toothWidth;
    float valleyArcAngle = anglePerTooth * (1.0f - params.toothWidth);
    int toothSegments = params.segments * 2;
    for (int i = 0; i < toothSegments; ++i) {
        float t = (float)i / toothSegments;
        float angle = toothStartAngle - toothArcAngle/2 + t * toothArcAngle;
        
        float radiusAdjust = params.outerRadius - (params.outerRadius - params.pitchRadius) * 0.1f * std::sin(t * M_PI);
        
        profile.push_back({
            radiusAdjust * std::cos(angle),
            radiusAdjust * std::sin(angle),
            0
        });
    }
    
    float valleyStart = toothStartAngle + toothArcAngle/2;
    int valleySegments = params.segments;
    for (int i = 0; i < valleySegments; ++i) {
        float t = (float)i / valleySegments;
        float angle = valleyStart + t * valleyArcAngle;
        profile.push_back({
            params.innerRadius * std::cos(angle),
            params.innerRadius * std::sin(angle),
            0
        });
    }
    
    return profile;
}

std::vector<Point3D> generateRoundedToothProfile(GearParams params, float toothStartAngle) {
    std::vector<Point3D> profile;
    float anglePerTooth = 2.0f * M_PI / params.numTeeth;
    float toothArcAngle = anglePerTooth * params.toothWidth;
    float valleyArcAngle = anglePerTooth * (1.0f - params.toothWidth);
    float roundRad = params.roundingRadius;
   
    int toothSegments = params.segments * 2;
    for (int i = 0; i < toothSegments; ++i) {
        float t = (float)i / toothSegments;
        float angle = toothStartAngle - toothArcAngle/2 + t * toothArcAngle;
        float roundingFactor = 1.0f;
        if (t < 0.2f) {
            roundingFactor = 0.9f + 0.1f * (t / 0.2f);
        } else if (t > 0.8f) {
            roundingFactor = 0.9f + 0.1f * ((1.0f - t) / 0.2f);
        }
        
        float radius = params.outerRadius * roundingFactor;
        
        profile.push_back({
            radius * std::cos(angle),
            radius * std::sin(angle),
            0
        });
    }
    
    float valleyStart = toothStartAngle + toothArcAngle/2;
    int valleySegments = params.segments;
    for (int i = 0; i < valleySegments; ++i) {
        float t = (float)i / valleySegments;
        float angle = valleyStart + t * valleyArcAngle;
        float filletFactor = 1.0f;
        if (t < 0.3f || t > 0.7f) {
            filletFactor = 1.02f;
        }
        
        profile.push_back({
            params.innerRadius * filletFactor * std::cos(angle),
            params.innerRadius * filletFactor * std::sin(angle),
            0
        });
    }
    
    return profile;
}

std::vector<Point3D> generateWormProfile(GearParams params) {
    std::vector<Point3D> profile;
    int segments = params.segments * 4;
    
    // Create a trapezoidal thread profile
    float threadDepth = params.outerRadius - params.innerRadius;
    
    for (int i = 0; i < segments; ++i) {
        float t = (float)i / segments;
        float angle = t * 2.0f * M_PI;
        
        // Trapezoidal wave for thread
        float radius;
        float phase = std::fmod(t * params.numTeeth, 1.0f);
        
        if (phase < 0.3f) {
            // Rising edge
            radius = params.innerRadius + (phase / 0.3f) * threadDepth;
        } else if (phase < 0.7f) {
            // Top of thread
            radius = params.outerRadius;
        } else {
            // Falling edge
            radius = params.outerRadius - ((phase - 0.7f) / 0.3f) * threadDepth;
        }
        
        profile.push_back({
            radius * std::cos(angle),
            radius * std::sin(angle),
            0
        });
    }
    
    return profile;
}

std::vector<Point3D> generateGearProfile(GearParams params) {
    std::vector<Point3D> profile;
    
    // WORM gears use different profile generation
    if (params.type == WORM) {
        return generateWormProfile(params);
    }
    
    float anglePerTooth = 2.0f * M_PI / params.numTeeth;
    
    for (int tooth = 0; tooth < params.numTeeth; ++tooth) {
        float toothStartAngle = tooth * anglePerTooth;
        std::vector<Point3D> toothProfile;
        
        if (params.useInvolute) {
            toothProfile = generateInvoluteToothProfile(params, toothStartAngle);
        } else if (params.roundedTeeth) {
            toothProfile = generateRoundedToothProfile(params, toothStartAngle);
        } else {
            // Basic rectangular teeth
            float toothArcAngle = anglePerTooth * params.toothWidth;
            float valleyArcAngle = anglePerTooth * (1.0f - params.toothWidth);
            
            for (int i = 0; i < params.segments; ++i) {
                float angle = toothStartAngle + (toothArcAngle * i) / params.segments;
                profile.push_back({
                    params.outerRadius * std::cos(angle),
                    params.outerRadius * std::sin(angle),
                    0
                });
            }
            
            float valleyStartAngle = toothStartAngle + toothArcAngle;
            for (int i = 0; i < params.segments; ++i) {
                float angle = valleyStartAngle + (valleyArcAngle * i) / params.segments;
                profile.push_back({
                    params.innerRadius * std::cos(angle),
                    params.innerRadius * std::sin(angle),
                    0
                });
            }
            continue;
        }
        
        profile.insert(profile.end(), toothProfile.begin(), toothProfile.end());
    }
    
    return profile;
}

// NEW: Unified gear generation with optional hole
void generateGear(std::ofstream &file, GearParams params) {
    std::vector<Point3D> bottomProfile = generateGearProfile(params);
    std::vector<Point3D> topProfile;
    
    // Create top profile based on gear type
    if (params.type == HELICAL) {
        for (const auto& point : bottomProfile) {
            Point3D topPoint = {point.x, point.y, params.thickness};
            topPoint = applyHelixAngle(topPoint, params.helixAngle, params.thickness);
            topProfile.push_back(topPoint);
        }
    } else if (params.type == BEVEL) {
        float taperRatio = 0.7f;
        for (const auto& point : bottomProfile) {
            topProfile.push_back({
                point.x * taperRatio,
                point.y * taperRatio,
                params.thickness
            });
        }
    } else if (params.type == WORM) {
        // Helical twist for worm
        float totalRotation = params.thickness / (params.outerRadius * 2) * 2.0f * M_PI * params.numTeeth;
        for (size_t i = 0; i < bottomProfile.size(); ++i) {
            const auto& point = bottomProfile[i];
            float t = (float)i / bottomProfile.size();
            float twist = t * totalRotation;
            Point3D topPoint = {point.x, point.y, params.thickness};
            topPoint = rotateZ(topPoint, twist);
            topProfile.push_back(topPoint);
        }
    } else {
        for (const auto& point : bottomProfile) {
            topProfile.push_back({point.x, point.y, params.thickness});
        }
    }
    
    int numPoints = bottomProfile.size();
    
    // Side faces
    for (int i = 0; i < numPoints; ++i) {
        int next = (i + 1) % numPoints;
        
        Point3D b1 = bottomProfile[i];
        Point3D b2 = bottomProfile[next];
        Point3D t1 = topProfile[i];
        Point3D t2 = topProfile[next];
        
        Point3D normal1 = calculateNormal(b1, b2, t1);
        writeTriangle(file, normal1.x, normal1.y, normal1.z,
                      b1.x, b1.y, b1.z, b2.x, b2.y, b2.z, t1.x, t1.y, t1.z);
        
        Point3D normal2 = calculateNormal(t1, b2, t2);
        writeTriangle(file, normal2.x, normal2.y, normal2.z,
                      t1.x, t1.y, t1.z, b2.x, b2.y, b2.z, t2.x, t2.y, t2.z);
    }
    
    // Handle hole or solid caps
    if (params.holeRadius > 0) {
        // Generate hole
        int holeSegments = 32;
        std::vector<Point3D> bottomHole = generateCirclePoints(params.holeRadius, holeSegments, 0);
        std::vector<Point3D> topHole = generateCirclePoints(params.holeRadius, holeSegments, params.thickness);
        
        // Hole inner walls (inverted normals)
        for (int i = 0; i < holeSegments; ++i) {
            int next = (i + 1) % holeSegments;
            Point3D b1 = bottomHole[i];
            Point3D b2 = bottomHole[next];
            Point3D t1 = topHole[i];
            Point3D t2 = topHole[next];
            
            Point3D normal1 = calculateNormal(b1, t1, b2);
            writeTriangle(file, normal1.x, normal1.y, normal1.z,
                          b1.x, b1.y, b1.z, t1.x, t1.y, t1.z, b2.x, b2.y, b2.z);
            
            Point3D normal2 = calculateNormal(t1, t2, b2);
            writeTriangle(file, normal2.x, normal2.y, normal2.z,
                          t1.x, t1.y, t1.z, t2.x, t2.y, t2.z, b2.x, b2.y, b2.z);
        }
        
        // Bottom ring
        for (int i = 0; i < numPoints; ++i) {
            int next = (i + 1) % numPoints;
            int holeIdx = (i * holeSegments) / numPoints;
            int holeNext = (holeIdx + 1) % holeSegments;
            
            Point3D outerP1 = bottomProfile[i];
            Point3D outerP2 = bottomProfile[next];
            Point3D innerP1 = bottomHole[holeIdx];
            Point3D innerP2 = bottomHole[holeNext];
            
            writeTriangle(file, 0, 0, -1,
                          outerP1.x, outerP1.y, outerP1.z,
                          innerP1.x, innerP1.y, innerP1.z,
                          outerP2.x, outerP2.y, outerP2.z);
            
            writeTriangle(file, 0, 0, -1,
                          outerP2.x, outerP2.y, outerP2.z,
                          innerP1.x, innerP1.y, innerP1.z,
                          innerP2.x, innerP2.y, innerP2.z);
        }
        
        // Top ring
        for (int i = 0; i < numPoints; ++i) {
            int next = (i + 1) % numPoints;
            int holeIdx = (i * holeSegments) / numPoints;
            int holeNext = (holeIdx + 1) % holeSegments;
            
            Point3D outerP1 = topProfile[i];
            Point3D outerP2 = topProfile[next];
            Point3D innerP1 = topHole[holeIdx];
            Point3D innerP2 = topHole[holeNext];
            
            writeTriangle(file, 0, 0, 1,
                          outerP1.x, outerP1.y, outerP1.z,
                          outerP2.x, outerP2.y, outerP2.z,
                          innerP1.x, innerP1.y, innerP1.z);
            
            writeTriangle(file, 0, 0, 1,
                          outerP2.x, outerP2.y, outerP2.z,
                          innerP2.x, innerP2.y, innerP2.z,
                          innerP1.x, innerP1.y, innerP1.z);
        }
    } else {
        // Solid caps (no hole)
        Point3D bottomCenter = {0, 0, 0};
        for (int i = 0; i < numPoints; ++i) {
            int next = (i + 1) % numPoints;
            writeTriangle(file, 0, 0, -1,
                          bottomCenter.x, bottomCenter.y, bottomCenter.z,
                          bottomProfile[next].x, bottomProfile[next].y, bottomProfile[next].z,
                          bottomProfile[i].x, bottomProfile[i].y, bottomProfile[i].z);
        }
        
        Point3D topCenter = {0, 0, params.thickness};
        if (params.type == BEVEL) {
            topCenter.x *= 0.7f;
            topCenter.y *= 0.7f;
        }
        
        for (int i = 0; i < numPoints; ++i) {
            int next = (i + 1) % numPoints;
            writeTriangle(file, 0, 0, 1,
                          topCenter.x, topCenter.y, topCenter.z,
                          topProfile[i].x, topProfile[i].y, topProfile[i].z,
                          topProfile[next].x, topProfile[next].y, topProfile[next].z);
        }
    }
}

int main() {
    // 1. Basic SPUR gear
    std::ofstream spurFile("gear_spur2.stl");
    spurFile << "solid spur_gear\n";
    GearParams spur = {20, 35.0f, 45.0f, 50.0f, 10.0f, 0.4f, 3, SPUR, 0, 20.0f, false, false, 2.0f, 0};
    generateGear(spurFile, spur);
    spurFile << "endsolid spur_gear\n";
    spurFile.close();
    std::cout << "✓ gear_spur.stl\n";
    
    // 2. SPUR with hole
    std::ofstream spurHoleFile("gear_spur_hole2.stl");
    spurHoleFile << "solid spur_gear_hole\n";
    spur.holeRadius = 10.0f;
    generateGear(spurHoleFile, spur);
    spurHoleFile << "endsolid spur_gear_hole\n";
    spurHoleFile.close();
    std::cout << "✓ gear_spur_hole.stl\n";
    
    // 3. Rounded teeth
    std::ofstream roundedFile("gear_rounded2.stl");
    roundedFile << "solid rounded_gear\n";
    GearParams rounded = {16, 30.0f, 40.0f, 45.0f, 12.0f, 0.45f, 4, SPUR, 0, 20.0f, false, true, 1.5f, 8.0f};
    generateGear(roundedFile, rounded);
    roundedFile << "endsolid rounded_gear\n";
    roundedFile.close();
    std::cout << "✓ gear_rounded.stl\n";
    
    // 4. Involute teeth
    std::ofstream involuteFile("gear_involute2.stl");
    involuteFile << "solid involute_gear\n";
    GearParams involute = {18, 32.0f, 42.0f, 47.0f, 10.0f, 0.5f, 4, SPUR, 0, 20.0f, true, false, 2.0f, 10.0f};
    generateGear(involuteFile, involute);
    involuteFile << "endsolid involute_gear\n";
    involuteFile.close();
    std::cout << "✓ gear_involute.stl\n";
    
    // 5. Helical gear
    std::ofstream helicalFile("gear_helical2.stl");
    helicalFile << "solid helical_gear\n";
    GearParams helical = {24, 38.0f, 48.0f, 53.0f, 20.0f, 0.4f, 3, HELICAL, 30.0f, 20.0f, false, false, 2.0f, 12.0f};
    generateGear(helicalFile, helical);
    helicalFile << "endsolid helical_gear\n";
    helicalFile.close();
    std::cout << "✓ gear_helical.stl\n";
    
    // 6. Bevel gear
    std::ofstream bevelFile("gear_bevel2.stl");
    bevelFile << "solid bevel_gear\n";
    GearParams bevel = {20, 35.0f, 45.0f, 50.0f, 15.0f, 0.4f, 3, BEVEL, 0, 20.0f, false, false, 2.0f, 0};
    generateGear(bevelFile, bevel);
    bevelFile << "endsolid bevel_gear\n";
    bevelFile.close();
    std::cout << "✓ gear_bevel.stl\n";
    
    // 7. Worm gear (FIXED)
    std::ofstream wormFile("gear_worm2.stl");
    wormFile << "solid worm_gear\n";
    GearParams worm = {3, 8.0f, 12.0f, 15.0f, 80.0f, 0.3f, 8, WORM, 0, 20.0f, false, false, 2.0f, 5.0f};
    generateGear(wormFile, worm);
    wormFile << "endsolid worm_gear\n";
    wormFile.close();
    std::cout << "✓ gear_worm.stl\n";
    
    std::cout << "\n✅ All gears generated!\n";
    
    return 0;
}