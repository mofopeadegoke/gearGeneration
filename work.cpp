#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#define M_PI 3.14159265358979323846

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
    float holeRadius;
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
    float pressureAngleRad = params.pressureAngle * M_PI / 180.0f;
    
    // Base circle radius (where involute curve begins)
    float baseRadius = params.pitchRadius * std::cos(pressureAngleRad);
    
    // Involute function: inv(α) = tan(α) - α
    auto involuteFunc = [](float angle) {
        return std::tan(angle) - angle;
    };
    
    // Tooth thickness at pitch circle (in radians)
    float pitchCircumference = 2.0f * M_PI * params.pitchRadius;
    float toothArcAtPitch = (pitchCircumference / params.numTeeth) * params.toothWidth;
    float toothAngleAtPitch = toothArcAtPitch / params.pitchRadius;
    
    // Involute parameter at pitch radius
    float alphaPitch = pressureAngleRad;
    float invAlphaPitch = involuteFunc(alphaPitch);
    
    // Function to get involute angle at any radius
    auto getInvoluteAngle = [&](float radius, float side) {
        if (radius < baseRadius) {
            // Below base circle: straight radial line
            return toothStartAngle + side * toothAngleAtPitch / 2.0f;
        } else {
            // Involute curve
            float alpha = std::acos(baseRadius / radius);
            float invAlpha = involuteFunc(alpha);
            // The key: tooth gets narrower as we go up
            float toothAngleAtRadius = toothAngleAtPitch * (params.pitchRadius / radius);
            return toothStartAngle + side * (toothAngleAtRadius / 2.0f + invAlphaPitch - invAlpha);
        }
    };
    
    int toothSegments = params.segments * 4;
    
    // Generate LEFT side of tooth (from root up to tip)
    for (int i = 0; i <= toothSegments; ++i) {
        float t = (float)i / toothSegments;
        float r = params.innerRadius + t * (params.outerRadius - params.innerRadius);
        float theta = getInvoluteAngle(r, -1.0f);  // Left side
        profile.push_back({r * std::cos(theta), r * std::sin(theta), 0});
    }
    
    // Generate RIGHT side of tooth (from tip down to root)
    for (int i = toothSegments; i >= 0; --i) {
        float t = (float)i / toothSegments;
        float r = params.innerRadius + t * (params.outerRadius - params.innerRadius);
        float theta = getInvoluteAngle(r, 1.0f);  // Right side
        profile.push_back({r * std::cos(theta), r * std::sin(theta), 0});
    }
    
    // Generate VALLEY (the gap/fillet at the root)
    float rightRootAngle = getInvoluteAngle(params.innerRadius, 1.0f);
    float leftNextToothAngle = getInvoluteAngle(params.innerRadius, -1.0f) + anglePerTooth;
    
    float gapAngle = leftNextToothAngle - rightRootAngle;
    int valleySegments = params.segments * 2;
    
    for (int i = 1; i < valleySegments; ++i) {
        float t = (float)i / valleySegments;
        float theta = rightRootAngle + t * gapAngle;
        profile.push_back({
            params.innerRadius * std::cos(theta),
            params.innerRadius * std::sin(theta),
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
    float threadDepth = params.outerRadius - params.innerRadius;
    
    for (int i = 0; i < segments; ++i) {
        float t = (float)i / segments;
        float angle = t * 2.0f * M_PI;
        
        float radius;
        float phase = std::fmod(t * params.numTeeth, 1.0f);
        
        if (phase < 0.3f) {
            radius = params.innerRadius + (phase / 0.3f) * threadDepth;
        } else if (phase < 0.7f) {
            radius = params.outerRadius;
        } else {
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

void generateGear(std::ofstream &file, GearParams params) {
    std::vector<Point3D> bottomProfile = generateGearProfile(params);
    std::vector<Point3D> topProfile;
    
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
    
    if (params.holeRadius > 0) {
        int holeSegments = 32;
        std::vector<Point3D> bottomHole = generateCirclePoints(params.holeRadius, holeSegments, 0);
        std::vector<Point3D> topHole = generateCirclePoints(params.holeRadius, holeSegments, params.thickness);
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
    std::cout << "Generating gears with different pressure angles...\n\n";
    
    // 1. Standard 20° pressure angle
    std::ofstream gear20File("gear_pressure_20.stl");
    gear20File << "solid gear_20_degrees\n";
    GearParams gear20 = {20, 35.0f, 45.0f, 50.0f, 10.0f, 0.45f, 4, SPUR, 0, 20.0f, true, false, 2.0f, 8.0f};
    generateGear(gear20File, gear20);
    gear20File << "endsolid gear_20_degrees\n";
    gear20File.close();
    std::cout << "✓ gear_pressure_20.stl (Standard - 20° pressure angle)\n";
    
    // 2. Low 14.5° pressure angle
    std::ofstream gear14File("gear_pressure_14.stl");
    gear14File << "solid gear_14_degrees\n";
    GearParams gear14 = {20, 35.0f, 45.0f, 50.0f, 10.0f, 0.45f, 4, SPUR, 0, 14.5f, true, false, 2.0f, 8.0f};
    generateGear(gear14File, gear14);
    gear14File << "endsolid gear_14_degrees\n";
    gear14File.close();
    std::cout << "✓ gear_pressure_14.stl (Low - 14.5° pressure angle)\n";
    
    // 3. High 25° pressure angle
    std::ofstream gear25File("gear_pressure_25.stl");
    gear25File << "solid gear_25_degrees\n";
    GearParams gear25 = {20, 35.0f, 45.0f, 50.0f, 10.0f, 0.45f, 4, SPUR, 0, 25.0f, true, false, 2.0f, 8.0f};
    generateGear(gear25File, gear25);
    gear25File << "endsolid gear_25_degrees\n";
    gear25File.close();
    std::cout << "✓ gear_pressure_25.stl (High - 25° pressure angle)\n";

       // 4. Basic SPUR gear
    std::ofstream spurFile("gear_spur3.stl");
    spurFile << "solid spur_gear\n";
    GearParams spur = {20, 35.0f, 45.0f, 50.0f, 10.0f, 0.4f, 3, SPUR, 0, 20.0f, false, false, 2.0f, 0};
    generateGear(spurFile, spur);
    spurFile << "endsolid spur_gear\n";
    spurFile.close();
    std::cout << "✓ gear_spur.stl (Basic rectangular teeth)\n";
    
    // 5. Rounded teeth
    std::ofstream roundedFile("gear_rounded3.stl");
    roundedFile << "solid rounded_gear\n";
    GearParams rounded = {16, 30.0f, 40.0f, 45.0f, 12.0f, 0.45f, 4, SPUR, 0, 20.0f, false, true, 1.5f, 8.0f};
    generateGear(roundedFile, rounded);
    roundedFile << "endsolid rounded_gear\n";
    roundedFile.close();
    std::cout << "✓ gear_rounded.stl (Rounded teeth)\n";
    
    // 6. Helical gear
    std::ofstream helicalFile("gear_helical3.stl");
    helicalFile << "solid helical_gear\n";
    GearParams helical = {24, 38.0f, 48.0f, 53.0f, 20.0f, 0.45f, 4, HELICAL, 30.0f, 20.0f, true, false, 2.0f, 12.0f};
    generateGear(helicalFile, helical);
    helicalFile << "endsolid helical_gear\n";
    helicalFile.close();
    std::cout << "✓ gear_helical.stl (Helical with 20° pressure angle)\n";
    
    // 7. Bevel gear
    std::ofstream bevelFile("gear_bevel3.stl");
    bevelFile << "solid bevel_gear\n";
    GearParams bevel = {20, 35.0f, 45.0f, 50.0f, 15.0f, 0.4f, 3, BEVEL, 0, 20.0f, false, false, 2.0f, 0};
    generateGear(bevelFile, bevel);
    bevelFile << "endsolid bevel_gear\n";
    bevelFile.close();
    std::cout << "✓ gear_bevel.stl (Conical gear)\n";
    
    // 8. Worm gear
    std::ofstream wormFile("gear_worm3.stl");
    wormFile << "solid worm_gear\n";
    GearParams worm = {3, 8.0f, 12.0f, 15.0f, 80.0f, 0.3f, 8, WORM, 0, 20.0f, false, false, 2.0f, 5.0f};
    generateGear(wormFile, worm);
    wormFile << "endsolid worm_gear\n";
    wormFile.close();
    std::cout << " gear_worm.stl (Screw-type gear)\n";
    std::cout << "\nPressure Angle Comparison:\n";
    std::cout << "- gear_pressure_14.stl: 14.5° (smoother engagement, wider teeth)\n";
    std::cout << "- gear_pressure_20.stl: 20°   (standard, balanced) \n";
    std::cout << "- gear_pressure_25.stl: 25°   (stronger, narrower teeth)\n";
    
    return 0;
}