#pragma once
#include "JuceHeader.h"
#include <vector>

struct Vec3 { float x, y, z; };

struct STLTriangle {
    Vec3 normal;
    Vec3 v[3];
};

class STLParser {
public:
    static std::vector<STLTriangle> parseSTL(const juce::File& file) {
        std::vector<STLTriangle> triangles;
        if (!file.existsAsFile()) return triangles;
        
        juce::MemoryBlock block;
        if (!file.loadFileAsData(block)) return triangles;
        
        if (block.getSize() < 84) return triangles;
        
        const char* data = (const char*)block.getData();
        // Check if binary by counting triangles
        uint32_t numTriangles = *(uint32_t*)(data + 80);
        size_t expectedSize = 84 + numTriangles * 50;
        
        if (block.getSize() >= expectedSize && numTriangles < 10000000) {
            // Binary STL
            size_t offset = 84;
            for (uint32_t i = 0; i < numTriangles; ++i) {
                if (offset + 50 > block.getSize()) break;
                float* floats = (float*)(data + offset);
                STLTriangle t;
                t.normal = {floats[0], floats[1], floats[2]};
                t.v[0] = {floats[3], floats[4], floats[5]};
                t.v[1] = {floats[6], floats[7], floats[8]};
                t.v[2] = {floats[9], floats[10], floats[11]};
                triangles.push_back(t);
                offset += 50;
            }
        } else {
            // Very naive ASCII STL parse
            juce::StringArray lines;
            lines.addLines(file.loadFileAsString());
            STLTriangle t;
            int vertexCount = 0;
            for (int i = 0; i < lines.size(); ++i) {
                juce::String line = lines[i].trim();
                if (line.startsWithIgnoreCase("vertex")) {
                    auto tokens = juce::StringArray::fromTokens(line, " \t", "");
                    if (tokens.size() >= 4 && vertexCount < 3) {
                        t.v[vertexCount] = {tokens[1].getFloatValue(), tokens[2].getFloatValue(), tokens[3].getFloatValue()};
                        vertexCount++;
                        if (vertexCount == 3) {
                            triangles.push_back(t);
                            vertexCount = 0;
                        }
                    }
                }
            }
        }
        // Invert X and Y axes
        for (auto& t : triangles) {
            t.normal.x = -t.normal.x;
            t.normal.y = -t.normal.y;
            for (int i = 0; i < 3; ++i) {
                t.v[i].x = -t.v[i].x;
                t.v[i].y = -t.v[i].y;
            }
        }
        
        return triangles;
    }
};
