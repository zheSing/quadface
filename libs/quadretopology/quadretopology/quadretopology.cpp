/***************************************************************************/
/* Copyright(C) 2021


The authors of

Reliable Feature-Line Driven Quad-Remeshing
Siggraph 2021


 All rights reserved.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
****************************************************************************/

#include "quadretopology.h"

#include "includes/qr_convert.h"
#include "includes/qr_utils.h"
#include "includes/qr_patterns.h"
#include "includes/qr_mapping.h"
#include <map>

#include <vcg/complex/algorithms/polygonal_algorithms.h>

#ifdef QUADRETOPOLOGY_DEBUG_SAVE_MESHES
#include <igl/writeOBJ.h>
#endif

namespace QuadRetopology {

template<class TriangleMesh>
ChartData computeChartData(
        TriangleMesh& mesh,
        const std::vector<std::vector<size_t>>& meshPartitions,
        const std::vector<std::vector<size_t>>& meshCorners)
{
    std::vector<int> faceLabel(mesh.face.size(), -1);
    for (size_t pId = 0; pId < meshPartitions.size(); pId++) {
        for (const size_t& fId : meshPartitions[pId]) {
            assert(faceLabel[fId] == -1);
            faceLabel[fId] = static_cast<int>(pId);
        }
    }

    ChartData chartData = computeChartData(mesh, faceLabel, meshCorners);
    return chartData;
}

//It works just on triangle meshes
template<class TriangleMeshType>
ChartData computeChartData(
        TriangleMeshType& mesh,
        const std::vector<int>& faceLabel,
        const std::vector<std::vector<size_t>>& corners)
{
    typedef std::map<std::pair<size_t, size_t>, std::pair<int, int>> EdgeLabelMap;
    typedef std::map<std::pair<size_t, size_t>, int> EdgeSubSideMap;

    ChartData chartData;

    if (mesh.face.size() == 0)
        return chartData;

    vcg::tri::UpdateTopology<TriangleMeshType>::FaceFace(mesh);

    //Region growing algorithm for getting charts
    internal::findChartFacesAndBorderFaces(mesh, faceLabel, chartData);

    //TODO SPLIT IN FUNCTIONS
    EdgeSubSideMap edgeSubSideMap;
    for (const int& pId : chartData.labels) {
        Chart& chart = chartData.charts[pId];

        assert(chart.label == pId);

        if (chart.faces.size() == 0)
            continue;

        std::unordered_set<size_t> cornerSet(corners[pId].begin(), corners[pId].end());

#ifndef NDEBUG
        if (cornerSet.size() < 3 || cornerSet.size() > 6) {
            std::cout << "Warning 3: Given as input for " << pId << ": " << cornerSet.size() << " sides." << std::endl;
        }
#endif

        EdgeLabelMap edgeLabelMap;
        std::vector<std::vector<size_t>> vertexNextMap(mesh.vert.size());

        std::set<size_t> remainingVertices;

        //Fill edge map and next vertex map
        for (const size_t& fId : chart.borderFaces) {
            typename TriangleMeshType::FaceType* currentFacePointer = &mesh.face[fId];
            vcg::face::Pos<typename TriangleMeshType::FaceType> pos(currentFacePointer, 0);

            for (int k = 0; k < currentFacePointer->VN(); k++) {
                pos.FlipF();
                size_t adjFace = vcg::tri::Index(mesh, pos.F());
                int adjLabel = faceLabel[adjFace];

                bool isBorderEdge = false;
                int adjChartLabel = -2;

                if (currentFacePointer == pos.F()) {
                    adjChartLabel = -1;
                    isBorderEdge = true;
                }
                else if (adjLabel != chart.label) {
                    adjChartLabel = adjLabel;
                    isBorderEdge = true;
                }
                pos.FlipF();

                //For each border edge
                if (isBorderEdge) {
                    assert(adjChartLabel > -2);

                    typename TriangleMeshType::VertexType* vStart = pos.V();
                    pos.FlipV();
                    typename TriangleMeshType::VertexType* vEnd = pos.V();
                    pos.FlipV();

                    size_t vStartId = vcg::tri::Index(mesh, vStart);
                    size_t vEndId = vcg::tri::Index(mesh, vEnd);

                    std::pair<size_t, size_t> edge(vStartId, vEndId);
                    if (edge.first > edge.second) {
                        std::swap(edge.first, edge.second);
                    }

                    edgeLabelMap.insert(std::make_pair(edge, std::make_pair(chart.label, adjChartLabel)));
                    vertexNextMap[vStartId].push_back(vEndId);

                    remainingVertices.insert(vStartId);
                    remainingVertices.insert(vEndId);
                }

                pos.FlipV();
                pos.FlipE();
            }
        }

        do {
            //Find first label
            size_t vStartId;
            size_t vCurrentId;
            size_t vNextId;

            //Corner detection variables
            typename TriangleMeshType::CoordType lastEdgeVec;
            bool isCorner = false;

            vCurrentId = *remainingVertices.begin();

            std::vector<size_t> nextConfiguration = internal::findVertexChainPath(vCurrentId, vertexNextMap);
            vNextId = vertexNextMap[vCurrentId][nextConfiguration[vCurrentId]];

            //Get last edge vector
            lastEdgeVec = mesh.vert[vNextId].P() - mesh.vert[vCurrentId].P();
            lastEdgeVec.Normalize();

//            std::pair<size_t, size_t> startEdge(vCurrentId, vNextId);
//            if (startEdge.first > startEdge.second) {
//                std::swap(startEdge.first, startEdge.second);
//            }

            int currentLabel;

            //Iterate in the borders to get the first corner
            vStartId = vCurrentId;
            size_t firstCornerIterations = 0;
            do {
                //Next border edge
                vCurrentId = vertexNextMap[vCurrentId][nextConfiguration[vCurrentId]];
                vNextId = vertexNextMap[vCurrentId][nextConfiguration[vCurrentId]];

                typename TriangleMeshType::CoordType currentEdgeVec = mesh.vert[vNextId].P() - mesh.vert[vCurrentId].P();
                currentEdgeVec.Normalize();

                //Check if it is a corner
                isCorner = cornerSet.find(vCurrentId) != cornerSet.end();

                lastEdgeVec = currentEdgeVec;

                firstCornerIterations++;
            } while (!isCorner && vCurrentId != vStartId && firstCornerIterations < MAXITERATIONS);

#ifndef NDEBUG
            if (firstCornerIterations >= MAXITERATIONS) {
                std::cout << "Error: error iterating! Cannot find the first corner or get back to the start vertex." << std::endl;
            }
#endif

#ifndef NDEBUG
            if (vCurrentId == vStartId) {
                std::cout << "Warning 1: input mesh is not well-defined: no corners!" << std::endl;
            }
#endif
            vStartId = vCurrentId;
            vNextId = vertexNextMap[vCurrentId][nextConfiguration[vCurrentId]];

            ChartSide currentSide;
            size_t chartSideId = 0;
            currentSide.length = 0;
            currentSide.size = 0;


            int adjChartLabel;
            do {
                size_t subsideId = chartData.subsides.size();
                ChartSubside currentSubSide;

                //Get edge
                std::pair<size_t, size_t> edge(vCurrentId, vNextId);
                if (edge.first > edge.second) {
                    std::swap(edge.first, edge.second);
                }
                //Get current label on the other side
                const std::pair<int,int>& currentEdgeLabels = edgeLabelMap.at(edge);
                assert(currentEdgeLabels.first == chart.label || currentEdgeLabels.second == chart.label);
                adjChartLabel = currentEdgeLabels.first == chart.label ? currentEdgeLabels.second : currentEdgeLabels.first;

                std::unordered_set<size_t> cornerSetAdj;
                if (adjChartLabel >= 0)
                    cornerSetAdj.insert(corners[adjChartLabel].begin(), corners[adjChartLabel].end());

                double length = 0;

                bool newSubSide = false;

                bool firstIteration = true;
                isCorner = false;
                bool isAdjCorner = false;
                size_t vSubSideStartId = vCurrentId;
                size_t iterations = 0;
                do {
                    typename TriangleMeshType::CoordType currentEdgeVec = mesh.vert[vNextId].P() - mesh.vert[vCurrentId].P();
                    currentEdgeVec.Normalize();

                    std::pair<size_t, size_t> edge(vCurrentId, vNextId);
                    if (edge.first > edge.second) {
                        std::swap(edge.first, edge.second);
                    }

                    //Check if it is a corner
                    if (!firstIteration) {
                        isCorner = cornerSet.find(vCurrentId) != cornerSet.end();
                        isAdjCorner = cornerSetAdj.find(vCurrentId) != cornerSetAdj.end();
                    }

                    //Get current label on the other subside
                    const std::pair<int,int>& currentEdgeLabels = edgeLabelMap.at(edge);
                    assert(currentEdgeLabels.first == chart.label || currentEdgeLabels.second == chart.label);
                    currentLabel = currentEdgeLabels.first == chart.label ? currentEdgeLabels.second : currentEdgeLabels.first;

                    if (!isCorner && !isAdjCorner && currentLabel == adjChartLabel) {
                        EdgeSubSideMap::iterator findIt = edgeSubSideMap.find(edge);

                        //If the subside has already been processed
                        if (findIt == edgeSubSideMap.end()) {
                            currentSubSide.vertices.push_back(vCurrentId);

                            length += (mesh.vert[vNextId].P() - mesh.vert[vCurrentId].P()).Norm();

                            edgeSubSideMap.insert(std::make_pair(edge, subsideId));

                            newSubSide = true;
                        }
                        else if (firstIteration) {
                            subsideId = findIt->second;
                        }
                        firstIteration = false;

                        remainingVertices.erase(vCurrentId);

                        //Next border edge
                        vCurrentId = vertexNextMap[vCurrentId][nextConfiguration[vCurrentId]];
                        vNextId = vertexNextMap[vCurrentId][nextConfiguration[vCurrentId]];

                        lastEdgeVec = currentEdgeVec;
                    }
                } while (!isCorner && !isAdjCorner && currentLabel == adjChartLabel && vCurrentId != vSubSideStartId && iterations < MAXITERATIONS);
#ifndef NDEBUG
                if (iterations >= MAXITERATIONS) {
                    std::cout << "Error: error iterating! Cannot find a corner or get back to the start vertex." << std::endl;
                }
#endif
#ifndef NDEBUG
                if (vCurrentId == vSubSideStartId) {
                    std::cout << "Warning 2: input mesh is not well-defined: single border chart with no corners!" << std::endl;
                }
#endif

                //True if the subside is reversed (from the last to the first vertex)
                bool reversed;

                if (newSubSide) {
                    //Add last vertex
                    currentSubSide.vertices.push_back(vCurrentId);

                    //Create new side
                    chart.chartSubsides.push_back(subsideId);

                    currentSubSide.incidentCharts[0] = chart.label;
                    currentSubSide.incidentCharts[1] = adjChartLabel;

                    currentSubSide.length = length;

                    currentSubSide.incidentChartSubsideId[0] = subsideId;
                    currentSubSide.incidentChartSideId[0] = chartSideId;

                    if (adjChartLabel >= 0) {
                        currentSubSide.isOnBorder = false;

                        chart.adjacentCharts.push_back(adjChartLabel);
                    }
                    else {
                        currentSubSide.isOnBorder = true;
                        currentSubSide.incidentChartSubsideId[1] = -1;
                        currentSubSide.incidentChartSideId[1] = -1;
                    }

                    assert(currentSubSide.vertices.size() >= 2);
                    currentSubSide.size = currentSubSide.vertices.size() - 1;

                    chartData.subsides.push_back(currentSubSide);


                    //Pop last vertex
                    if (currentSide.vertices.size() > 0) {
                        assert(currentSide.vertices.back() == currentSubSide.vertices.front());
                        currentSide.vertices.pop_back();
                    }
                    currentSide.vertices.insert(
                                currentSide.vertices.end(),
                                currentSubSide.vertices.begin(),
                                currentSubSide.vertices.end());

                    reversed = false;
                }
                else {
                    assert(currentSubSide.vertices.size() == 0);

                    //Add the side to other chart
                    if (adjChartLabel >= 0) {
                        chart.chartSubsides.push_back(subsideId);

                        assert(chartData.subsides[subsideId].incidentCharts[1] == chart.label);
                        chartData.subsides[subsideId].incidentChartSubsideId[1] = subsideId;
                        chartData.subsides[subsideId].incidentChartSideId[1] = chartSideId;

                        chart.adjacentCharts.push_back(adjChartLabel);

                        //Pop last vertex
                        if (currentSide.vertices.size() > 0) {
                            assert(currentSide.vertices.back() == chartData.subsides[subsideId].vertices.back());
                            currentSide.vertices.pop_back();
                        }

                        //Add side vertices
                        currentSide.vertices.insert(
                                    currentSide.vertices.end(),
                                    chartData.subsides[subsideId].vertices.rbegin(),
                                    chartData.subsides[subsideId].vertices.rend());
                    }

                    reversed = true;
                }

                currentSide.subsides.push_back(subsideId);
                currentSide.reversedSubside.push_back(reversed);

                currentSide.length += chartData.subsides[subsideId].length;
                currentSide.size += chartData.subsides[subsideId].size;

                if (isCorner) {
                    chart.chartSides.push_back(currentSide);
                    currentSide = ChartSide();
                    chartSideId++;
                }

            } while (vCurrentId != vStartId);

#ifndef NDEBUG
            if (!isCorner) {
                std::cout << "Warning 4: Chart has no final corner!" << std::endl;
            }
#endif

        } while (!remainingVertices.empty());

#ifndef NDEBUG
        if (chart.chartSides.size() < 3 || chart.chartSides.size() > 6) {
            std::cout << "Warning 3: Chart " << pId << " has " << chart.chartSides.size() << " sides." << std::endl;
        }
#endif
    }

    return chartData;
}

inline std::vector<double> computeChartEdgeLength(
        const ChartData& chartData,
        const size_t& iterations,
        const std::vector<int>& ilpResults,
        const double& weight)
{
    std::vector<double> avgLengths(chartData.charts.size() , -1);


    std::vector<bool> isFixed(chartData.subsides.size(), false);
    std::vector<bool> isComputable(chartData.subsides.size(), true);
    for (size_t subsideId = 0; subsideId < chartData.subsides.size(); ++subsideId) {
        if (ilpResults[subsideId] == ILP_IGNORE) {
            isComputable[subsideId] = false;
        }
        else if (ilpResults[subsideId] >= 0) {
            isFixed[subsideId] = true;
        }
    }

    //Fill charts with a border
    for (size_t i = 0; i < chartData.charts.size(); i++) {
        const Chart& chart = chartData.charts[i];
        if (chart.faces.size() > 0) {
            double currentQuadLength = 0;
            int numSides = 0;

            for (size_t sId : chart.chartSubsides) {
                const ChartSubside& subside = chartData.subsides[sId];
                if (isComputable[sId] && isFixed[sId]) {
                    currentQuadLength += subside.length / subside.size;
                    numSides++;
                }
            }

            if (numSides > 0) {
                currentQuadLength /= numSides;
                avgLengths[i] = currentQuadLength;
            }
        }
    }

    //Fill charts with no borders
    bool done;
    do {
        done = true;
        for (size_t i = 0; i < chartData.charts.size(); i++) {
            const Chart& chart = chartData.charts[i];
            if (chart.faces.size() > 0 && avgLengths[i] < 0) {
                double currentLength = 0;
                size_t numAdjacentCharts = 0;

                for (size_t adjId : chart.adjacentCharts) {
                    if (avgLengths[adjId] > 0) {
                        currentLength += avgLengths[adjId];
                        numAdjacentCharts++;
                        done = false;
                    }
                }

                if (currentLength > 0) {
                    currentLength /= numAdjacentCharts;
                    avgLengths[i] = currentLength;
                }
            }
        }
    } while (!done);


    //Smoothing
    for (size_t k = 0; k < iterations; k++) {
        std::vector<double> lastAvgLengths = avgLengths;

        for (size_t i = 0; i < chartData.charts.size(); i++) {
            const Chart& chart = chartData.charts[i];
            if (chart.faces.size() > 0) {
                assert(lastAvgLengths[i] > 0);

                double adjValue = 0.0;
                size_t numAdjacentCharts = 0;

                for (size_t adjId : chart.adjacentCharts) {
                    if (avgLengths[adjId] > 0) {
                        adjValue += lastAvgLengths[adjId];
                        numAdjacentCharts++;
                    }
                }

                if (adjValue > 0.0) {
                    adjValue = weight * lastAvgLengths[i] + (1.0 - weight) * adjValue;
                }
            }
        }
    }

    return avgLengths;
}

inline void findSubdivisions(
        const ChartData& chartData,
        const std::vector<double>& chartEdgeLength,
        const Parameters& parameters,
        double& gap,
        std::vector<int>& ilpResults)
{
    return findSubdivisions(
        chartData,
        chartEdgeLength,
        parameters.isAdapt,
        parameters.subsideEdgeLength,
        parameters.ilpMethod,
        parameters.alpha,
        parameters.isometry,
        parameters.regularityQuadrilaterals,
        parameters.regularityNonQuadrilaterals,
        parameters.regularityNonQuadrilateralsWeight,
        parameters.alignSingularities,
        parameters.alignSingularitiesWeight,
        parameters.repeatLosingConstraintsIterations,
        parameters.repeatLosingConstraintsQuads,
        parameters.repeatLosingConstraintsNonQuads,
        parameters.repeatLosingConstraintsAlign,
        parameters.feasibilityFix,
        parameters.hardParityConstraint,
        parameters.timeLimit,
        parameters.gapLimit,
        parameters.callbackTimeLimit,
        parameters.callbackGapLimit,
        parameters.minimumGap,
        gap,
        ilpResults);
}

inline void findSubdivisions(
        const ChartData& chartData,
        const std::vector<double>& chartEdgeLength,
        const bool isAdapt,
        const std::vector<double>& subsideEdgeLength,
        const ILPMethod& method,
        const double alpha,
        const bool isometry,
        const bool regularityQuadrilaterals,
        const bool regularityNonQuadrilaterals,
        const double regularityNonQuadrilateralsWeight,
        const bool alignSingularities,
        const double alignSingularitiesWeight,
        const int repeatLosingConstraintsIterations,
        const bool repeatLosingConstraintsQuads,
        const bool repeatLosingConstraintsNonQuads,
        const bool repeatLosingConstraintsAlign,
        const bool feasibilityFix,
        const bool hardParityConstraint,
        const double timeLimit,
        const double gapLimit,
        const std::vector<float>& callbackTimeLimit,
        const std::vector<float>& callbackGapLimit,
        const double minimumGap,
        double& gap,
        std::vector<int>& ilpResults)
{
    if (chartData.charts.size() <= 0)
        return;

    ILPStatus status;

    //Solve ILP to find the best patches
    std::vector<int> result = ilpResults;
    internal::solveILP(
        chartData,
        chartEdgeLength,
        isAdapt, 
        subsideEdgeLength,
        method,
        alpha,
        isometry,
        regularityQuadrilaterals,
        regularityNonQuadrilaterals,
        regularityNonQuadrilateralsWeight,
        alignSingularities,
        alignSingularitiesWeight,
        repeatLosingConstraintsIterations,
        repeatLosingConstraintsQuads,
        repeatLosingConstraintsNonQuads,
        repeatLosingConstraintsAlign,
        feasibilityFix,
        hardParityConstraint,
        timeLimit,
        gapLimit,
        callbackTimeLimit,
        callbackGapLimit,
        gap,
        status,
        result);

    if (status == ILPStatus::SOLUTIONFOUND && gap < minimumGap) {
        std::cout << "Solution found! Gap: " << gap << std::endl;
        ilpResults = result;
    }
    else if (status == ILPStatus::SOLUTIONWRONG && !hardParityConstraint) {
        std::cout << std::endl << " >>>>>> Solution wrong! Trying with hard constraints for parity. It should not happen..." << std::endl << std::endl;

        return findSubdivisions(
            chartData,
            chartEdgeLength,
            isAdapt,
            subsideEdgeLength,
            method,
            alpha,
            isometry,
            regularityQuadrilaterals,
            regularityNonQuadrilaterals,
            regularityNonQuadrilateralsWeight,
            alignSingularities,
            alignSingularitiesWeight,
            repeatLosingConstraintsIterations,
            repeatLosingConstraintsQuads,
            repeatLosingConstraintsNonQuads,
            repeatLosingConstraintsAlign,
            feasibilityFix,
            true,
            timeLimit,
            gapLimit,
            callbackTimeLimit,
            callbackGapLimit,
            minimumGap,
            gap,
            ilpResults);
    }
    else if (method != ILPMethod::ABS ||
             alignSingularities ||
             repeatLosingConstraintsIterations > 0 ||
             repeatLosingConstraintsQuads ||
             repeatLosingConstraintsNonQuads ||
             repeatLosingConstraintsAlign)
    {
        std::cout << std::endl << " >>>>>> Minimum gap has been not reached. Trying with ABS (linear optimization method), singularity align disabled, minimum gap 1.0 and timeLimit x10. Gap was: " << gap << std::endl << std::endl;

        return findSubdivisions(
            chartData,
            chartEdgeLength,
            isAdapt,
            subsideEdgeLength,
            ILPMethod::ABS,
            alpha,
            isometry,
            regularityQuadrilaterals,
            regularityNonQuadrilaterals,
            regularityNonQuadrilateralsWeight,
            false,
            alignSingularitiesWeight,
            0,
            false,
            false,
            false,
            feasibilityFix,
            hardParityConstraint,
            timeLimit*10,
            gapLimit,
            callbackTimeLimit,
            callbackGapLimit,
            1.0,
            gap,
            ilpResults);
    }
}

template<class TriangleMeshType, class PolyMeshType>
void quadrangulate(
        TriangleMeshType& newSurface,
        const ChartData& chartData,
        const std::vector<size_t> fixedPositionSubsides,
        const std::vector<int>& ilpResult,
        const Parameters& parameters,
        PolyMeshType& quadrangulation,
        std::vector<int>& quadrangulationFaceLabel,
        std::vector<std::vector<size_t>>& quadrangulationPartitions,
        std::vector<std::vector<size_t>>& quadrangulationCorners)
{
    return QuadRetopology::quadrangulate(
            newSurface,
            chartData,
            fixedPositionSubsides,
            ilpResult,
            parameters.chartSmoothingIterations,
            parameters.quadrangulationFixedSmoothingIterations,
            parameters.quadrangulationNonFixedSmoothingIterations,
            parameters.doubletRemoval,
            quadrangulation,
            quadrangulationFaceLabel,
            quadrangulationPartitions,
            quadrangulationCorners);
}

template<class TriangleMeshType, class PolyMeshType>
void quadrangulate(
        TriangleMeshType& newSurface,
        const ChartData& chartData,
        const std::vector<size_t> fixedPositionSubsides,
        const std::vector<int>& ilpResult,
        const int chartSmoothingIterations,
        const int quadrangulationFixedSmoothingIterations,
        const int quadrangulationNonFixedSmoothingIterations,
        const bool doubletRemoval,
        PolyMeshType& quadrangulation,
        std::vector<int>& quadrangulationFaceLabel,
        std::vector<std::vector<size_t>>& quadrangulationPartitions,
        std::vector<std::vector<size_t>>& quadrangulationCorners)
{
    if (newSurface.face.size() <= 0)
        return;
    if (ilpResult.size() == 0)
        return;

    std::vector<std::vector<size_t>> subsideVertexMap(chartData.subsides.size());
    std::vector<int> cornerVertices(newSurface.vert.size(), -1);

    quadrangulationPartitions.resize(chartData.charts.size());
    quadrangulationCorners.resize(chartData.charts.size());

    std::vector<bool> isFixed(chartData.subsides.size(), false);
    for (int sId : fixedPositionSubsides) {
        isFixed[sId] = true;
    }

    //Fill fixed vertices (subsides corners)
    for (const ChartSubside& subside : chartData.subsides) {
        size_t vStart = subside.vertices[0];
        size_t vEnd = subside.vertices[subside.vertices.size() - 1];

        if (cornerVertices[vStart] == -1) {
            cornerVertices[vStart] = quadrangulation.vert.size();
            vcg::tri::Allocator<PolyMeshType>::AddVertex(
                        quadrangulation,
                        newSurface.vert[vStart].P());
        }

        if (cornerVertices[vEnd] == -1) {
            cornerVertices[vEnd] = quadrangulation.vert.size();
            vcg::tri::Allocator<PolyMeshType>::AddVertex(
                        quadrangulation,
                        newSurface.vert[vEnd].P());
        }
    }

    //Fill subside map for fixed borders
    std::set<size_t> fixedVerticesSet;
    for (size_t subsideId = 0; subsideId < chartData.subsides.size(); subsideId++) {
        const ChartSubside& subside = chartData.subsides[subsideId];
        if (isFixed[subsideId]) {
            for (size_t k = 0; k < subside.vertices.size(); k++) {
                const size_t& vId = subside.vertices[k];

                size_t newVertexId;

                if (cornerVertices[vId] == -1) {
                    assert(k > 0 && k < subside.vertices.size() - 1);

                    newVertexId = quadrangulation.vert.size();
                    vcg::tri::Allocator<PolyMeshType>::AddVertex(
                                quadrangulation,
                                newSurface.vert[vId].P());
                }
                else {
                    newVertexId = cornerVertices[vId];
                    assert(newVertexId >= 0);
                }

                fixedVerticesSet.insert(newVertexId);
                subsideVertexMap[subsideId].push_back(newVertexId);
            }

            if (ilpResult[subsideId] > subside.size) {
                int vToSplit = -1;
                double maxLength = 0.0;
                for (size_t k = 0; k < subside.vertices.size() - 1; k++) {
                    const size_t& vId1 = subside.vertices[k];
                    const size_t& vId2 = subside.vertices[k + 1];
                    double length = (newSurface.vert[vId2].P() - newSurface.vert[vId1].P()).Norm();
                    if (length >= maxLength) {
                        vToSplit = k;
                    }
                }

                if (vToSplit >= 0) {
                    const size_t& vId1 = subside.vertices[vToSplit];
                    const size_t& vId2 = subside.vertices[vToSplit + 1];
                    size_t splitVertexId = quadrangulation.vert.size();
                    vcg::tri::Allocator<PolyMeshType>::AddVertex(
                                quadrangulation,
                                (newSurface.vert[vId1].P() + newSurface.vert[vId2].P()) / 2.0);
                    vcg::tri::Allocator<PolyMeshType>::AddFace(
                                quadrangulation, subsideVertexMap[subsideId][vToSplit], subsideVertexMap[subsideId][vToSplit + 1], splitVertexId);

                    subsideVertexMap[subsideId].insert(subsideVertexMap[subsideId].begin() + vToSplit + 1, splitVertexId);

                    std::cout << "Triangle added in subside " << subsideId << ": +1!" << std::endl;
                }
                else {
                    std::cout << "ERROR: impossible to augment the subside " << subsideId << ": +1! Target vertex not found." << std::endl;
                }
            }
            else if (ilpResult[subsideId] < subside.size) {
                if (subside.size >= 2) {
                    int vToSkip = -1;
                    double minLength = std::numeric_limits<double>::max();
                    for (size_t k = 0; k < subside.vertices.size() - 2; k++) {
                        const size_t& vId1 = subside.vertices[k];
                        const size_t& vId2 = subside.vertices[k + 1];
                        const size_t& vId3 = subside.vertices[k + 2];
                        double length = (newSurface.vert[vId2].P() - newSurface.vert[vId1].P()).Norm() + (newSurface.vert[vId3].P() - newSurface.vert[vId2].P()).Norm();
                        if (length <= minLength) {
                            vToSkip = k;
                        }
                    }

                    if (vToSkip >= 0) {
                        vcg::tri::Allocator<PolyMeshType>::AddFace(
                                    quadrangulation, subsideVertexMap[subsideId][vToSkip], subsideVertexMap[subsideId][vToSkip + 1], subsideVertexMap[subsideId][vToSkip + 2]);

                        subsideVertexMap[subsideId].erase(subsideVertexMap[subsideId].begin() + vToSkip + 1);

                        std::cout << "Triangle added in subside " << subsideId << ": -1!" << std::endl;
                    }
                    else {
                        std::cout << "ERROR: impossible to reduce the subside " << subsideId << ": -1! Target vertex not found." << std::endl;
                    }
                }
                else {
                    std::cout << "ERROR: impossible to reduce the subside " << subsideId << ": -1! Subside is less than 2." << std::endl;
                }
            }
        }
    }


    //For each chart
    for (size_t cId = 0; cId < chartData.charts.size(); cId++) {
        const Chart& chart = chartData.charts[cId];

        if (chart.faces.size() == 0)
            continue;

        const std::vector<ChartSide>& chartSides = chart.chartSides;
        if (chartSides.size() < 3 || chartSides.size() > 6) {
            std::cout << "Chart " << cId << " with corners less than 3 or greater than 6!" << std::endl;
            continue;
        }

        bool ilpSolvedForAll = true;
        for (size_t sId : chart.chartSubsides) {
            if (ilpResult[sId] < 0)
                ilpSolvedForAll = false;
        }

        if (!ilpSolvedForAll) {
            std::cout << "Chart " << cId << " not computed. ILP was not solved." << std::endl;
            continue;
        }

        //Input mesh
        Eigen::MatrixXd chartV;
        Eigen::MatrixXi chartF;
        std::vector<std::vector<double>> chartUV;
        vcg::tri::UpdateFlags<TriangleMeshType>::FaceClearS(newSurface);
        vcg::tri::UpdateFlags<TriangleMeshType>::VertexClearS(newSurface);
        for (const size_t& fId : chart.faces) {
            newSurface.face[fId].SetS();
            for (int k = 0; k < newSurface.face[fId].VN(); k++) {
                newSurface.face[fId].V(k)->SetS();
            }
        }
        std::vector<int> vMap, fMap;
        QuadRetopology::internal::VCGToEigen(newSurface, chartV, chartF, vMap, fMap, true, 3);
        QuadRetopology::internal::VCGToEigenUV(newSurface, chartV, chartF, chartUV, true);
        
        // printf("DEBUG: 2\n");

        //Input subdivisions
        Eigen::VectorXi l(chartSides.size());

        std::vector<std::vector<double>> chartSideLength(chartSides.size());
        std::vector<std::vector<std::vector<size_t>>> chartSideVertices(chartSides.size());
        std::vector<std::vector<size_t>> chartSideSubdivision(chartSides.size());

        for (size_t i = 0; i < chartSides.size(); i++) {
            const ChartSide& chartSide = chartSides[i];

            chartSideLength[i].resize(chartSide.subsides.size());
            chartSideVertices[i].resize(chartSide.subsides.size());
            chartSideSubdivision[i].resize(chartSide.subsides.size());

            size_t targetSideSubdivision = 0;
            for (size_t j = 0; j < chartSide.subsides.size(); j++) {
                const size_t& subSideId = chartSides[i].subsides[j];
                const ChartSubside& subSide = chartData.subsides[subSideId];

                if (ilpResult[subSideId] < 0) {
                    std::cout << "Error: ILP not valid" << std::endl;
                    return;
                }

                targetSideSubdivision += ilpResult[subSideId];

                chartSideLength[i][j] = subSide.length;
                chartSideVertices[i][j] = subSide.vertices;
                chartSideSubdivision[i][j] = ilpResult[subSideId];

                if (chartSide.reversedSubside[j]) {
                    std::reverse(chartSideVertices[i][j].begin(), chartSideVertices[i][j].end());
                }

                for (size_t k = 0; k < chartSideVertices[i][j].size(); k++) {
                    size_t vId = chartSideVertices[i][j][k];
                    assert(vMap[vId] >= 0);
                    chartSideVertices[i][j][k] = vMap[vId];
                }

            }

            l(static_cast<int>(i)) = targetSideSubdivision;
        }
        
        // printf("DEBUG: 3\n");

        //Pattern quadrangulation
        Eigen::MatrixXd patchV;
        Eigen::MatrixXi patchF;
        std::vector<size_t> patchBorders;
        std::vector<size_t> patchCorners;
        PolyMeshType patchMesh;
        std::vector<std::vector<size_t>> patchSides;
        QuadRetopology::internal::computePattern(l, patchV, patchF, patchMesh, patchBorders, patchCorners, patchSides);

        // printf("DEBUG: 4\n");

#ifdef QUADRETOPOLOGY_DEBUG_SAVE_MESHES
        igl::writeOBJ(std::string("results/") + std::to_string(cId) + std::string("_patch.obj"), patchV, patchF);
#endif

        assert(chartSides.size() == patchCorners.size());
        assert(chartSides.size() == patchSides.size());

#ifdef QUADRETOPOLOGY_DEBUG_SAVE_MESHES
        igl::writeOBJ(std::string("results/") + std::to_string(cId) + std::string("_chart.obj"), chartV, chartF);
#endif

        //Compute quadrangulation
        Eigen::MatrixXd uvMapV;
        Eigen::MatrixXi uvMapF;
        Eigen::MatrixXd quadrangulationV;
        Eigen::MatrixXi quadrangulationF;
        Eigen::MatrixXd quadrangulationUV;
        // QuadRetopology::internal::computeQuadrangulation(chartV, chartF, patchV, patchF, chartSideVertices, chartSideLength, chartSideSubdivision, patchSides, uvMapV, uvMapF, quadrangulationV, quadrangulationF);
        QuadRetopology::internal::computeQuadrangulation(chartV, chartF, chartUV, patchV, patchF, chartSideVertices, chartSideLength, chartSideSubdivision, patchSides, uvMapV, uvMapF, quadrangulationV, quadrangulationF, quadrangulationUV);


#ifdef QUADRETOPOLOGY_DEBUG_SAVE_MESHES
        Eigen::MatrixXd uvMesh(uvMapV.rows(), 3);
        for (int i = 0; i < uvMapV.rows(); i++) {
            uvMesh(i, 0) = uvMapV(i, 0);
            uvMesh(i, 1) = uvMapV(i, 1);
            uvMesh(i, 2) = 0;
        }

        std::string uvFile = std::string("results/") + std::to_string(cId) + std::string("_uv.obj");
        igl::writeOBJ(uvFile, uvMesh, uvMapF);
#endif
        assert(chartV.rows() == uvMapV.rows());

        //Get polymesh
        PolyMeshType quadrangulatedChartMesh;
        QuadRetopology::internal::eigenToVCG(quadrangulationV, quadrangulationF, quadrangulatedChartMesh, 4);
        QuadRetopology::internal::eigenUVToVCG(quadrangulationV, quadrangulationF, quadrangulationUV, quadrangulatedChartMesh, 4);


#ifdef QUADRETOPOLOGY_DEBUG_SAVE_MESHES
        igl::writeOBJ(std::string("results/") + std::to_string(cId) + std::string("_quadrangulation.obj"), quadrangulationV, quadrangulationF);
#endif

        //Smoothing
        // todo: not finish yet.
        if (chartSmoothingIterations > 0) {
            vcg::tri::UpdateSelection<PolyMeshType>::VertexAll(quadrangulatedChartMesh);
            for (size_t vId : patchBorders) {
                quadrangulatedChartMesh.vert[vId].ClearS();
            }
            vcg::PolygonalAlgorithm<PolyMeshType>::LaplacianReproject(quadrangulatedChartMesh, chartSmoothingIterations, 0.5, true);
        }

        std::vector<int> currentVertexMap(quadrangulatedChartMesh.vert.size(), -1);

        //Map subsides on the vertices of the current mesh (create if necessary)
        for (size_t i = 0; i < chartSides.size(); i++) {
            const ChartSide& side = chartSides[i];
            const std::vector<size_t>& patchSide = patchSides[i];

            size_t currentPatchSideVertex = 0;

            for (size_t j = 0; j < side.subsides.size(); j++) {
                const size_t& subsideId = side.subsides[j];
                const bool& reversed = side.reversedSubside[j];
                const ChartSubside& subside = chartData.subsides[subsideId];

                //Create new vertices of the subsides
                if (subsideVertexMap[subsideId].empty()) {
                    assert(!isFixed[subsideId]);

                    //Get fixed corners of the subside
                    size_t vStart = subside.vertices[0];
                    size_t vEnd = subside.vertices[subside.vertices.size() - 1];
                    assert(cornerVertices[vStart] >= 0 && cornerVertices[vEnd] >= 0);

                    currentVertexMap[patchSide[currentPatchSideVertex]] = cornerVertices[vStart];
                    currentVertexMap[patchSide[currentPatchSideVertex + ilpResult[subsideId]]] = cornerVertices[vEnd];

                    for (int k = 0; k <= ilpResult[subsideId]; k++) {
                        size_t patchSideVId = patchSide[currentPatchSideVertex];

                        if (currentVertexMap[patchSideVId] == -1) {
                            assert(k > 0 && k < ilpResult[subsideId]);

                            //Add new vertex
                            size_t newVertexId = quadrangulation.vert.size();

                            const typename PolyMeshType::CoordType& coord = quadrangulatedChartMesh.vert[patchSideVId].P();
                            vcg::tri::Allocator<PolyMeshType>::AddVertex(quadrangulation, coord);

                            currentVertexMap[patchSideVId] = newVertexId;

                            subsideVertexMap[subsideId].push_back(newVertexId);
                        }
                        else {
                            //Use the existing vertex
                            int existingVertexId = currentVertexMap[patchSideVId];
                            assert(existingVertexId >= 0);
                            subsideVertexMap[subsideId].push_back(existingVertexId);
                        }

                        currentPatchSideVertex++;
                    }

                    if (reversed) {
                        std::reverse(subsideVertexMap[subsideId].begin(), subsideVertexMap[subsideId].end());
                    }
                }
                //Set the existing vertices
                else {
                    assert(subsideVertexMap[subsideId].size() == ilpResult[subsideId] + 1);

                    for (int k = 0; k <= ilpResult[subsideId]; k++) {
                        int patchSideVId = patchSide[currentPatchSideVertex];

                        size_t subSideVertexIndex = reversed ? ilpResult[subsideId] - k : k;

                        currentVertexMap[patchSideVId] = subsideVertexMap[subsideId][subSideVertexIndex];

                        size_t existingVertexId = currentVertexMap[patchSideVId];

                        //If it is not a corner or if it is not on border
                        if (!isFixed[subsideId] && k > 0 && k < ilpResult[subsideId]) {
                            //Average
                            const typename PolyMeshType::CoordType& coord = quadrangulatedChartMesh.vert[patchSideVId].P();
                            quadrangulation.vert[existingVertexId].P() =
                                    (coord + quadrangulation.vert[existingVertexId].P())/2;
                        }

                        currentPatchSideVertex++;
                    }
                }

                currentPatchSideVertex--;
            }

            assert(currentPatchSideVertex+1 == patchSide.size());
        }

        //Internal vertices
        for (size_t i = 0; i < quadrangulatedChartMesh.vert.size(); i++) {
            if (currentVertexMap[i] == -1) {
                size_t newId = quadrangulation.vert.size();

                const typename PolyMeshType::CoordType& coord = quadrangulatedChartMesh.vert[i].P();
                vcg::tri::Allocator<PolyMeshType>::AddVertex(quadrangulation, coord);

                currentVertexMap[i] = newId;
            }
        }

        //Set faces
        for (size_t i = 0; i < quadrangulatedChartMesh.face.size(); i++) {
            assert(quadrangulatedChartMesh.face[i].VN() == 4);

            size_t newFaceId = quadrangulation.face.size();

            vcg::tri::Allocator<PolyMeshType>::AddFaces(quadrangulation, 1);

            quadrangulation.face[newFaceId].Alloc(quadrangulatedChartMesh.face[i].VN());
            
            for (int j = 0; j < quadrangulatedChartMesh.face[i].VN(); j++) {
                int vId = currentVertexMap[vcg::tri::Index(quadrangulatedChartMesh, quadrangulatedChartMesh.face[i].V(j))];
                assert(vId >= 0);

                quadrangulation.face[newFaceId].V(j) = &quadrangulation.vert[vId];

                // for uv 
                quadrangulation.face[newFaceId].WT(j).P() = quadrangulatedChartMesh.face[i].WT(j).P();
            }

            quadrangulationFaceLabel.push_back(chart.label);
            quadrangulationPartitions[chart.label].push_back(newFaceId);
            
        }

        //Fill corners vertices
        for (size_t i = 0; i < chartSides.size(); i++) {
            const ChartSide& side = chartSides[i];
            size_t vStart = side.vertices[0];

            quadrangulationCorners[chart.label].push_back(cornerVertices.at(vStart));
        }
    }

#ifdef QUADRETOPOLOGY_DEBUG_SAVE_MESHES
    vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(quadrangulation, "results/quadrangulation_1_original.obj", vcg::tri::io::Mask::IOM_NONE);
#endif

    //Duplicate vertices
    vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(quadrangulation);
    int numDuplicateVertices = QuadRetopology::internal::removeDuplicateVertices(quadrangulation, false);
    if (numDuplicateVertices > 0) {
        std::cout << "Warning: removed " << numDuplicateVertices << " duplicate vertices in quadrangulation." << std::endl;
    }

    //Degenerate faces
    int numDegenerateFaces = QuadRetopology::internal::removeDegenerateFaces(quadrangulation, false, true);
    if (numDegenerateFaces > 0) {
        std::cout << "Warning: removed " << numDegenerateFaces << " degenerate faces in quadrangulation." << std::endl;
        vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(quadrangulation);
    }

    //Update attributes and re-orient faces
    vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(quadrangulation);
    vcg::PolygonalAlgorithm<PolyMeshType>::UpdateFaceNormalByFitting(quadrangulation);
    QuadRetopology::internal::OrientFaces<PolyMeshType>::AutoOrientFaces(quadrangulation);

    //Reupdate attributes
    vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(quadrangulation);
    vcg::PolygonalAlgorithm<PolyMeshType>::UpdateFaceNormalByFitting(quadrangulation);

    //Doublet removal
    if (doubletRemoval) {
        int numDoublets = QuadRetopology::internal::removeDoubletFaces(quadrangulation, false, true);
        if (numDoublets > 0) {
            std::cout << "Removed " << numDoublets << " doublets in quadrangulation." << std::endl;
            vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(quadrangulation);
        }
    }

    //Unreferenced vertices
    int numUnreferencedVertices = QuadRetopology::internal::removeUnreferencedVertices(quadrangulation, false);
    if (numUnreferencedVertices > 0) {
        std::cout << "Warning: removed " << numUnreferencedVertices << " unreferenced vertices in quadrangulation." << std::endl;
    }


#ifdef QUADRETOPOLOGY_DEBUG_SAVE_MESHES
    vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(quadrangulation, "results/quadrangulation_2_cleaned.obj", vcg::tri::io::Mask::IOM_NONE);
#endif

    //Set birth info in quality
    vcg::tri::UpdateQuality<PolyMeshType>::VertexConstant(quadrangulation, -1);
    vcg::tri::UpdateQuality<PolyMeshType>::FaceConstant(quadrangulation, -1);
    for (size_t i=0; i < quadrangulation.vert.size(); i++) {
        if (quadrangulation.vert[i].IsD())
            continue;

        quadrangulation.vert[i].Q() = i;
    }
    for (size_t i = 0;i < quadrangulation.face.size(); i++) {
        if (quadrangulation.face[i].IsD())
            continue;

        quadrangulation.face[i].Q() = i;
    }

    //Compact
    PolyMeshType tmpMesh;
    vcg::tri::Append<PolyMeshType, PolyMeshType>::Mesh(tmpMesh, quadrangulation);

    //Maps for remapping infos
    std::vector<int> tmpToQuadrangulationVertex(tmpMesh.vert.size(), -1);
    std::vector<int> quadrangulationToTmpVertex(quadrangulation.vert.size(), -1);
    std::vector<int> tmpToQuadrangulationFace(tmpMesh.face.size(), -1);
    std::vector<int> quadrangulationToTmpFace(quadrangulation.face.size(), -1);
    for (size_t i = 0; i < tmpMesh.vert.size(); i++) {
        if (tmpMesh.vert[i].IsD())
            continue;

        int id = tmpMesh.vert[i].Q();
        assert(id >= 0);

        tmpToQuadrangulationVertex[i] = id;
        quadrangulationToTmpVertex[id] = i;
    }
    for (size_t i = 0; i < tmpMesh.face.size(); i++) {
        if (tmpMesh.face[i].IsD())
            continue;

        int id = tmpMesh.face[i].Q();
        assert(id >= 0);

        tmpToQuadrangulationFace[i] = id;
        quadrangulationToTmpFace[id] = i;
    }

    //Remap all infos
    std::vector<int> newFaceLabel(tmpMesh.face.size(), -1);
    for (size_t i = 0; i < tmpMesh.face.size(); i++) {
        if (tmpMesh.face[i].IsD())
            continue;

        assert(tmpToQuadrangulationFace[i] >= 0);
        newFaceLabel[i] = quadrangulationFaceLabel[tmpToQuadrangulationFace[i]];
    }
    quadrangulationFaceLabel = newFaceLabel;

    for (std::vector<size_t>& corners : quadrangulationCorners) {
        std::vector<size_t> newCorners;
        for (const size_t& vId : corners) {
            if (quadrangulation.vert[vId].IsD())
                continue;

            assert(quadrangulationToTmpVertex[vId] >= 0);
            newCorners.push_back(quadrangulationToTmpVertex[vId]);
        }
        corners = newCorners;
    }

    for (std::vector<size_t>& partition : quadrangulationPartitions) {
        std::vector<size_t> newPartition;
        for (const size_t& fId : partition) {
            if (quadrangulation.face[fId].IsD())
                continue;

            assert(quadrangulationToTmpFace[fId] >= 0);
            newPartition.push_back(quadrangulationToTmpFace[fId]);
        }
        partition = newPartition;
    }


    std::vector<size_t> quadrangulationVerticesBetweenPatch;
    for (const size_t& vId : fixedVerticesSet) {
        if (quadrangulation.vert[vId].IsD())
            continue;

        assert(quadrangulationToTmpFace[vId] >= 0);
        quadrangulationVerticesBetweenPatch.push_back(quadrangulationToTmpVertex[vId]);
    }

    //Recopy in quadrangulation mesh
    quadrangulation.Clear();
    vcg::tri::Append<PolyMeshType, PolyMeshType>::Mesh(quadrangulation, tmpMesh);
    vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(quadrangulation);
    vcg::tri::UpdateFlags<PolyMeshType>::FaceBorderFromFF(quadrangulation);
    vcg::tri::UpdateFlags<PolyMeshType>::VertexBorderFromFaceAdj(quadrangulation);

#ifdef QUADRETOPOLOGY_DEBUG_SAVE_MESHES
    vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(quadrangulation, "results/quadrangulation_3_recompacted.obj", vcg::tri::io::Mask::IOM_NONE);
#endif

    vcg::tri::UpdateNormal<TriangleMeshType>::PerFaceNormalized(newSurface);
    vcg::tri::UpdateNormal<TriangleMeshType>::PerVertexNormalized(newSurface);
    vcg::tri::UpdateBounding<TriangleMeshType>::Box(newSurface);

    vcg::GridStaticPtr<typename TriangleMeshType::FaceType,typename TriangleMeshType::FaceType::ScalarType> Grid;
    Grid.Set(newSurface.face.begin(),newSurface.face.end());

    //Reproject
    vcg::tri::UpdateBounding<PolyMeshType>::Box(quadrangulation);
    typename TriangleMeshType::ScalarType maxD=quadrangulation.bbox.Diag();
    typename TriangleMeshType::ScalarType minD=0;

    for (size_t i=0;i<quadrangulation.vert.size();i++) {
        if (quadrangulation.vert[i].IsD())
            continue;

        typename TriangleMeshType::CoordType closestPT;
        typename TriangleMeshType::FaceType *f=
                vcg::tri::GetClosestFaceBase<TriangleMeshType>(
                    newSurface,
                    Grid,
                    quadrangulation.vert[i].P(),
                    maxD,minD,
                    closestPT);

        quadrangulation.vert[i].P()=closestPT;
    }

#ifdef QUADRETOPOLOGY_DEBUG_SAVE_MESHES
    vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(quadrangulation, "results/quadrangulation_4_reprojected.obj", vcg::tri::io::Mask::IOM_NONE);
#endif

    // todo: not finished yet.
    if (quadrangulationFixedSmoothingIterations > 0) {
        vcg::tri::UpdateSelection<PolyMeshType>::VertexAll(quadrangulation);
        for (const size_t& fixedVertexId : quadrangulationVerticesBetweenPatch) {
            if (quadrangulation.vert[fixedVertexId].IsD())
                continue;

            quadrangulation.vert[fixedVertexId].ClearS();
        }

        vcg::PolygonalAlgorithm<PolyMeshType>::template LaplacianReproject<TriangleMeshType>(quadrangulation, newSurface, quadrangulationFixedSmoothingIterations, 0.7, 0.7, true);
    }

#ifdef QUADRETOPOLOGY_DEBUG_SAVE_MESHES
    vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(quadrangulation, "results/quadrangulation_5_smoothed_fixed.obj", vcg::tri::io::Mask::IOM_NONE);
#endif

    // todo: not finished yet.
    if (quadrangulationNonFixedSmoothingIterations > 0) {
        vcg::tri::UpdateSelection<PolyMeshType>::VertexAll(quadrangulation);
        for (size_t i=0;i<quadrangulation.vert.size();i++) {
            if (quadrangulation.vert[i].IsD())
                continue;
            if (quadrangulation.vert[i].IsB()) {
                quadrangulation.vert[i].ClearS();
            }
        }

        vcg::PolygonalAlgorithm<PolyMeshType>::template LaplacianReproject<TriangleMeshType>(quadrangulation, newSurface, quadrangulationNonFixedSmoothingIterations, 0.7, 0.7, true);
    }

    vcg::PolygonalAlgorithm<PolyMeshType>::UpdateFaceNormalByFitting(quadrangulation);
    vcg::tri::UpdateNormal<PolyMeshType>::PerVertexNormalized(quadrangulation);
    vcg::tri::UpdateBounding<PolyMeshType>::Box(quadrangulation);

#ifdef QUADRETOPOLOGY_DEBUG_SAVE_MESHES
    vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(quadrangulation, "results/quadrangulation_6_final.obj", vcg::tri::io::Mask::IOM_NONE);
#endif
}

}
