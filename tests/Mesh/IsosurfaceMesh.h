#pragma once
#include "FiniteElementMesh.h"
#include <math.h>
#include <Eigen/Dense>

template<class T>
struct IsosurfaceMesh : public FiniteElementMesh<T>
{
	using Base = FiniteElementMesh<T>;

	// from AnimatedMesh
	using Base::m_meshElements;
	using Base::m_particleX;
	using Base::initializeUSD;
	using Base::initializeTopology;
	using Base::initializeParticles;
	using Vector3 = typename Base::Vector3;

	// from FiniteElementMesh
	using Base::m_surfaceParticles;
	using Base::m_particleV;
	using Base::m_particleMass;
	using Base::initializeUndeformedConfiguration;
	using Base::m_stepEndTime;

	std::array<int, 3> m_cellSize; // dimensions in grid cells
	T m_gridDX;


	std::vector<Vector3> m_particleUndeformedX;

	IsosurfaceMesh()
		:Base(1.e0, 5.e1, 2.e2, .002)
	{

	}
	struct LatticePoint {
		std::array<int, 3> gridLoc;
		Vector3 spaceLoc;
		float dist;
		int index;
		LatticePoint(std::array<int, 3> loc, Vector3 x, float d) {
			gridLoc = loc;
			spaceLoc = x;
			dist = d;
			index = -1;
		}
	};
	struct CutPoint {
		std::array<int, 3> outside;
		Vector3 spaceLoc;
		std::array<int, 3> inside;
		float alpha;
		int index;
		bool isBlack;
		CutPoint(std::array<int, 3> pos, Vector3 x, std::array<int, 3> neg, bool black) {
			outside = pos;
			spaceLoc = x;
			inside = neg;
			index = -1;
			isBlack = black;
			if (black)
				alpha = 0.24999;
			else
				alpha = 0.40173;
		}
	};
	void initialize(const float viewBound)
	{
		//generate all the grid points that'll be within our mesh and all the cut points
		std::vector<LatticePoint> candidatePoints;
		std::vector<CutPoint> cutPoints;
		for (int del = 0; del < 2; del++) {
			for (int cell_i = del; cell_i < m_cellSize[0]; cell_i += 2) {
				for (int cell_j = del; cell_j < m_cellSize[1]; cell_j += 2) {
					for (int cell_k = del; cell_k < m_cellSize[2]; cell_k += 2) {
						std::array<int, 3> cand = std::array<int, 3>{cell_i, cell_j, cell_k};
						Vector3 location = convertGridPoint(cand);
						float distance = distanceFunction(location);
						if (distance <= 0.) {
							candidatePoints.push_back(LatticePoint(cand, location, distance));
						}
						else {
							//our point is outside the surface.
							bool added = false;
							std::vector<std::array<int, 3>> connections = getConnections(cand);
							for (int index = 0; index < connections.size(); index++) {
								std::array<int, 3> conn = connections[index];
								if (distanceFunction(convertGridPoint(conn)) < 0.) {
									if (!added) {
										candidatePoints.push_back(LatticePoint(cand, location, distance));
										added = true;
									}
									if (conn[0] >= 0 && conn[0] < m_cellSize[0]
										&& conn[1] >= 0 && conn[1] < m_cellSize[1]
										&& conn[2] >= 0 && conn[2] < m_cellSize[2]) {
										std::array<int, 3> pos = cand;
										std::array<int, 3> neg = conn;
										Vector3 positive = convertGridPoint(pos);
										Vector3 negative = convertGridPoint(neg);
										Vector3 cutPoint = 0.5 * (positive + negative);
										float dist = distanceFunction(cutPoint);
										while (abs(dist) > 0.00001) {
											if (dist > 0)
												positive = cutPoint;
											else
												negative = cutPoint;
											cutPoint = 0.5 * (positive + negative);
											dist = distanceFunction(cutPoint);
										}
										if (index < 6)
											cutPoints.push_back(CutPoint(pos, cutPoint, neg, true));
										else
											cutPoints.push_back(CutPoint(pos, cutPoint, neg, false));
									}
								}
							}
						}
					}
				}
			}
		}

		//create the map between a grid location and an index in the lattice list.
		std::map<std::array<int, 3>, int> pointIndices;
		for (int index = 0; index < candidatePoints.size(); index++) {
			std::array<int, 3> point = candidatePoints[index].gridLoc;
			auto search = pointIndices.find(point);
			if (search != pointIndices.end())
				throw std::logic_error("duplicate point");
			else
				pointIndices.insert({ point, index });
		}

		//create the map between a pos,neg grid index pair and the location of the cut point
		std::map<std::array<int, 2>, CutPoint> cutMap;
		for (const auto& elem : cutPoints) {
			std::array<int, 2> pair = std::array<int, 2>{
				pointIndices.find(elem.outside)->second, 
				pointIndices.find(elem.inside)->second, };
			cutMap.insert({ pair, elem });
		}

		//snap lattice points to nearby cut points
		for (auto& elem : candidatePoints) {
			int index = pointIndices.find(elem.gridLoc)->second;
			bool shouldSnap = false;
			Vector3 snapTo;
			float distance;
			for (std::array<int, 3> conn : getConnections(elem.gridLoc)) {
				int cind = -1; if (pointIndices.find(conn) != pointIndices.end()) { cind = pointIndices.find(conn)->second; }
				if (cind == -1)
					continue;
				const auto& elec = candidatePoints[cind];
				Vector3 ed = -1 * elec.spaceLoc + elem.spaceLoc;
				float edgeLength = sqrt(ed[0] * ed[0] + ed[1] * ed[1] + ed[2] * ed[2]);
				if (elem.dist > 0 && elec.dist < 0) {
					CutPoint cutPointFound = cutMap.find(std::array<int, 2>{index, cind})->second;
					Vector3 snapCan = cutPointFound.spaceLoc;
					Vector3 d = -1 * snapCan + elem.spaceLoc;
					float cutLength = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
					if (cutLength / edgeLength < cutPointFound.alpha) {
						if (shouldSnap && cutLength < distance) {
							distance = cutLength;
							snapTo = snapCan;
						}
						else {
							shouldSnap = true;
							distance = cutLength;
							snapTo = snapCan;
						}
					}
				}
				else if (elem.dist < 0 && elec.dist > 0) {
					CutPoint cutPointFound = cutMap.find(std::array<int, 2>{cind, index})->second;
					Vector3 snapCan = cutPointFound.spaceLoc;
					Vector3 d = -1 * snapCan + elem.spaceLoc;
					float cutLength = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
					if (cutLength / edgeLength < cutPointFound.alpha) {
						if (shouldSnap && cutLength < distance) {
							distance = cutLength;
							snapTo = snapCan;
						}
						else {
							shouldSnap = true;
							distance = cutLength;
							snapTo = snapCan;
						}
					}
				}
			}

			if (shouldSnap) {
				elem.spaceLoc = snapTo;
				elem.dist = 0;
				for (std::array<int, 3> conn : getConnections(elem.gridLoc)) {
					int cind = -1; if (pointIndices.find(conn) != pointIndices.end()) { cind = pointIndices.find(conn)->second; }
					if (cind == -1)
						continue;
					const auto& elec = candidatePoints[cind];
					if (elem.dist > 0 && elec.dist < 0) {
						cutMap.erase(std::array<int, 2>{index, cind});
					}else if (elem.dist < 0 && elec.dist > 0) {
						cutMap.erase(std::array<int, 2>{cind, index});
					}
				}
			}
		}

		//generate the mesh!
		std::map<std::array<int, 3>, int>::iterator iter;
		for (iter = pointIndices.begin(); iter != pointIndices.end(); iter++) {
			int cell_i = iter->first[0];
			int cell_j = iter->first[1];
			int cell_k = iter->first[2];
			if (cell_i % 2 != 1 || cell_j % 2 != 1 || cell_k %2 != 1) {
				continue;
			}
			int centI = iter->second;
			std::array<int, 3 > left = std::array<int, 3>{cell_i - 2, cell_j, cell_k};
			int leftI = -1; if (pointIndices.find(left) != pointIndices.end()) { leftI = pointIndices.find(left)->second; }
			std::array<int, 3 > down = std::array<int, 3>{cell_i, cell_j - 2, cell_k};
			int downI = -1; if (pointIndices.find(down) != pointIndices.end()) { downI = pointIndices.find(down)->second; }
			std::array<int, 3 > outw = std::array<int, 3>{cell_i, cell_j, cell_k - 2};
			int outwI = -1; if (pointIndices.find(outw) != pointIndices.end()) { outwI = pointIndices.find(outw)->second; }
			std::array<int, 3 > clll = std::array<int, 3>{cell_i - 1, cell_j - 1, cell_k - 1};
			int clllI = -1; if (pointIndices.find(clll) != pointIndices.end()) { clllI = pointIndices.find(clll)->second; }
			std::array<int, 3 > cllu = std::array<int, 3>{cell_i - 1, cell_j - 1, cell_k + 1};
			int clluI = -1; if (pointIndices.find(cllu) != pointIndices.end()) { clluI = pointIndices.find(cllu)->second; }
			std::array<int, 3 > clul = std::array<int, 3>{cell_i - 1, cell_j + 1, cell_k - 1};
			int clulI = -1; if (pointIndices.find(clul) != pointIndices.end()) { clulI = pointIndices.find(clul)->second; }
			std::array<int, 3 > cluu = std::array<int, 3>{cell_i - 1, cell_j + 1, cell_k + 1};
			int cluuI = -1; if (pointIndices.find(cluu) != pointIndices.end()) { cluuI = pointIndices.find(cluu)->second; }
			std::array<int, 3 > cull = std::array<int, 3>{cell_i + 1, cell_j - 1, cell_k - 1};
			int cullI = -1; if (pointIndices.find(cull) != pointIndices.end()) { cullI = pointIndices.find(cull)->second; }
			std::array<int, 3 > culu = std::array<int, 3>{cell_i + 1, cell_j - 1, cell_k + 1};
			int culuI = -1; if (pointIndices.find(culu) != pointIndices.end()) { culuI = pointIndices.find(culu)->second; }
			std::array<int, 3 > cuul = std::array<int, 3>{cell_i + 1, cell_j + 1, cell_k - 1};
			int cuulI = -1; if (pointIndices.find(cuul) != pointIndices.end()) { cuulI = pointIndices.find(cuul)->second; }

			if (leftI != -1) {
				if (clllI != -1 && clluI != -1)
					addMeshElement(centI, leftI, clllI, clluI, candidatePoints, cutMap);
				if (clluI != -1 && cluuI != -1)
					addMeshElement(centI, leftI, clluI, cluuI, candidatePoints, cutMap);
				if (cluuI != -1 && clulI != -1)
					addMeshElement(centI, leftI, cluuI, clulI, candidatePoints, cutMap);
				if (clulI != -1 && clllI != -1)
					addMeshElement(centI, leftI, clulI, clllI, candidatePoints, cutMap);
			}
			if (downI != -1) {
				if (clluI != -1 && clllI != -1)
					addMeshElement(centI, downI, clluI, clllI, candidatePoints, cutMap);
				if (culuI != -1 && clluI != -1)
					addMeshElement(centI, downI, culuI, clluI, candidatePoints, cutMap);
				if (cullI != -1 && culuI != -1)
					addMeshElement(centI, downI, cullI, culuI, candidatePoints, cutMap);
				if (clllI != -1 && cullI != -1)
					addMeshElement(centI, downI, clllI, cullI, candidatePoints, cutMap);
			}
			if (outwI != -1) {
				if (clllI != -1 && clulI != -1)
					addMeshElement(centI, outwI, clllI, clulI, candidatePoints, cutMap);
				if (clulI != -1 && cuulI != -1)
					addMeshElement(centI, outwI, clulI, cuulI, candidatePoints, cutMap);
				if (cuulI != -1 && cullI != -1)
					addMeshElement(centI, outwI, cuulI, cullI, candidatePoints, cutMap);
				if (cullI != -1 && clllI != -1)
					addMeshElement(centI, outwI, cullI, clllI, candidatePoints, cutMap);
			}
		}

		initializeTopology(viewBound);
		initializeParticles();

		// Also resize the velocities to match
		m_particleV.resize(m_particleX.size(), Vector3::Zero());

		// Initialize rest shape matrices and particle mass
		initializeUndeformedConfiguration();

		// Also record rest shape
		m_particleUndeformedX = m_particleX;
	}
	void addMeshElement(int one, int two, int thr, int fou, 
		std::vector<LatticePoint>& candidatePoints, std::map<std::array<int, 2>, CutPoint>& cutMap) {
		int numPos = 0;
		int numZip = 0;
		int numNeg = 0;
		if (candidatePoints[one].dist > 0)
			numPos++;
		else if (candidatePoints[one].dist == 0)
			numZip++;
		else
			numNeg++;
		if (candidatePoints[two].dist > 0)
			numPos++;
		else if (candidatePoints[two].dist == 0)
			numZip++;
		else
			numNeg++;
		if (candidatePoints[thr].dist > 0)
			numPos++;
		else if (candidatePoints[thr].dist == 0)
			numZip++;
		else
			numNeg++;
		if (candidatePoints[fou].dist > 0)
			numPos++;
		else if (candidatePoints[fou].dist == 0)
			numZip++;
		else
			numNeg++;

		if (numPos == 0) {
			addStencilFull(one, two, thr, fou, candidatePoints, cutMap);
		}
		else if (numNeg == 1) {
			if (candidatePoints[one].dist < 0) addStencilMini(one, two, thr, fou, candidatePoints, cutMap, numZip);
			else if (candidatePoints[two].dist < 0) addStencilMini(two, thr, one, fou, candidatePoints, cutMap, numZip);
			else if (candidatePoints[thr].dist < 0) addStencilMini(thr, one, two, fou, candidatePoints, cutMap, numZip);
			else if (candidatePoints[fou].dist < 0) addStencilMini(fou, one, thr, two, candidatePoints, cutMap, numZip);
		}
		else if (numZip == 1 && numPos == 1) {
			if (candidatePoints[one].dist == 0) addStencilSliced(one, two, thr, fou, candidatePoints, cutMap);
			else if (candidatePoints[two].dist == 0) addStencilSliced(two, thr, one, fou, candidatePoints, cutMap);
			else if (candidatePoints[thr].dist == 0) addStencilSliced(thr, one, two, fou, candidatePoints, cutMap);
			else if (candidatePoints[fou].dist == 0) addStencilSliced(fou, one, thr, two, candidatePoints, cutMap);
		}
		else if (numPos == 2 && numNeg == 2) {
			if (candidatePoints[one].dist > 0) {
				if (candidatePoints[two].dist > 0) addStencilChopped(one, two, thr, fou, candidatePoints, cutMap);
				else if (candidatePoints[thr].dist > 0) addStencilChopped(one, thr, fou, two, candidatePoints, cutMap);
				else if (candidatePoints[fou].dist > 0) addStencilChopped(one, fou, two, thr, candidatePoints, cutMap);
			}
			else if (candidatePoints[two].dist > 0) {
				if (candidatePoints[thr].dist > 0) addStencilChopped(two, thr, one, fou, candidatePoints, cutMap);
				else if (candidatePoints[fou].dist > 0) addStencilChopped(two, fou, thr, one, candidatePoints, cutMap);
			}
			else if (candidatePoints[thr].dist > 0) {
				addStencilChopped(thr, fou, one, two, candidatePoints, cutMap);
			}
		}
		else if (numPos == 1 && numNeg == 3) {
			if (candidatePoints[one].dist > 0) addStencilClipped(one, two, thr, fou, candidatePoints, cutMap);
			else if (candidatePoints[two].dist > 0) addStencilClipped(two, thr, one, fou, candidatePoints, cutMap);
			else if (candidatePoints[thr].dist > 0) addStencilClipped(thr, one, two, fou, candidatePoints, cutMap);
			else if (candidatePoints[fou].dist > 0) addStencilClipped(fou, one, thr, two, candidatePoints, cutMap);
		}
	}
	void addStencilClipped(int a, int b, int c, int d, std::vector<LatticePoint>& lat, std::map<std::array<int, 2>, CutPoint>& cut) {
		CutPoint cutOne = cut.find(std::array<int, 2>{a, b})->second;
		if (cutOne.index == -1) {
			cutOne.index = m_particleX.size();
			m_particleX.emplace_back(cutOne.spaceLoc);
			m_surfaceParticles.emplace_back(cutOne.index);
			cut.insert_or_assign(std::array<int, 2>{a, b}, cutOne);
		}
		CutPoint cutTwo = cut.find(std::array<int, 2>{a, c})->second;
		if (cutTwo.index == -1) {
			cutTwo.index = m_particleX.size();
			m_particleX.emplace_back(cutTwo.spaceLoc);
			m_surfaceParticles.emplace_back(cutTwo.index);
			cut.insert_or_assign(std::array<int, 2>{a, c}, cutTwo);
		}
		CutPoint cutThr = cut.find(std::array<int, 2>{a, d})->second;
		if (cutThr.index == -1) {
			cutThr.index = m_particleX.size();
			m_particleX.emplace_back(cutThr.spaceLoc);
			m_surfaceParticles.emplace_back(cutThr.index);
			cut.insert_or_assign(std::array<int, 2>{a, d}, cutThr);
		}
		if (lat[b].index == -1) {
			lat[b].index = m_particleX.size();
			m_particleX.emplace_back(lat[b].spaceLoc);
		}
		if (lat[c].index == -1) {
			lat[c].index = m_particleX.size();
			m_particleX.emplace_back(lat[c].spaceLoc);
		}
		if (lat[d].index == -1) {
			lat[d].index = m_particleX.size();
			m_particleX.emplace_back(lat[d].spaceLoc);
		}

		if (cutOne.isBlack || true) {
			m_meshElements.push_back(std::array<int, 4>{cutOne.index, lat[b].index, lat[c].index, lat[d].index});
			int parity = 0;
			if (lat[d].spaceLoc[0] > cutTwo.spaceLoc[0]) parity++;
			if (lat[d].spaceLoc[1] > cutTwo.spaceLoc[1]) parity++;
			if (lat[d].spaceLoc[2] > cutTwo.spaceLoc[2]) parity++;
			if ((parity % 2 == 1) == ((lat[d].gridLoc[0] + lat[d].gridLoc[1] + lat[d].gridLoc[2]) % 2 == 0) || true) {
				m_meshElements.push_back(std::array<int, 4>{cutThr.index, cutOne.index, cutTwo.index, lat[d].index});
				m_meshElements.push_back(std::array<int, 4>{cutTwo.index, cutOne.index, lat[c].index, lat[d].index});
			}
			else {
				m_meshElements.push_back(std::array<int, 4>{cutThr.index, cutOne.index, cutTwo.index, lat[c].index});
				m_meshElements.push_back(std::array<int, 4>{cutThr.index, cutOne.index, lat[c].index, lat[d].index});
			}
		}
		else if (cutTwo.isBlack) {
			m_meshElements.push_back(std::array<int, 4>{cutTwo.index, lat[c].index, lat[d].index, lat[b].index});
			int parity = 0;
			if (lat[b].spaceLoc[0] > cutThr.spaceLoc[0]) parity++;
			if (lat[b].spaceLoc[1] > cutThr.spaceLoc[1]) parity++;
			if (lat[b].spaceLoc[2] > cutThr.spaceLoc[2]) parity++;
			if ((parity % 2 == 1) == ((lat[b].gridLoc[0] + lat[b].gridLoc[1] + lat[b].gridLoc[2]) % 2 == 0)) {
				m_meshElements.push_back(std::array<int, 4>{cutOne.index, cutTwo.index, cutThr.index, lat[b].index});
				m_meshElements.push_back(std::array<int, 4>{cutThr.index, cutTwo.index, lat[d].index, lat[b].index});
			}
			else {
				m_meshElements.push_back(std::array<int, 4>{cutOne.index, cutTwo.index, cutThr.index, lat[d].index});
				m_meshElements.push_back(std::array<int, 4>{cutOne.index, cutTwo.index, lat[d].index, lat[b].index});
			}
		}
		else if (cutThr.isBlack) {
			m_meshElements.push_back(std::array<int, 4>{cutThr.index, lat[d].index, lat[b].index, lat[c].index});
			int parity = 0;
			if (lat[c].spaceLoc[0] > cutOne.spaceLoc[0]) parity++;
			if (lat[c].spaceLoc[1] > cutOne.spaceLoc[1]) parity++;
			if (lat[c].spaceLoc[2] > cutOne.spaceLoc[2]) parity++;
			if ((parity % 2 == 1) == ((lat[c].gridLoc[0] + lat[c].gridLoc[1] + lat[c].gridLoc[2]) % 2 == 0)) {
				m_meshElements.push_back(std::array<int, 4>{cutTwo.index, cutThr.index, cutOne.index, lat[c].index});
				m_meshElements.push_back(std::array<int, 4>{cutOne.index, cutThr.index, lat[b].index, lat[c].index});
			}
			else {
				m_meshElements.push_back(std::array<int, 4>{cutTwo.index, cutThr.index, cutOne.index, lat[b].index});
				m_meshElements.push_back(std::array<int, 4>{cutTwo.index, cutThr.index, lat[b].index, lat[c].index});
			}
		}
	}
	void addStencilChopped(int a, int b, int c, int d, std::vector<LatticePoint>& lat, std::map<std::array<int, 2>, CutPoint>& cut) {
		CutPoint cutOne = cut.find(std::array<int, 2>{a, c})->second;
		if (cutOne.index == -1) {
			cutOne.index = m_particleX.size();
			m_particleX.emplace_back(cutOne.spaceLoc);
			m_surfaceParticles.emplace_back(cutOne.index);
			cut.insert_or_assign(std::array<int, 2>{a, c}, cutOne);
		}
		CutPoint cutTwo = cut.find(std::array<int, 2>{a, d})->second;
		if (cutTwo.index == -1) {
			cutTwo.index = m_particleX.size();
			m_particleX.emplace_back(cutTwo.spaceLoc);
			m_surfaceParticles.emplace_back(cutTwo.index);
			cut.insert_or_assign(std::array<int, 2>{a, d}, cutTwo);
		}
		CutPoint cutThr = cut.find(std::array<int, 2>{b, c})->second;
		if (cutThr.index == -1) {
			cutThr.index = m_particleX.size();
			m_particleX.emplace_back(cutThr.spaceLoc);
			m_surfaceParticles.emplace_back(cutThr.index);
			cut.insert_or_assign(std::array<int, 2>{b, c}, cutThr);
		}
		CutPoint cutFou = cut.find(std::array<int, 2>{b, d})->second;
		if (cutFou.index == -1) {
			cutFou.index = m_particleX.size();
			m_particleX.emplace_back(cutFou.spaceLoc);
			m_surfaceParticles.emplace_back(cutFou.index);
			cut.insert_or_assign(std::array<int, 2>{b, d}, cutFou);
		}
		if (lat[c].index == -1) {
			lat[c].index = m_particleX.size();
			m_particleX.emplace_back(lat[c].spaceLoc);
		}
		if (lat[d].index == -1) {
			lat[d].index = m_particleX.size();
			m_particleX.emplace_back(lat[d].spaceLoc);
		}
		if (!cutOne.isBlack && !cutTwo.isBlack && !cutThr.isBlack && !cutFou.isBlack) {

			int parity = 0;
			if (lat[d].spaceLoc[0] > cutOne.spaceLoc[0]) parity++;
			if (lat[d].spaceLoc[1] > cutOne.spaceLoc[1]) parity++;
			if (lat[d].spaceLoc[2] > cutOne.spaceLoc[2]) parity++;
			if ((parity % 2 == 1) == ((lat[d].gridLoc[0] + lat[d].gridLoc[1] + lat[d].gridLoc[2]) % 2 == 0)) {
				m_meshElements.push_back(std::array<int, 4>{cutOne.index, cutFou.index, lat[c].index, lat[d].index});
				m_meshElements.push_back(std::array<int, 4>{cutTwo.index, cutFou.index, cutOne.index, lat[d].index});
				m_meshElements.push_back(std::array<int, 4>{cutOne.index, cutThr.index, lat[c].index, cutFou.index});
			}
			else {
				m_meshElements.push_back(std::array<int, 4>{cutTwo.index, cutThr.index, lat[c].index, lat[d].index});
				m_meshElements.push_back(std::array<int, 4>{cutFou.index, cutThr.index, cutTwo.index, lat[d].index});
				m_meshElements.push_back(std::array<int, 4>{cutTwo.index, cutOne.index, lat[c].index, cutThr.index});
			}
			
		}
		else {
			if (cutTwo.isBlack && cutThr.isBlack) {
				m_meshElements.push_back(std::array<int, 4>{cutOne.index, cutThr.index, lat[c].index, cutTwo.index});
				m_meshElements.push_back(std::array<int, 4>{cutTwo.index, cutThr.index, lat[c].index, lat[d].index});
				m_meshElements.push_back(std::array<int, 4>{cutTwo.index, cutFou.index, cutThr.index, lat[d].index});
			}
			else if (cutOne.isBlack && cutFou.isBlack) {
				m_meshElements.push_back(std::array<int, 4>{cutTwo.index, cutFou.index, cutOne.index, lat[d].index});
				m_meshElements.push_back(std::array<int, 4>{cutOne.index, cutFou.index, lat[c].index, lat[d].index});
				m_meshElements.push_back(std::array<int, 4>{cutOne.index, cutThr.index, lat[c].index, cutFou.index});
			}
			
		}
	}
	void addStencilSliced(int a, int b, int c, int d, std::vector<LatticePoint>& lat, std::map<std::array<int, 2>, CutPoint>& cut) {
		if (lat[b].dist > 0) addStencilFou(a, b, c, d, lat, cut);
		else if (lat[c].dist > 0) addStencilFou(a, c, d, b, lat, cut);
		else if(lat[d].dist > 0) addStencilFou(a, d, b, c, lat, cut);
	}
	void addStencilFou(int a, int b, int c, int d, std::vector<LatticePoint>& lat, std::map<std::array<int, 2>, CutPoint>& cut) {
		if (lat[a].index == -1) {
			lat[a].index = m_particleX.size();
			m_particleX.emplace_back(lat[a].spaceLoc);
			m_surfaceParticles.emplace_back(lat[a].index);
		}
		CutPoint cutOne = cut.find(std::array<int, 2>{b, c})->second;
		if (cutOne.index == -1) {
			cutOne.index = m_particleX.size();
			m_particleX.emplace_back(cutOne.spaceLoc);
			m_surfaceParticles.emplace_back(cutOne.index);
			cut.insert_or_assign(std::array<int, 2>{b, c}, cutOne);
		}
		if (lat[c].index == -1) {
			lat[c].index = m_particleX.size();
			m_particleX.emplace_back(lat[c].spaceLoc);
		}
		CutPoint cutTwo = cut.find(std::array<int, 2>{b, d})->second;
		if (cutTwo.index == -1) {
			cutTwo.index = m_particleX.size();
			m_particleX.emplace_back(cutTwo.spaceLoc);
			m_surfaceParticles.emplace_back(cutTwo.index);
			cut.insert_or_assign(std::array<int, 2>{b, d}, cutTwo);
		}
		if (lat[d].index == -1) {
			lat[d].index = m_particleX.size();
			m_particleX.emplace_back(lat[d].spaceLoc);
		}
		if (cutOne.isBlack == cutTwo.isBlack) {
			int parity = 0;
			if (lat[d].spaceLoc[0] > cutOne.spaceLoc[0]) parity++;
			if (lat[d].spaceLoc[1] > cutOne.spaceLoc[1]) parity++;
			if (lat[d].spaceLoc[2] > cutOne.spaceLoc[2]) parity++;
			if ((parity % 2 == 1) == ((lat[d].gridLoc[0] + lat[d].gridLoc[1] + lat[d].gridLoc[2]) % 2 == 0)) {
				m_meshElements.push_back(std::array<int, 4>{lat[a].index, cutTwo.index, cutOne.index, lat[d].index});
				m_meshElements.push_back(std::array<int, 4>{lat[a].index, cutOne.index, lat[c].index, lat[d].index});
			}
			else {
				m_meshElements.push_back(std::array<int, 4>{lat[a].index, cutTwo.index, cutOne.index, lat[c].index});
				m_meshElements.push_back(std::array<int, 4>{lat[a].index, cutTwo.index, lat[c].index, lat[d].index});
			}
		}
		else if (cutOne.isBlack) {
			m_meshElements.push_back(std::array<int, 4>{lat[a].index, cutTwo.index, cutOne.index, lat[d].index});
			m_meshElements.push_back(std::array<int, 4>{lat[a].index, cutOne.index, lat[c].index, lat[d].index});
		}
		else if(cutTwo.isBlack){
			m_meshElements.push_back(std::array<int, 4>{lat[a].index, cutTwo.index, cutOne.index, lat[c].index});
			m_meshElements.push_back(std::array<int, 4>{lat[a].index, cutTwo.index, lat[c].index, lat[d].index});
		}
	}
	void addStencilMini(int a, int b, int c, int d, std::vector<LatticePoint>& lat, std::map<std::array<int, 2>, CutPoint>& cut, int numZip) {
		if (numZip == 0) {
			addStencilThr(a, b, c, d, lat, cut);
		}
		else if (numZip == 1) {
			if (lat[b].dist == 0) addStencilTwo(a, b, c, d, lat, cut);
			else if (lat[c].dist == 0) addStencilTwo(a, c, d, b, lat, cut);
			else if (lat[d].dist == 0) addStencilTwo(a, d, b, c, lat, cut);
		}
		else if (numZip == 2) {
			if (lat[b].dist > 0) addStencilOne(a, b, c, d, lat, cut);
			else if (lat[c].dist > 0) addStencilOne(a, c, d, b, lat, cut);
			else if (lat[d].dist > 0) addStencilOne(a, d, b, c, lat, cut);
		}
	}
	void addStencilTwo(int a, int b, int c, int d, std::vector<LatticePoint>& lat, std::map<std::array<int, 2>, CutPoint>& cut) {
		if (lat[a].index == -1) {
			lat[a].index = m_particleX.size();
			m_particleX.emplace_back(lat[a].spaceLoc);
			if (lat[a].dist == 0) m_surfaceParticles.emplace_back(lat[a].index);
		}
		if (lat[b].index == -1) {
			lat[b].index = m_particleX.size();
			m_particleX.emplace_back(lat[b].spaceLoc);
			if (lat[b].dist == 0) m_surfaceParticles.emplace_back(lat[b].index);
		}
		CutPoint cutOne = cut.find(std::array<int, 2>{c, a})->second;
		if (cutOne.index == -1) {
			cutOne.index = m_particleX.size();
			m_particleX.emplace_back(cutOne.spaceLoc);
			m_surfaceParticles.emplace_back(cutOne.index);
			cut.insert_or_assign(std::array<int, 2>{c, a}, cutOne);
		}
		CutPoint cutTwo = cut.find(std::array<int, 2>{d, a})->second;
		if (cutTwo.index == -1) {
			cutTwo.index = m_particleX.size();
			m_particleX.emplace_back(cutTwo.spaceLoc);
			m_surfaceParticles.emplace_back(cutTwo.index);
			cut.insert_or_assign(std::array<int, 2>{d, a}, cutTwo);
		}
		m_meshElements.push_back(std::array<int, 4>{lat[a].index, lat[b].index, cutOne.index, cutTwo.index});
	}
	void addStencilOne(int a, int b, int c, int d, std::vector<LatticePoint>& lat, std::map<std::array<int, 2>, CutPoint>& cut) {
		if (lat[a].index == -1) {
			lat[a].index = m_particleX.size();
			m_particleX.emplace_back(lat[a].spaceLoc);
			if (lat[a].dist == 0) m_surfaceParticles.emplace_back(lat[a].index);
		}
		CutPoint cutOne = cut.find(std::array<int, 2>{b, a})->second;
		if (cutOne.index == -1) {
			cutOne.index = m_particleX.size();
			m_particleX.emplace_back(cutOne.spaceLoc);
			m_surfaceParticles.emplace_back(cutOne.index);
			cut.insert_or_assign(std::array<int, 2>{b, a}, cutOne);
		}
		if (lat[c].index == -1) {
			lat[c].index = m_particleX.size();
			m_particleX.emplace_back(lat[c].spaceLoc);
			if (lat[c].dist == 0) m_surfaceParticles.emplace_back(lat[c].index);
		}
		if (lat[d].index == -1) {
			lat[d].index = m_particleX.size();
			m_particleX.emplace_back(lat[d].spaceLoc);
			if (lat[d].dist == 0) m_surfaceParticles.emplace_back(lat[d].index);
		}
		m_meshElements.push_back(std::array<int, 4>{lat[a].index, cutOne.index, lat[c].index, lat[d].index});
	}
	void addStencilThr(int a, int b, int c, int d, std::vector<LatticePoint>& lat, std::map<std::array<int, 2>, CutPoint>& cut) {
		if (lat[a].index == -1) {
			lat[a].index = m_particleX.size();
			m_particleX.emplace_back(lat[a].spaceLoc);
			if (lat[a].dist == 0) m_surfaceParticles.emplace_back(lat[a].index);
		}
		CutPoint cutOne = cut.find(std::array<int, 2>{b, a})->second;
		CutPoint cutTwo = cut.find(std::array<int, 2>{c, a})->second;
		CutPoint cutThr = cut.find(std::array<int, 2>{d, a})->second;
		if (cutOne.index == -1) {
			cutOne.index = m_particleX.size();
			m_particleX.emplace_back(cutOne.spaceLoc);
			m_surfaceParticles.emplace_back(cutOne.index);
			cut.insert_or_assign(std::array<int, 2>{b, a}, cutOne);
		}
		if (cutTwo.index == -1) {
			cutTwo.index = m_particleX.size();
			m_particleX.emplace_back(cutTwo.spaceLoc);
			m_surfaceParticles.emplace_back(cutTwo.index);
			cut.insert_or_assign(std::array<int, 2>{c, a}, cutTwo);
		}
		if (cutThr.index == -1) {
			cutThr.index = m_particleX.size();
			m_particleX.emplace_back(cutThr.spaceLoc);
			m_surfaceParticles.emplace_back(cutThr.index);
			cut.insert_or_assign(std::array<int, 2>{d, a}, cutThr);
		}
		m_meshElements.push_back(std::array<int, 4>{lat[a].index, cutOne.index, cutTwo.index, cutThr.index});
	}
	void addStencilFull(int a, int b, int c, int d, std::vector<LatticePoint>& lat, std::map<std::array<int, 2>, CutPoint>& cut) {
		if (lat[a].index == -1) {
			lat[a].index = m_particleX.size();
			m_particleX.emplace_back(lat[a].spaceLoc);
			if (lat[a].dist == 0) m_surfaceParticles.emplace_back(lat[a].index);
		}
		if (lat[b].index == -1) {
			lat[b].index = m_particleX.size();
			m_particleX.emplace_back(lat[b].spaceLoc);
			if (lat[b].dist == 0) m_surfaceParticles.emplace_back(lat[b].index);
		}
		if (lat[c].index == -1) {
			lat[c].index = m_particleX.size();
			m_particleX.emplace_back(lat[c].spaceLoc);
			if (lat[c].dist == 0) m_surfaceParticles.emplace_back(lat[c].index);
		}
		if (lat[d].index == -1) {
			lat[d].index = m_particleX.size();
			m_particleX.emplace_back(lat[d].spaceLoc);
			if (lat[d].dist == 0) m_surfaceParticles.emplace_back(lat[d].index);
		}
		m_meshElements.push_back(std::array<int, 4>{lat[a].index, lat[b].index, lat[c].index, lat[d].index});
	}
	Vector3 convertGridPoint(std::array<int, 3 > loc) {
		return Vector3(m_gridDX * (T)loc[0], m_gridDX * (T)loc[1], m_gridDX * (T)loc[2]);
	}
	std::vector<std::array<int, 3>> getConnections(std::array<int, 3 > loc) {
		std::vector<std::array<int, 3>> connections;
		connections.push_back(std::array<int, 3>{loc[0] + 2, loc[1], loc[2]});
		connections.push_back(std::array<int, 3>{loc[0], loc[1] + 2, loc[2]});
		connections.push_back(std::array<int, 3>{loc[0], loc[1], loc[2] + 2});
		connections.push_back(std::array<int, 3>{loc[0], loc[1], loc[2] - 2});
		connections.push_back(std::array<int, 3>{loc[0], loc[1] - 2, loc[2]});
		connections.push_back(std::array<int, 3>{loc[0] - 2, loc[1], loc[2]});
		connections.push_back(std::array<int, 3>{loc[0] + 1, loc[1] + 1, loc[2] + 1});
		connections.push_back(std::array<int, 3>{loc[0] + 1, loc[1] - 1, loc[2] + 1});
		connections.push_back(std::array<int, 3>{loc[0] - 1, loc[1] + 1, loc[2] + 1});
		connections.push_back(std::array<int, 3>{loc[0] - 1, loc[1] - 1, loc[2] + 1});
		connections.push_back(std::array<int, 3>{loc[0] + 1, loc[1] + 1, loc[2] - 1});
		connections.push_back(std::array<int, 3>{loc[0] + 1, loc[1] - 1, loc[2] - 1});
		connections.push_back(std::array<int, 3>{loc[0] - 1, loc[1] + 1, loc[2] - 1});
		connections.push_back(std::array<int, 3>{loc[0] - 1, loc[1] - 1, loc[2] - 1});
		return connections;
	}

	void initializeDeformation()
	{
		// No need to apply any deformation; this example is driven by moving handles
	}

	virtual void addExternalForce(std::vector<Vector3>& f)
	{
		Vector3 aGravity(0., -9.81, 0);

		for (int p = 0; p < m_particleX.size(); p++)
			f[p] += m_particleMass[p] * aGravity;
	}

	void clearConstrainedParticles(std::vector<Vector3>& x) override
	{

	}

	void setBoundaryConditions() override
	{

	}
	virtual float distanceFunction(const Vector3 x) { return 0.; }

private:
	inline int gridToParticleID(const int i, const int j, const int k) const {
		return i * (m_cellSize[1] + 1) * (m_cellSize[2] + 1) + j * (m_cellSize[2] + 1) + k;
	}
};