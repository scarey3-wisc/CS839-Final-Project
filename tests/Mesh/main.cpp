#include "IsosurfaceMesh.h"
#include <math.h>
#include <Eigen/Dense>
#include <random>

template<class T>
struct Shape
{
	using Vector3 = Eigen::Matrix<T, 3, 1>;
	virtual float distanceFunction(const Vector3 x) = 0;
};

template<class T>
struct Subtract : public Shape<T>
{
	using Base = Shape<T>;
	using Vector3 = typename Base::Vector3;

	Shape<T>* shapeOne;
	Shape<T>* shapeTwo;

	Subtract(Shape<T>* one, Shape<T>* two) {
		shapeOne = one;
		shapeTwo = two;
	}

	virtual float distanceFunction(const Vector3 x) override
	{
		return std::max(shapeOne->distanceFunction(x), -1 * shapeTwo->distanceFunction(x));
	}
};

template<class T>
struct Union : public Shape<T>
{
	using Base = Shape<T>;
	using Vector3 = typename Base::Vector3;

	Shape<T>* shapeOne;
	Shape<T>* shapeTwo;
	
	Union(Shape<T>* one, Shape<T>* two) {
		shapeOne = one;
		shapeTwo = two;
	}

	virtual float distanceFunction(const Vector3 x) override
	{
		return std::min(shapeOne->distanceFunction(x), shapeTwo->distanceFunction(x));
	}
};

template<class T>
struct Intersect : public Shape<T>
{
	using Base = Shape<T>;
	using Vector3 = typename Base::Vector3;

	Shape<T>* shapeOne;
	Shape<T>* shapeTwo;

	Intersect(Shape<T>* one, Shape<T>* two) {
		shapeOne = one;
		shapeTwo = two;
	}

	virtual float distanceFunction(const Vector3 x) override
	{
		return std::max(shapeOne->distanceFunction(x), shapeTwo->distanceFunction(x));
	}
};

template<class T>
struct Sphere : public Shape<T>
{
	using Base = Shape<T>;
	using Vector3 = typename Base::Vector3;
	
	T m_radius;
	Vector3 m_center;
	Sphere(const T x, const T y, const T z, const T r) {
		m_radius = r;
		m_center = Vector3(x, y, z);
	}

	virtual float distanceFunction(const Vector3 x) override
	{
		Vector3 d = -1 * x + m_center;
		T distance = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
		return distance - m_radius;
	}
};
template<class T>
struct Rect : public Shape<T>
{
	using Base = Shape<T>;
	using Vector3 = typename Base::Vector3;

	T mx;
	T my;
	T mz;
	T mw;
	T mh;
	T ml;
	T rounding;
	Rect(const T x, const T y, const T z, const T w, const T h, const T l, const T r) {
		mx = x;
		my = y;
		mz = z;
		mw = w;
		mh = h;
		ml = l;
		rounding = r;
	}

	virtual float distanceFunction(const Vector3 x) override
	{
		T dx = abs(x[0] - mx) - mw;
		T dy = abs(x[1] - my) - mh;
		T dz = abs(x[2] - mz) - ml;
		T zi = 0.;
		T ox = std::max(dx, zi);
		T oy = std::max(dy, zi);
		T oz = std::max(dz, zi);
		T le = sqrt(ox * ox + oy * oy + oz * oz);
		return le + std::min(std::max(dx, std::max(dy, dz)), zi) - rounding;
	}
};

template<class T>
struct Cylinder : public Shape<T>
{
	using Base = Shape<T>;
	using Vector3 = typename Base::Vector3;

	T m_radius;
	T m_height;
	Vector3 m_center;
	Cylinder(const T x, const T y, const T z, const T r, const T h) {
		m_radius = r;
		m_height = h;
		m_center = Vector3(x, y, z);
	}

	virtual float distanceFunction(const Vector3 x) override
	{
		Vector3 d = -1 * x + m_center;
		T radi = sqrt(d[0] * d[0] + d[2] * d[2]);
		T radDis = radi - m_radius;
		T heiDis = abs(d[1]) - m_height;
		if (heiDis > 0 && radDis < 0) {
			return heiDis;
		}
		else if (heiDis < 0 && radDis > 0) {
			return radDis;
		}
		else if (heiDis < 0 && radDis < 0) {
			return std::max(heiDis, radDis);
		}
		else if (heiDis > 0 && radDis > 0) {
			return sqrt(heiDis * heiDis + radDis * radDis);
		}
		else {
			return 0;
		}
	}
};

template<class T>
struct ShapeMesh : public IsosurfaceMesh<T>
{
	using Base = IsosurfaceMesh<T>;
	using Vector3 = typename Base::Vector3;

	Shape<T> * shape;

	ShapeMesh(Shape<T> * object)
	{
		shape = object;
	}
	

	float distanceFunction(const Vector3 x) override
	{
		return shape->distanceFunction(x);
	}
};

int main(int argc, char* argv[])
{

	Cylinder<float> l1(0.25, 0.25, 0.25, 0.1, 0.2);
	Cylinder<float> l2(0.25, 0.25, 0.75, 0.1, 0.2);
	Cylinder<float> l3(0.75, 0.25, 0.25, 0.1, 0.2);
	Cylinder<float> l4(0.75, 0.25, 0.75, 0.1, 0.2);
	Union<float> u1(&l1, &l2);
	Union<float> u2(&l3, &l4);
	Union<float> u3(&u1, &u2);
	Rect<float> counter(0.5, 0.5, 0.5, 0.4, 0.05, 0.4, 0.02);
	Rect<float> box(0.25, 0.6, 0.45, 0.15, 0.1, 0.15, 0.01);
	Sphere<float> inner(0.25, 0.65, 0.45, 0.12);
	Subtract<float> cont(&box, &inner);
	Sphere<float> fid(0.7, 0.72, 0.3, 0.2);
	Sphere<float> fed(0.5, 0.58, 0.75, 0.1);
	Union<float> u4(&fid, &fed);
	Union<float> u5(&counter, &cont);
	Union<float> u6(&u4, &u5);

	ShapeMesh<float> simulationMesh(&Union<float>(&u3, &u6));
	simulationMesh.m_cellSize = { 50, 50, 50 };
	simulationMesh.m_gridDX = 0.02;
	simulationMesh.m_nFrames = 0;
	simulationMesh.m_subSteps = 1;
	simulationMesh.m_frameDt = 0.02;

	// Initialize the simulation example
	simulationMesh.initializeUSD("Testing.usda");
	simulationMesh.initialize(1.);
	simulationMesh.initializeDeformation();

	// Output the initial shape of the mesh
	simulationMesh.writeFrame(0);

	// Perform the animation, output results at each frame
	for (int frame = 1; frame <= simulationMesh.m_nFrames; frame++) {
		simulationMesh.simulateFrame(frame);
		simulationMesh.writeFrame(frame);
	}

	// Write the entire timeline to USD
	simulationMesh.writeUSD();

	return 0;
}