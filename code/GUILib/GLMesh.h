#pragma once

#include <MathLib/MathLib.h>
#include <Utils/Utils.h>
#include <MathLib/P3D.h>
#include <MathLib/V3D.h>
#include <MathLib/Quaternion.h>
#include "GLShaderMaterial.h"
#include <assert.h>

#include <vector>

/*=====================================================================================================================================================================*
 * This file implements the classes needed for the storage and display of simple 3d meshes. The classes are implemented in a way that promotes the use of vertex       *
 * arrays for increased efficiency.                                                                                                                                    *
 * At some point a better mesh data structure suitable for editing operations should be used (e.g. the half-Edge structure - http://www.holmes3d.net/graphics/dcel/)   *
 *                                                                                                                                                                     *
 * The following classes are implemented:                                                                                                                              *
 *                                                                                                                                                                     *
 * CLASS VertexNeighbourInfo: For every vertex in every polygon, this class will keep track of the neighbours - this will be used for fast normal computations.                 *
 *                                                                                                                                                                     *
 * CLASS SharedVertexInfo: This class is used to store, for each vertex, a collection of VertexNeighbourInfo classes - this way, when computing the normal at each vertex,      *
 *					we only need to iterate through this collection in order to get the normals from each polygon that contains this vertex, and then we average them  *
 *                                                                                                                                                                     *
 * CLASS GLIndexedPoly: This class contains a list of vertex indices that correspond to the vertices that make up the polygon. This class should be populated by the   *
 *					class that reads in the mesh file, and then passes it to the GLMeshObject.                                                                         *
 *                                                                                                                                                                     *
 * CLASS GLPolyIndexList: This class contains an array of tightly packed indeces of vertices that belong to polygons that have the same number of vertices. As an      *
 *					example the indices of all vertices that belong to traingles in the mesh, would be tightly packed together in one instance of this class. There    *
 *					will be one of these objects for all the triangles in the mesh, one for all the quads, and so on.                                                  *
 *                                                                                                                                                                     *
 * CLASS GLPolyCategory: The last in the set of helper classes, this one acts as a container of GLPolyIndexList, one for triangles, one for quads, and so on.          *
 *					Everytime a new polygon is added, this class is responsible with choosing the correct PolyIndexList that its vertices should be added to.          *
 *                                                                                                                                                                     *
 * CLASS GLMesh:	This is the class that has lists of vertex coordinates, texture coordinates and normals that will be used when handling the object. It implements  *
 *					methods for printing the model.                                                                                                                    *
 *                                                                                                                                                                     *
 *=====================================================================================================================================================================*/

/**
	Store, for a given vertex, and for a given polygon, the two neighbouring vertices. Used to quickly generate the normals for each vertex. 
*/
class VertexNeighbourInfo  {
public:
	int n1Index, n2Index;		//these are the indices of the two neighbouring vertices (the one before and the one after in the polygon)
	VertexNeighbourInfo(int n1Index, int n2Index){this->n1Index = n1Index; this->n2Index = n2Index;}
	VertexNeighbourInfo(const VertexNeighbourInfo& other){
		this->n1Index = other.n1Index;
		this->n2Index = other.n2Index;
	}
	~VertexNeighbourInfo(){
	}
};


/**
	This class is used to store, for each vertex, a collection of VertexNeighbourInfo classes - this way, when computing the normal at each vertex,
	we only need to iterate through this collection in order to get the normals from each polygon that contains this vertex
*/
class SharedVertexInfo  {
friend class GLMesh;
private:
	//the VertexNeighbourInfo class only has two ints, so we can use the POD DynamicArray
	DynamicArray<VertexNeighbourInfo> vertexInstances;
	/**
		Default constructor
	*/
	SharedVertexInfo(){
		vertexInstances = DynamicArray<VertexNeighbourInfo>();
	}
	SharedVertexInfo(const SharedVertexInfo &other){
		for(uint i=0;i<other.vertexInstances.size();i++)
			vertexInstances.push_back(VertexNeighbourInfo(other.vertexInstances[i].n1Index,other.vertexInstances[i].n2Index));
	}

	/**
		Destructor.
	*/
	~SharedVertexInfo(){
	}

	/**
		Adding a new VertexNeighbourInfo object.
	*/
	void addVertexInfo(const VertexNeighbourInfo &vi){
		vertexInstances.push_back(vi);
	}
};


/**
	This class contains a list of vertex indices that correspond to the vertices that make up the polygon. This class should be populated by the
	class that reads in the mesh file, and then passes it to the GLMeshObject.
*/
class GLIndexedPoly  {
	friend class GLPolyIndexList;
	friend class GLPolyCategory;
	friend class GLMesh;
	friend class OBJReader;
//private:
public:
	DynamicArray<int> indexes;
	
	/**
		clear the indices
	*/
	void clear(){
		indexes.clear();
	}
		
	/**
		adding the index of a new vertex
	*/
	void addVertexIndex(int index){
		indexes.push_back(index);
	}

	/**
		a copy constructor.
	*/

	GLIndexedPoly& operator = (const GLIndexedPoly& other){
		this->indexes = other.indexes;
		return *this;
	}

	/**
		default constructor.
	*/
	GLIndexedPoly(void){
		indexes = DynamicArray<int>();
	}
	/**
		destructor.
	*/
	~GLIndexedPoly(void){
	}
};

class GLIndexedQuad : public GLIndexedPoly {
public:

	/**
		default constructor.
	*/
	GLIndexedQuad(int i1, int i2, int i3, int i4, bool flipNormal = false){
		if (flipNormal == false){
			addVertexIndex(i1); addVertexIndex(i2); addVertexIndex(i3); addVertexIndex(i4);
		}else{
			addVertexIndex(i1); addVertexIndex(i4); addVertexIndex(i3); addVertexIndex(i2);
		}
	}
	/**
		destructor.
	*/
	~GLIndexedQuad(void){
	}
};

class GLIndexedTriangle : public GLIndexedPoly {
public:

	/**
		default constructor.
	*/
	GLIndexedTriangle(int i1, int i2, int i3, bool flipNormal = false){
		if (flipNormal == false){
			addVertexIndex(i1); addVertexIndex(i2); addVertexIndex(i3);
		}else{
			addVertexIndex(i1); addVertexIndex(i3); addVertexIndex(i2);
		}
	}
	/**
		destructor.
	*/
	~GLIndexedTriangle(void){
	}
};

/**
	This class contains an array of tightly packed indeces of vertices that belong to polygons that have the same number of vertices. As an
	example the indices of all vertices that belong to traingles in the mesh, would be tightly packed together in one instance of this class. 
	There will be one of these objects for all the triangles in the mesh, one for all the quads, and so on. 
*/
class GLPolyIndexList  {
	friend class ParticleSystemLoader;
	friend class GLPolyCategory;
	friend class GLMesh;
	friend class CDMA_MechanicalComponent;
private:
	DynamicArray<unsigned int> indexList;		//a growing list of indexes
	unsigned int polyVertexCount;				//that's how many vertices belong to each polygon... list has listSize/polyVertexCount polygons
public:
	GLPolyIndexList(unsigned int vertCount){indexList = DynamicArray<unsigned int>(); this->polyVertexCount = vertCount;}
	GLPolyIndexList(const GLPolyIndexList &other): indexList(other.indexList), polyVertexCount(other.polyVertexCount){};
	~GLPolyIndexList(){
	}


	/**
		This method adds the indices of the vertices of the polygon p to indexList.
	*/
	int addPoly(GLIndexedPoly &p){
		if (p.indexes.size()!=this->polyVertexCount)
			return -1;
		for (uint i=0;i<p.indexes.size();i++)
			indexList.push_back(p.indexes[i]);
		return (indexList.size()%this->polyVertexCount)-1;
	}
};


/**
	finally, this class acts as a container of GLIndexLists, and everytime a new poly (GLIndexedPoly) is added,
	it checks to see if there is a list of the corect poly indexed count, and adds the data to that one, or it creates a new
	category.
*/
class GLPolyCategory  {
friend class ParticleSystemLoader;
friend class GLMesh;
friend class CDMA_MechanicalComponent;
private:
	DynamicArray<GLPolyIndexList*> categories;

	/**
		This method adds a new polygon to the category that has the same vertex count. If none exists then one is created.
	*/
	int addPoly(GLIndexedPoly &p){
		for (uint i=0;i<categories.size();i++){
			if (categories[i]->polyVertexCount == p.indexes.size())
				return categories[i]->addPoly(p);
		}
		categories.push_back(new GLPolyIndexList((unsigned int)p.indexes.size()));
		return categories[categories.size()-1]->addPoly(p);
	}

	/**
		constructor...
	*/
	GLPolyCategory(){categories = DynamicArray<GLPolyIndexList*>();};
	
	GLPolyCategory(const GLPolyCategory& other){
		categories.resize(other.categories.size());
		for(uint i=0;i<this->categories.size();i++)
			categories[i] = new GLPolyIndexList(*(other.categories[i]));
	}
	/**
		destructor - make sure we free the memory
	*/
	~GLPolyCategory(){
		for (uint i=0;i<categories.size();i++)
			delete categories[i];
	}

	/**
		This method returns the PolyIndexList whose polygons all have vertexCount vertices.
	*/
	GLPolyIndexList* getIndexList(int vertexCount){
		for (size_t i=0;i<categories.size();i++){
			if (categories[i]->polyVertexCount == vertexCount)
				return categories[i];
		}
		return NULL;
	}
};


/**
	This class is responsible with the storage and drawing of 3d static meshes. It is designed to work with OpenGL
	array lists in order to improve performance.
*/
class GLMesh   {
	friend class ParticleSystem;
	friend class ParticleSystemApp;
	friend class ParticleSystemLoader;
	friend class CDMA_MechanicalComponent;
private:
	//this is the array of vertex coordinates.
	DynamicArray<double> vertexList;
	//this array stores the normals for each vertex - NOTE: there can only be one for each vertex
	DynamicArray<double> normalList;
	//this array holds the texture coordinates - NOTE: there can only be one set of texture coordinates for each vertex
	DynamicArray<double> texCoordList;
	//this array holds the tangents for each vertex - NOTE: there can only be one for each vertex
	DynamicArray<double> tangentList;
	//this variable is used to store information related to which polys share each vertex. Only used for computation of normals
	DynamicArray<SharedVertexInfo*> sharedVertices;	

	//keep a list of the triangles in this mesh... they might need to be accessed later
	DynamicArray<GLIndexedTriangle> triangles;

	//this is the total number of vertices in the mesh
	int vertexCount;

	//this PolyCategory object holds the list of triangles, quads, etc in the mesh.
	GLPolyCategory *polygons;

	//variable that indicates whether or not the normals are to be used.
	bool useNormals;

	//this is the total number of polygons in the mesh.
	int nrPolys;

	//material
	GLShaderMaterial material;

	/**
	draws all the polygons of the mesh
	*/
	void drawMeshElements();

public:
	/**
		this is the default constructor
	*/
	GLMesh(void);

	GLMesh* clone();

	/**
		this is the destructor
	*/
	~GLMesh(void);

	/**
		this method removes all vertices and polygons from the mesh
	*/
	void clear();
	/**
		this method is used to add a new vertex to the mesh
	*/
	void addVertex(P3D &coords);

	/**
		this method also adds texture coordinates for the vertex that is to be added.
	*/
	void addVertex(P3D &coords, P3D &texCoords);

	/**
		this method add a prism to the current mesh.
	*/
	void addPrism(std::vector<P3D> &firstPoints, std::vector<P3D> &secondPoints, bool invertNormals = false);

	/**
		this method add a cylinder to the current mesh.
	*/
	void addCylinder(const P3D &p1, const P3D &p2, double radius, int resolution = 16);

	/**
		this method add a sphere to the current mesh.
	*/
	void addSphere(const P3D &center, double radius, int resolution = 16);

	/**
		this method add a cone to the current mesh.
	*/
	void addCone(P3D &basePoint, P3D &topPoint, double radius, int resolution = 16);

    /**
        this method increase the mesh density by split the face into smaller ones
    */
    GLMesh increaseMeshDensity();
	/**
		this method sets the coordinates for an existing vertex
	*/
	void setVertexCoordinates(int index, P3D &coords);

	/**
		This method sets the texture coordinates for an existing vertex.
	*/
	void setVertexTexCoordinates(int index, P3D &texCoords);
	
	/**
		This method sets the normal for an existing vertex.
	*/
	void setVertexNormal(int index, const V3D &normal);

	/**
		This method is used to compute the normals. The modifier passed in as a parameter should be either 1 or -1, and it
		should indicate if the polygons in this mesh have their vertices expressed in clockwise or anticlockwise order.
	*/
	void computeNormals(double modifier = 1);

	/**
		This method is used to compute the tangents from the UVs and vertex coordinates.
	*/
	void computeTangents();
	
	/**
		This is the method that adds new polygons to the mesh. The polygons have to be populated by the class that reads in the
		mesh from a file.
	*/
	void addPoly(GLIndexedPoly &p, bool duplicateVertices = false);

	/**
		returns true if the vertex is already indexed by any mesh poly, false otherwise
	*/
	bool vertexIsInUse(int vIndex);

	/**
		This method draws the model.
	*/
	void drawMesh();

	/**
		This method prints out the normals of the model - for testing purposes.
	*/
	void drawNormals();

	// Renders the mesh by appending it to the OBJ file that is passed in. The vertices are transformed using the provided orientation and translation	
	// vertexIdxOffset indicates the index of the first vertex for this object, making it possible to render multiple meshes to the same OBJ file
	// Returns the number of vertices written to the file
	uint renderToObjFile( FILE* fp, uint vertexIdxOffset, const Quaternion& rotation, const P3D& translation );


	/**
		This method returns the number of polygons in the current mesh.
	*/
	int getPolyCount() const {return nrPolys;}
    /**
    This method returns the number of triangles in the current mesh.
    */
    int getTriangleCount() const { return triangles.size(); }
	/**
		This method returns the number of vertices in the current mesh.
	*/
	int getVertexCount() const {
		assert(vertexList.size() % 3 == 0);
		return (int)vertexList.size()/3;
	}

	/**
		This method returns a reference to the dynamic array that stores the vertex positions.
	*/
	double* getVertexArray(){
		assert(!vertexList.empty());
		return &vertexList[0];
	}

	GLShaderMaterial &getMaterial() { return material; }
	const GLShaderMaterial &getMaterial() const { return material; }

	void setMaterial(const GLShaderMaterial &material) {
		this->material = material;
	}

	P3D getVertex(int i){
		return P3D(vertexList[3*i+0], vertexList[3*i+1], vertexList[3*i+2]);
	}

	P3D getTexCoordinates(int i){
		assert(texCoordList.size() > 0);
		assert(texCoordList.size() % 3 == 0);

		return P3D(texCoordList[3*i+0], texCoordList[3*i+1], texCoordList[3*i+2]);
	}

	V3D getNormal(int i) {
		assert(normalList.size() > 0);
		assert(normalList.size() % 3 == 0);

		return V3D(normalList[3*i+0], normalList[3*i+1], normalList[3*i+2]);
	}

	void setVertex(int i, P3D p){
		vertexList[3*i+0] = p[0];
		vertexList[3*i+1] = p[1];
		vertexList[3*i+2] = p[2];
	}

	P3D getCenterOfMass()
	{
		P3D ans(0, 0, 0);
		for (int k = 0; k < 3; ++k)
			for (int i = 0; i < getVertexCount(); ++i)
			{
				ans[k] += vertexList[3 * i + k] / getVertexCount();
			}
		return ans;
	}

	void translate(const P3D &offset)
	{
		for (int k = 0; k < 3; ++k)
			for (int i = 0; i < getVertexCount(); ++i)
			{
				vertexList[3 * i + k] += offset[k];
			}
	}

	void rotate(const Quaternion &q, const P3D &center)
	{
		for (int i = 0; i < getVertexCount(); ++i)
		{
			V3D tmp = q.rotate(V3D(vertexList[3 * i + 0] - center[0], 
				vertexList[3 * i + 1] - center[1], vertexList[3 * i + 2] - center[2]));
			for (int k = 0; k < 3; ++k)
				vertexList[3 * i + k] = tmp[k] + center[k];
		}
	}

	void scale(double scale, const P3D &center)
	{
		for (int i = 0; i < getVertexCount(); ++i)
			for (int k = 0; k < 3; ++k)
			{
				vertexList[3 * i + k] = (vertexList[3 * i + k] - center[k]) * scale + center[k];
			}
	}

	GLIndexedTriangle &getTriangle(int ind)
	{
		//return polygons[0].getIndexList
		return triangles[ind];
	}

//	void applyTransformationMatrix(const TransformationMatrix &Rt);

};
