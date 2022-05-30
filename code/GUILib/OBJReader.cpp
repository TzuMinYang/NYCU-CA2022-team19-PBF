/**
	TODO: test this code with texture coordinates. The part that duplicates vertices has not been very well tested yet.
*/

#include <GUILib/OBJReader.h>
#include <Utils/Utils.h>
#include <Utils/Logger.h>

OBJReader::OBJReader(void){
}

OBJReader::~OBJReader(void){
}

/**
	This method checks a line of input from an ASCII OBJ file, and determines its type.
*/
int OBJReader::getLineType(char* line){
	if (line == NULL || strlen(line)<3)
		return NOT_IMPORTANT;
	if (line[0] == 'v' && line[1] == ' ')
		return VERTEX_INFO;
	if (line[0] == 'v' && line[1] == 't' && line[2] == ' ')
		return TEXTURE_INFO;
	if (line[0] == 'v' && line[1] == 'n' && line[2] == ' ')
		return NORMAL_INFO;
	if (line[0] == 'f' && line[1] == ' ')
		return FACE_INFO;
	return NOT_IMPORTANT;
}

/**
	This method reads three double values separated by white space and returns a point3d populated with the information. It is assumed that the values
	are at the very beginning of the character array.
*/
P3D OBJReader::readCoordinates(char* line){
	P3D p;
	// int n = sscanf(line, "%lf %lf %lf", &p.x, &p.y, &p.z);
	sscanf(line, "%lf %lf %lf", &p[0], &p[1], &p[2]);
	return p;
}

/**
	This method takes a character array and does the following:
	if a number is encountered, that number is read - this will be the vertex index. If a '/' is encountered then the number after it is also read (if any).
	The method then returns the pointer to the first white space after the last / encountered from the beginning of the character array. If there were
	no valid indices to read, this method returns NULL. The flag that is passed in is set to reflect which indices were read.
*/
char* OBJReader::getNextIndex(char* line, int &vertexIndex, int &texcoordIndex, int &normalIndex, int &flags){
	flags = 0;
	//skip the white spaces...
	while (*line == ' ')
		line++;
	char* result = line;
	if (sscanf(line, "%d",&vertexIndex)!=1)
		return NULL;
	flags |= READ_VERTEX_INDEX;
	//check for the '/'
	while (*result != '\0' && *result != ' ' && *result != '/')
		result++;
	//did we find a '/'?
	if (*result == '/')
	{
		if (*(result+1) != '/' && sscanf(result+1, "%d",&texcoordIndex)==1){
			flags |= READ_TEXTCOORD_INDEX;
			result++;
		}
		else if(*(result+1) == '/')
		{
			result++;
		}

		//check for another '/'
		while (*result != '\0' && *result != ' ' && *result != '/')
			result++;
		if (*result == '/')
		{
			if (sscanf(result+1, "%d", &normalIndex) == 1) {
				flags |= READ_NORMAL_INDEX;
				result++;
			}
		}
	}

	//and now move result to the next available index to be read.
	while (*result!='\0' && *result!=' ')
		result++;

	return result;
}


/**
	This static method reads an obj file, whose name is sent in as a parameter, and returns a pointer to a GLMesh object that it created based on the file information.
	This method throws errors if the file doesn't exist, is not an obj file, etc.
*/
GLMesh* OBJReader::loadOBJFile(const char* fileName){
	if (fileName == NULL)
		throwError("fileName is NULL.");
	
//	Logger::out()<< "Loading mesh: " << fileName <<std::endl;

	//have a temporary buffer used to read the file line by line...
	char buffer[200];

	//an array of vertex positions - need to store it initially because vertices might have to be inserted/duplicated, and that may mess up the indices of polygons defined afterwards
	DynamicArray<P3D> vertexCoordinates;
	//and this is an array of texture coordinates
	DynamicArray<P3D> texCoordinates;
	//this is an array of normals
	DynamicArray<V3D> normals;
	//this is to store whether we assigned a set of uvs to vertex i - if yes, and we come across another set of UVs, we will need to
	//duplicate this vertex.
	DynamicArray<int> texCoordsOnVertex;
	//... and the same for normals
	DynamicArray<int> normalOnVertex;

	bool readNormalsFromObj = false;
	bool readTextureCoordinatesFromObj = false;

	GLMesh* result = new GLMesh();

	//pass 1: read all vertices, texture coordinates and normals...
	FILE* f = fopen(fileName, "r");
	if (f == NULL)
		throwError("Cannot open file \'%s\'.", fileName);

	//this is where it happens.
	while (!feof(f)){
		//get a line from the file...
		fgets(buffer, 200, f);
		//see what line it is...
		int lineType = getLineType(buffer);
		if (lineType == VERTEX_INFO){
			//we need to read in the three coordinates - skip over the v
			vertexCoordinates.push_back(readCoordinates(buffer + 1));
			result->addVertex(vertexCoordinates[vertexCoordinates.size()-1]);
			texCoordsOnVertex.push_back(-1);
			normalOnVertex.push_back(-1);
		}

		if (lineType == TEXTURE_INFO){
			texCoordinates.push_back(readCoordinates(buffer + 2));
			readTextureCoordinatesFromObj = true;
		}

		if (lineType == NORMAL_INFO){
			normals.push_back(V3D(readCoordinates(buffer + 2)));
			readNormalsFromObj = true;
		}

		if (lineType == FACE_INFO) {} //skip it this time around... we need to collect all vertices first so that we don't mess up indexing order...
	}
	fclose(f);

	//pass 2: read all faces and create polys
	f = fopen(fileName, "r");
	if (f == NULL)
		throwError("Cannot open file (again) \'%s\'.", fileName);

	//this variable will keep getting populated with face information
	GLIndexedPoly temporaryPolygon;

	//this is where it happens.
	while (!feof(f)) {
		//get a line from the file...
		fgets(buffer, 200, f);
		//see what line it is...
		int lineType = getLineType(buffer);
		if (lineType == VERTEX_INFO) {} //skip, already processed...

		if (lineType == TEXTURE_INFO) {} //skip, already processed...

		if (lineType == NORMAL_INFO) {} //skip, already processed...

		if (lineType == FACE_INFO) {
			temporaryPolygon.indexes.clear();
			int vIndex, tIndex, nIndex;
			int flag;
			char* tmpPointer = buffer + 1;
			while (tmpPointer = getNextIndex(tmpPointer, vIndex, tIndex, nIndex, flag)) {
				if (flag & READ_TEXTCOORD_INDEX) {
					if (tIndex < 0)
						tIndex = (int)texCoordinates.size() + tIndex;

					if (texCoordsOnVertex[vIndex - 1] >= 0 && texCoordsOnVertex[vIndex - 1] != tIndex - 1) {
						// We have an additional set of texcoords for this vertex - will duplicate it 
						vertexCoordinates.push_back(vertexCoordinates[vIndex-1]);
						result->addVertex(vertexCoordinates[vertexCoordinates.size() - 1]);
						// copy over the old normal, as it might have been set properly already...
						result->setVertexNormal(vertexCoordinates.size() - 1, result->getNormal(vIndex - 1));
						// have a way of tracking the normal/texture coords associated with this new vertex...
						normalOnVertex.push_back(normalOnVertex[vIndex - 1]);
						texCoordsOnVertex.push_back(texCoordsOnVertex[vIndex - 1]);
						// and update the index of this new vertex, keeping the convention that was there before - vIndex is 1-relative
						vIndex = vertexCoordinates.size();
					}
					result->setVertexTexCoordinates(vIndex - 1, texCoordinates[tIndex - 1]);
					texCoordsOnVertex[vIndex - 1] = tIndex - 1;
				}

				if (flag & READ_NORMAL_INDEX) {
					if (nIndex < 0)
						nIndex = (int)normals.size() + nIndex;

					if (normalOnVertex[vIndex - 1] >= 0 && normalOnVertex[vIndex - 1] != nIndex - 1) {
						// We have an additional set of texcoords for this vertex - will duplicate it 
						vertexCoordinates.push_back(vertexCoordinates[vIndex - 1]);
						result->addVertex(vertexCoordinates[vertexCoordinates.size() - 1]);
						// copy over the old texture info, as it might have been set properly already...
						P3D texCoords = result->getTexCoordinates(vIndex - 1);
						result->setVertexTexCoordinates(vertexCoordinates.size() - 1, texCoords);
						// have a way of tracking the normal/texture coords associated with this new vertex...
						normalOnVertex.push_back(normalOnVertex[vIndex - 1]);
						texCoordsOnVertex.push_back(texCoordsOnVertex[vIndex - 1]);
						// and update the index of this new vertex, keeping the convention that was there before - vIndex is 1-relative
						vIndex = vertexCoordinates.size();
					}
					result->setVertexNormal(vIndex - 1, normals[nIndex - 1]);
					normalOnVertex[vIndex - 1] = nIndex - 1;
				}

				temporaryPolygon.indexes.push_back(vIndex - 1);
			}
			if (temporaryPolygon.indexes.size() == 0)
				Logger::consolePrint("Found a polygon with zero vertices.\n");
			else
				result->addPoly(temporaryPolygon);
		}

	}

	fclose(f);

	if (readNormalsFromObj == false)
		result->computeNormals();
	if (readTextureCoordinatesFromObj == true)
		result->computeTangents();

 	return result;
}




