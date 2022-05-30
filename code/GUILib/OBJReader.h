#pragma once

#include <GUILib/GLMesh.h>
#include <MathLib/P3D.h>

/*======================================================================================================================================================================*
 | This class implements the routines that are needed to load a mesh from an OBJ file. This class loads a GLMesh with the polygonal mesh that is stored in a file.      |
 | Only vertex coordinates, vertex texture coordinates and connectivity information are loaded from the file.                                                           |
 *======================================================================================================================================================================*/
class OBJReader{
public:
	OBJReader(void);
	~OBJReader(void);

	/**
		This static method reads an obj file, whose name is sent in as a parameter, and returns a pointer to a GLMesh object that it created based on the file information.
		This method throws errors if the file doesn't exist, is not an obj file, etc.
	*/
	static GLMesh* loadOBJFile(const char* fileName);

protected:
	static char* getNextIndex(char* line, int &vertexIndex, int &texcoordIndex, int &normalIndex, int &flags);
	static P3D readCoordinates(char* line);
	static int getLineType(char* line);
	
	static const int VERTEX_INFO = 0;
	static const int FACE_INFO = 1;
	static const int TEXTURE_INFO = 2;
	static const int NORMAL_INFO = 3;
	static const int NOT_IMPORTANT = 5;
	
	static const int READ_VERTEX_INDEX = 0x01;
	static const int READ_TEXTCOORD_INDEX = 0x02;
	static const int READ_NORMAL_INDEX = 0x04;
};
