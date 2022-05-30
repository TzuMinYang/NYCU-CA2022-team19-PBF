#pragma once

/*==================================================================================================================================================*
 | This abstract class presents the methods that need to be implemented by the extending classes - classes which will be used for imageIO using     |
 | different image formats.                                                                                                                         |
 *==================================================================================================================================================*/

#define RGB_MODEL 0
#define RGBA_MODEL 1


#include "Image.h"
#include <string>

class ImageIO  {
protected:
	//the file name of the image to be written or read
	std::string m_fileName;

public:
	ImageIO(const char* fileName){
		this->m_fileName = fileName;
	}
	virtual ~ImageIO(void){
	}

	/**
		This method should be implemented by the extending classes - it will load an image from the given file.
		The imageModel indicates whether the image will be an RGB or RGBA image.
	*/
	virtual Image* loadFromFile(int imageModel = RGB_MODEL) = 0;
	/**
		This method should be implemented by the extending classes - it will save an image in a specific file format.
	*/
	virtual void writeToFile(Image* img) = 0;
};
